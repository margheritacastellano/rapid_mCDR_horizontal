# La cascade : pourquoi le champ partagé ne peut pas marcher

Ce document explique **la seule vraie différence de conception** entre le modèle d'Océane
(`rapidmCDR_ALK_3D.ipynb`) et le nôtre (`test.ipynb`), et pourquoi cette différence est
nécessaire — pas cosmétique.

> **En une phrase** : Océane stocke l'alcalinité dans un **champ partagé** et doit la déplacer
> « à la main » quand une particule bouge ; nous l'attachons **à la particule**, et elle suit
> toute seule. Le champ partagé n'a aucune notion de **propriétaire**, et c'est fatal dès qu'il
> y a plus d'une particule.

---

## 1. Les deux conceptions, côte à côte

### Chez Océane : l'alcalinité vit dans le champ

```
   ALK0[profondeur, y, x]          ← UN grand tableau partagé (408 Mo)
        ▲
        │  la particule n'est qu'un MARQUEUR : elle ne porte rien
        ●  P

   Quand P change de maille, il faut ALLER DÉPLACER la colonne à la main :

        ALK0[:, nouvelle] += ALK0[:, ancienne]      ← HorizontalAdvConcentration
        ALK0[:, ancienne]  = 0
```

D'où ses **3 kernels** :

| kernel | rôle |
|---|---|
| `VerticalConcentration` | advection verticale sur la colonne (dans le champ) |
| `HorizontalAdvConcentration` | **recoller** la colonne à la particule quand elle bouge |
| `HorizontalAdvParticle` | déplacer la particule (RK4) |

### Chez nous : l'alcalinité vit sur la particule

```
   pcol[p]  = [50 valeurs]         ← LA colonne de la particule p (1 550 nombres au total)
        ▲
        │  la position de la particule EST la position de la colonne
        ●  P  (pidx = p)

   Quand P bouge... il n'y a RIEN À FAIRE. La colonne suit.
```

D'où nos **2 kernels** :

| kernel | rôle |
|---|---|
| `VerticalAdvDiffImplicit` | advection **+ diffusion** verticales, sur `pcol[pidx]` |
| `parcels.AdvectionRK4` | déplacer la particule (le 2D) |

Le kernel `HorizontalAdvConcentration` **n'a plus d'objet**. Ce n'est pas qu'on l'a « optimisé » :
il n'existe plus parce que le problème qu'il résolvait n'existe plus.

---

## 2. Le problème de fond : le champ n'a pas de propriétaire

Regardons **la** ligne critique d'Océane :

```python
if xi != particle.cell_x or yi != particle.cell_y:
    ALK0[:, yi, xi] += ALK0[:, particle.cell_y, particle.cell_x]
    ALK0[:, particle.cell_y, particle.cell_x] -= ALK0[:, particle.cell_y, particle.cell_x]
```

Traduction : *« je bouge, donc j'emporte **tout ce qui se trouve** dans mon ancienne maille. »*

Et c'est là le drame :

> ### 🔑 Une maille contient **un nombre**. Pas « la dose de P0 » et « la dose de P1 ».
> ### Juste un nombre. Impossible de savoir à qui il appartient.

Donc quand une particule dit « j'emporte ma colonne », elle emporte en réalité **la colonne de
tout le monde**. Elle ne peut pas faire autrement : l'information « à qui appartient cette
alcalinité ? » **n'existe nulle part** dans le champ.

---

## 3. Le mécanisme, pas à pas

Le point crucial : les particules avancent d'environ **0,2°/jour**, alors qu'une maille fait
environ **0,33°**. Une particule ne change donc de maille que tous les **~1,5 jours**.
Conséquence : **deux particules se retrouvent régulièrement dans la même maille.**

Déroulons. P0 est partie hier, P1 vient de naître. Elles suivent le **même** courant.

```
   JOUR n
   ────────────────────────────────────────────────────
   maille :    3        4        5        6
   contenu:    0      dose1    dose0      0
   position:          ● P1     ● P0
```

**JOUR n+1** — P1 avance dans la maille 5, mais P0 **n'a pas encore changé de maille** (elle
attend son 1,5ᵉ jour) :

```
   P0 : xi = 5 = cell_x  →  elle NE BOUGE PAS  →  rien ne se passe
   P1 : cell_x = 4, xi = 5  →  ELLE BOUGE :

        ALK0[:, 5] += ALK0[:, 4]        ← P1 emporte "sa" dose vers la maille 5...
        ALK0[:, 4]  = 0                   ...mais la maille 5 contenait DÉJÀ dose0 !
```

```
   maille :    3        4        5             6
   contenu:    0        0    dose0+dose1       0        ← 💥 FUSIONNÉES
   position:                 ● P1  ● P0
```

**Et c'est irréversible.** Les deux doses ne sont plus qu'**un seul nombre**. La prochaine fois
que P0 **ou** P1 bouge, elle emportera **les deux**. Elles ne se sépareront **jamais plus**.

Répétez sur 31 jours, avec 31 particules alignées sur la même trajectoire : elles se rattrapent
sans cesse, fusionnent à chaque fois, et l'ensemble finit **empilé sur la maille de tête**.

> C'est ce qu'on appelle la **cascade** : les doses sont relayées vers l'avant, de particule en
> particule, et s'agglomèrent au bout de la traînée.

---

## 4. La preuve numérique

J'ai reproduit **la logique exacte** d'Océane sur un modèle jouet (particules alignées, 1 dose
chacune, elles avancent), en comparant les deux conceptions :

```
--- 1 particule ---
   champ PARTAGÉ (Océane)  : [0 0 0 0 0 0 0 1 0 0 0 0]     max = 1 dose / maille
   colonnes PAR PARTICULE  : [0 0 0 0 0 0 0 1 0 0 0 0]     max = 1 dose / maille
   →  IDENTIQUES ✅

--- 2 particules ---
   champ PARTAGÉ           : [0 0 0 0 0 0 0 2 0 0 0 0]     max = 2
   colonnes PAR PARTICULE  : [0 0 0 0 0 0 0 2 0 0 0 0]     max = 2
   →  IDENTIQUES ✅

--- 8 particules ---
   champ PARTAGÉ           : [0 0 0 0 0 0 0 8 0 0 0 0]     max = 8   ← 💥 TOUT sur 1 maille
   colonnes PAR PARTICULE  : [0 0 0 0 2 2 2 2 0 0 0 0]     max = 2   ← ✅ réparties
   →  DIVERGENCE
```

Et sur le **vrai run** (31 jours, données ACC réelles) :

| | champ partagé | colonnes par particule |
|---|---|---|
| mailles touchées | 18 | 18 |
| doses par maille | 17 mailles à **1 dose**… **+ 1 maille à 14 doses** 💥 | max **2 doses** (physique : 2 particules peuvent vraiment coexister) |
| masse totale | 31 × dose ✅ | 31 × dose ✅ |

> ⚠️ **Le piège** : la **masse est conservée dans les deux cas** ! Rien ne « disparaît ». C'est
> pour ça que le bug passe sous le radar : la vérification de conservation, elle, est **verte**.
> L'alcalinité est juste **au mauvais endroit**.

---

## 5. Pourquoi c'est très facile à rater 🎯

Regarde le tableau du §4 : **avec 1 ou 2 particules, les deux conceptions donnent exactement
le même résultat.**

Le bug n'apparaît qu'à partir du moment où :
1. il y a **plusieurs** particules (⟹ l'injection **quotidienne**), **et**
2. elles se **rattrapent** (⟹ elles suivent la même trajectoire — ce qui est exactement le cas :
   même source, mêmes courants)

Autrement dit : on développe et on valide avec **une** particule → tout marche. On active
l'injection quotidienne → le résultat est **silencieusement** faux, sans erreur, sans NaN, sans
perte de masse. Juste une traînée qui s'agglomère au mauvais endroit.

**Ce n'est pas une faute d'inattention.** C'est une limite intrinsèque de la conception
« champ partagé + marqueurs », qui n'apparaît que dans un régime précis.

---

## 6. Notre solution

On donne à chaque particule **sa propre colonne** :

```python
fieldset.pcol = np.zeros((last_time, 50))    # 31 lignes × 50 niveaux
pidx = parcels.Variable('pidx', dtype=np.int32)   # ← le "numéro de casier" de la particule
```

```
   pcol
   ┌──────────────────────────────┐
   │ ligne 0  → [50 niveaux] │ ← la dose du jour 0,  portée par P0 (pidx=0)
   │ ligne 1  → [50 niveaux] │ ← la dose du jour 1,  portée par P1
   │   ...                    │
   │ ligne 30 → [50 niveaux] │ ← la dose du jour 30
   └──────────────────────────────┘
```

Maintenant l'information **« à qui appartient cette alcalinité ? »** existe : c'est l'indice de
ligne. Deux particules dans la même maille ? Aucun problème : leurs colonnes sont dans **deux
lignes différentes**. Elles se croisent, puis se séparent — comme dans la réalité.

Le champ `ALK0` n'est plus qu'un **diagnostic**, reconstruit à la demande pour les tracés :

```python
fieldset.ALK0.data[:] = 0.0
for p in pset:
    fieldset.ALK0.data[0, :, p.cell_y, p.cell_x] += fieldset.pcol[p.pidx]
```

C'est une **addition** (`+=`) : si deux particules sont dans la même maille, leurs doses
s'additionnent **à l'affichage** — ce qui est correct — sans jamais fusionner **dans le modèle**.

---

## 7. Le renversement conceptuel

| | Océane | Nous |
|---|---|---|
| **l'état du modèle** | le champ `ALK0` (408 Mo) | `pcol` (**1 550 nombres**) |
| **la particule** | un marqueur (ne porte rien) | **porte sa colonne** |
| **`ALK0`** | *la vérité* | un **diagnostic** reconstruit |
| **transport horizontal de la concentration** | un kernel dédié qui déplace les colonnes | **rien à faire** (elle suit la particule) |
| **notion de propriétaire** | ❌ inexistante | ✅ `pidx` |
| **nb de kernels** | 3 | 2 |

> **L'idée en une image** : Océane peint l'alcalinité **sur le sol** et doit la re-peindre à
> chaque fois que quelqu'un se déplace. Nous la mettons **dans le sac à dos** de chacun — le sac
> suit son propriétaire, sans effort, et deux personnes peuvent se croiser sans mélanger leurs sacs.

---

## 8. Pour être juste

- Sa conception est **le réflexe naturel** : on veut un champ eulérien en sortie, donc on stocke
  dans le champ. C'est ce que fait n'importe quel modèle eulérien.
- Le bug est **invisible** en dessous de ~3 particules, **ne casse rien** (pas de NaN, pas
  d'explosion) et **conserve la masse**. Toutes les vérifications évidentes passent.
- Son `VerticalConcentration` (advection verticale explicite) est **correct** en tant que tel —
  on l'a étendu (diffusion) et rendu implicite pour des raisons de **stabilité** (D ≈ 7,7),
  ce qui est un problème **indépendant** de la cascade.
- Enfin : **trois autres bugs** empêchaient de toute façon les particules d'avancer (lâcher à 0 m
  hors-bornes, `lon_nextloop` non initialisé, particules ajoutées sans horodatage). Avec des
  particules figées, la cascade ne peut même pas se manifester — tout s'empile simplement à la
  source. Les deux problèmes se masquaient mutuellement.

---

## 9. Ce qu'il faut retenir

1. **La cascade n'est pas un bug d'implémentation, c'est une limite de conception.** On ne peut
   pas la corriger en réparant `HorizontalAdvConcentration` : l'information manquante
   (le propriétaire) n'existe nulle part dans un champ partagé.
2. **La conservation de la masse ne prouve rien.** Elle est verte dans les deux cas. Il faut
   regarder **la répartition** (doses par maille), pas seulement le total.
3. **Le lagrangien pris au sérieux** : si on choisit de suivre des particules, alors la quantité
   transportée doit **vivre sur la particule**. À moitié lagrangien (particules) et à moitié
   eulérien (le champ), on hérite des inconvénients des deux.

---

*Voir `SCHEMA_NUMERIQUE.md` pour les maths du schéma vertical, et `ALGO_3a_3b.md` pour le détail
des cellules de simulation.*
