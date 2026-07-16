# Le schéma numérique : maths & algorithme

Documentation du kernel `VerticalAdvDiffImplicit` (cellule **[12]** de `test.ipynb`) :
**advection + diffusion verticales**, en **volumes finis 1D**, **Euler implicite**.

> Notation : tout est écrit en Unicode (pas de LaTeX), lisible en brut comme en aperçu.

---

## 1. Le problème physique

On transporte **c = la perturbation d'alcalinité** (`ALK0`, en µM), qui part de zéro partout.
Ce n'est ni du CO₂ ni du DIC : c'est un **traceur passif**.

```
   ∂c/∂t  +  ∂(w·c)/∂z  −  ∂( K(z)·∂c/∂z )/∂z  =  f_source
            └─ advection ┘   └──── diffusion ────┘
```

| symbole | sens | dans le code |
|---|---|---|
| c | perturbation d'alcalinité (µM) | `C` = `fieldset.pcol[pidx]` |
| w | vitesse verticale (m/s, **> 0 vers le HAUT**) | `fieldset.W` |
| K(z) | diffusivité verticale (m²/s) | `K_diff_vertical_v` |
| f_source | injection en surface | `alkalinity_forcing` |

⚠️ La diapo 8 de la présentation écrit `∂(K·c)/∂z` — il **manque un ∂/∂z**. La bonne forme
(celle de l'Algorithme 2 et du rapport) est `∂(K·∂c/∂z)/∂z`.

---

## 2. Le maillage vertical (grille Arakawa C)

Les données fournissent **deux** grilles verticales distinctes :

```
      surface
  ═══════════════  ← interface k       (face,   RF[k])
        • cₖ                            (centre, RC[k])      ┐ hₖ
  ───────────────  ← interface k+1  wₖ₊₁ (face, RF[k+1])     ┘
        • cₖ₊₁                          (centre, RC[k+1])
  ───────────────  ← interface k+2
         ⋮
  ═══════════════  ← interface Nz      (fond)
```

- **Traceurs** (c) → aux **centres** : `RC = ALK.grid.depth` = 5, 15, 25, …
- **Vitesse verticale** (w) → aux **faces** : `RF = W.grid.depth` = 0, 10, 20, …

**Convention** : k croît vers le **bas**. L'**interface k** est la face du **haut** de la maille k
(elle sépare la maille k−1 au-dessus de la maille k en dessous). Donc `wₖ` s'utilise
**directement**, sans moyenne.

Les **deux** longueurs dont les volumes finis ont besoin :

```
   hₖ = 2·(RC[k] − RF[k])        épaisseur de la maille k     ← le "volume", pour DIVISER le bilan
   δₖ = RC[k] − RC[k−1]          distance centre-à-centre     ← pour le GRADIENT
```
```python
h_fv     = 2.0 * (RC - RF)   # (50,)
delta_fv = np.diff(RC)       # (49,)  -> 50 mailles => 49 interfaces intérieures
```

> ⚠️ **hₖ ≠ δₖ** dès que le maillage est irrégulier — et il l'est (10 m en surface → 450 m au fond).

---

## 3. Les flux aux interfaces

Le bilan sur la maille k : **variation = ce qui entre − ce qui sort**.

### 3.a Advection — décentrée amont (*upwind*)

Comme **w > 0 pointe vers le haut**, on sépare les deux sens :

```
   a↓ᵢ = max(−wᵢ, 0)     vitesse DESCENDANTE à l'interface i
   a↑ᵢ = max(+wᵢ, 0)     vitesse MONTANTE
```

**La règle upwind** : à l'interface i, on transporte toujours la maille **en amont du courant**.
- ça **descend** (a↓ᵢ > 0) → on prend `cᵢ₋₁` (celle du **dessus**)
- ça **monte** (a↑ᵢ > 0) → on prend `cᵢ` (celle du **dessous**)

### 3.b Diffusion — loi de Fick

Le gradient à l'interface i = **différence à 2 points entre les centres**, divisée par δᵢ :

```
   (∂c/∂z)|ᵢ  ≈  (cᵢ − cᵢ₋₁) / δᵢ
```

et le K de la maille est ramené à l'interface par une **moyenne** :

```
   Kᵢ(interface) = (K[i−1] + K[i]) / 2
```

d'où le flux diffusif :

```
   F_dif(i) = −Kᵢ · (∂c/∂z)|ᵢ  =  (K[i−1]+K[i]) / (2·δᵢ) · (cᵢ₋₁ − cᵢ)
```

---

## 4. ⭐ Le schéma, forme développée

En assemblant les deux interfaces de la maille k (celle du **haut** = k, celle du **bas** = k+1) :

```
  hₖ · (cₖⁿ⁺¹ − cₖⁿ) / Δt  =

      cₖ₋₁ⁿ⁺¹ · [  a↓ₖ  +  (Kₖ₋₁+Kₖ)/(2·δₖ)  ]                        ← ce qui ENTRE par le haut

    + cₖⁿ⁺¹   · [ −a↑ₖ   −  (Kₖ₋₁+Kₖ)/(2·δₖ)                          ← ce qui SORT par le haut
                  −a↓ₖ₊₁ −  (Kₖ+Kₖ₊₁)/(2·δₖ₊₁) ]                      ← ce qui SORT par le bas

    + cₖ₊₁ⁿ⁺¹ · [  a↑ₖ₊₁ +  (Kₖ+Kₖ₊₁)/(2·δₖ₊₁) ]                      ← ce qui ENTRE par le bas
```

**Lecture ligne par ligne** — chaque `a` apparaît **une fois en + (entrée)** et **une fois en − (sortie)** :

| terme | interface | sens | effet sur la maille k |
|---|---|---|---|
| `a↓ₖ` | haut (k) | ça descend | **entrée** depuis cₖ₋₁ → **+** |
| `a↑ₖ` | haut (k) | ça monte | **sortie** de cₖ vers le haut → **−** |
| `a↓ₖ₊₁` | bas (k+1) | ça descend | **sortie** de cₖ vers le bas → **−** |
| `a↑ₖ₊₁` | bas (k+1) | ça monte | **entrée** depuis cₖ₊₁ → **+** |

👉 C'est **exactement** ça qui donne la **conservation** : ce qui sort de k entre dans k±1.

### La règle de signe (vérification en 10 secondes) 🎯

```
   coefficients des VOISINS (cₖ₋₁, cₖ₊₁)  :  ≥ 0     ← on ne peut que RECEVOIR d'un voisin
   coefficient de cₖ (la diagonale)        :  ≤ 0     ← on ne peut que PERDRE de chez soi
```

**Test physique** : K = 0 et subsidence partout (a↓ = |w| > 0, a↑ = 0) :
```
   hₖ dcₖ/dt = a↓ₖ·cₖ₋₁ − a↓ₖ₊₁·cₖ
```
→ la maille reçoit du **dessus**, perd vers le **dessous** → le traceur **descend** ✅

---

## 5. Pourquoi **implicite** ?

Un schéma **explicite** est stable seulement si le **nombre de diffusion** vérifie

```
   D = K·Δt / Δz²  ≤  1/2
```

Or près de la surface :

```
   D = (8,91e−3 × 86400) / 10²  ≈  7,7   >>>  0,5
```

👉 **L'explicite diverge** à Δt = 1 jour. Deux issues : réduire Δt à ~1,5 h (16× plus de pas,
et ça casse la logique d'injection quotidienne), ou **passer en implicite** ← le choix fait ici.

En **Euler rétrograde**, on évalue tous les flux à l'instant **n+1** (l'inconnue) → **stable
quel que soit Δt**.

---

## 6. Du schéma à la matrice tridiagonale

Posons les **4 coefficients** de la maille k :

```
   Aₖ = a↓ₖ   + (Kₖ₋₁+Kₖ)/(2·δₖ)        ← coefficient de cₖ₋₁
   Bₖ = a↑ₖ   + (Kₖ₋₁+Kₖ)/(2·δₖ)        ← perte par le haut
   Cₖ = a↓ₖ₊₁ + (Kₖ+Kₖ₊₁)/(2·δₖ₊₁)      ← perte par le bas
   Dₖ = a↑ₖ₊₁ + (Kₖ+Kₖ₊₁)/(2·δₖ₊₁)      ← coefficient de cₖ₊₁
```

Le schéma s'écrit alors de façon compacte :

```
   hₖ·(cₖⁿ⁺¹ − cₖⁿ)/Δt  =  Aₖ·cₖ₋₁ⁿ⁺¹  −  (Bₖ + Cₖ)·cₖⁿ⁺¹  +  Dₖ·cₖ₊₁ⁿ⁺¹
```

On multiplie par **rₖ = Δt/hₖ** et on passe les inconnues à gauche :

```
   −rₖ·Aₖ·cₖ₋₁ⁿ⁺¹  +  [1 + rₖ·(Bₖ + Cₖ)]·cₖⁿ⁺¹  −  rₖ·Dₖ·cₖ₊₁ⁿ⁺¹  =  cₖⁿ
   └────┬────┘        └──────────┬──────────┘       └────┬────┘        └┬┘
   sous-diagonale         diagonale                sur-diagonale     2ᵉ membre
```

Chaque ligne ne touche que cₖ₋₁, cₖ, cₖ₊₁ → la matrice **A** est **tridiagonale**, et on résout

```
   A · cⁿ⁺¹ = cⁿ
```

### Conditions aux limites : **flux nul**
- **Surface** : rien ne traverse l'interface 0 → `A₀ = B₀ = 0`
- **Fond** : rien ne traverse l'interface Nz → `C_{Nz−1} = D_{Nz−1} = 0`

### Le raccourci du code (L et U)
Le code factorise en **deux** tableaux indexés par **interface** plutôt que par maille :

```
   Lᵢ = a↓ᵢ + Kᵢ/δᵢ        ← coefficient de la maille du DESSUS
   Uᵢ = a↑ᵢ + Kᵢ/δᵢ        ← coefficient de la maille du DESSOUS
```
La correspondance est directe :
```
   Aₖ = Lₖ        Bₖ = Uₖ        Cₖ = Lₖ₊₁        Dₖ = Uₖ₊₁
```
C'est la même chose, écrite une seule fois par interface (au lieu de deux fois par maille).

---

## 7. ⭐ L'algorithme

Pour **une particule**, sur **un pas de temps** (Δt = 1 jour) :

```
ENTRÉES : cⁿ = pcol[pidx]   (la colonne portée par cette particule)
          la position de la particule
SORTIE  : cⁿ⁺¹ → réécrit dans pcol[pidx]

 1. LOCALISER la particule            → maille (xi, yi)         [l. 26–39]
 2. LIRE la colonne de vitesse w      → W = fieldset.W[:, yi, xi]  [l. 43]
 3. NETTOYER w  (NaN / valeurs de remplissage sous le plancher → 0)   [l. 49–50]
 4. INJECTER  si c'est le jour du lâcher :  c₀ += f·Δt ; released = 1  [l. 53–55]

 5. POUR chaque interface i = 1 … Nz−1 :                          [l. 58–66]
        a↓ᵢ = max(−wᵢ, 0)
        a↑ᵢ = max(+wᵢ, 0)
        Kdᵢ = (K[i−1] + K[i]) / (2·δᵢ)
        Lᵢ  = a↓ᵢ + Kdᵢ
        Uᵢ  = a↑ᵢ + Kdᵢ

 6. POUR chaque maille k = 0 … Nz−1 : construire la ligne k       [l. 69–78]
        rₖ        = Δt / hₖ
        sous-diag = −rₖ · Lₖ                (0 si k = 0)
        diag      = 1 + rₖ · (Uₖ + Lₖ₊₁)    (Uₖ = 0 si k = 0 ; Lₖ₊₁ = 0 si k = Nz−1)
        sur-diag  = −rₖ · Uₖ₊₁              (0 si k = Nz−1)

 7. RÉSOUDRE  A · cⁿ⁺¹ = cⁿ           (solve_banded, direct)      [l. 80]
 8. ÉCRIRE    pcol[pidx] = cⁿ⁺¹                                   [l. 80]
```

### Deux points d'ordre qui comptent
1. **L'injection (4) est AVANT le solve (7)** → la dose du jour est **transportée dès le jour
   même**. Le 2ᵉ membre du système est bien `cⁿ + dose`.
2. **w est lu à la maille COURANTE (2)** → chaque colonne subit la vitesse verticale de l'endroit
   où elle se trouve, et cet endroit change au fil de la dérive.

---

## 8. Les 3 propriétés obtenues

### 8.a Conservation exacte
La **forme flux** garantit que ce qui sort de la maille k par l'interface k+1 **entre exactement**
dans la maille k+1 (même nombre, signe opposé). Avec des flux nuls aux bords, tout se télescope :

```
   Σₖ hₖ·cₖⁿ⁺¹  =  Σₖ hₖ·cₖⁿ
```

> ⚠️ L'invariant est **Σₖ hₖ·cₖ** (pondéré par l'épaisseur), **pas** Σₖ cₖ — les mailles n'ont pas
> la même taille. *Vérifié numériquement : erreur relative ~1e−15.*

### 8.b Positivité
Grâce à l'upwind, Lᵢ ≥ 0 et Uᵢ ≥ 0 **toujours**. Donc A a une diagonale > 0, des extra-diagonaux
≤ 0, et est à diagonale dominante : c'est une **M-matrice** → `cⁿ ≥ 0 ⟹ cⁿ⁺¹ ≥ 0`.
**Jamais de concentration négative.**

> Un schéma **centré** (cᵢ ≈ (cᵢ₋₁+cᵢ)/2 à l'interface) donnerait Lᵢ = −wᵢ/2 + Kᵢ/δᵢ, qui peut
> devenir **négatif** → perte de la M-matrice → oscillations. Il n'est sûr que si le Péclet de
> maille `Pe = |w|·δ/K ≤ 2`, ce qui casse en profondeur où K tombe à 1e−5.

### 8.c Stabilité inconditionnelle
Aucune contrainte sur Δt — c'est tout l'intérêt face à D ≈ 7,7.

---

## 9. Le format `solve_banded` (le piège d'implémentation)

`scipy.linalg.solve_banded((1,1), ab, b)` veut la matrice **empilée en bandes**, `ab` de forme
(3, Nz), avec la convention :

```
   ab[0, j] = A[j−1, j]    (sur-diagonale)
   ab[1, j] = A[j,   j]    (diagonale)
   ab[2, j] = A[j+1, j]    (sous-diagonale)
```

D'où le décalage d'indices, qui surprend à la lecture :

```python
ab[1, k]     = 1.0 + r * (Uk + Lk1)      # A[k,k]    -> colonne k
ab[2, k - 1] = -(dt / h[k]) * L[k]       # A[k,k-1]  -> rangé en colonne k-1 !
ab[0, k + 1] = -(dt / h[k]) * U[k + 1]   # A[k,k+1]  -> rangé en colonne k+1 !
```

---

## 10. Le couplage 2D + 1D

Le kernel ci-dessus, c'est le **1D vertical**. Le **2D horizontal** est ailleurs :

```python
kernels = pset.Kernel(VerticalAdvDiffImplicit) + pset.Kernel(parcels.AdvectionRK4)
```

- `VerticalAdvDiffImplicit` → fait évoluer **pcol[pidx]** (la colonne)
- `parcels.AdvectionRK4` → déplace **la particule** (RK4 sur u, v)

**La colonne suit la particule automatiquement**, puisqu'elle lui est attachée par `pidx`.
Il n'y a **aucun échange entre mailles** — c'est ce qui supprime la *cascade* (avec un champ
**partagé**, les échanges séquentiels entre particules relayaient une colonne vers l'avant et
empilaient l'alcalinité sur la maille de tête).

Le champ gridé n'est **reconstruit** que pour les sorties :

```
   ALK0[:, y, x]  =  Σ  pcol[p]        (somme sur les particules p présentes dans la maille (x,y))
```
```python
fieldset.ALK0.data[:] = 0.0
for p in pset:
    fieldset.ALK0.data[0, :, p.cell_y, p.cell_x] += fieldset.pcol[p.pidx]
```

> 🎯 **À retenir** : l'état réel du modèle, c'est **pcol (31 × 50 = 1 550 nombres)**.
> ALK0 (408 Mo) n'est qu'un **diagnostic** reconstruit à la demande.

---

## 11. Ce que ça donne (31 jours, ACC)

| Diagnostic | Résultat |
|---|---|
| Conservation Σₖ hₖ·cₖ | = 31 × la dose, erreur relative ~1e−15 |
| Concentrations négatives | aucune |
| Pénétration **advection seule** (K = 0) | ~5 m (tout reste piégé en surface : w est négligeable) |
| Pénétration **advection + diffusion** | **~75–85 m** |
| Dérive horizontale | ~0,2 °/jour vers l'est le long de l'ACC |

👉 C'est **la diffusion** qui fait descendre l'alcalinité, pas l'advection verticale.

---

## 12. Les limites actuelles (honnêtes)

1. **Pas de chimie du carbonate** : ni DIC, ni pCO₂, ni flux air–mer. Le modèle dit *où va
   l'alcalinité*, pas *combien de CO₂ est absorbé*. (= perspective n°3 de la présentation)
2. **Courants figés** : un seul pas de temps (moyenne mensuelle) réutilisé pour les 31 jours
   (`allow_time_extrapolation=True`).
3. **Hypothèse forte** : mêmes courants horizontaux à toutes les profondeurs (toute la colonne
   suit la trajectoire de surface).
4. **Pas de diffusion horizontale** (supposée négligeable devant l'advection).
5. L'alcalinité **de fond** (`fieldset.ALK`) est chargée mais **jamais lue** — le transport de la
   perturbation est linéaire. Elle ne servira qu'avec la chimie.
6. **Moyenne arithmétique** pour K aux interfaces. La moyenne **harmonique** serait plus rigoureuse
   vu le saut 9e−3 → 1e−5 entre les niveaux 1 et 7. À tester.
