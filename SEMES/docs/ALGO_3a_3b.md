# Les cellules 3a et 3b : la simulation

Explication des deux sections qui **lancent** le modèle :
- **3a** — préparer les particules (*qui* va simuler, et *avec quoi*)
- **3b** — la boucle de simulation (*comment* on avance dans le temps)

---

# PARTIE 1 — L'IDÉE

## L'idée en une image

On veut simuler une injection d'alcalinité **tous les jours pendant 31 jours**, au même endroit.

L'astuce : **1 jour = 1 particule**. Chaque particule représente **le devenir de la dose de son
jour**. Elle porte sa colonne d'alcalinité dans son « sac à dos » (`pcol[pidx]`), et dérive avec
le courant en la diffusant vers le bas.

```
   jour 0 :   ●P0                                 1 particule
   jour 1 :   ●P1 ●P0                             2 particules   (P0 a dérivé)
   jour 2 :   ●P2 ●P1 ●P0                         3 particules
    ...
   jour 30:   ●P30 ●P29 ... ●P1 ●P0              31 particules
              └─source──┘        └─la plus vieille, la plus loin
```

À la fin, la **traînée** qu'on observe, ce sont ces 31 doses à 31 endroits différents : chacune
est partie un jour différent, donc a dérivé plus ou moins loin.

## Le déroulé d'une journée

Chaque jour, on fait **exactement 3 choses** :

```
   1. NAÎTRE     →  on crée 1 nouvelle particule à la source
   2. AVANCER    →  on fait avancer TOUT LE MONDE d'exactement 1 jour
                      ├─ VerticalAdvDiffImplicit : chaque colonne diffuse/advecte verticalement
                      └─ parcels.AdvectionRK4    : chaque particule dérive horizontalement
   3. PHOTOGRAPHIER → (certains jours) on reconstruit le champ ALK0 et on l'écrit sur disque
```

## Les 3 pièges que ces cellules doivent éviter

Ces cellules ont l'air anodines, mais **trois détails** y sont vitaux — chacun a été un bug réel
qui figeait totalement le modèle :

| détail | si on l'oublie |
|---|---|
| lâcher à **5 m** (pas 0 m) | U/V sont définis à partir de 5 m → à 0 m on échantillonne **hors-bornes** → la particule **ne bouge jamais** |
| **horodater** chaque nouvelle particule (`time=day*86400`) | elle naît à t=0 alors que le set est au jour N → parcels la considère « dans le passé » → **elle n'avance jamais** |
| donner un **`pidx` unique** à chaque particule | deux particules partageraient la même ligne de `pcol` → leurs doses fusionneraient |

---

# PARTIE 2 — LIGNE PAR LIGNE

## 3a — Préparer les particules

### Cellule 14 : les réglages

```python
position = (-2, -49)   # Initial position of the particle (longitude, latitude)
last_time = 31         # Number of days to run the simulation
```

| ligne | explication |
|---|---|
| `position` | le point de lâcher : **2°W, 49°S**, en plein courant circumpolaire antarctique (ACC) |
| `last_time` | **31 jours** de simulation ⟹ **31 particules** (1/jour) ⟹ `pcol` aura **31 lignes** |

> 💡 `last_time` sert **deux** choses : la durée **et** le nombre de particules. C'est la même
> valeur parce que la règle est « 1 particule par jour ».

### Cellule 15 : la classe Particle + les colonnes

```python
class Particle(parcels.ScipyParticle):
    forcing_alk = parcels.Variable('f_alk', dtype=np.float64, initial=fieldset.alkalinity_forcing)
    cell_x = parcels.Variable('cell_x', dtype=np.int32, initial=-1)
    cell_y = parcels.Variable('cell_y', dtype=np.int32, initial=-1)
    released = parcels.Variable('released', dtype=np.int32, initial=0)
    pidx = parcels.Variable('pidx', dtype=np.int32, initial=0)
```

On **hérite** de `parcels.ScipyParticle` : on récupère gratuitement `lon`, `lat`, `depth`, `time`,
`dt`… et on **ajoute nos propres variables** :

| variable | type | rôle | qui l'écrit ? |
|---|---|---|---|
| `f_alk` | float | la dose portée — **purement diagnostique**, écrite dans le `.zarr` et tracée plus loin | personne (constante) |
| `cell_x`, `cell_y` | int | **la maille courante** de la particule. Initialisée à `-1` = « pas encore localisée » | le kernel vertical, à chaque pas |
| `released` | int | `0` tant que la dose n'a pas été déposée, `1` après | le kernel vertical, **une seule fois** |
| **`pidx`** | int | ⭐ **le numéro de ligne de `pcol` qui appartient à cette particule** | la boucle 3b, à la naissance |

> ⚠️ **Contrainte parcels** : une variable de particule doit être un **scalaire**. Impossible d'y
> accrocher un tableau de 50 valeurs. C'est **toute la raison d'être de `pidx`** : la colonne vit
> ailleurs (dans `pcol`), et la particule ne porte que le **numéro** qui permet de la retrouver.

```python
fieldset.pcol = np.zeros((last_time, fieldset.ALK0.data.shape[1]))
```

| | |
|---|---|
| **forme** | `(31, 50)` = (nb de particules, nb de niveaux verticaux) |
| **`pcol[p, j]`** | concentration d'ALK₀ (µM) au niveau `j` de la colonne portée par la particule `p` |
| **initialisation** | des **zéros** : on transporte une **perturbation**, qui part de rien |
| **où ?** | attaché au `fieldset` — parce que les kernels reçoivent `fieldset` en argument, c'est le seul moyen d'y accéder depuis un kernel |
| **`.data.shape[1]`** | on lit `50` depuis le champ plutôt que de l'écrire en dur |

> 🔑 **C'est ici que vit l'état réel du modèle.** 31 × 50 = **1 550 nombres**. Le champ `ALK0`
> (408 Mo) n'est qu'un diagnostic reconstruit. Voir `POURQUOI_LA_CASCADE.md`.

```python
release_depth = 5.0
```

**Le piège n°1.** Les vitesses U et V sont définies aux **centres** de maille, dont le plus haut
est à **5 m**. Une particule lâchée à **0 m** échantillonnerait les vitesses **au-dessus** de leur
premier niveau → erreur hors-bornes → l'advection échoue silencieusement → **la particule reste
figée à la source pour toujours**.

```python
pset = parcels.ParticleSet(fieldset=fieldset, pclass=Particle,
                           lon=[position[0]], lat=[position[1]],
                           depth=[release_depth], time=[0.0])
```

La **première** particule (celle du jour 0) :

| argument | valeur | pourquoi |
|---|---|---|
| `pclass=Particle` | notre classe | pour avoir `pidx`, `released`, `cell_x/y` |
| `lon`, `lat` | `(-2, -49)` | la source |
| `depth=[5.0]` | 5 m | piège n°1 |
| `time=[0.0]` | t = 0 | elle naît au début |

Son `pidx` vaut **0** — c'est la valeur `initial=0` de la classe. Elle prend donc la **ligne 0**
de `pcol`. (Les suivantes recevront `pidx = day` explicitement en 3b.)

---

## 3b — La boucle de simulation

### Le décor

```python
folder = "Test"
filename_save = f"{folder}/DailyRapidmCDR_LLC270_ACC_{date}"
os.makedirs(folder, exist_ok=True)
```
Le dossier de sortie (chemin **relatif** → créé là où tu lances le notebook).

```python
output_file = pset.ParticleFile(name=f"{filename_save}.zarr", outputdt=timedelta(days=1))
```
Le fichier des **trajectoires** (positions des particules), écrit **tous les jours** par parcels.
⚠️ Il ne contient **pas** le champ ALK₀ — celui-là, on l'écrit nous-mêmes (voir plus bas).

```python
kernels = pset.Kernel(VerticalAdvDiffImplicit) + pset.Kernel(parcels.AdvectionRK4)
```

⭐ **Le couplage 2D+1D tient dans cette ligne** :

```
   pset.Kernel(VerticalAdvDiffImplicit)  →  le 1D vertical (advection + diffusion sur pcol[pidx])
                    +
   pset.Kernel(parcels.AdvectionRK4)     →  le 2D horizontal (RK4 sur u, v)
```

Le `+` **enchaîne** les kernels : pour chaque particule, à chaque pas, parcels exécute le premier
puis le second. **La colonne suit la particule automatiquement** (elle lui est attachée par
`pidx`) — d'où l'absence d'un 3ᵉ kernel de transport horizontal.

### La sauvegarde

```python
snapshot_days = [1, 5, 13, last_time]

def save_ALK0_snapshot(day_index):
    if day_index not in snapshot_days:
        return
    fieldset.ALK0.data[:] = 0.0
    for p in pset:
        if p.cell_x >= 0:
            fieldset.ALK0.data[0, :, p.cell_y, p.cell_x] += fieldset.pcol[p.pidx]
    fieldset.ALK0.write(f"{filename_save}_{day_index:04d}")
```

C'est **la reconstruction** : on passe du lagrangien (`pcol`) à l'eulérien (`ALK0`), uniquement
pour les sorties.

| ligne | explication |
|---|---|
| `snapshot_days` | chaque fichier pèse ~0,4 Go → on n'en écrit que 4 (au lieu de 31 = 13 Go) |
| `if day_index not in snapshot_days: return` | on sort tout de suite les autres jours (aucun calcul) |
| `fieldset.ALK0.data[:] = 0.0` | ⚠️ **on repart de zéro** : le champ est un *cache d'affichage*, pas un état à conserver |
| `for p in pset` | on parcourt **toutes** les particules vivantes |
| `if p.cell_x >= 0` | garde-fou : une particule jamais localisée aurait encore `-1` |
| `ALK0[0, :, p.cell_y, p.cell_x] += pcol[p.pidx]` | ⭐ on **dépose** la colonne de `p` à **sa** maille. Le `+=` fait que deux particules dans la même maille **s'additionnent à l'affichage** — correct — **sans fusionner** dans `pcol` |
| `.write(f"..._{day_index:04d}")` | produit `..._0005ALK0.nc` (parcels ajoute `ALK0.nc`) — le nom que lisent les cellules de tracé |

### Le jour 0

```python
pset.execute(kernels, dt=timedelta(days=1), runtime=timedelta(days=1), output_file=output_file)
save_ALK0_snapshot(1)
```

| argument | rôle |
|---|---|
| `dt=1 jour` | le **pas de temps** du schéma |
| `runtime=1 jour` | **la durée de cet appel** → exactement **1 pas** |

> 🔑 `dt == runtime` ⟹ **1 seul pas de temps par `execute`**. C'est ce qui rend la boucle
> lisible : *un tour de boucle = un jour*. (La version d'origine faisait
> `runtime=last_time` puis re-`execute` avec des durées décroissantes — c'était la source
> d'une belle confusion temporelle.)

### La boucle quotidienne

```python
for day in range(1, last_time):
    new_particle = parcels.ParticleSet(
        fieldset=fieldset, pclass=Particle,
        lon=[position[0]], lat=[position[1]],
        depth=[release_depth], time=[day * 86400.0],
    )
    new_particle.pidx[:] = day
    pset.add(new_particle)

    pset.execute(kernels, dt=timedelta(days=1), runtime=timedelta(days=1), output_file=output_file)
    save_ALK0_snapshot(day + 1)
```

Déroulé, ligne par ligne :

| ligne | explication |
|---|---|
| `for day in range(1, last_time)` | jours **1 → 30** (le jour 0 est déjà fait) ⟹ 30 + 1 = **31 particules** |
| `parcels.ParticleSet(...)` | on crée la particule du jour, **à la source** |
| `depth=[release_depth]` | **piège n°1** : 5 m, sinon elle ne bougera jamais |
| `time=[day * 86400.0]` | ⭐ **piège n°2** : on l'**horodate** à l'instant courant. Sans ça elle naît à `t=0` alors que le set est au jour `N` → parcels la voit « dans le passé » et **ne l'avance jamais** (bug réel : toutes les particules restaient figées) |
| `new_particle.pidx[:] = day` | ⭐ **piège n°3** : on lui donne **sa** ligne de `pcol`. La particule du jour 5 → ligne 5. Sans ça, toutes utiliseraient la ligne 0 et **leurs doses fusionneraient** |
| `pset.add(new_particle)` | elle rejoint le troupeau |
| `pset.execute(..., runtime=1 jour)` | **tout le monde** avance d'un jour (l'ancienne **et** la nouvelle) |
| `save_ALK0_snapshot(day + 1)` | photo (si le jour est dans `snapshot_days`) |

> 💡 **Où est l'injection de la dose ?** Nulle part ici ! Elle est **dans le kernel vertical** :
> ```python
> if particle.released == 0:
>     c[0] += fieldset.alkalinity_forcing * dt
>     particle.released = 1
> ```
> Chaque particule dépose **sa** dose toute seule, **le jour de sa naissance** (`released` passe
> à 1 et ne redescend jamais). C'est ce qui garantit **exactement 1 dose par particule**, donc
> **1 dose par jour**.

---

## Le bilan chiffré

Après 31 jours :

| | |
|---|---|
| particules créées | **31** (1/jour) |
| doses injectées | **31** (1/particule, via `released`) |
| appels du kernel | 1 + 2 + 3 + … + 31 = **496** |
| masse finale `Σ hₖ·cₖ` | **31 × la dose**, à 10⁻¹⁵ près ✅ |
| concentrations négatives | **aucune** ✅ |
| mailles touchées | **18** |
| dérive | ~0,2 °/jour vers l'est |
| pénétration verticale | ~**75–85 m** (grâce à la diffusion ; ~5 m sans elle) |

---

## Le récapitulatif visuel

```
   3a  ─────────────────────────────────────────
        position, last_time          les réglages
        class Particle               pidx, released, cell_x/y
        fieldset.pcol (31×50)        ⭐ l'état du modèle
        release_depth = 5.0          ⚠️ piège n°1
        pset = ParticleSet(...)      la particule du jour 0

   3b  ─────────────────────────────────────────
        kernels = Vertical + RK4     ⭐ le couplage 2D+1D
        save_ALK0_snapshot()         lagrangien → eulérien (sorties)

        execute(1 jour) ; snapshot   ← jour 0

        POUR day = 1 … 30 :
            créer la particule       depth=5 ⚠️, time=day*86400 ⚠️, pidx=day ⚠️
            pset.add()
            execute(1 jour)          ← tout le monde avance
            snapshot                 ← si day+1 ∈ snapshot_days
```

---

*Voir `SCHEMA_NUMERIQUE.md` (les maths du kernel vertical) et `POURQUOI_LA_CASCADE.md`
(pourquoi les colonnes sont par particule).*
