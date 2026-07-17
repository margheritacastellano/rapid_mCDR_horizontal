# Partie 4 : charger et visualiser les résultats

La simulation (partie 3) a produit des fichiers sur le disque. La partie 4 les **relit** et
en fait **trois** choses : les trajectoires, la carte de concentration, et la vérification
de conservation.

> ⚠️ La partie 4 **ne calcule rien**. Elle ne fait que relire `Test/`. Tu peux donc la
> relancer autant de fois que tu veux (changer de jour, de profondeur…) **sans refaire la
> simulation** — qui prend ~1 min.

---

## 1. Ce que la simulation a laissé sur le disque

```
Test/
├── DailyRapidmCDR_LLC270_ACC_1995-01.zarr          ← les TRAJECTOIRES (positions des particules)
├── DailyRapidmCDR_LLC270_ACC_1995-01_0001ALK0.nc   ┐
├── ...                                _0005ALK0.nc  │ les CHAMPS ALK0 complets
├── ...                                _0013ALK0.nc  │ (~0,4 Go chacun, 50 niveaux)
└── ...                                _0031ALK0.nc  ┘
```

Les 4 jours sauvegardés viennent de `snapshot_days = [1, 5, 13, last_time]` (cellule de la
partie 3b). **Ce sont les seuls jours traçables.**

---

## 2. Les trois cellules

### 📈 Les trajectoires

```python
filename  = f"Test/DailyRapidmCDR_LLC270_ACC_{date}"
pset_traj = xr.open_zarr(f"{filename}.zarr")
...
plt.plot(pset_traj["lon"].T, pset_traj["lat"].T, ".-")
plt.savefig("figures/trajectoires.png", ...)
```

Chaque courbe = **le chemin d'une particule** = le trajet d'une dose quotidienne.
Elles partent toutes du même point (2°W, 49°S) et suivent le courant circumpolaire vers l'est
(~0,2 °/jour). Les plus longues sont les plus vieilles (lâchées le jour 0).

### 🗺️ La carte de concentration (2 panneaux)

C'est **la** figure du modèle. Les réglages sont **en haut de la cellule** :

```python
var_name = "ALK0"
dayFirst = "0001"     # panneau de GAUCHE
dayLast  = "0013"     # panneau de DROITE
depth    = 1          # le niveau vertical
```

Elle affiche **deux instants côte à côte**, à une profondeur donnée :

| panneau | contenu |
|---|---|
| **gauche** (`dayFirst`) | le champ au jour 1 — la dose initiale à la source |
| **droite** (`dayLast`) | le champ au jour N — **masqué** là où rien n'a changé depuis le jour 1 (fond blanc) |

### ✅ La vérification de conservation

```python
mass_model = np.nansum(h_fv[:, None, None] * field_last)     # Σ hₖ·cₖ
```

Doit afficher :
```
=> effective number of surface injections = 31.000  (expected ~ last_time = 31)
```
Si ce n'est pas ≈ `last_time`, quelque chose cloche dans la simulation.

> 🔑 On somme **`Σ hₖ·cₖ`** (pondéré par l'épaisseur de maille), **pas** `Σ cₖ`. Les mailles
> font de 10 m (surface) à 450 m (fond) : une somme brute n'aurait aucun sens physique.

---

## 3. Comment choisir le **jour**

```python
dayLast = "0013"
```

⚠️ **Deux contraintes** :

1. **4 chiffres, entre guillemets** : `"0005"` et non `5`
2. **le jour doit être dans `snapshot_days`**, sinon → `FileNotFoundError`

Avec `snapshot_days = [1, 5, 13, 31]`, les seules valeurs valides sont :

```
   "0001"     "0005"     "0013"     "0031"
```

**Pour tracer un autre jour** (ex. le 20), il faut **relancer la simulation** après avoir ajouté
le jour à la liste :
```python
snapshot_days = [1, 5, 13, 20, last_time]    # partie 3b
```
*(chaque jour ajouté = +0,4 Go sur le disque)*

---

## 4. Comment choisir la **profondeur** ⚠️

```python
depth = 1
```

> ### 🚨 `depth` est un **INDICE DE NIVEAU**, pas une profondeur en mètres !
> `depth = 1` ne veut **pas** dire « 1 mètre » — c'est le **niveau n°1**, soit **15 m**.

La correspondance (les centres de maille de la grille LLC270) :

| `depth` | profondeur réelle | |
|---|---|---|
| `0` | **5 m** | la surface — c'est là qu'on injecte |
| `1` | 15 m | |
| `2` | 25 m | |
| `3` | 35 m | |
| `4` | 45 m | |
| `5` | 55 m | |
| `6` | 65 m | |
| `7` | 75 m | |
| `8` | 85 m | ← ~la limite du signal |
| `9` et + | 95 m et + | **quasiment vide** |

💡 L'alcalinité descend jusqu'à ~**85 m** grâce à la diffusion verticale (contre ~5 m avec la
seule advection). Au-delà de `depth = 8`, la carte sera blanche.

---

## 5. Les combinaisons intéressantes

```python
# carte "deux jours, une profondeur"
depth = 0 ;  dayLast = "0031"     # ⭐ LA figure : la traînée en surface, à la fin
depth = 0 ;  dayLast = "0005"     # la traînée à ses débuts

# carte "un jour, deux profondeurs"  (echelle_commune = True)
jour = "0031" ; niveau_A = 0 ; niveau_B = 4    # ⭐ montre le gradient vertical à la source
jour = "0005" ; niveau_A = 0 ; niveau_B = 4    # le même, quand la traînée est jeune
```

### 🔬 Le résultat physique : la colonne s'homogénéise en vieillissant

Compare `niveau_A = 0` et `niveau_B = 4` au jour 31 (avec la cellule « deux profondeurs ») :

| longitude | 5 m | 45 m | rapport |
|---|---|---|---|
| **−1,83 (la source)** | **0,750** | **0,093** | **×8** ← énorme gradient |
| −0,83 | 0,432 | 0,307 | ×1,4 |
| 1,50 | 0,373 | 0,333 | ×1,1 |
| **3,50 (le plus loin)** | **0,352** | **0,328** | **×1,07** ← quasi identique |

> **Le seul endroit avec un vrai gradient vertical, c'est la SOURCE.**
> La dose fraîche du jour est encore collée en surface (0,75 en haut contre 0,09 à 45 m).
> En aval, les doses plus vieilles ont eu le temps de se mélanger : la colonne est
> **quasi uniforme** sur les 50 premiers mètres (0,35 contre 0,33).

C'est **ça**, la signature de la diffusion verticale : elle homogénéise les ~50 premiers mètres
en quelques jours. Sur la traînée, on voit donc le **vieillissement** de chaque dose : plus elle
est loin, plus sa colonne est mélangée.

### 🚨 Le piège des colorbars

Sur la carte, **chaque panneau auto-ajuste sa propre échelle**. Au jour 31 :
```
   niveau 0 :  colorbar 0,20 → 0,75
   niveau 4 :  colorbar 0,09 → 0,33
```
Donc **le « jaune » à 45 m (0,33) est plus FONCÉ que le « bleu » à 5 m (0,35)**.

👉 Comparer les **couleurs** entre deux panneaux auto-échelonnés **n'a aucun sens** : les deux
cartes *paraissent* se ressembler alors que les valeurs diffèrent d'un facteur 8 à la source.

C'est pour ça que la cellule « deux profondeurs » utilise **`echelle_commune = True`** par
défaut. Avec une échelle commune, l'inversion à la source saute aux yeux : ⭐ jaune vif à 5 m
vs ⭐ bleu foncé à 45 m, au même point.

---

## 6. Où sont sauvegardées les figures

```
figures/
├── trajectoires.png
├── ALK0_jour0001_vs_jour0013_niveau1.png    ← deux JOURS, une profondeur
└── ALK0_jour0031_niveaux0-4.png             ← un jour, deux PROFONDEURS
```

Le nom de la carte **encode les paramètres** : changer `depth` ou `dayLast` produit un
**nouveau** fichier au lieu d'écraser le précédent. Tu peux donc accumuler tes comparaisons.

*(`figures/` est gitignoré : ce sont des sorties régénérables.)*

---

## 7. Ce qui diffère du notebook d'Océane

| | Océane | Nous |
|---|---|---|
| **snapshots** | `_initial` + `_final` → **2** (début/fin seulement) | `snapshot_days` → **n'importe quel jour**, configurable |
| **la carte** | initial **vs** final | **jour A vs jour B** → on voit l'évolution |
| **`depth_levels`** | recopiés **en dur** dans le notebook | lus depuis la **vraie grille** |
| **vérif de masse** | `Σ cₖ` brut ❌ | **`Σ hₖ·cₖ`** ✅ |
| **sauvegarde des figures** | ❌ (à la main) | ✅ `savefig` automatique |
| **plot 3D voxel** | ✅ présent | ❌ **absent** |

⚠️ **La dernière ligne est une vraie régression** : Océane a un plot 3D en voxels
(`%matplotlib tk`), pas nous. Il en existe une version générée hors notebook
(`depth_plots/plume_3d.png`), à réintégrer si besoin.

---

## 8. Les pièges

1. **`dayLast` doit être dans `snapshot_days`** — sinon `FileNotFoundError`.
2. **`depth` est un indice, pas des mètres** — `depth=1` → 15 m.
3. **Relancer la partie 3b efface `Test/`** (le nettoyage est au début de la cellule de
   simulation). Les figures déjà dans `figures/`, elles, ne sont pas touchées.
4. **Le panneau de droite est masqué** là où rien n'a changé depuis `dayFirst` → un fond blanc
   est **normal**, ce n'est pas un bug.

---

*Voir aussi : `SCHEMA_NUMERIQUE.md` (les maths), `ALGO_3a_3b.md` (la simulation),
`POURQUOI_LA_CASCADE.md` (la conception par particule).*
