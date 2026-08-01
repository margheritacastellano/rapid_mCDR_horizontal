# Extraire soi-même les données pour la chimie CO₂

Comment produire `sensitivity_pco2_1995-01.nc` (le fichier qui active la chimie) **sans passer
par Océane**. Trois étapes : installer 2 paquets, télécharger 3 fichiers, lancer un script.

> Résultat : `data/sensitivity_pco2_1995-01.nc` → `test.ipynb` passe automatiquement en
> `CHIMIE = True`. Voir `CHIMIE_CO2.md` pour ce que fait la chimie une fois activée.

---

## Étape 1 — installer les paquets (2 min)

```bash
/opt/anaconda3/bin/pip install PyCO2SYS tqdm
```
*(déjà fait sur cette machine ; à refaire si tu changes d'environnement.)*

## Étape 2 — télécharger 3 binaires ECCO-Darwin

Portail **public** (pas de compte nécessaire), organisé par variable. Le pas de temps
**`0000081144`** correspond à **janvier 1995** (le même que pour ALK/U/V/W).

Les mettre dans **`SEMES/data/ECCO-Darwin_data/`** :

```bash
mkdir -p /Users/margheritacastellano/Documents/These-Recherche/SEMES/data/ECCO-Darwin_data
cd       /Users/margheritacastellano/Documents/These-Recherche/SEMES/data/ECCO-Darwin_data

BASE=https://data.nas.nasa.gov/ecco/llc_270/ecco_darwin_v5/output/monthly
wget $BASE/DIC/DIC.0000081144.data            # ~181 Mo
wget $BASE/SALTanom/SALTanom.0000081144.data  # ~181 Mo
wget $BASE/SST/SST.0000081144.data            # ~4 Mo
wget $BASE/SIarea/SIarea.0000081144.data      # ~4 Mo  (optionnel : fraction de glace réelle)
```
*(si `wget` n'est pas installé : `curl -O <url>` marche aussi, ou clic droit → enregistrer
depuis le navigateur.)*

| fichier | ce que c'est | taille |
|---|---|---|
| `DIC.0000081144.data` | carbone inorganique dissous (fond) | ~181 Mo |
| `SALTanom.0000081144.data` | anomalie de salinité | ~181 Mo |
| `SST.0000081144.data` | température de surface | ~4 Mo |
| `SIarea.0000081144.data` | fraction de glace *(optionnel)* | ~4 Mo |

> **L'ALK de fond n'est pas à télécharger** : le script réutilise `data/ALK_1995-01.nc`.
> **La grille non plus** : elle est déjà dans le dépôt (`flat_XC_inter.npy`, `flat_YC_inter.npy`).
> On ne lit que la **surface** de chaque binaire, donc pas besoin de charger les 181 Mo en mémoire.

## Étape 3 — lancer le script

```bash
/opt/anaconda3/bin/python \
  /Users/margheritacastellano/Documents/These-Recherche/SEMES/rapid_mCDR_horizontal/DatasLLC270/extraire_sensibilite_pco2.py
```

Il lit DIC/SALT/SST (+ ALK existant), calcule les sensibilités avec **PyCO2SYS**, et écrit
`data/sensitivity_pco2_1995-01.nc`. Contrôles affichés (doivent être vrais) :
```
dpco2/dDIC : median > 0
dpco2/dALK : median < 0
```

---

## Ce qu'il reste après (petits réglages)

Une fois le fichier en place, `test.ipynb` détecte `CHIMIE = True`. Deux constantes restent
à ajuster dans la cellule de setup (aujourd'hui des placeholders) :

- **`siarea`** : le script affiche la valeur réelle à la source (≈ 0 à l'ACC, donc le placeholder
  `0` est déjà bon).
- **`k_surf`** : coefficient d'échange air-mer, placeholder `3e-5`. Valeur exacte : soit celle
  d'Océane (données Suselj), soit calculée depuis `wspeed` + SST/SALT (formule de Wanninkhof).
  Un ordre de grandeur suffit pour des résultats **relatifs**.

⚠️ **Caveat d'unités** (voir `CHIMIE_CO2.md`, §4) : notre `ALK0` n'est pas divisé par l'épaisseur
de maille à l'injection, il est donc ~10× l'échelle µM attendue par ces sensibilités. Le **signe
et la mécanique** du flux sont corrects ; la **valeur absolue** devra être recalée (le plus propre :
aligner l'injection sur `c[0] += forcing/h[0]*dt`).

---

## Détails techniques (pour référence)

- Format des binaires : MDS MITgcm, **big-endian float32**, dimensions brutes
  `(nz=50, ny=13×270, nx=270)`. Le script lit le niveau 0 (`np.fromfile(dtype=">f4", count=ny*nx)`)
  puis remet en grille globale `(945, 1080)` via `flat()` de `data_functions.py`.
- La grille de sortie est **identique** à celle de nos champs (vérifié : écart 0.0 avec
  `ALK_1995-01.nc`), donc `dpco2[yi,xi]` s'aligne exactement sur `ALK0[yi,xi]`.
- Le calcul de sensibilité : `compute_dpco2_sensitivity` (`data_functions.py`) appelle
  `PyCO2SYS.sys` avec une perturbation de ±10 µmol/kg sur DIC puis ALK.
- Autres mois : changer `FILE` (le pas de temps) et `DTIME` dans le script. Le mapping
  date → pas de temps est dans `read-create_LLC270data.ipynb` (cellule des réglages).

---

*Voir aussi : `CHIMIE_CO2.md` (ce que fait la chimie), `SCHEMA_NUMERIQUE.md` (le schéma vertical).*
