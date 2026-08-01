# La chimie du carbonate : DIC + flux de CO₂

Comment le modèle passe de *« où va l'alcalinité ? »* à *« combien de CO₂ est absorbé ? »*.

> ⚠️ **État actuel : le code est prêt, mais les DONNÉES manquent.** Sans le fichier de
> sensibilité pCO₂, la chimie est **désactivée** (`CHIMIE = False`) et le notebook ne simule que
> le transport de l'alcalinité, exactement comme avant. Ce document explique la physique et
> **comment obtenir les données** pour l'activer.

---

## 1. La chaîne physique de l'OAE

```
   + Alcalinité  →  ALK0 ↑  →  pCO₂ (océan) ↓  →  flux CO₂ atmosphère→océan  →  DIC0 ↑
   (ce qu'on injecte)          (déséquilibre chimique)   (le CO₂ ENTRE)        (carbone stocké)
```

Le modèle simulait déjà les 2 premières flèches (transport de l'ALK). La chimie ajoute les 2
dernières : le **flux de CO₂** (`f_co2`) et le **carbone dissous** (`DIC0`).

## 2. La formule (linéarisée, méthode rapid-mCDR de Suselj et al.)

Pas de chimie du carbonate résolue au run-time. On **linéarise** la réponse de pCO₂ autour de
l'état de fond, avec deux **sensibilités** pré-calculées :

```
   f_co2 = -k_surf · (1 - siarea) · ( dpco2_ddic · DIC0  +  dpco2_dalk · ALK0 )

   DIC0[surface] += (f_co2 / h₀) · Δt          # le CO₂ absorbé entre dans la maille de surface
```

| symbole | sens | signe |
|---|---|---|
| `dpco2_dalk` | ∂pCO₂/∂ALK : combien pCO₂ baisse si on ajoute de l'ALK | **< 0** |
| `dpco2_ddic` | ∂pCO₂/∂DIC : combien pCO₂ monte si on ajoute du DIC | **> 0** |
| `k_surf` | coefficient d'échange air-mer (vitesse de piston × solubilité) | > 0 |
| `siarea` | fraction de glace de mer (0 = pas de glace) | 0–1 |

**Le signe** : `ALK0 > 0` et `dpco2_dalk < 0` → `f_co2 > 0` → `DIC0` augmente. Le CO₂ **entre**.
Puis, à mesure que `DIC0` s'accumule, le terme `dpco2_ddic · DIC0 > 0` **réduit** le flux :
c'est la **rétroaction** (saturation) qui ralentit l'absorption. ✓ *(vérifié en test synthétique.)*

## 3. Comment c'est intégré chez nous

- **Un 2ᵉ traceur par particule** : `fieldset.pcol_dic` (à côté de `fieldset.pcol` pour l'ALK).
- **Dans le kernel** : après le transport de l'ALK, si `chimie == 1`, on calcule `f_co2` en surface
  puis on transporte `DIC0` par **le même schéma vertical implicite** (advection + diffusion).
- **Grâce à la linéarité**, chaque particule calcule sa propre contribution → pas de cascade.
- **Bonus vs Océane** : son modèle est 2D (DIC en surface seulement). Le nôtre est **3D** : le CO₂
  absorbé **diffuse en profondeur**, comme l'alcalinité.
- **Sorties** : champ `DIC0` reconstruit (`..._{NNNN}DIC0.nc`) + flux `f_co2` par particule (dans
  le `.zarr`) → carte du flux (`figures/flux_co2.png`).

## 4. ✅ Les unités de l'injection (vérifié contre le code source d'Océane et le papier Suselj et al. 2025)

Aucun recalage nécessaire : `alkalinity_forcing = 0.000012409` est **déjà** un taux volumique
(meq/(m³·s), donc en µM ALK par seconde), pas un flux de surface. Preuve, remontée jusqu'à la
source :

- Dans le **modèle 2D d'origine** d'Océane (`RapidmCDR/LLC270_OceanParcels_rapidmCDR.ipynb`), le
  flux **brut** chargé depuis les données de Suselj vaut `f_alk = alk_forcing/area = **1,2409e-4**
  meq/(m²·s)` (vérifié par son propre print de contrôle). Son kernel fait alors
  `ALK0[0] += (f_alk / dz) * dt` — c'est cette **division par `dz` (10 m)** qui produit
  `1,2409e-4 / 10 = **1,2409e-5**`.
- C'est cette valeur **déjà divisée** qui a été copiée dans les notebooks 3D — la nôtre et
  l'intermédiaire `rapidmCDR_ALK_3D.ipynb`, dont le kernel confirme :
  `ALK_new[0] += fieldset.alkalinity_forcing * dt` — **sans** division supplémentaire.
- Le papier (Suselj et al. 2025, éq. 10 et §2.2.3) confirme la structure attendue : le forçage de
  surface doit être divisé par l'épaisseur de la couche supérieure (`Δz₁ = 10 m`, explicitement
  citée dans le papier) avant d'être ajouté comme terme source volumique — exactement ce que fait
  la division ci-dessus, **une seule fois**, déjà faite en amont.

👉 **Notre code (`c[0] += forcing·Δt`) est donc correct tel quel.** Diviser une seconde fois
(par erreur) reproduirait le bug — désormais identifié — qui affecte le kernel `VerticalConcentration`
d'origine si on le comparait à la version 2D (celui-ci ne divise pas non plus, cohérent).

*(Ancienne hypothèse, invalidée par cette vérification : on avait d'abord cru l'inverse — qu'il
manquait une division par `h₀`. C'était faux : la division existe déjà, en amont, dans la valeur
numérique elle-même.)*

## 5. 📥 Obtenir les données (la vraie étape)

### Ce qu'il faut
1. **PyCO2SYS** installé : `pip install PyCO2SYS`
2. Quatre champs de fond **ECCO-Darwin**, région ACC, 1995-01 (comme ALK/U/V/W) :
   `DIC`, `ALK`, `SALT`, `SST` — à télécharger et mettre au format NetCDF grille C
   (via `DatasLLC270/read-create_LLC270data.ipynb`, comme pour les autres champs).
3. Le **forçage mensuel** (`k_surf`, `siarea`) — fichier `monthly_mean_forcing_data_RapidmCDR`
   du modèle 1D de Suselj (variables : `k_diff`, `k_surf`, `siarea`, `area`, `alk_forcing`, …).

### Générer le fichier de sensibilité
Lancer `DatasLLC270/create_forcing_data.ipynb`, qui appelle
`compute_dpco2_sensitivity(DIC, ALK, SALT, SST)` (`data_functions.py`, ligne ~563 — un appel à
**PyCO2SYS** avec une perturbation de ±10 µmol/kg). Il produit :

```
sensitivity_pco2_1995-01.nc
├── dpco2_ov_ddic   (year_month, y, x)   atm / (mol C/kg)
├── dpco2_ov_dalk   (year_month, y, x)   atm / (mol ALK/kg)
├── Latitude, Longitude
```

### Activer la chimie
Déposer `sensitivity_pco2_1995-01.nc` dans `SEMES/data/`. Le notebook le détecte automatiquement
(`CHIMIE = os.path.exists(...)`) et passe en `CHIMIE = True`.
⚠️ Renseigner aussi les vraies valeurs de `k_surf` et `siarea` (aujourd'hui : placeholders
`3e-5` et `0.0` dans la cellule de setup).

## 6. Vérifier que ça marche

En mode chimie, à la fin :
```
Chimie du carbonate ACTIVEE (fichier : sensitivity_pco2_1995-01.nc)
...
flux CO2 : min=... max=...  (positif = absorption)
figure sauvegardee -> figures/flux_co2.png
```
Contrôles physiques (validés en synthétique) : `f_co2 > 0` (absorption), `DIC0 > 0` en surface et
qui diffuse vers le bas, `ALK0` toujours conservé, aucune concentration négative.

---

*Voir aussi : `SCHEMA_NUMERIQUE.md` (le schéma vertical, commun aux 2 traceurs),
`ALGO_3a_3b.md` (la simulation), `PARTIE_4_LES_PLOTS.md` (les tracés).*
