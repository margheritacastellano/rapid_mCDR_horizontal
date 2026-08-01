# Documentation du modèle 2D+1D

Les quatre documents à lire pour comprendre `test.ipynb`, dans l'ordre conseillé.

| # | document | ce qu'il explique |
|---|---|---|
| 1 | **[POURQUOI_LA_CASCADE.md](POURQUOI_LA_CASCADE.md)** | 🧠 **La conception.** Pourquoi chaque particule porte sa propre colonne, et pourquoi le champ partagé du modèle d'origine ne pouvait pas marcher (l'effet « cascade »). *À lire en premier : tout le reste en découle.* |
| 2 | **[SCHEMA_NUMERIQUE.md](SCHEMA_NUMERIQUE.md)** | 📐 **Les maths.** Les volumes finis, le maillage Arakawa C, pourquoi le schéma doit être implicite (D ≈ 7,7), la matrice tridiagonale `A·cⁿ⁺¹ = cⁿ`, et le kernel ligne par ligne. |
| 3 | **[ALGO_3a_3b.md](ALGO_3a_3b.md)** | ⚙️ **La simulation** (parties 3a et 3b). L'idée (1 jour = 1 particule), puis chaque ligne : la classe `Particle`, `pcol`, `pidx`, la boucle quotidienne, et les 3 pièges qui figeaient le modèle. |
| 4 | **[PARTIE_4_LES_PLOTS.md](PARTIE_4_LES_PLOTS.md)** | 📊 **Les tracés** (partie 4). Comment choisir le jour et la profondeur, où sont sauvées les figures, et comment lire les résultats. |
| 5 | **[CHIMIE_CO2.md](CHIMIE_CO2.md)** | 🧪 **La chimie du carbonate** (DIC + flux CO₂). La physique OAE→pCO₂→flux→DIC, comment c'est intégré, et **comment obtenir les données** pour l'activer. |
| 6 | **[EXTRAIRE_DONNEES_CO2.md](EXTRAIRE_DONNEES_CO2.md)** | 📥 **Extraire soi-même les données** de sensibilité pCO₂ depuis ECCO-Darwin (3 fichiers + un script), sans dépendre d'Océane. |

**[demo_cascade.py](demo_cascade.py)** — un modèle jouet (sans parcels, ~40 lignes) qui démontre
l'effet cascade en quelques secondes :
```bash
/opt/anaconda3/bin/python docs/demo_cascade.py
```

---

## Le modèle en trois phrases

1. On injecte une dose d'alcalinité **chaque jour** au même point (2°W, 49°S, dans l'ACC).
2. Chaque dose est portée par **une particule** qui dérive avec le courant (le **2D**), tout en
   diffusant vers le bas dans sa colonne d'eau (le **1D**).
3. Après 31 jours, on obtient une **traînée** de 31 doses, chacune plus ou moins loin et plus ou
   moins profonde selon son âge.

## Les trois chiffres à retenir

| | |
|---|---|
| **D ≈ 7,7** | le nombre de diffusion explicite → **>> 0,5** → le schéma **doit** être implicite |
| **~85 m** | la pénétration verticale avec diffusion (contre **~5 m** sans) |
| **1 550** | les nombres qui constituent l'état réel du modèle (`pcol` = 31 × 50), face aux 408 Mo du champ `ALK0` qui n'est qu'un **diagnostic** |

## Ce que le modèle ne fait PAS encore

La **chimie du carbonate** (DIC + flux CO₂) est **codée mais désactivée** faute de données de
sensibilité pCO₂ : voir **[CHIMIE_CO2.md](CHIMIE_CO2.md)**. Une fois les données téléchargées, le
modèle répondra à *« combien de CO₂ est absorbé ? »* (perspective n°3 de la présentation).
