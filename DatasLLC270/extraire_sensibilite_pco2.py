"""
Extraction des sensibilités pCO2 (dpco2/dDIC, dpco2/dALK) pour activer la chimie CO2
dans test.ipynb -- SANS retélécharger toute la grille ni MITgcmutils.

Ce qu'il faut AVANT de lancer :
  1) pip install PyCO2SYS tqdm
  2) Télécharger 3 binaires ECCO-Darwin pour le pas de temps 0000081144 (= 1995-01)
     et les déposer dans le dossier DOSSIER_BINAIRES ci-dessous :
        - DIC.0000081144.data       (~181 Mo)
        - SALTanom.0000081144.data  (~181 Mo)
        - SST.0000081144.data       (~4 Mo)
     (optionnel) SIarea.0000081144.data (~4 Mo)  -> pour la fraction de glace réelle
     Voir docs/EXTRAIRE_DONNEES_CO2.md pour les URLs exactes.

Ce que le script produit :
  data/sensitivity_pco2_1995-01.nc   (dpco2_ov_ddic, dpco2_ov_dalk, Latitude, Longitude)
  -> déposé directement au bon endroit ; test.ipynb passe alors en CHIMIE = True.

L'ALK de fond (absolue) est prise dans notre fichier existant data/ALK_1995-01.nc
(pas besoin de la retélécharger).
"""
import os, sys, numpy as np, xarray as xr

# ---- chemins ----
ICI = os.path.dirname(os.path.abspath(__file__))                 # .../DatasLLC270
DATA = os.path.abspath(os.path.join(ICI, "..", "..", "data"))    # .../SEMES/data
DOSSIER_BINAIRES = os.path.join(DATA, "ECCO-Darwin_data")        # où mettre les .data téléchargés
FILE = "0000081144"      # pas de temps ECCO-Darwin pour 1995-01
DTIME = "1995-01"
NX, NY, NZ = 270, 13 * 270, 50   # dimensions brutes LLC270 (tuiles empilées)
FILL = -1e10             # seuil terre / valeurs de remplissage -> NaN

sys.path.insert(0, ICI)
try:
    from data_functions import flat, compute_dpco2_sensitivity
except Exception as e:
    sys.exit(f"[ERREUR] Impossible d'importer data_functions ({e}).\n"
             f"  -> pip install PyCO2SYS tqdm  (data_functions importe PyCO2SYS)")


def lire_surface(nom):
    """Lit UNIQUEMENT le niveau de surface d'un binaire MDS big-endian float32,
    puis le remet en grille globale (945, 1080) via flat()."""
    chemin = os.path.join(DOSSIER_BINAIRES, f"{nom}.{FILE}.data")
    if not os.path.exists(chemin):
        sys.exit(f"[ERREUR] Fichier manquant : {chemin}\n"
                 f"  -> à télécharger (voir docs/EXTRAIRE_DONNEES_CO2.md).")
    # le niveau 0 = les NY*NX premiers flottants (les niveaux sont empilés)
    brut = np.fromfile(chemin, dtype=">f4", count=NY * NX).reshape(NY, NX)
    champ = np.asarray(flat(brut), dtype=np.float64)   # -> (945, 1080)
    champ[champ < FILL] = np.nan                        # masque terre
    return champ


def main():
    print("Lecture des champs de surface (DIC, SALT, SST) ...")
    DIC = lire_surface("DIC")          # µM C
    SALT = lire_surface("SALTanom")    # anomalie de salinité (compute ajoute +35)
    SST = lire_surface("SST")          # °C

    # ALK de fond (absolue) depuis notre fichier existant, niveau de surface
    alk_ds = xr.open_dataset(os.path.join(DATA, f"ALK_{DTIME}.nc"))
    ALK = alk_ds["ALK"].values[0, 0, :, :].astype(np.float64)   # (945, 1080), µM ALK
    ALK[ALK < FILL] = np.nan
    lon = alk_ds["Longitude"].values
    lat = alk_ds["Latitude"].values
    alk_ds.close()

    for nom, ch in [("DIC", DIC), ("SALT", SALT), ("SST", SST), ("ALK", ALK)]:
        print(f"  {nom:4s} shape={ch.shape}  min={np.nanmin(ch):.3g}  max={np.nanmax(ch):.3g}  "
              f"NaN={np.isnan(ch).mean()*100:.0f}%")

    print("Calcul des sensibilités pCO2 (PyCO2SYS, perturbation +-10 µmol/kg) ...")
    dpco2_ddic, dpco2_dalk = compute_dpco2_sensitivity(DIC, ALK, SALT, SST)

    # sanity: dpco2/dALK doit etre NEGATIF, dpco2/dDIC POSITIF (sur l'océan)
    print(f"  dpco2/dDIC : median={np.nanmedian(dpco2_ddic):.3g}  (attendu > 0)")
    print(f"  dpco2/dALK : median={np.nanmedian(dpco2_dalk):.3g}  (attendu < 0)")

    # ---- écriture au format attendu par create_forcing_data / test.ipynb ----
    ds = xr.Dataset(
        {
            "dpco2_ov_ddic": (("year_month", "y", "x"), dpco2_ddic[None, :, :],
                              {"long_name": "d pCO2/d DIC", "units": "atm / (mol C/kg)"}),
            "dpco2_ov_dalk": (("year_month", "y", "x"), dpco2_dalk[None, :, :],
                              {"long_name": "d pCO2/d ALK", "units": "atm / (mol ALK/kg)"}),
            "Latitude": (("y", "x"), lat),
            "Longitude": (("y", "x"), lon),
        },
        coords={"year_month": [DTIME]},
    )
    sortie = os.path.join(DATA, f"sensitivity_pco2_{DTIME}.nc")
    ds.to_netcdf(sortie)
    print(f"\nOK -> {sortie}")

    # ---- SIarea (optionnel) : fraction de glace réelle ----
    si = os.path.join(DOSSIER_BINAIRES, f"SIarea.{FILE}.data")
    if os.path.exists(si):
        siarea = np.asarray(flat(np.fromfile(si, dtype=">f4", count=NY * NX).reshape(NY, NX)), dtype=np.float64)
        siarea[siarea < FILL] = np.nan
        # valeur à la source (2°W, 49°S)
        iy, ix = np.unravel_index(np.nanargmin((lat + 49) ** 2 + (lon + 2) ** 2), lat.shape)
        print(f"SIarea à la source (ACC) = {siarea[iy, ix]:.3f}  (à reporter dans test.ipynb si != 0)")
    else:
        print("SIarea non téléchargé : siarea = 0 (raisonnable à l'ACC).")

    print("\nDépose ce fichier est déjà fait. Relance test.ipynb -> il détecte CHIMIE = True.")


if __name__ == "__main__":
    main()
