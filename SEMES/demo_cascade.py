"""Demonstration minimale de la cascade, SANS parcels : juste la logique du champ partage."""
import numpy as np

def simule(n_part, n_jours, partage=True, cellules=12):
    """Particules alignees sur une ligne, avancant de 1 cellule tous les 2 jours
    (comme dans le vrai run : ~0.2 deg/jour, mailles ~0.33 deg)."""
    champ = np.zeros(cellules)          # le champ ALK0 PARTAGE (1 valeur par maille)
    pcol  = np.zeros(n_part)            # la colonne PROPRE a chaque particule
    cell  = np.full(n_part, -1)         # maille courante de chaque particule
    pos   = np.full(n_part, np.nan)     # position continue
    ne    = np.zeros(n_part, bool)      # particule nee ?
    for jour in range(n_jours):
        if jour < n_part:               # 1 nouvelle particule par jour, a la source
            ne[jour] = True; pos[jour] = 0.0
        for p in range(n_part):         # parcels traite les particules dans l'ordre de creation
            if not ne[p]: continue
            new = int(pos[p] // 2)      # sa maille
            if partage:
                if cell[p] == -1:  champ[new] += 1.0          # injection de sa dose
                elif new != cell[p]:
                    champ[new] += champ[cell[p]]              # <-- HorizontalAdvConcentration
                    champ[cell[p]] = 0.0
            else:
                if cell[p] == -1:  pcol[p] += 1.0             # la dose reste SUR la particule
            cell[p] = new
            pos[p] += 1.0               # avance
    if partage: return champ
    champ2 = np.zeros(cellules)
    for p in range(n_part):
        if ne[p]: champ2[cell[p]] += pcol[p]                  # reconstruction
    return champ2

for n in [1, 2, 8]:
    a = simule(n, 16, partage=True)
    b = simule(n, 16, partage=False)
    print(f"--- {n} particule(s) ---")
    print(f"   champ PARTAGE (Oceane)   : {np.round(a,2)}   total={a.sum():.0f}")
    print(f"   colonnes PAR PARTICULE   : {np.round(b,2)}   total={b.sum():.0f}")
    print(f"   -> max doses sur 1 maille : partage={a.max():.0f}  vs  par-particule={b.max():.0f}")
    print()
