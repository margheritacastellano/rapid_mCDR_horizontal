# Le schéma numérique : maths & code

Documentation du kernel `VerticalAdvDiffImplicit` (cellule **[12]** de `test.ipynb`) :
**advection + diffusion verticales**, résolues par **volumes finis 1D** en **Euler implicite**.

---

## 1. Le problème physique

On transporte $c$ = **la perturbation d'alcalinité** (`ALK0`, en µM), qui part de zéro partout.
Ce n'est ni du CO₂ ni du DIC : c'est un **traceur passif**.

$$\boxed{\;\partial_t c \;+\; \underbrace{\partial_z(w\,c)}_{\text{advection}} \;-\; \underbrace{\partial_z\!\left(K(z)\,\partial_z c\right)}_{\text{diffusion}} \;=\; f_{source}\;}$$

| symbole | sens | dans le code |
|---|---|---|
| $c$ | perturbation d'alcalinité (µM) | `C` = `fieldset.pcol[pidx]` |
| $w$ | vitesse verticale (m/s, **>0 vers le haut**) | `fieldset.W` |
| $K(z)$ | diffusivité verticale turbulente (m²/s) | `Kint_fv` |
| $f_{source}$ | injection en surface | `alkalinity_forcing` |

---

## 2. Le maillage vertical (grille Arakawa C)

C'est **la** clé du schéma : les données fournissent **deux** grilles verticales distinctes.

```
      surface
  ═══════════════  ← interface 0        (face,   RF[0] = 0 m)
        • c_0                            (centre, RC[0] = 5 m)     ┐ h_0
  ───────────────  ← interface 1  w_1    (face,   RF[1] = 10 m)    ┘
        • c_1                            (centre, RC[1] = 15 m)
  ───────────────  ← interface 2  w_2
        • c_2
         ⋮
  ═══════════════  ← interface N_z       (fond)
```

- **Traceurs** ($c$) → aux **centres** : `RC = ALK.grid.depth` = $[5, 15, 25, \dots]$
- **Vitesse verticale** ($w$) → aux **faces** : `RF = W.grid.depth` = $[0, 10, 20, \dots]$

On en déduit les **deux** longueurs dont les volumes finis ont besoin :

$$h_k = 2\,(RC_k - RF_k) \qquad\text{(épaisseur de la maille } k\text{, le « volume »)}$$
$$\delta_i = RC_i - RC_{i-1} \qquad\text{(distance centre-à-centre à l'interface } i\text{)}$$

```python
h_fv     = 2.0 * (RC - RF)   # (50,)  -> le "volume" de chaque maille
delta_fv = np.diff(RC)       # (49,)  -> 50 mailles => 49 interfaces intérieures
```

> ⚠️ Ne pas confondre : $h_k$ sert à **diviser** (bilan de volume), $\delta_i$ sert au **gradient** de diffusion. Ce sont deux distances différentes dès que le maillage est irrégulier — et il l'est (de 10 m en surface à ~450 m au fond).

**Convention d'indices** : $k$ croît vers le **bas**. L'**interface $i$** est la face **supérieure** de la maille $i$ : elle sépare la maille $i-1$ (au-dessus) de la maille $i$ (en-dessous). Donc $w_i$ = `W[i]` s'utilise **directement**, sans moyenne.

---

## 3. La discrétisation en volumes finis

Le principe : on intègre l'équation sur chaque maille. Le taux de variation du contenu = **ce qui entre − ce qui sort**.

$$h_k \,\frac{c_k^{n+1}-c_k^{n}}{\Delta t} \;=\; F_k \;-\; F_{k+1}$$

où $F_i$ est le **flux à travers l'interface $i$, compté positif vers le bas** :
- $F_k$ : ce qui **entre** par le haut de la maille $k$
- $F_{k+1}$ : ce qui **sort** par le bas

### 3.a Flux advectif (décentré amont / *upwind*)

Comme $w>0$ pointe vers le **haut**, la vitesse **descendante** vaut $-w_i$. On sépare les deux sens :

$$a^{\downarrow}_i = \max(-w_i,\,0), \qquad a^{\uparrow}_i = \max(w_i,\,0)$$

$$F^{adv}_i = a^{\downarrow}_i\, c_{i-1} \;-\; a^{\uparrow}_i\, c_i$$

*Lecture* : si ça **descend**, on transporte la maille **du dessus** ($c_{i-1}$) vers le bas. Si ça **monte**, on transporte celle **du dessous** ($c_i$) vers le haut → contribution **négative** au flux descendant. On prend toujours la valeur **en amont** du courant : c'est ce qui rend le schéma stable et positif.

### 3.b Flux diffusif (loi de Fick)

La diffusion va du **fort** vers le **faible** :

$$F^{dif}_i = \frac{K_i}{\delta_i}\left(c_{i-1} - c_i\right), \qquad K_i = \tfrac12\left(K_{i-1}+K_i\right)$$

### 3.c Forme compacte

En regroupant, le flux total à l'interface $i$ est **linéaire** en $c$ :

$$\boxed{\;F_i = L_i\, c_{i-1} \;-\; U_i\, c_i\;}
\qquad\text{avec}\qquad
\begin{cases}
L_i = a^{\downarrow}_i + \dfrac{K_i}{\delta_i} & \text{(coeff. de la maille du dessus)}\\[2mm]
U_i = a^{\uparrow}_i + \dfrac{K_i}{\delta_i} & \text{(coeff. de la maille du dessous)}
\end{cases}$$

```python
for i in range(1, Nz):
    a_down[i] = max(-W[i], 0.0)
    a_up[i]   = max( W[i], 0.0)
    Kd[i]     = Kint[i] / delta[i - 1]
L = a_down + Kd     # L_i
U = a_up   + Kd     # U_i
```

> 💡 `delta[i-1]` et non `delta[i]` : `delta` a 49 cases (0…48) pour les interfaces 1…49.

---

## 4. Pourquoi **implicite** ? (le cœur du problème)

Un schéma **explicite** (flux évalués à $n$) est stable seulement si le **nombre de diffusion** vérifie

$$D = \frac{K\,\Delta t}{\Delta z^{2}} \;\le\; \frac12$$

Or ici, près de la surface :

$$D = \frac{8{,}91\times10^{-3} \times 86400}{10^{2}} \;\approx\; \boxed{7{,}7} \;\ggg\; 0{,}5$$

👉 **L'explicite diverge** à $\Delta t = 1$ jour. Deux issues :
- réduire $\Delta t$ à ~1,5 h (≈ 16× plus de pas… et ça casse la logique d'injection quotidienne) ;
- **passer en implicite** ← le choix fait ici.

En **Euler rétrograde**, on évalue les flux à l'instant $n+1$ (l'inconnue) :

$$h_k \,\frac{c_k^{n+1}-c_k^{n}}{\Delta t} = \Big(L_k c_{k-1}^{n+1} - U_k c_k^{n+1}\Big) - \Big(L_{k+1} c_k^{n+1} - U_{k+1} c_{k+1}^{n+1}\Big)$$

C'est **inconditionnellement stable** : n'importe quel $\Delta t$ passe.

---

## 5. Le système tridiagonal

On réarrange en posant $r_k = \dfrac{\Delta t}{h_k}$ et en mettant les inconnues à gauche :

$$\boxed{\;-\,r_k L_k\; c_{k-1}^{n+1} \;+\; \Big[1 + r_k\big(U_k + L_{k+1}\big)\Big]\, c_k^{n+1} \;-\; r_k U_{k+1}\; c_{k+1}^{n+1} \;=\; c_k^{n}\;}$$

Chaque ligne ne touche que $c_{k-1}, c_k, c_{k+1}$ → la matrice $A$ est **tridiagonale** :

$$A_{k,k-1} = -r_k L_k, \qquad A_{k,k} = 1 + r_k(U_k + L_{k+1}), \qquad A_{k,k+1} = -r_k U_{k+1}$$

et on résout simplement

$$A\, c^{n+1} = c^{n}$$

### Conditions aux limites : **flux nul**

- **Surface** : rien ne traverse l'interface 0 → $L_0 = U_0 = 0$
- **Fond** : rien ne traverse l'interface $N_z$ → $L_{N_z} = U_{N_z} = 0$

```python
Uk  = U[k]     if k >= 1      else 0.0   # interface du haut  -> 0 en surface
Lk1 = L[k + 1] if k <= Nz - 2 else 0.0   # interface du bas   -> 0 au fond
```

---

## 6. Les 3 propriétés qu'on obtient gratuitement

### 6.a Conservation exacte de la masse
La **forme flux** garantit que ce qui sort de la maille $k$ par l'interface $k+1$ **entre exactement** dans la maille $k+1$ (même $F_{k+1}$, signe opposé). Avec des flux nuls aux bords, tout se télescope :

$$\sum_k h_k\, c_k^{n+1} \;=\; \sum_k h_k\, c_k^{n}$$

> ⚠️ L'invariant est $\sum_k h_k c_k$ (**pondéré par l'épaisseur**), **pas** $\sum_k c_k$ — puisque les mailles n'ont pas la même taille. *Vérifié numériquement : erreur relative ~$10^{-15}$.*

### 6.b Positivité
$A$ est une **M-matrice** : diagonale $>0$, extra-diagonaux $\le 0$, à diagonale dominante (car $L_i, U_i \ge 0$). Donc $c^n \ge 0 \Rightarrow c^{n+1} \ge 0$ : **jamais de concentration négative**.

C'est l'*upwind* qui offre ça — un schéma centré ne le garantirait pas.

### 6.c Stabilité inconditionnelle
Aucune contrainte sur $\Delta t$. C'est tout l'intérêt face à $D \approx 7{,}7$.

---

## 7. Le format `solve_banded` (le piège d'implémentation)

`scipy.linalg.solve_banded((1,1), ab, b)` veut la matrice **empilée en bandes**, `ab` de forme (3, Nz), avec la convention :

$$\texttt{ab[0, j]} = A_{j-1,\,j} \quad\text{(sur-diagonale)}, \qquad
\texttt{ab[1, j]} = A_{j,\,j} \quad\text{(diagonale)}, \qquad
\texttt{ab[2, j]} = A_{j+1,\,j} \quad\text{(sous-diagonale)}$$

D'où le décalage d'indices, qui surprend à la lecture :

```python
ab[1, k]     = 1.0 + r * (Uk + Lk1)      # A[k,k]     -> colonne k
ab[2, k - 1] = -(dt / h[k]) * L[k]       # A[k,k-1]   -> rangé en colonne k-1 !
ab[0, k + 1] = -(dt / h[k]) * U[k + 1]   # A[k,k+1]   -> rangé en colonne k+1 !
```

---

## 8. Le code, ligne par ligne

```python
def VerticalAdvDiffImplicit(particle, fieldset, time):
```

| lignes | ce que ça fait |
|---|---|
| **25–26** | `dt = particle.dt` ; `_ = fieldset.UV[particle]` → force parcels à localiser la particule (met à jour ses indices de maille) |
| **28–37** | récupère la maille `(xi, yi)` et **corrige** l'indice si une maille voisine est plus proche (grille curviligne) |
| **38–39** | mémorise `cell_x/cell_y` → servira à **reconstruire** le champ `ALK0` |
| **41** | `Nz = fieldset.ALK0.data.shape[1]` → **seul** usage d'`ALK0` : récupérer 50 |
| **42** | ⭐ `C = fieldset.pcol[particle.pidx].copy()` → **la colonne propre à cette particule** |
| **43** | `W` = la colonne de vitesse verticale **à la position courante** |
| **49–50** | garde-fous : `NaN` / valeurs de remplissage sous le plancher → $w=0$ (donc flux nul) |
| **53–55** | 💧 **l'injection** : le jour du lâcher, on ajoute la dose en surface, puis `released = 1` (⇒ **une** dose par particule) |
| **58–66** | les coefficients d'interface $a^{\downarrow}, a^{\uparrow}, K/\delta$ → $L_i, U_i$ |
| **69–78** | l'assemblage de la **matrice tridiagonale** (+ conditions de flux nul) |
| **80** | ⭐ `fieldset.pcol[pidx] = solve_banded((1,1), ab, C)` → on résout et on **réécrit la colonne** |

### Un détail d'ordre qui compte
L'injection (l. 53) se fait **avant** le solve (l. 80) : la dose du jour est donc **transportée dès le jour même**. Le second membre du système est bien $c^n + \text{dose}$.

---

## 9. Le couplage 2D + 1D

Le kernel ci-dessus, c'est le **1D vertical**. Le **2D horizontal** est ailleurs :

```python
kernels = pset.Kernel(VerticalAdvDiffImplicit) + pset.Kernel(parcels.AdvectionRK4)
```

- `VerticalAdvDiffImplicit` → fait évoluer **`pcol[pidx]`** (la colonne)
- `parcels.AdvectionRK4` → déplace **la particule** (RK4 sur $u, v$)

**La colonne suit la particule automatiquement**, puisqu'elle lui est attachée par `pidx`. Il n'y a **aucun échange entre mailles** — c'est précisément ce qui supprime la *cascade* (avec un champ **partagé**, les échanges séquentiels entre particules relayaient une colonne vers l'avant et empilaient l'alcalinité sur la maille de tête).

Le champ gridé n'est **reconstruit** que pour les sorties :

$$\text{ALK0}[:,\,y,\,x] \;=\!\!\sum_{p\ \text{dans la maille}\ (x,y)}\!\! \text{pcol}[p]$$

```python
fieldset.ALK0.data[:] = 0.0
for p in pset:
    fieldset.ALK0.data[0, :, p.cell_y, p.cell_x] += fieldset.pcol[p.pidx]
```

> 🎯 **À retenir** : l'état réel du modèle, c'est **`pcol` (31 × 50 = 1 550 nombres)**. `ALK0` (408 Mo) n'est qu'un **diagnostic** reconstruit à la demande.

---

## 10. Ce que ça donne (31 jours, ACC)

| Diagnostic | Résultat |
|---|---|
| Conservation $\sum_k h_k c_k$ | = $31 \times$ la dose, erreur relative ~$10^{-15}$ |
| Concentrations négatives | aucune |
| Pénétration **advection seule** ($K=0$) | ~5 m (tout reste piégé en surface : $w$ est négligeable) |
| Pénétration **advection + diffusion** | **~75–85 m** |
| Dérive horizontale | ~0,2 °/jour vers l'est le long de l'ACC |

👉 C'est **la diffusion** qui fait descendre l'alcalinité, pas l'advection verticale.

---

## 11. Les limites actuelles (honnêtes)

1. **Pas de chimie du carbonate** : ni DIC, ni pCO₂, ni flux air–mer. Le modèle dit *où va l'alcalinité*, pas *combien de CO₂ est absorbé*. (= perspective n°3 de la présentation)
2. **Courants figés** : un seul pas de temps (moyenne mensuelle) réutilisé pour les 31 jours (`allow_time_extrapolation=True`).
3. **Hypothèse forte** : les mêmes courants horizontaux à toutes les profondeurs (toute la colonne suit la trajectoire de surface).
4. **Pas de diffusion horizontale** (supposée négligeable devant l'advection).
5. L'alcalinité **de fond** (`fieldset.ALK`) est chargée mais **jamais lue** — le transport de la perturbation est linéaire. Elle ne servira qu'avec la chimie.
