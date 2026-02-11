# Analytical Additive Release Extension  
Branch: `feature/analytical-additive-release`

This branch introduces an **optional additive release module** to FRAGMENT-MNP based on the analytical diffusion framework proposed by Ronan’s team.

The core polymer fragmentation–dissolution–mineralisation model remains **unchanged**.

---

# 1. Design Philosophy

The additive extension was implemented with the following principles:

• Preserve the existing ODE system  
• Avoid destabilising the solver  
• Maintain full backward compatibility  
• Keep additive physics modular and reviewable  

To achieve this, additive behaviour is implemented as a **post-processing operator** applied after the polymer ODE solution is computed.

---

# 2. What Was Added

## 2.1 Additive State Variables

When enabled, the model tracks:

- `A_part (N × T)`  
  Additive mass remaining in each particle size class  

- `A_aq (T)`  
  Additive mass released into the aqueous phase  

If additive inputs are not provided, the model behaves exactly as before.

---

## 2.2 Additive Transport During Fragmentation

At each timestep:

1. Polymer mass fragments according to the existing FSD matrix.
2. Additive mass is redistributed proportionally, assuming uniform mixing within each size class.

This guarantees strict mass conservation.

---

## 2.3 Analytical Diffusion Model

Additive release from particles is computed using spherical radial diffusion with Biot-number regime selection.

For particle radius \( r \):

Fourier number:

\[
Fo = \frac{D_p t}{r^2}
\]

Mass transfer and Biot number:

\[
K_m = \frac{K_{pw} D_w}{r}
\]

\[
Bi = \frac{K_m r}{D_p} = \frac{K_{pw} D_w}{D_p}
\]

Remaining fraction \( F = M/M_0 \) is computed by regime:

### Case I — External control (Bi ≤ 1)

\[
F = \exp(-Bi \cdot Fo)
\]

### Case II — Intermediate (1 < Bi < 100)

Smooth bridge approximation:

\[
F = \exp\left(-Fo \frac{Bi}{1+Bi}\right)
\]

### Case III — Internal diffusion control (Bi ≥ 100)

Classical diffusion-in-sphere series:

\[
F = \sum_{n=1}^{\infty}
\frac{6}{n^2\pi^2}
\exp(-n^2\pi^2 Fo)
\]

The series is truncated (default 50 terms).

Released fraction per timestep:

\[
f_{rel} = 1 - F
\]

---

# 3. New Inputs (Optional)

Additive tracking is activated only if both fields are provided:

```python
"initial_additive_concs": [ ... ],  # length = n_size_classes

"additive_release": {
    "model": "analytical",
    "params": {
        "D_p": 1e-16,
        "D_w": 1e-9,
        "K_pw": 1e4,
        "n_terms": 50
    }
}

# 4. Validation and Testing

Added tests verify:

Additive outputs are created when enabled

Total additive mass conservation:

Correct analytical limits:

Bi → 0

Bi → ∞

All existing polymer tests continue to pass.

# Important Notes

Polymer ODE system is unchanged.

Additive release is implemented via operator splitting.

The model assumes spherical particles.

Additive concentration is assumed uniform within each size class per timestep.

Future work could include:

Fully coupled additive ODE integration

Non-uniform internal concentration profiles

Time-dependent diffusivities
