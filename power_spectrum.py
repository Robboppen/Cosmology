import numpy as np
import matplotlib.pyplot as plt
import camb
from classy import Class

# Parámetros Planck 2018
h, ombh2, omch2, As, ns = 0.6736, 0.02237, 0.12, 2.1e-9, 0.9649
redshifts = [0, 1, 2]

# ---------- CAMB ----------
p = camb.set_params(H0=100*h, ombh2=ombh2, omch2=omch2, As=As, ns=ns)
p.set_matter_power(redshifts=redshifts, kmax=10)
r = camb.get_results(p)
kh, z, pk = r.get_matter_power_spectrum(minkh=1e-4, maxkh=10, npoints=400)
# kh en h/Mpc, pk[i] en (Mpc/h)^3 para z[i]

# ---------- CLASS ----------
c = Class()
c.set({'h': h, 'omega_b': ombh2, 'omega_cdm': omch2, 'A_s': As, 'n_s': ns,
       'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06,
       'output': 'mPk', 'P_k_max_h/Mpc': 12, 'z_max_pk': 2})
c.compute()

k_inv_Mpc = kh * h            # h/Mpc -> 1/Mpc
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for i, zi in enumerate(z):
    # CLASS da P en Mpc^3 -> multiplicar por h^3 para (Mpc/h)^3
    Pk_class = np.array([c.pk_lin(k, zi) for k in k_inv_Mpc]) * h**3
    ax[0].loglog(kh, pk[i], label=f"CAMB z={zi:g}")
    ax[0].loglog(kh, Pk_class, 'k:', lw=1)
    ax[1].semilogx(kh, Pk_class/pk[i] - 1, label=f"z={zi:g}")

print("sigma8 CAMB :", r.get_sigma8_0())
print("sigma8 CLASS:", c.sigma8())
print("D(z=1)/D(0) :", c.scale_independent_growth_factor(1.0))

ax[0].set(xlabel=r"$k\,[h/\mathrm{Mpc}]$", ylabel=r"$P_m\,[(\mathrm{Mpc}/h)^3]$",
          title="CAMB (color) vs CLASS (punteado)")
ax[1].set(xlabel=r"$k\,[h/\mathrm{Mpc}]$", ylabel="CLASS/CAMB − 1", ylim=(-0.01, 0.01))
ax[0].legend(); ax[1].legend()
plt.tight_layout(); plt.show()

c.struct_cleanup()