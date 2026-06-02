"""Canonical fitted optical-response parameters extracted from the project's existing fit outputs.

These define the SiC-like and Si-like 'truth' systems used to generate synthetic spectra in the
Phase-3 studies, so the numerical experiments are anchored to the actual measured systems rather than
arbitrary toy values. Provenance is recorded per block.

QA NOTE: these are read-only fixtures copied from committed JSON fit outputs; they are inputs to studies,
not themselves manuscript claims. The study scripts re-fit / re-derive where a claim is made.
"""

# --- SiC epilayer on SiC substrate (4H), MDF/L-STC + Drude ---
# Provenance: a project fitted-parameter reference (hybrid two-beam fit, ~10deg basis).
SIC_4H = dict(
    eps_inf=6.56,
    wT_cm=798.0,           # TO phonon (4H-SiC), from the reference SiC model fit_stack_enhanced
    wL_cm=970.0,           # LO phonon (4H-SiC)
    # film
    log10_N_film=17.41905524575606,
    mu_film_cm2_Vs=799.9995828006565,
    GT_film_cm=3.9999999566557145,
    GL_film_cm=4.959960656896943,
    # substrate
    log10_N_sub=18.579150765844922,
    mu_sub_cm2_Vs=67.94410789501784,
    GT_sub_cm=3.9999943007300383,
    GL_sub_cm=8.283688785008353,
    # nominal thickness + reported uncertainty (for cross-check, NOT a manuscript claim here)
    d_um=7.426834835410501,
    d_um_std_error=0.014297631973016136,   # Jacobian std error from the source fit
    d_um_ci95=(7.3947880841078915, 7.458881586713111),
)

# --- Si epilayer stack: air / SiO2 / Si-epi / Si-substrate, Drude + Lorentz-SiO2 + empirical HF correction ---
# Provenance: a project joint-fit reference (joint 10deg/15deg NLLS fit).
SI_STACK = dict(
    d_um=4.336448662187574,
    t_ox_nm=5.999999999999999,
    N_epi=5.2395641826102675e17,
    mu_epi=599.9999999999999,
    N_sub=5.421882686423222e19,
    mu_sub=70.68243993678816,
    eps_inf_e=11.828660255994379,
    eps_inf_s=13.299999999999999,
    # SiO2 double-Lorentz
    S1=1.4999999999999998, w01=1066.1744070575442, g1=45.03211367879835,
    S2=0.5999999999999999, w02=814.4138700011705, g2=20.000000000000004,
    # empirical HF correction / dispersion (the reference implements the HF term in Re(epsilon))
    A_disp=2.3175950429433616e-27,
    B_damp=0.49999999999999994,
    B_width=582.6750663227855,
    sigma_theta_deg=0.020000000000000004,
    d_jitter_um=0.12018086822232724,
    scales=(0.9556948008761972, 1.0583880609632075),
    d_um_ci95=(4.299878683949063, 4.373018640426085),
)
