# Quadratic Scalar-Field Dark Matter in CLASS v3.3.4

This repository contains a port of a scalar-field dark-matter (SFDM)
implementation originally developed with CLASS v2.6.3 to CLASS v3.3.4.

The currently validated model is one real scalar field with a quadratic
potential,

\[
V(\phi)=\frac{1}{2}m_\phi^2\phi^2.
\]

The repository also contains Cobaya input files for DESI DR2 BAO and
Pantheon supernova analyses.

## Project status

The code modification and tested CLASS–Cobaya integration are complete
for the one-field quadratic-SFDM configuration documented here.

Validated items:

- Native CLASS compilation
- SFDM parameter parsing
- Background evolution
- Linear perturbations and matter power spectrum
- Unlensed and lensed CMB spectra
- Comparison with the original CLASS v2.6.3 SFDM implementation
- Python `classy` interface
- Cobaya fixed-point evaluation
- DESI DR2 BAO and Pantheon likelihood evaluation
- MCMC smoke sampling
- Four-process MPI execution

The following scientific stages are not yet complete:

- Long production MCMC chains
- Final convergence assessment
- Marginalized posterior constraints
- Best-fit analysis and publication plots

Therefore, the current repository is a validated implementation and
sampling pipeline, not a final cosmological-constraint result.

## Validated Git revision

The final legacy-cutoff code revision is:

```text
Branch: sfdm-port-v3.3.4
Commit: b9235b48
Subject: Restore legacy SFDM theta cutoff
```

The immediate parent used when testing the final cutoff change was:

```text
f30cb9112376bf0e721162555bd4dd4cd1c3c13e
```

## Scalar-field model

For a homogeneous canonical scalar field,

\[
\ddot{\phi}+3H\dot{\phi}
+\frac{dV}{d\phi}=0.
\]

The energy density and pressure are

\[
\rho_\phi=
\frac{1}{2}\dot{\phi}^{\,2}+V(\phi),
\]

\[
p_\phi=
\frac{1}{2}\dot{\phi}^{\,2}-V(\phi),
\]

and its equation-of-state parameter is

\[
w_\phi=\frac{p_\phi}{\rho_\phi}.
\]

The SFDM implementation evolves polar-like background variables,
including `alpha_sfdm_1`, `theta_sfdm_1`, and `y1_sfdm_1`.

In the implemented variables, the SFDM pressure is evaluated as

\[
p_{\rm sfdm}
=
-C(\theta)\cos(\theta)\rho_{\rm sfdm},
\]

and consequently

\[
w_{\rm sfdm}(\theta)
=
-C(\theta)\cos(\theta).
\]

The angular evolution contains

\[
\frac{d\theta}{d\ln a}
=
-3C(\theta)\sin(\theta)+y_1.
\]

## Legacy oscillation cutoff

The original CLASS v2.6.3 SFDM implementation used

\[
C(\theta)
=
\frac{1}{2}
\left[
1-\tanh
\left(
\theta_{\rm tol}
\left(
\theta^2-\theta_{\rm thresh}^2
\right)
\right)
\right],
\]

with

\[
\theta_{\rm thresh}=100,
\qquad
\theta_{\rm tol}=1.
\]

Since \(\theta_{\rm tol}=1\), this becomes

\[
C(\theta)
=
\frac{1}{2}
\left[
1-\tanh(\theta^2-100^2)
\right].
\]

The modified helper functions are

\[
{\tt cos\_sfdm}(\theta)
=
C(\theta)\cos(\theta),
\]

\[
{\tt sin\_sfdm}(\theta)
=
C(\theta)\sin(\theta).
\]

For \(|\theta|\ll100\), \(C(\theta)\simeq1\), and the explicit
oscillations are retained.

Near \(|\theta|=100\), the cutoff changes smoothly.

For \(|\theta|\gg100\), \(C(\theta)\simeq0\), suppressing rapid
trigonometric oscillations.

## Derivative of the cutoff

Define

\[
T(\theta)=\tanh(\theta^2-100^2).
\]

Then

\[
C(\theta)=\frac{1}{2}[1-T(\theta)]
\]

and

\[
\frac{dC}{d\theta}
=
-\theta\left[1-T^2(\theta)\right].
\]

Because

\[
w(\theta)=-C(\theta)\cos(\theta),
\]

the required derivative is

\[
\frac{dw}{d\theta}
=
\theta
\left[
1-\tanh^2(\theta^2-100^2)
\right]
\cos(\theta)
+
C(\theta)\sin(\theta).
\]

The second term is implemented through `sin_sfdm()`.

The pressure derivative used in the background module is

\[
\frac{dp_{\rm sfdm}}{d\ln a}
=
\rho_{\rm sfdm}
\left[
\frac{dw}{d\theta}
\frac{d\theta}{d\ln a}
-
3w(1+w)
\right].
\]

## Final source implementation

The final modification is in:

```text
source/background.c
```

The cutoff used in `cos_sfdm()` and `sin_sfdm()` is:

```c
const double theta_thresh = 1.e2;

cutoff =
  0.5 * (
    1. -
    tanh(
      theta_sfdm * theta_sfdm
      - theta_thresh * theta_thresh
    )
  );
```

The derivative contribution is:

```c
cutoff_tanh_sfdm_1 =
  tanh(
    theta_sfdm_1 * theta_sfdm_1
    - 100. * 100.
  );

dw_dtheta_sfdm_1 =
  theta_sfdm_1
  * (
      1.
      - cutoff_tanh_sfdm_1
        * cutoff_tanh_sfdm_1
    )
  * cos(theta_sfdm_1)
  + sin_sfdm(pba, theta_sfdm_1);
```

This replaced the temporary v3.3.4-port cutoff

```c
0.5 * (1. - tanh(theta_sfdm - 30.*_PI_));
```

with the legacy squared-angle cutoff.

## CLASS input parameters

The validated quadratic-SFDM parameter block is:

```ini
Omega_sfdm_1 = 0.261205693250012
attractor_ic_sfdm_1 = yes
sfdm_parameters_1 = -22.0, 0.0, 1.e-2
sfdm_tuning_index_1 = 2
```

For the quadratic configuration used here:

- The first component of `sfdm_parameters_1` is
  \(\log_{10}(m_\phi/{\rm eV})\).
- The second component is fixed to `0.0`.
- The third component is `1.e-2`.
- `sfdm_tuning_index_1 = 2` identifies the shooting parameter.
- `attractor_ic_sfdm_1 = yes` activates the attractor initial
  conditions.

## Numerical mass domain

The earlier MontePython prior extended to

\[
-24\leq\log_{10}(m_\phi/{\rm eV})\leq-7.
\]

That complete range is not numerically valid for the current
initial-condition prescription at the tested cosmological point.

The boundary scan found:

```text
Largest tested passing value:  -8.90
Smallest tested failing value: -8.88
```

The failing region produced the CLASS initial-condition rejection that
the quadratic SFDM starts too late relative to the required initial
mass-to-Hubble criterion.

The final production prior is therefore

\[
-24\leq\log_{10}(m_\phi/{\rm eV})\leq-17.
\]

The upper limit `-17` is a conservative numerical-domain selection. It
is not an observational constraint.

## Native CLASS validation

The modified code successfully completed:

1. Background-only test
2. Perturbation and matter-power-spectrum test
3. CMB and lensing test

CLASS read all required SFDM parameters, and no SFDM input parameter was
reported as unused.

Representative outputs were:

```text
Age             = 13.770598 Gyr
Omega_sfdm      = 0.261206
sigma8          = approximately 0.82499
100 theta_s     = approximately 1.041785
```

## Old-to-new port comparison

The physical SFDM response was compared through

\[
R_{\rm old}(k)
=
\frac{P_{\rm SFDM}^{v2.6.3}(k)}
     {P_{\Lambda{\rm CDM}}^{v2.6.3}(k)},
\]

\[
R_{\rm new}(k)
=
\frac{P_{\rm SFDM}^{v3.3.4}(k)}
     {P_{\Lambda{\rm CDM}}^{v3.3.4}(k)}.
\]

The port residual was

\[
\Delta_{\rm port}(k)
=
100
\left[
\frac{R_{\rm new}(k)}
     {R_{\rm old}(k)}
-1
\right].
\]

For the high-precision comparison:

| Common wavenumber interval | Maximum absolute residual |
|---|---:|
| \(k<0.1\,h\,{\rm Mpc}^{-1}\) | 0.00317% |
| \(0.1\leq k\leq1\,h\,{\rm Mpc}^{-1}\) | 0.00525% |
| \(1\leq k\leq5.228\,h\,{\rm Mpc}^{-1}\) | 0.00627% |

Representative lensed-CMB port residuals were:

| Spectrum | Maximum absolute residual |
|---|---:|
| TT | approximately 0.0730% |
| EE | approximately 0.1233% |
| lensing \(\phi\phi\) | approximately 0.0469% |

The standard-precision high-\(k\) discrepancy decreased strongly after
increasing precision. This identifies the larger standard-setting
difference as a numerical-convergence effect rather than evidence of a
missing cutoff term.

## Python environment used for validation

```text
Python:   3.10.20
NumPy:    1.26.4
Cython:   3.2.9
classy:   3.3.4.0
Cobaya:   3.6.2
mpi4py:   4.1.2
MPI:      Open MPI 4.1.2
```

## Creating the Conda environment

```bash
conda env create -f environment_sfdm_v334.yml
conda activate sfdm_v334
```

Alternatively:

```bash
conda create -n sfdm_v334 python=3.10 -y
conda activate sfdm_v334

python -m pip install \
  "numpy==1.26.4" \
  "Cython==3.2.9" \
  "cobaya==3.6.2"
```

## Compiling CLASS

```bash
cd class_sfdm_v334

make clean
make -j2
```

## Installing the modified classy wrapper

The wrapper must be compiled from this repository:

```bash
export CC=gcc

python -m pip install \
  --no-build-isolation \
  --no-cache-dir \
  --force-reinstall \
  --no-deps \
  .
```

Verify the installation:

```bash
python - <<'PY'
import sys
import numpy
import classy
from classy import Class
from importlib.metadata import version

print("Python:", sys.executable)
print("NumPy:", numpy.__version__)
print("classy version:", version("classy"))
print("classy location:", classy.__file__)
print("Class interface:", Class)
PY
```

The displayed `classy` path must belong to the active SFDM environment.

## Installing Cobaya likelihood data

Choose a local likelihood-data directory:

```bash
export COBAYA_PACKAGES_PATH=/path/to/cobaya_packages
```

Install and test the likelihoods:

```bash
cobaya-install \
  bao.desi_dr2 \
  sn.pantheon \
  --packages-path "$COBAYA_PACKAGES_PATH" \
  --test
```

## Cobaya input files

The repository provides:

```text
cobaya_inputs/
├── sfdm_bao_dr2_pantheon_parameterized_evaluate.yaml
├── sfdm_bao_dr2_pantheon_mcmc_smoke.yaml
└── sfdm_bao_dr2_pantheon_production.yaml
```

The theory configuration uses

```yaml
theory:
  classy:
    path: global
```

This forces Cobaya to import the modified `classy` installed in the
active environment.

The parameter vector is dynamically assembled with

```yaml
sfdm_parameters_1:
  value: "lambda log10_m_sfdm: '%g, 0.0, 1.e-2' % log10_m_sfdm"
  derived: false
```

## Fixed-point likelihood test

```bash
conda activate sfdm_v334

export COBAYA_PACKAGES_PATH=/path/to/cobaya_packages

cobaya-run \
  cobaya_inputs/sfdm_bao_dr2_pantheon_parameterized_evaluate.yaml \
  -f
```

The validated reference evaluation returned approximately:

```text
chi2_bao.desi_dr2 = 20.1219
chi2_sn.pantheon  = 1035.05
Omega_m           = 0.311384
Omega_Lambda      = 0.688537
rs_drag           = 147.041 Mpc
age               = 13.7526 Gyr
```

These values are a pipeline test at one fixed cosmological point. They
are not a best fit or a posterior constraint.

## MCMC smoke test

```bash
cobaya-run \
  cobaya_inputs/sfdm_bao_dr2_pantheon_mcmc_smoke.yaml \
  -f
```

The validated smoke run produced:

```text
Accepted samples: 200
Acceptance rate:  approximately 0.178
Stored rows:      200
```

All four sampled parameters were finite, remained within their priors,
and moved away from the starting values.

The reported \(R-1\) value from this short run is not a convergence
result. The smoke test is only intended to verify sampling, parameter
movement, likelihood calls, and output writing.

## Production MPI run

After the evaluation and smoke test pass:

```bash
mpirun -np 4 cobaya-run \
  cobaya_inputs/sfdm_bao_dr2_pantheon_production.yaml \
  2>&1 |
tee sfdm_bao_dr2_pantheon_production.log
```

The production run should not be interpreted until all chains have
finished and convergence has been checked.

The current stopping requirements are:

```yaml
Rminus1_stop: 0.01
Rminus1_cl_stop: 0.2
Rminus1_cl_level: 0.95
```

## Rebuilding after moving the repository

CLASS embeds build-dependent paths. If the repository is moved or
extracted into a different directory, rebuild and reinstall it:

```bash
make clean
make -j2

python -m pip install \
  --no-build-isolation \
  --force-reinstall \
  --no-deps \
  .
```

## Important limitations

- The validated implementation covers the tested one-field quadratic
  model.
- It does not validate every scalar-field potential or two-field
  combination included in the historical code.
- CLASS v2.6.3 and CLASS v3.3.4 are not expected to produce bitwise
  identical results.
- Fixed-point likelihood values are not cosmological constraints.
- The 200-point chain is only a smoke test.
- Production posterior results are still pending.
- The mass upper prior of `-17` is a numerical-domain choice, not an
  observational bound.

## Upstream CLASS and citation

This repository is a modified research fork of CLASS and is not the
official upstream CLASS distribution.

Official CLASS repository:

https://github.com/lesgourg/class_public

CLASS home page:

https://class-code.net/

Users should cite the relevant CLASS papers, including:

D. Blas, J. Lesgourgues and T. Tram,  
*The Cosmic Linear Anisotropy Solving System (CLASS). Part II:
Approximation schemes*,  
JCAP 07 (2011) 034,  
https://arxiv.org/abs/1104.2933

Likelihood citations can be obtained with:

```bash
cobaya-bib bao.desi_dr2 sn.pantheon
```

The original upstream CLASS documentation has been preserved in:

```text
README_CLASS_UPSTREAM.md
```

## Handoff statement

The legacy-cutoff restoration, tested CLASS v3.3.4 port, Python wrapper,
Cobaya interface, BAO+Pantheon likelihood evaluation, MCMC smoke test,
and MPI import test are complete.

Converged production chains and the resulting scientific constraints
remain to be completed.
