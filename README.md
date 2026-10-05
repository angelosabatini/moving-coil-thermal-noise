# Thermal fluctuations and electronic readout using a moving-coil movement

MATLAB code accompanying the paper:

**Thermal fluctuations and electronic readout using a moving-coil movement**

The scripts reproduce the main numerical calculations presented in the paper. 
No external data files are required.

## Files

### `monte_carlo_deflection_acf.m`

Reproduces the Monte Carlo analysis and deflection-angle autocorrelation function 
shown in Fig. 2.

The script:

- reconstructs the parameters of the numerical moving-coil movement from its geometry;
- computes the critical-damping resistance and the relevant mechanical and electrical time scales;
- constructs the stationary covariance predicted by equipartition;
- propagates the complete third-order stochastic model using exact Van Loan discretization;
- initializes each realization directly from the equilibrium covariance, without a burn-in interval;
- compares the Monte Carlo deflection-angle autocorrelation function with the analytical result for the reduced critically damped model.

The random-number generator is initialized with a fixed seed for reproducibility.

### `active_readout_noise_budget.m`

Reproduces the numerical noise-budget calculations for the active voltage-to-current 
readout discussed in Sec. IV.

The script evaluates:

- the intrinsic thermal noise floor of the moving-coil movement;
- the white and 1/f noise contributions of the electronic readout;
- the equivalent noise bandwidth of the critically damped mechanical response;
- the input-referred noise for the three sense-resistor values considered in the paper;
- the corresponding full-scale-to-noise dynamic range and noise in parts per million of full scale.

## Numerical movement

The parameters used in the numerical example are derived from the geometry and 
material properties of the moving-coil movement wherever possible. In particular, 
the electromechanical coupling coefficient, moment of inertia, coil resistance, 
and coil inductance are reconstructed from the geometrical parameters specified 
in `monte_carlo_deflection_acf.m`.

The remaining mechanical parameters and the ambient temperature are specified 
explicitly in the script. The electrical termination required for critical damping 
is then calculated from the deterministic model rather than entered as a fitted 
parameter.

## Requirements

The scripts were developed in MATLAB.

`monte_carlo_deflection_acf.m` uses `xcorr` for the autocorrelation calculation and 
therefore requires the Signal Processing Toolbox.

No additional MATLAB toolboxes or external data files are required by 
`active_readout_noise_budget.m`.

## Usage

Run either script directly from MATLAB:

```matlab
monte_carlo_deflection_acf
```

or

```matlab
active_readout_noise_budget
```

The Monte Carlo script generates the autocorrelation figure and prints the numerical
checks reported in the paper. The active-readout script prints the quantities used
in the readout noise-budget table.

## Reproducibility

The Monte Carlo simulation uses a fixed random-number seed. The continuous-time
linear stochastic system is propagated between stored samples using exact Van Loan
discretization, and each realization is initialized directly from the analytical
equilibrium covariance.

The scripts require no measured or externally supplied data.

## Citation

If you use this code, please cite the accompanying paper:

> **A. M. Sabatini**,  
> **Thermal fluctuations and electronic readout using a moving-coil movement**,  
> *European Journal of Physics*, submitted.  