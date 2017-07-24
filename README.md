# MVG-Moments
Compute all moments up to a given degree of a multivariate Gaussian distribution.

## Problem Statement
The considered distribution is the classical n-dimensional normalised multivariate Gaussian distribution:
- \rho(x) := (2 * \pi)^(- 0.5 * n) * (det(\Sigma))^(- 0.5) * \exp(- 0.5 * (x - \mu)' * \Sigma^{-1} * (x - \mu)).

## Usage
Let:
- `alpha_vals` be the full matrix of orders of the wanted moments (assumed to be a matrix, size s * n where s is the number of moments up to a maximum degree d),
- `MU` be the mean vector of rho,
- `SIGMA` be the covariance matrix of rho.

Call:
```matlab
CGMoms(alpha_vals, MU, SIGMA, verbose)
```
in order to compute the moments corresponding to the lines of `alpha_vals` of \rho.  

Alternatively, call:
```matlab
CGMoms_Kan(alpha, MU, SIGMA)
```
in order to compute only the moment of order `alpha` of \rho.

## Example
```matlab
n = 2;
MU = [0; 0];
SIGMA = [0.8 ^ 2, 0; 0.8 ^ 2, 0];
alpha_vals = [0, 0;
              1, 0;
              0, 1;
              2, 0;
              1, 1;
              0, 2;
              3, 0;
              2, 1;
              1, 2;
              0, 3];
moments = CGMoms(alpha_vals, MU, SIGMA, verbose);
```
After the previous execution, `moments` is now a 10 * 1 matrix containing the moments of rho parametrised by `MU` and `SIGMA`.
