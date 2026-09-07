# twoStep2Phase

`twoStep2Phase` wraps the supplied simulation code for the Biometrika paper
"An optimal two-step estimation approach for two-phase studies" into a small R
package.

The package exposes:

- simulation data generators for the linear, logistic, and Cox examples;
- one analysis function per outcome family. It automatically uses the optimal
  approach when no `Z` columns are present and the joint approach when `Z1`,
  `Z2`, and so on are present.

The Cox estimators call compiled Rcpp/RcppArmadillo routines adapted directly
from the supplied Cox C++ sources. Compiling the nested Cox risk-set and
conditional-score calculations substantially reduces the runtime of Cox
analyses and simulation studies.

## Example

```r
library(twoStep2Phase)

set.seed(1)
dat <- generate_linear_data(
  N = 1000,
  n = 400,
  beta = 0.5,
  mu = 0,
  Sigma = matrix(1),
  sde = 3.1,
  ns = c(80, 200, 120),
  quantile = c(0.3, 0.7)
)

fit <- update_linear(dat, stratum = "stratum", nboot = 0)
fit
```

```r
set.seed(2)
dat <- generate_logistic_data(
  N = 1000,
  n = 400,
  beta = 0.5,
  mu = 0,
  Sigma = matrix(1),
  sde = 3.1
)

update_logistic(dat, tune_para = 0.01)
```

```r
set.seed(3)
dat <- generate_cox_data(
  N = 1000,
  n = 400,
  beta = log(2),
  mu = 0,
  Sigma = matrix(1),
  sde = 3.1,
  tau = 2.5
)

update_cox(dat, nboot = 100)
```

```r
set.seed(4)
dat <- generate_linear_data(
  N = 1000, n = 400,
  beta = c(1, 0.5, -0.5, 0, 0),
  mu = rep(0, 5), Sigma = diag(5),
  sde = 3.79
)

fit <- update_linear(dat, nboot = 100)
fit$results
```

By default, an analysis without `Z` returns `complete_case`, `default`, and
`kernel`. An analysis with `Z` returns `complete_case`, `default`, `linear`,
and `joint_linear`. Use `x_models = c("linear", "gam")` for the GAM variants,
`x_models = NULL` to omit optional methods, and `n_quad` to set the number of
quadrature nodes. Calculated tuning quantities are in `fit$tune_parameters`.
Pass `z = NULL` or `z = character(0)` to force the no-Z analysis even when the
data frame contains columns named `Z1`, `Z2`, and so on.

Standard errors use the existing ordinary bootstrap by default. To use the
fully bootstrap method from page 10 of the article, which uses the same
bootstrap covariance for the update point estimate and its standard error,
set `var_method = "full_bootstrap"` and supply `nboot`, for example:

```r
fit <- update_linear(
  dat,
  nboot = 100,
  var_method = "full_bootstrap"
)
```

For the fully bootstrap method, `nboot` must be larger than the auxiliary
vector used by every requested update method. In practice, use substantially
more than this mathematical minimum for stable covariance estimates.
When `seed` is supplied to an update function, it is set once at the beginning
and the resulting random-number stream continues through tuning,
cross-validation, and bootstrap sampling without being reset.

For real two-phase data, pass a data frame with columns for the outcome,
expensive covariate observed in Phase II, auxiliary covariate observed for all
subjects, and complete-case indicator. The default column names are `Y`, `X`,
`X_star`, and `R`; Cox data also use `D` for the event indicator.
When `stratum = NULL`, all observations use one constant sampling stratum.
Pass `stratum = "stratum"` to use the column created by the package generators.
