## Possible outcomes
O <- c('d', 't', 'g')

## External covariates CR
dim_ext_cr <- 3L
dim_covariates_cr <- 2L + dim_ext_cr

set.seed(123)
beta_covariates_d <- rnorm(dim_covariates_cr, 0, 1)
beta_covariates_t <- rnorm(dim_covariates_cr, 0, 1)
beta_covariates_g <- rnorm(dim_covariates_cr, 0, 1)

beta_covariates <- c(beta_covariates_d, beta_covariates_t, beta_covariates_g)

## Number of years
number_years_before <- 6L
number_years_after <- 2L

## Year-related intercept
set.seed(124)

# increasing intercepts for `d` and `t`
beta0_d <- rnorm(number_years_before, 0, 1)
beta0_t <- rnorm(number_years_before, 0, 1)
beta0_g <- rnorm(number_years_after, 0, 1)

# increasing intercepts for `g` ?
beta0_after <- rnorm(number_years_after, 0, 1) |> sort()

## parameters CR
beta <- c(beta_covariates, beta0_d, beta0_t, beta0_g)

extract_params_CR(
  THETA_CR = beta,
  DIM_EXT = dim_ext_cr,
  NYB = number_years_before,
  NYA = number_years_after
)

extract_params_cr(
  THETA_CR = beta,
  DIM_EXT = dim_ext_cr,
  NYB = number_years_before,
  NYA = number_years_after,
  OPTION = 4
)
beta0_d
### pars IRT ####
grades <- 4L
exams <- 3L

alpha <- rnorm(exams, 0, 1)
beta_exams <- list()
for(ex in 1:exams){
  beta_exams[[ex]] <- rnorm(grades, 0, 1)
}

theta_irt <- c(alpha, unlist(beta_exams))
extract_params_irt(
  THETA_IRT = theta_irt,
  N_GRADES = grades,
  N_EXAMS = exams,
  OPTION = 2,
  EXAM = 1
)
alpha
beta_exams[[1]]

extract_params(
  THETA = c(beta, theta_irt),
  DIM_EXT = dim_ext_cr,
  NYB = number_years_before,
  NYA = number_years_after,
  N_GRADES = grades,
  N_EXAMS = exams,
  PAR_TYPE = 2,
  OUTCOME = 'g',
  EXAM = NA
)
beta0_g
beta_covariates_g

set.seed(123)
x <- rnorm(dim_ext_cr+2, 0, 1)
year <- 7
hazard(
  OUTCOME = 2,
  YEAR = year,
  THETA_CR = beta,
  COVARIATES = x,
  NYB = number_years_before,
  NYA = number_years_after
)
expeta_d <- exp(beta0_d[year] + t(beta_covariates_d)%*%x)
expeta_t <- exp(beta0_t[year] + t(beta_covariates_t)%*%x)
expeta_d/(1+expeta_d+expeta_t); expeta_t/(1+expeta_d+expeta_t)

expeta_g <- exp(beta0_g[year] + t(beta_covariates_g)%*%x)
expeta_g/(1+expeta_g)





