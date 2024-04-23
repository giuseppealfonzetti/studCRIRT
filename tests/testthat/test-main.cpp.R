## change the seed to try different draws
seed <- 123

## Setup dimensions
external_covariates <- 1L     # number of external covariates
years_before <- 6L            # number of possible years without graduation risk
years_after <- 2L             # number of possible years  with graduation risk
grades <- 4L                  # number of grades
exams <- 10L                  # number of exams

# dimension of the parameter vector related to the competing risk model
dim_cr <- 3*(external_covariates+2) + 2*(years_before) + years_after

# dimension of the parameter vector related to the irt model
dim_irt <- 3*exams + exams*grades

# dimension of the parameter vector related to latent variables
dim_lat <- 2

# draw a random parameter vector of the suitable length
set.seed(seed)
theta <- rnorm(dim_cr+dim_irt+dim_lat)

# To inspect the parameters we can convert theta to a named list
params_list <- paramsVec2list(THETA = theta,
                              DIM_EXT = external_covariates,
                              NYB = years_before,
                              NYA = years_after,
                              N_GRADES = grades,
                              N_EXAMS = exams)

# here you can manually change parameters (chek the list names)
# I do that just to simulate reasonable data
params_list[['LAT']][['Corr']] <- .5
params_list[['LAT']][['Speed_variability']] <- 1.5
params_list[['IRT']][['Exams_average_time']] <- sapply(sort(rep(1:5, (round(exams/6)+1)))[1:exams], function(x) log(x*180))
params_list[['IRT']][['Exams_variability_time']] <- runif(exams, .5,2)
params_list$IRT$Exams_slopes <- runif(exams, 0.5,1)
params_list$CR$Graduation$Slope_ability <- 2
params_list$CR$Graduation$Slope_speed <- 2

# and project them back on a new theta
theta <- paramsList2vec(PARAMS_LIST = params_list,
                        DIM_EXT = external_covariates,
                        NYB = years_before,
                        NYA = years_after,
                        N_GRADES = grades,
                        N_EXAMS = exams)

## read latent parameters from list
rho <- params_list[['LAT']][['Corr']]
sig <- params_list[['LAT']][['Speed_variability']]

## latent covariance matrix
S <- matrix(c(1,rho*sig,rho*sig,sig^2), 2, 2)

## draw ability, speed and external covariates value
set.seed(seed+2)
lat <- mvtnorm::rmvnorm(1, sigma = S) # lat[1] = ability, lat[2] = speed
# lat <- matrix(c(.5,.5), ncol = 2)
x <- rnorm(external_covariates)       # external covariates

## draw grades and times for all exams
grades_vec <- rngGrades(
  THETA_IRT = theta[(dim_cr+1):(dim_cr+dim_irt)],
  N_GRADES = grades,
  N_EXAMS = exams,
  ABILITY = lat[1],
  SEED = seed
)$grades
times_vec <- rngTimes(
  THETA_IRT = theta[(dim_cr+1):(dim_cr+dim_irt)],
  N_GRADES = grades,
  N_EXAMS = exams,
  SPEED = lat[2],
  SEED = seed
)

## Read the year of the last exam
yle <- round(max(times_vec)/365,0)
yle

# grid for ability and speed divided by outcome and year
dt <- expand.grid(
  ability = seq(-1,1, by=1),
  speed = seq(-1,1, by=1),
  outcome = c(0,1,2,3),
  year = 1:(yle+1)
)

for (r in 1:nrow(dt)) {
  # use only exams observed up to `year`
  maxTime <- 365*dt$year[r]
  obs_vec <- grades_vec>0&times_vec<maxTime

  test_that('check complete data likelihood', {
    v1 <- complete_likelihood(
      THETA = theta,
      EXTCOVARIATES = x,
      EXAMS_GRADES = grades_vec,
      EXAMS_DAYS = times_vec,
      EXAMS_OBSFLAG = obs_vec,
      OUTCOME = dt$outcome[r],
      YEAR = dt$year[r],
      N_GRADES = grades,
      N_EXAMS = exams,
      NYB = years_before,
      NYA = years_after,
      ABILITY = dt$ability[r],
      SPEED = dt$speed[r],
      YEAR_LAST_EXAM = yle,
      LOGFLAG = TRUE
    )

    v2 <- complete_likelihood(
      THETA = theta,
      EXTCOVARIATES = x,
      EXAMS_GRADES = grades_vec,
      EXAMS_DAYS = times_vec,
      EXAMS_OBSFLAG = obs_vec,
      OUTCOME = dt$outcome[r],
      YEAR = dt$year[r],
      N_GRADES = grades,
      N_EXAMS = exams,
      NYB = years_before,
      NYA = years_after,
      ABILITY = dt$ability[r],
      SPEED = dt$speed[r],
      YEAR_LAST_EXAM = yle,
      LOGFLAG = FALSE
    )

    expect_equal(exp(v1), v2)
  })
}

# r <- 1
# maxTime <- 365*dt$year[r]
# obs_vec <- grades_vec>0&times_vec<maxTime
# Rfun <- function(LAT){
#   complete_likelihood(
#     THETA = theta,
#     EXTCOVARIATES = x,
#     EXAMS_GRADES = grades_vec,
#     EXAMS_DAYS = times_vec,
#     EXAMS_OBSFLAG = obs_vec,
#     OUTCOME = dt$outcome[r],
#     YEAR = dt$year[r],
#     N_GRADES = grades,
#     N_EXAMS = exams,
#     NYB = years_before,
#     NYA = years_after,
#     ABILITY = LAT[1],
#     SPEED = LAT[2],
#     YEAR_LAST_EXAM = yle,
#     LOGFLAG = FALSE
#   )
# }
# marginal(
#   THETA = theta,
#   EXTCOVARIATES = x,
#   EXAMS_GRADES = grades_vec,
#   EXAMS_DAYS = times_vec,
#   EXAMS_OBSFLAG = obs_vec,
#   OUTCOME = dt$outcome[r],
#   YEAR = dt$year[r],
#   N_GRADES = grades,
#   N_EXAMS = exams,
#   NYB = years_before,
#   NYA = years_after,
#   YEAR_LAST_EXAM = yle,
#   LOGFLAG = FALSE
# )



# cubature::cubintegrate(f= Rfun, lower = c(-Inf, -Inf), upper = c(Inf, Inf))
# marginal_likelihood(
#   THETA = theta,
#   EXTCOVARIATES = x,
#   EXAMS_GRADES = grades_vec,
#   EXAMS_DAYS = times_vec,
#   EXAMS_OBSFLAG = obs_vec,
#   OUTCOME = dt$outcome[r],
#   YEAR = dt$year[r],
#   N_GRADES = grades,
#   N_EXAMS = exams,
#   NYB = years_before,
#   NYA = years_after,
#   YEAR_LAST_EXAM = yle,
#   LOGFLAG = FALSE,
#   CUBATURE_METHOD = 'hcubature'
# )
