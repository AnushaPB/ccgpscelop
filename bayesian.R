install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"))
library(INLA)
library(sf)
library(sp)
model_df <- read_csv(here("analysis/genetic_diversity/outputs/model_df.csv"))

# Convert to sf and sp objects (sp needed for semivariograms)
model_sf <- st_as_sf(model_df, coords = c("x", "y"), crs = 4326) %>% st_transform(3310)

# Convert sf to SpatialPoints
sp_data <- as(model_sf, "Spatial")

# https://www.paulamoraga.com/book-spatial/bayesian-spatial-models.html

# Create neighbors list based on range value from semivariogram
coords <- st_coordinates(model_sf)
nb <- dnearneigh(coords, d1 = 0, d2 = range_value)
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")

model_sf$re_u <- 1
model_sf$re_v <- 1
formula <- Ho ~ cur_lgm +
  f(re_u, model = "besag", graph = g, scale.model = TRUE) +
  f(re_v, model = "iid")

formula <- Ho ~ cur_lgm + f(re_u, model = "bym2", graph = g)
res <- inla(formula, family = "gaussian", data = model_sf,
control.predictor = list(compute = TRUE),
control.compute = list(return.marginals.predictor = TRUE))

res$summary.fixed

library(spBayes)
library(sf)
library(sp)
library(spdep)
library(readr)
library(here)

# Load dataset
model_df <- read_csv(here("analysis/genetic_diversity/outputs/model_df.csv"))

# Convert to sf object and reproject
model_sf <- st_as_sf(model_df, coords = c("x", "y"), crs = 4326) %>% 
  st_transform(3310)

# Extract coordinates
coords <- st_coordinates(model_sf)

# Define spatial neighbors using semivariogram range
nb <- dnearneigh(coords, d1 = 0, d2 = range_value)

# Convert neighbors to spatial weights matrix
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Define priors
priors <- list(
  beta = list(rep(0, 2), diag(1000, 2)),  # Prior for regression coefficients
  phi.Unif = c(0.01, 1),                  # Uniform prior for spatial decay parameter
  sigma.sq.IG = c(2, 1)                   # Inverse-Gamma prior for variance
)

# Define starting values
starting <- list(
  phi = 0.1,       # Initial value for spatial decay
  sigma.sq = 1     # Initial value for variance
)

# Fit Bayesian Spatial Model
bayes_error_model <- spLM(Ho ~ lgm_lig + bio1 + tmean_dif, 
                          data = as.data.frame(model_sf), 
                          coords = coords, 
                          starting = starting, 
                          priors = priors, 
                          cov.model = "exponential",  # Use exponential spatial covariance
                          tuning = list(phi = 0.05, sigma.sq = 0.05), 
                          n.samples = 5000)

# Print summary
summary(bayes_error_model)

# Recover MCMC samples for beta (fixed effects) and spatial random effects
bayes_recovery <- spRecover(bayes_error_model, start = 1000, thin = 5, verbose = FALSE)

# Extract posterior samples for regression coefficients
beta_samples <- bayes_recovery$p.beta.recover.samples

# Summarize posterior distribution of beta coefficients
beta_summary <- apply(beta_samples, 2, function(x) {
  c(mean = mean(x), 
    sd = sd(x), 
    l95 = quantile(x, 0.025),  # 2.5% quantile (lower bound of 95% CI)
    u95 = quantile(x, 0.975),  # 97.5% quantile (upper bound of 95% CI)
    l90 = quantile(x, 0.05),  # 5% quantile (lower bound of 90% CI)
    u90 = quantile(x, 0.95))  # 95% quantile (upper bound of 90% CI)
})

# Convert to a dataframe for better readability
beta_summary_df <- as.data.frame(t(beta_summary))

# Print results
print(beta_summary_df)

hist(beta_samples[, "lgm_lig"], main = "Posterior Distribution of cur_lgm", xlab = "cur_lgm", col = "lightblue", breaks = 30)
abline(v = beta_summary_df["lgm_lig", "l95"], col = "red", lwd = 2, lty = 2)  # Lower CI
abline(v = beta_summary_df["lgm_lig", "u95"], col = "red", lwd = 2, lty = 2)  # Upper CI
abline(v = beta_summary_df["lgm_lig", "mean"], col = "blue", lwd = 2)         # Posterior mean
