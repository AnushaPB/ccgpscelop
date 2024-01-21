library(here)
library(tidyverse)
library(vcfR)
library(wingen)
library(terra)
library(viridis)
library(tmap)


get_gsd <- function(filepath){
  gsd <- read.csv(filepath)
  #assign IDs to rownames
  rownames(gsd) <- gsd$idx
  #extract env and z values
  gsd <-
    gsd %>%
    tidyr::extract(z, into = "z1", regex = "\\[(.*)\\]", remove = FALSE) %>%
    tidyr::extract(e, into = c("env0", "klyr", "envlyr"), regex = "\\[(.*), (.*), (.*)\\]", remove = FALSE) %>%
    dplyr::mutate(across(c(z1, env0, klyr, envlyr), as.numeric))
  
  if (all(is.na(gsd$z1))) gsd <- gsd %>% mutate(z1 = parse_number(z))
  
  #correct coordinates
  gsd$y <- -gsd$y

  return(gsd)
}

# DISTINCT SIM 1
filepath <- here("gnx", "GNX_mod-distinct_sim", "it-0", "spp-spp_0")
gsd <- get_gsd(here(filepath, "mod-distinct_sim_it-0_t-1000_spp-spp_0.csv"))
vcf <- read.vcfR(here(filepath, "mod-distinct_sim_it-0_t-1000_spp-spp_0.vcf"))
dos <- vcf_to_dosage(vcf)

df <- 
bind_cols(gsd, dos[,1:8]) %>%
pivot_longer(starts_with("0_"))

ggplot(gsd) +
  geom_point(aes(x = x, y = y, col = z1)) +
  coord_fixed()

ggplot(gsd) +
  geom_point(aes(x = x, y = y, col = factor(envlyr), pch = factor(klyr))) +
  coord_fixed()

ggplot(df) +
  geom_point(aes(x = x, y = y, col = value)) +
  facet_wrap(~name, nrow = 2) +
  coord_fixed()

lyr <- rast(matrix(rep(1,100*100),nrow = 100, ncol = 100))
ext(lyr) <- c(0, 100, -100, 0)
s <- sample(nrow(gsd), 1000)

#pstk <- window_p2(dos[s,1:2], coords = gsd[s,c("x", "y")], lyr = lyr, wdim = 3)
# Convert to SpatVector
points <- vect(bind_cols(gsd, dos[,1:8]), geom=c("x", "y"))
# Rasterize the points
rasterized <- rasterize(points, lyr, field=paste0("0_", 0:7), fun=mean)
freqs <- aggregate(rasterized, 2, FUN = mean, na.rm = TRUE)
dlstk <- divloss_p(freqs, prop = TRUE)

plot(freqs)
par(mfrow = c(2,2))
plot(mean(scale(dlstk), na.rm = TRUE), col = inferno(100, direction = -1))
plot(mean(dlstk, na.rm = TRUE), col = inferno(100, direction = -1))
plot(mean(freqs, na.rm = TRUE), col = inferno(100))

# DISTINCT SIM 2
filepath <- here("gnx", "GNX_mod-distinct_gradient_sim", "it-0", "spp-spp_0")
gsd <- get_gsd(here(filepath, "mod-distinct_gradient_sim_it-0_t-1000_spp-spp_0.csv"))
vcf <- read.vcfR(here(filepath, "mod-distinct_gradient_sim_it-0_t-1000_spp-spp_0.vcf"))
dos <- vcf_to_dosage(vcf)

df <- 
bind_cols(gsd, dos[,1:8]) %>%
pivot_longer(starts_with("0_"))

ggplot(gsd) +
  geom_point(aes(x = x, y = y, col = z1)) +
  coord_fixed()

ggplot(gsd) +
  geom_point(aes(x = x, y = y, col = envlyr, pch = factor(klyr))) +
  coord_fixed()

ggplot(df) +
  geom_point(aes(x = x, y = y, col = value)) +
  facet_wrap(~name, nrow = 2) +
  coord_fixed()

lyr <- rast(matrix(rep(1,100*100),nrow = 100, ncol = 100))
ext(lyr) <- c(0, 100, -100, 0)
s <- sample(nrow(gsd), 1000)

#pstk <- window_p2(dos[s,1:2], coords = gsd[s,c("x", "y")], lyr = lyr, wdim = 3)
# Convert to SpatVector
points <- vect(bind_cols(gsd, dos[,1:8]), geom=c("x", "y"))
# Rasterize the points
rasterized <- rasterize(points, lyr, field=paste0("0_", 0:7), fun=mean)
freqs <- aggregate(rasterized, 2, FUN = mean, na.rm = TRUE)
dlstk <- divloss_p(freqs, prop = TRUE)
plot(freqs)


a <- rast(matrix(c(1,1,1,1,0,0), ncol = 3, nrow = 2, byrow = TRUE))
b <- rast(matrix(c(0,0,1,0,0,1), ncol = 3, nrow = 2, byrow = TRUE))
d <- divloss_p(c(a,b), prop = TRUE)
par(mfrow = c(1,3))
plot(mean(a,b),  col = inferno(100, direction = -1))
plot(mean(d),  col = inferno(100, direction = 1))
plot((a*mean(1 - values(a)) + b*mean(1 - values(b)))/2, col = inferno(100, direction = -1))


ALL1 <- matrix(c(1, 1, 1, 1, 0, 0), nrow = 2, byrow = TRUE)
ALL2 <- matrix(c(0, 0, 1, 0, 0, 1), nrow = 2, byrow = TRUE)

# Function to calculate distinctness
calc_distinct <- function(x) {
  # Create weighted allele frequency map
  wx <- 
    map(terra::as.list(x), ~{
      p0 <- as.numeric(global(.x, fun = "sum", na.rm = TRUE))
      return(.x/p0)
    }) %>%
    rast()

  # calculate average of weighted allele frequency maps
  dmap <- mean(wx, na.rm = TRUE)
  
  return(dmap)
}

ALL1 <- rast(matrix(c(1,1,1,1,0,0), ncol = 3, nrow = 2, byrow = TRUE))
ALL2 <- rast(matrix(c(0,0,1,0,0,1), ncol = 3, nrow = 2, byrow = TRUE))
calc_distinct(c(ALL1, ALL2))


a <- mean(rast("outputs/pstk_pc1.tif"), na.rm = TRUE)
b <- rast("outputs/dlavg_pstk1.tif")
newb <- calc_distinct(a)
par(mfrow = c(1,2))
plot(b, col = inferno(100, direction = -1))
plot(-1*newb, col = inferno(100, direction = -1))
cor(values(b), values(-1*newb), use = "pairwise.complete.obs")

ALL1 <- rast(matrix(c(1,1,1,1,0,0), ncol = 3, nrow = 2, byrow = TRUE))
ALL2 <- rast(matrix(c(0,0,1,0,0,1), ncol = 3, nrow = 2, byrow = TRUE))
opt1 <- distinct_map(c(ALL1, ALL2), rare_allele = FALSE)
opt2 <- distinct_map(c(ALL1, ALL2), rare_allele = TRUE)
opt3 <- distinct_map(c(ALL1, ALL2, 1-ALL1, 1-ALL2), rare_allele = FALSE)

par(mfrow = c(1,3))
plot(opt1[[1]], col = inferno(100, direction = -1), main = "Original")
plot(opt2[[1]], col = inferno(100, direction = -1), main = "Rare Allele Transformation")
plot(opt3[[1]], col = inferno(100, direction = -1), main = "All Alleles")
plot()


a <- mean(rast("outputs/pstk_pc1.tif"), na.rm = TRUE)
opt1 <- distinct_map(a, rare_allele = FALSE)
opt2 <- distinct_map(a, rare_allele = TRUE)
opt3 <- distinct_map(c(a, 1-a), rare_allele = FALSE)

par(mfrow = c(1,3))
plot(opt1[[1]], col = inferno(100, direction = -1), main = "Original")
plot(opt2[[1]], col = inferno(100, direction = -1), main = "Rare Allele Transformation")
plot(opt3[[1]], col = inferno(100, direction = -1), main = "All Alleles")
plot()