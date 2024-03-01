Genetic diversity
================

  - [Heterozygosity](#heterozygosity)
  - [ROH](#roh)
  - [Models](#models)

# Heterozygosity

``` r
# get het data
het <- get_het()

# Combine coords and het data
coords <- left_join(coords, het, by = "SampleID") 

ggplot() + 
  geom_sf(data = ca) +
  geom_sf(data = coords, aes(col = Ho), cex = 3) +
  scale_color_viridis_c(option = "mako") +
  theme_void()
```

![](genetic_diversity_files/figure-gfm/heterozygosity-1.png)<!-- -->

``` r
coords_proj <- st_transform(coords, 3310)
ca_proj <- st_transform(ca, 3310)
lyr <- coords_to_raster(coords_proj, buffer = 1, res = 10000)
window_Ho <- window_general(coords_proj$Ho, coords = coords_proj, lyr = lyr, wdim = 9, stat = mean, na.rm = TRUE)
ggplot_gd(window_Ho, bkg = ca_proj, col = magma(100))
krig_Ho <- krig_gd(window_Ho, grd = lyr, disagg_grd = 2)
mask_Ho <- mask_gd(krig_Ho, ca_proj, minval = 2)
ggplot_gd(mask_Ho, bkg = ca_proj, col = magma(100))
```

# ROH

``` r
roh_df <- get_roh()
dem <- get_dem()

ggplot(roh_df) +
 geom_sf(data = ca) +
 geom_sf(aes(geometry = geometry, fill = froh), col = "black", pch = 21, cex = 2.5) +
 scale_fill_viridis_c(option = "inferno", labels = scales::label_number(scale = 1, accuracy = 0.001)) +
 theme_void()
```

![](genetic_diversity_files/figure-gfm/ROH-1.png)<!-- -->

``` r
ggplot(roh_df) +
 geom_sf(data = ca) +
 geom_sf(aes(geometry = geometry, fill = froh), col = "black", pch = 21, cex = 2.5) +
 labs(fill = "froh\n(log color scale)") +
 scale_fill_viridis_c(option = "inferno",  trans = "log", labels = scales::label_number(scale = 1, accuracy = 0.001)) +
 theme_void()
```

![](genetic_diversity_files/figure-gfm/ROH-2.png)<!-- -->

``` r
ggplot(roh_df) +
 geom_sf(data = ca) +
 geom_sf(aes(geometry = geometry, col = pop), cex = 2.5) +
 theme_void()
```

![](genetic_diversity_files/figure-gfm/ROH-3.png)<!-- -->

``` r
ggplot(roh_df) +
 geom_sf(data = ca) +
 geom_raster(data = dem, aes(x = x, y = y, fill = elev)) + 
 geom_sf(aes(geometry = geometry, col = froh), pch = 16, cex = 2.5) +
 scale_color_viridis_c(option = "inferno", trans = "log", labels = scales::label_number(scale = 1, accuracy = 0.001)) +
 scale_fill_gradient(low = "black", high = "white", na.value = NA) +
 theme_void()
```

![](genetic_diversity_files/figure-gfm/ROH-4.png)<!-- -->

# Models

``` r
envdata <- get_env()
mod_df <- 
  left_join(roh_df, envdata) %>%
  st_drop_geometry() 

# Climate model:
mod_clim <- lm(froh0 ~ CA_rPCA1 + CA_rPCA2 + CA_rPCA3, data = mod_df)
summary(mod_clim)
```

    ## 
    ## Call:
    ## lm(formula = froh0 ~ CA_rPCA1 + CA_rPCA2 + CA_rPCA3, data = mod_df)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.024613 -0.006506 -0.003416 -0.000737  0.091316 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.0124348  0.0020489   6.069 2.21e-08 ***
    ## CA_rPCA1    -0.0016348  0.0009055  -1.805 0.073972 .  
    ## CA_rPCA2    -0.0028828  0.0007569  -3.808 0.000239 ***
    ## CA_rPCA3    -0.0018577  0.0010765  -1.726 0.087439 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01665 on 102 degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## Multiple R-squared:  0.1553, Adjusted R-squared:  0.1304 
    ## F-statistic:  6.25 on 3 and 102 DF,  p-value: 0.0006144

``` r
mod_clim2 <- lm(froh0 ~ CA_rPCA2 , data = mod_df)
summary(mod_clim2)
```

    ## 
    ## Call:
    ## lm(formula = froh0 ~ CA_rPCA2, data = mod_df)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.021028 -0.009425 -0.003462  0.000690  0.090924 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.0137804  0.0019782   6.966 3.05e-10 ***
    ## CA_rPCA2    -0.0024537  0.0006849  -3.583 0.000519 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01692 on 104 degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## Multiple R-squared:  0.1099, Adjusted R-squared:  0.1013 
    ## F-statistic: 12.84 on 1 and 104 DF,  p-value: 0.0005192

``` r
# Elevation model:
mod_elev <- lm(froh0 ~ elevation, data = mod_df)
summary(mod_elev)
```

    ## 
    ## Call:
    ## lm(formula = froh0 ~ elevation, data = mod_df)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.019146 -0.007176 -0.003883 -0.000020  0.096116 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 4.328e-03  2.211e-03   1.957 0.052848 .  
    ## elevation   8.389e-06  2.358e-06   3.558 0.000555 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01674 on 109 degrees of freedom
    ## Multiple R-squared:  0.104,  Adjusted R-squared:  0.09582 
    ## F-statistic: 12.66 on 1 and 109 DF,  p-value: 0.0005553

``` r
# Elevation and CA_rPCA2 are correlated
mat <- as.matrix(drop_na(mod_df[,c("CA_rPCA1", "CA_rPCA2", "CA_rPCA3", "elevation")]))
mod_clim_elev <- lm(froh0 ~ CA_rPCA1 + elevation + CA_rPCA3, data = mod_df)
summary(mod_clim_elev)
```

    ## 
    ## Call:
    ## lm(formula = froh0 ~ CA_rPCA1 + elevation + CA_rPCA3, data = mod_df)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.019777 -0.006584 -0.003581 -0.000530  0.096657 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)  3.143e-03  2.651e-03   1.186  0.23857   
    ## CA_rPCA1    -3.795e-04  1.042e-03  -0.364  0.71651   
    ## elevation    8.782e-06  2.887e-06   3.042  0.00299 **
    ## CA_rPCA3    -7.749e-04  1.017e-03  -0.762  0.44801   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01703 on 102 degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## Multiple R-squared:  0.1154, Adjusted R-squared:  0.08941 
    ## F-statistic: 4.436 on 3 and 102 DF,  p-value: 0.005664

``` r
# Development and imperviousness model
mod_imperv <- lm(froh0 ~ nlcd_imperviousness, data = mod_df)
summary(mod_imperv)
```

    ## 
    ## Call:
    ## lm(formula = froh0 ~ nlcd_imperviousness, data = mod_df)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.010941 -0.009505 -0.006339  0.000011  0.101348 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          1.094e-02  1.818e-03   6.017 2.42e-08 ***
    ## nlcd_imperviousness -1.375e-04  8.946e-05  -1.537    0.127    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01749 on 109 degrees of freedom
    ## Multiple R-squared:  0.02121,    Adjusted R-squared:  0.01223 
    ## F-statistic: 2.362 on 1 and 109 DF,  p-value: 0.1272

``` r
mod_devel <- lm(froh0 ~ LF_evh_developed_bin, data = mod_df)
summary(mod_devel)
```

    ## 
    ## Call:
    ## lm(formula = froh0 ~ LF_evh_developed_bin, data = mod_df)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.010497 -0.009021 -0.006613 -0.000106  0.104204 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               0.010497   0.001985   5.287 6.45e-07 ***
    ## LF_evh_developed_binTRUE -0.002413   0.003698  -0.652    0.515    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01765 on 109 degrees of freedom
    ## Multiple R-squared:  0.00389,    Adjusted R-squared:  -0.005248 
    ## F-statistic: 0.4257 on 1 and 109 DF,  p-value: 0.5155

``` r
# Comparison of AIC suggests that elevation alone is the best model
AIC(mod_clim_elev)
```

    ## [1] -556.6455

``` r
AIC(mod_clim)
```

    ## [1] -561.533

``` r
AIC(mod_elev)
```

    ## [1] -589.0453

``` r
AIC(mod_clim2)
```

    ## [1] -559.9809

``` r
# plot of elevation vs roh
# doesn't seem suggestive of a strong relationship
ggplot(mod_df) +
  geom_point(aes(x = elevation, y = froh0, col = pop)) +
  theme_classic()
```

![](genetic_diversity_files/figure-gfm/models-1.png)<!-- -->

``` r
ggplot(mod_df) +
  geom_point(aes(x = CA_rPCA2, y = froh0, col = pop)) +
  theme_classic()
```

![](genetic_diversity_files/figure-gfm/models-2.png)<!-- -->
