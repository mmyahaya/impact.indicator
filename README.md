Code to calculate the impact indicator
================
Mukhtar Yahaya, Sabrina Kumschick, Sandra MacFadyen, Pietro Landi, Cang
Hui
2024-11-10

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

## Description

This R Markdown demonstrates the computation of an impact indicator for
biological invasions using the `impact_indicator()`. The
`impact_indicator()` feeds in species occurrence cube from the
`b3gbi::process_cube()` using `taxaFun()` and processed Environmental
Impact Classification of Alien Taxa (EICAT) impact score of species
using `impact_cat()`.

``` r
# Load packages
library(b3gbi)       # Biodiversity indicators for data cubes
library(tidyverse)   # Data wrangling and visualisation
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
#> ✔ purrr     1.0.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(sf)          # Spatial features
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.3.1; sf_use_s2() is TRUE


# Source functions
source("R/impact_indicator.R")
source("R/taxa_cube.R")
source("R/impact_cat.R")
```

``` r
# load the shapefile for study region
countries_sf<-readRDS("Data/countries_shapefile.rds")
SA.sf<-filter(countries_sf,name=="South Africa") %>% select(name,geometry)
plot(SA.sf, main = "South African map")
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

# Process occurrence cube

``` r
# load the GBIF occurrence data for taxa
taxa_Acacia<-readRDS("Data/taxa_Acacia.rds")
acacia_cube<-taxaFun(taxa=taxa_Acacia,
                    region=SA.sf,
                    res=0.25,
                    first_year=2010)
#> Warning: attribute variables are assumed to be spatially constant throughout
#> all geometries

acacia_cube$cube
#> 
#> Simulated data cube for calculating biodiversity indicators
#> 
#> Date Range: 2010 - 2024 
#> Number of cells: 366 
#> Grid reference system: custom 
#> Coordinate range:
#> [1] "Coordinates not provided"
#> 
#> Total number of observations: 5323 
#> Number of species represented: 25 
#> Number of families represented: Data not present 
#> 
#> Kingdoms represented: Data not present 
#> 
#> First 10 rows of data (use n = to show more):
#> 
#> # A tibble: 5,323 × 6
#>    scientificName   taxonKey minCoordinateUncertaintyInMe…¹  year cellCode   obs
#>    <chr>               <dbl>                          <dbl> <dbl>    <int> <dbl>
#>  1 Acacia implexa    2979232                              1  2010      206     1
#>  2 Acacia cyclops    2980425                            122  2010      668     1
#>  3 Acacia saligna    2978552                              1  2010      206     1
#>  4 Acacia pycnantha  2978604                              1  2010      206     1
#>  5 Acacia mearnsii   2979775                            110  2010      215     1
#>  6 Acacia mearnsii   2979775                              1  2010      215     1
#>  7 Acacia mearnsii   2979775                              8  2010     1376     1
#>  8 Acacia saligna    2978552                              1  2011      206     1
#>  9 Acacia saligna    2978552                             15  2011     1312     1
#> 10 Acacia mearnsii   2979775                              1  2011      230     1
#> # ℹ 5,313 more rows
#> # ℹ abbreviated name: ¹​minCoordinateUncertaintyInMeters
head(acacia_cube$coords)
#>   siteID        X       Y
#> 1      1 16.60833 -34.697
#> 2      2 16.85833 -34.697
#> 3      3 17.10833 -34.697
#> 4      4 17.35833 -34.697
#> 5      5 17.60833 -34.697
#> 6      6 17.85833 -34.697
```

# Aggregate impact scores for each species

There are often several impact records per species in different
mechanisms and regions. There is need to aggregate these impact records
for each species. The `impact_cat()` aggregates impact using ***max***,
***mean*** and ***max_mech*** as metrics as used by different studies.

- ***max***: maximum impact score across all records for the species  
- ***mean***: mean impact score across all records  
- ***max_mech***: sum of the maximum impact per mechanisms

``` r
eicat_data<-readRDS("Data/eicat_data.rds")
full_species_list<-sort(unique(acacia_cube$cube$data$scientificName))
#display part of data
eicat_data %>% select(scientific_name,impact_region,impact_mechanism,
                      impact_category) %>% 
  head(10)
#> # A tibble: 10 × 4
#>    scientific_name     impact_region            impact_mechanism impact_category
#>    <chr>               <chr>                    <chr>            <chr>          
#>  1 Dryophytes cinereus <NA>                     <NA>             DD - Data defi…
#>  2 Acacia cyclops      Cape Province            (9) Chemical Im… DD - Data defi…
#>  3 Acacia cyclops      Melkbosstrand            (9) Chemical Im… DD - Data defi…
#>  4 Acacia cyclops      Agulhas Plain            (9) Chemical Im… MR - Major     
#>  5 Acacia cyclops      South-western part of t… (11) Structural… MO - Moderate  
#>  6 Acacia cyclops      Cape of Good Hope Natur… (11) Structural… MO - Moderate  
#>  7 Acacia cyclops      Millers Point            (9) Chemical Im… MC - Minimal c…
#>  8 Acacia cyclops      Cape Province            (10) Physical I… MN - Minor     
#>  9 Acacia cyclops      Western Cape             (11) Structural… MO - Moderate  
#> 10 Acacia dealbata     Ourense (NW of Spain)    (9) Chemical Im… MR - Major
agg_impact<-impact_cat(impact_data=eicat_data,
                     species_list=full_species_list,
                     col_category="impact_category",
                     col_species="scientific_name",
                     col_mechanism="impact_mechanism",
                     trans=1)
agg_impact
#>                       max     mean max_mech
#> Acacia acinacea        NA       NA       NA
#> Acacia adunca          NA       NA       NA
#> Acacia baileyana       NA       NA       NA
#> Acacia binervata       NA       NA       NA
#> Acacia crassiuscula    NA       NA       NA
#> Acacia cultriformis    NA       NA       NA
#> Acacia cyclops          3 1.500000        6
#> Acacia dealbata         3 1.812500       18
#> Acacia decurrens        3 1.500000        6
#> Acacia elata           NA       NA       NA
#> Acacia falciformis     NA       NA       NA
#> Acacia implexa         NA       NA       NA
#> Acacia longifolia       3 1.562500       20
#> Acacia mearnsii         3 1.769231       15
#> Acacia melanoxylon     NA       NA       NA
#> Acacia paradoxa        NA       NA       NA
#> Acacia piligera        NA       NA       NA
#> Acacia podalyriifolia  NA       NA       NA
#> Acacia provincialis    NA       NA       NA
#> Acacia pycnantha        3 3.000000        3
#> Acacia saligna         NA       NA       NA
#> Acacia schinoides      NA       NA       NA
#> Acacia stricta         NA       NA       NA
#> Acacia ulicifolia      NA       NA       NA
#> Acacia viscidula       NA       NA       NA
```

# Compute impact risk map

The impact risk map shows the impact score for each site, where multiple
species can be present. To compute the risk per site, aggregated scores
across species at each site are needed. The `impact_indicator()` uses
***max***, ***sum*** and ***mean*** metrics to aggregate impact scores
as proposed by Boulesnane-Guengant et al., (in preparation). The
combinations of aggregation metrics for each species and site leads to
five type of indicators, namely, ***precautionary***, ***precautionary
cumulative***, ***mean***, ***mean cumulative*** and ***cumulative***.  

- ***precautionary***: maximum score across species’ max in each site  
- ***precautionary cumulative***: cumulative score across species’ max
  in each site  
- ***mean***: mean score across species’ mean in each site  
- ***mean cumulative***: cumulative score across species’ mean in each
  site  
- ***cumulative***: cumulative score across species’ sum of maximum
  score per mechanism  

``` r
impact_value<-impact_indicator(cube=acacia_cube$cube,
                               impact_data = eicat_data,
                               col_category="impact_category",
                               col_species="scientific_name",
                               col_mechanism="impact_mechanism",
                               trans=1,
                               type = "mean cumulative",
                               coords=acacia_cube$coords)
```

### impact risk map

``` r
#impact risk map
#visualise last four years for readability
impact_value$sitedf%>% 
  gather(-c(siteID,X,Y),key="year",value="impact") %>% 
  na.omit() %>% 
  filter(year>=2021) %>% 
  ggplot() +
  geom_tile(
    aes(x=X,y=Y,fill=impact),color="black")+
  geom_sf(data = SA.sf, fill = NA, color = "black", alpha = 0.5)+
  scale_fill_gradient(low = "yellow",
                       high = "red")+
  theme_minimal() +
   labs(
    title = "impact risk map (mean cumulative)",
    y = "Latitude", x="Longitude"
  )+
  theme(text=element_text(size=14))+
  facet_wrap(~year)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# Compute impact indicators

To compute the impact indicator of alien taxa, we sum all the yearly
impact scores of each site of the study region. To correct for sampling
effort we divide the yearly impact scores by number of sites in the
study region with at least a single occurrence throughout the whole
year.  
$$I_i = \frac{\sum{S_i}}{N}$$ $I_i$ is impact score at year $i$.  
$S_i$ is the sum of risk map value, where $S=\{s_1,s_2,...,s_n\}$ and
$s_n$ is the site score for site $n$  
$N$ is number of sites occupied through out the study years of the
region.  
**Note**: This is the only method incorporated as at now. Other methods
will be considered later.  
**Note**: A function `impact_uncertainty()` is being developed to use
bootstrap method to compute confidence interval of the indicator rather
than using `geom_smooth()` used below.

``` r
#sum of impact risk map for each year

ggplot(data = impact_value$impact_values) +
  geom_line(aes(y = value, x = year),colour="red",
            stat="identity",
            linewidth=2)+
  geom_smooth(aes(y = value, x = year),linetype=2)+
  labs(
    title = "Impact indicator",
    y = "impact score"
  )+
  theme_minimal() +
  theme(text=element_text(size=14))
#> `geom_smooth()` using method = 'loess' and formula = 'y ~ x'
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# Impact indicator per species

We compute the impact indicator per species by summing the impact risk
map per species and correct for sampling effort by dividing by $N$.

``` r
#  impact indicator per species
impact_value$species_values %>%
  rownames_to_column("year") %>%
  mutate(year=as.numeric(year)) %>%
  gather(-year,key = "Alien_species", value = "impact_score") %>%
  ggplot(aes(x = year, y = impact_score)) +
  geom_line(aes(color = Alien_species),linewidth=1.5)+
  theme_minimal() +
  labs(
    title = "sum of species impact",
    y = "impact score"
  )+
  theme(text=element_text(size=14))
#> Warning: Removed 5 rows containing missing values or values outside the scale range
#> (`geom_line()`).
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Comparing type of indicators

To compare type of impact indicators for a case study, we provide a plot
which can be adapted by a user to compare a set of method.

``` r
# plot all type of impact indicators
types<-c("precautionary",
         "precautionary cumulative",
         "mean",
         "mean cumulative",
         "cumulative")

all_impact<-data.frame("year"=unique(acacia_cube$cube$data$year))
for(type in types){
  impact_value<-impact_indicator(cube=acacia_cube$cube,
                                 impact_data = eicat_data,
                                 col_category="impact_category",
                                 col_species="scientific_name",
                                 col_mechanism="impact_mechanism",
                                 trans=1,
                                 type = type,
                                 coords=acacia_cube$coords)
  
  all_impact[type]<-impact_value$impact_values$value
}

all_impact %>%
  gather(-year,key = "indicator_type", value = "impact_score") %>% 
  ggplot(aes(x = year, y = impact_score)) +
  geom_line(aes(color = indicator_type),linewidth=1.5)+
  theme_minimal() +
   labs(
    title = "Type of indicators",
    y = "impact score"
  )+
  theme(text=element_text(size=14))
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
