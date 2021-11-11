# UDMAD
Shiny App for Single Species Seasonal Spatio-temporal Models using Areal Data

This shiny application used areal data to fit spatio-temporal models of seasonal species counts. When a fixed domain is partitioned into a finite number of subregions from which results are aggregated, areal or lattice data are produced. The Besag-York-Mollié (BYM) model, besag, besagproper, and BYM2 (Besag, York, and Mollié, 1991) have taken into account that data may be spatially correlated and observations in neighbouring areas may be more comparable than observations in areas further away. This model incorporates a spatial random effect that smooths the data based on a neighbourhood structure, as well as an unstructured exchangeable component that mimics uncorrelated noise. Spatio-temporal models that account not only for spatial structure but also for temporal correlations and spatio-temporal interactions are utilised in spatio-temporal contexts where species counts are monitored throughout time. We determined the mean expected value of the counts for each fitted model (with 95% credible intervals). 

Integrated Nested Laplace Approximations (INLA) replaces long MCMC simulations with accurate, deterministic approximations to posterior marginal distributions for fitting these models, gaining speed. We have considered the Poisson, Negative Binomial, Zero-inflated Poisson, Zero-inflated Negative Binomial, Poisson Hurdle, or Negative Binomial Hurdle models to fit the INLA model. An interaction effect exists in regression when the influence of one or more independent variables on a dependent variable changes depending on the value(s) of one or more other independent variables. Users can add interaction terms to their regression models in this section.

### Layout
The application consists of two tabs. The sidebarPanel is shared by both tabs. The user must first upload the .CSV file and decide whether or not to normalise the numeric predictors. The user should next upload the shapefile (a combination of .shp,.dbf,.shx, and .prj files) to construct the adjacency matrix. The area name is a unique identifier for the area, and it should match the area names on the map. This shiny application used areal data to fit spatio-temporal models of seasonal species counts. This can be used to count data on a monthly or daily basis by substituting the values seasonID and Season below with monthly/daily values.

The data CSV should include: Trend - Detected Trend, Count - Species count, Area - Area Name, areaID - Area ID, Year - Detected Year, Season - Detected Season, seasonID - Season ID with or without numeric/factor predictor variables. These names are case sensitive. Data should be ordered according to factor levels as in sample “Sample Data.csv”. Factor variables should add as factor/character without numeric/integer. This can also be used for monthly or daily data. Simply replace seasonID and Season with the appropriate values. More details can be found in vignatte.

### Case study: 
#### Effect of biomes for Masked Lapwing in Australia

For this study we have excluded the four marine bioregions Coral Sea, Indian Tropical Islands, Pacific Subtropical Islands and Subantarctic Islands which have both terrestrial and marine areas. To use only 85 bioregions, add the following section to the "generate the adjacency matrix" section soon after reading the shapefile in app.R. The removed shapefile has included here. 

#### Interim Biogeographic Regionalisation for Australia (IBRA) version 7 Regions shapefile 
This can be accessed from 
https://www.environment.gov.au/fed/catalog/search/resource/details.page?uuid=%7B4A2321F0-DD57-454E-BE34-6FD4BDE64703%7D 

Interim Biogeographic Regionalisation for Australia (IBRA) version 7.0 is a landscape-based technique to classifying land surface in Australia. There are 89 biogeographic regions, each reflecting a unified set of key environmental forces that impact the occurrence of flora and fauna and their interaction with the physical environment across Australia and its foreign territories (excluding Antarctica). IBRA bioregions are a larger scale regional classification of homogeneous ecosystems, whilst subregions are more localised. 

```r
map <- subset(map, !REG_NAME_7 == "Coral Sea")
map <- subset(map, !REG_NAME_7 == "Pacific Subtropical Islands")
map <- subset(map, !REG_NAME_7 == "Subantarctic Islands")
map <- subset(map, !REG_NAME_7 == "Indian Tropical Islands")
```

```r
Udani Wijewardhana (udaniwijewardhana@gmail.com)
```

