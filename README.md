# UDMAD
Shiny App for Songle Species Seasonal Spatio-temporal Models using Areal Data

### Layout
The application consists of two tabs. The sidebarPanel is shared by both tabs. The user must first upload the .CSV file and decide whether or not to normalise the numeric predictors. The user should next upload the shapefile (a combination of .shp,.dbf,.shx, and .prj files) to construct the adjacency matrix. The area name is a unique identifier for the area, and it should match the area names on the map. This shiny application used areal data to fit spatio-temporal models of seasonal species counts. This can be used to count data on a monthly or daily basis by substituting the values seasonID and Season below with monthly/daily values.

The data CSV should include: Trend - Detected Trend, Count - Species count, Area - Area Name, areaID - Area ID, Year - Detected Year, Season - Detected Season, seasonID - Season ID with or without numeric/factor predictor variables. These names are case sensitive. Data should be ordered according to factor levels as in sample “Sample Data.csv”. Factor variables should add as factor/character without numeric/integer. This can also be used for monthly or daily data. Simply replace seasonID and Season with the appropriate values. More details can be found in vignatte.

### Case study: 
Effect of biomes for inland resident shorebirds in Australia. 

For this studt we have excluded the four marine bioregions Coral Sea, Indian Tropical Islands, Pacific Subtropical Islands and Subantarctic Islands which have both terrestrial and marine areas. To use only 85 bioregions, add the following section to the "generate the adjacency matrix" section soon after reading the shapefile in app.R.

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



