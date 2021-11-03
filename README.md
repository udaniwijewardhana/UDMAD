# UDMAD
Shiny App for Seasonal Abundance Models using Areal Data

### Interim Biogeographic Regionalisation for Australia (IBRA) version 7 Regions shapefile 
This can be accessed from 
https://www.environment.gov.au/fed/catalog/search/resource/details.page?uuid=%7B4A2321F0-DD57-454E-BE34-6FD4BDE64703%7D 

Interim Biogeographic Regionalisation for Australia (IBRA) version 7.0 is a landscape-based technique to classifying land surface in Australia. There are 89 biogeographic regions, each reflecting a unified set of key environmental forces that impact the occurrence of flora and fauna and their interaction with the physical environment across Australia and its foreign territories (excluding Antarctica). IBRA bioregions are a larger scale regional classification of homogeneous ecosystems, whilst subregions are more localised. We have excluded the four marine bioregions Coral Sea, Indian Tropical Islands, Pacific Subtropical Islands and Subantarctic Islands which have both terrestrial and marine areas. To use only 85 bioregions, add the following section to the "generate the adjacency matrix" section soon after reading the shapefile in app.R.

```r
map <- subset(map, !REG_NAME_7 == "Coral Sea")
map <- subset(map, !REG_NAME_7 == "Pacific Subtropical Islands")
map <- subset(map, !REG_NAME_7 == "Subantarctic Islands")
map <- subset(map, !REG_NAME_7 == "Indian Tropical Islands")
```



