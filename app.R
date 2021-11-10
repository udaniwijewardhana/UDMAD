### Load Packages ###
library(shiny); library(rgdal); library(mapview); library(spdep) 
library(sf); library(INLA); library(DT); library(plyr) 
library(forge); library(dplyr); library(tidyverse)
library(shinythemes); library(ggplot2); library(leaflet)

################################################################################################################
# Shiny App for Seasonal Abundance Models using Areal Data
################################################################################################################

### Shiny User Interface ###

# By default the file size limit is 5MB. Here limit is 70MB.
options(shiny.maxRequestSize = 70*1024^2)

# Increase memory limit
memory.size(max = FALSE)

ui <- fluidPage(
  
# define theme 
theme = shinytheme("lumen"),

# use custom css  
tags$head(tags$link(href = "style.css", rel = "stylesheet")),
  
titlePanel(strong("UDMAD - A shiny App for Spatio-temporal modelling of areal data", 
           titleWidth = 450, windowTitle = "UDMAD")),

tags$div(h5(style="text-align: justify;", "First the user needs to upload the data csv file into the application and 
         then select whether normalize the numerical predictors or not.
         The data file should include only (with or without numeric/factor predictor variables):",
strong("Trend - Detected Trend, Count - Species count, Area - Area Name,
         areaID - Area ID, Year - Detected Year, Season - Detected Season,
         seasonID - Season ID. 
         These names are case sensitive."),
         "Data should be ordered according to factor levels as in sample 'Sample Data.csv'.
         Factor variables should add as factor/character without numeric/integer. 
         This can also be used for monthly or daily data.
         Simply replace seasonID and Season with the appropriate values.
         A sample format of the data can be found in",
         em(tags$u("https://github.com/uwijewardhana/UDMAD")),".")),

sidebarPanel(div(style='height: 890px; overflow: scroll',
                  
    # Loading the data file
    fileInput("file", "Choose data CSV File", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
    selectInput("prednorm", "Numeric predictors normalization:", choices=c("rnorm", "stand", "none"), selected = "none"),
    hr(),   
    
    # Input the 4 files to generate the shapefile (.dbf, .prj, .shp and .shx)
    fileInput("filemap", "Upload all map files: shp, dbf, shx and prj", accept=c('.shp','.dbf','.sbn','.sbx','.shx','.prj'), multiple=TRUE),
    helpText("Area name is a unique identifier of the area which should be the same as the area names in the map."),
    hr(),     
    
    selectInput("distribution", "Distribution:", choices=c("Poisson", "Negative Binomial",
                                                           "Zeroinflated Poisson", "Zeroinflated Negative Binomial",
                                                           "Poisson Hurdle", "Negative Binomial Hurdle"), selected = "Negative Binomial"),
    selectInput("tempeffect", "temporal random effect model:", choices=c("'ar1'", "'iid'", "'rw1'", "'rw2'"), selected = "'ar1'"),
    selectInput("speffect", "spatial random effect model:", choices=c("'bym'", "'bym2'", "'besag'", "'besagproper'"), selected = "'besagproper'"),
    selectInput("factor", "Include factor variables in the model:", choices=c("No", "Yes"), selected = "No"),
    h5('Generate Interaction Variables Here (if applicable)'),
    uiOutput("independent"),
    uiOutput("makeInteract1"), uiOutput("makeInteract2"),
    uiOutput("uiAdded"), actionButton("actionBtnAdd", "Create Interaction Term"),
    hr(),
    actionButton("summary", "Summary"))),

mainPanel(
    tabsetPanel(
    tabPanel("Data",
             column(12, DT::dataTableOutput("contents"))),
    tabPanel("Results",
             fluidRow(h4("Shapefile and adjacency matrix:"), align = "left"),
             fluidRow(column(6, div(style='height: 400px', mapviewOutput("mapview"))),
                      column(6, div(style='height: 400px', plotOutput("adj")))),
             fluidRow(h4("INLA model results:"), align = "left"),
             fluidRow(h5("Fixed effects:"), align = "left"),
             fluidRow(column(12, verbatimTextOutput("summary"))),
             fluidRow(column(12, htmlOutput("val"))))))
)

### Shiny Server ###

server <- function(input, output, session){
  
# Generate shapefile map
observe({
    shpdf <- input$filemap
    if(is.null(shpdf)){return(NULL)}
    
    previouswd <- getwd()
    uploaddirectory <- dirname(shpdf$datapath[1])
    setwd(uploaddirectory)
    for(i in 1:nrow(shpdf)){
      file.rename(shpdf$datapath[i], shpdf$name[i])
    }
    setwd(previouswd)
    
    map <- readOGR(paste(uploaddirectory, shpdf$name[grep(pattern="*.shp$", shpdf$name)], sep="/")) # Delete_null_obj=TRUE)
    map <- spTransform(map, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    
    output$mapview <- renderLeaflet({mapview(map)@map})
})
  
# Generate the adjacency matrix 
p <- reactive({
    shpdf <- input$filemap
    if(is.null(shpdf)){return(NULL)}
    
    previouswd <- getwd()
    uploaddirectory <- dirname(shpdf$datapath[1])
    setwd(uploaddirectory)
    setwd(previouswd)
    
    map <- st_read(paste(uploaddirectory, shpdf$name[grep(pattern="*.shp$", shpdf$name)], sep="/")) # Delete_null_obj=TRUE)
    p <- poly2nb(map) # construct the neighbour list
    return(p)
})

# Plot the adjacency matrix  
output$adj <- renderPlot({
  if(!is.null(input$filemap)){
  td <- tempdir()
  nb2INLA(paste(td, "map.adj", sep = "/"), p())  # create adjacency matrix
  inla.setOption(scale.model.default = F)
  g = inla.read.graph(filename = paste(td, "map.adj", sep="/")) # convert the adjacency matrix into a file in the INLA format
  return(image(inla.graph2matrix(g), xlab = "", ylab = ""))}
  #return(plot(g))}
  })

# Read Data CSV file
filedata <- reactive({
  inFile <- input$file
  if (is.null(inFile)){return(NULL)}
  
  x <- as.data.frame(read.csv(inFile$datapath, fileEncoding="UTF-8-BOM"))
  x$Count <- as.character(x$Count)
  x$Count <- as.numeric(x$Count)
  
  fac = x %>% select_if(negate(is.numeric))
  if(ncol(fac)>1){
    fac = fac[ , !(names(fac) %in% c("Area"))]    
    for(i in 1:ncol(fac)){fac[,i] = as.factor(fac[,i])}
    z = dplyr::select_if(x, is.numeric)
    R <- cbind(z,fac, Area = x$Area)
  }else {
    R = x
  }
  
  num = dplyr::select_if(x, is.numeric)
  if(ncol(num)>5){
    q = subset(num, select = -c(Count, Trend))
    q <- unique(q)
    p = subset(q, select = -c(Year, seasonID, areaID))
  }else {
    p = NULL
  }
  
  if(!is.null(p)){
    for(i in 1:ncol(p)){
      if(input$prednorm == "rnorm"){p[,i] <- round(rnorm(p[,i]), digits = 6)
      } else if(input$prednorm == "stand"){p[,i] <- round(scale(p[,i]), digits = 6)
      } else {p[,i] <- round(p[,i], digits = 6)}}
  } 
  
  p = cbind(p, Year = q$Year, seasonID = q$seasonID, areaID = q$areaID)
  p$ID = paste(p$Year, p$seasonID, p$areaID, sep = "-", collapse = NULL)
  p = subset(p, select = -c(Year, seasonID, areaID))
  
  y = dplyr::select_if(R, is.factor)
  Q = cbind(y, subset(num, select = c(Year, seasonID, areaID, Count, Trend)))
  Q$ID =  paste(Q$Year, Q$seasonID, Q$areaID, sep = "-", collapse = NULL)
  
  Q <- left_join(Q,p, by  = "ID")
  Q = subset(Q, select = -c(ID))
  return(Q)
})

# Output of the data table
output$contents <- DT::renderDataTable({
  req(input$file)
  df <- filedata()
  return(DT::datatable(df, options = list(scrollX = TRUE)))
})

# distribution
distribution <- reactive({
  if(input$distribution == "Poisson"){distribution = "poisson"
  } else if(input$distribution == "Negative Binomial"){distribution = "nbinomial"
  } else if(input$distribution == "Zeroinflated Poisson") {distribution = "zeroinflatedpoisson1"
  } else if(input$distribution == "Zeroinflated Negative Binomial") {distribution = "zeroinflatednbinomial1"
  } else if(input$distribution == "Poisson Hurdle") {distribution = "zeroinflatedpoisson0"
  } else {distribution = "zeroinflatednbinomial0"}
  return(distribution)
})

# Rendering the list to the ui
output$uiAdded <- renderUI({checkboxGroupInput('added', 'List of combinations', choices = names(interacts))})

# The main named list that will be used in other tasks
interacts <- reactiveValues()
makeReactiveBinding("interacts")

observe({
  input$actionBtnAdd # Trigger Add actions
  isolate({
    a <- c(input$makeInteract1,input$makeInteract2)
    b <- a %>% paste(collapse = "*")
    if(b != ""){interacts[[b]] <- a}
  })
})

# Checkbox list of all numeric variables to use 
independent <- reactive({
  inFile <- input$file
  
  if(!is.null(inFile)){
    x <- as.data.frame(read.csv(inFile$datapath, fileEncoding="UTF-8-BOM"))
    df = x[ , !(names(x) %in% c("Count", "seasonID", "areaID", "Area"))]
    return(names(df))
  }
})

output$independent <- renderUI({checkboxGroupInput("independent", "Independent (Predictor) Variables:", independent())})

# Variables to Add to the List of Combinations  
makeInteract <- reactive({
  inFile <- input$file
  
  if(!is.null(inFile)){
    x <- as.data.frame(read.csv(inFile$datapath, fileEncoding="UTF-8-BOM"))
    df = x[ , !(names(x) %in% c("Count", "seasonID", "areaID", "Area"))]    
    return(names(df))
  }
})

output$makeInteract1 <- renderUI({selectInput("makeInteract1", "Variable1 For Interaction:", makeInteract())})
output$makeInteract2 <- renderUI({selectInput("makeInteract2", "Variable2 For Interaction:", makeInteract())})

# Create dataframe for regression with predictor variables (numeric/factor)
M <- reactive({
  df <- filedata()
  if(is.null(df)){return(NULL)}
  
  c = c(input$independent)
  if(!is.null(c)){
    n = length(c)
    if(n == 1){
      M = data.frame(df[ , (names(df) %in% c)])
      names(M)[1] <- c[1]
      if("Year" %in% colnames(M)){M = NULL}else {M = M}  
    }else {
      M = data.frame(df[ , (names(df) %in% c)])
      if("Year" %in% colnames(M)){M = subset(M, select = -c(Year))}else {M = M}
    }}
  return(M)
})

# Final dataframe to fit the model
Final <- reactive({
  df <- filedata()
  if(is.null(df)){return(NULL)}
  
  # Count for the specific distribution
  if(input$distribution == "Poisson Hurdle" | input$distribution == "Negative Binomial Hurdle"){
    df$Count[df$Count == 0] <- NA
  }else {
    df$Count = df$Count
  }
  
  # Dataframe for models include factor variables and models include only numerical variables
  if(input$factor == "No"){
    
    z = cbind(M(), areaID = df$areaID, Year = df$Year, 
        effect = df$Trend, Count = df$Count, seasonID = df$seasonID)
    z$ID = paste(z$Year, z$seasonID, z$areaID, sep = "-", collapse = NULL)
    zz = subset(z, select = -c(Count))
    Final <- unique(zz)
    
    p <- aggregate(Count ~ Year + seasonID + areaID, z, FUN = sum)
    p$ID = paste(p$Year, p$seasonID, p$areaID, sep = "-", collapse = NULL)
    p = subset(p, select = -c(Year, seasonID, areaID))
    Final <- left_join(Final, p, by = "ID")
    Final = subset(Final, select = -c(ID))
    return(Final)
    
  }else {
    
    w = dplyr::select_if(M(), is.factor)
    w = cbind(w, areaID = df$areaID, seasonID = df$seasonID, Year = df$Year)
    w = unite(w, ID, sep = "-", remove = FALSE, na.rm = FALSE)
    Final = unique(w)
    w$Count = df$Count
    w$effect = df$Trend
    
    v = w[ , c("ID","effect")]
    Final <- left_join(Final, v, by = "ID")
    ww <- aggregate(Count ~ ID, w, FUN = sum)
    Final <- left_join(Final, ww, by = "ID")
    
    c = dplyr::select_if(M(), is.numeric)
    c$ID = w$ID
    c = unique(c)
    
    Final = merge(x = Final, y = c, by = "ID", all.x = TRUE)
    Final = subset(Final, select = -c(ID))
    return(Final)
  }
  
})

# Fit SDM using R-INLA
fitsummary <- reactive({
  if(is.null(input$filemap) | is.null(Final())){return(NULL)}
  df <- Final()
  
  td <- tempdir()
  nb2INLA(paste(td, "map.adj", sep = "/"), p())  # create adjacency matrix
  g = inla.read.graph(filename = paste(td, "map.adj", sep="/")) # convert the adjacency matrix into a file in the INLA format
  
  # formula
  if(!is.null(input$added)){
    formula = paste("Count ~ 1 +", paste(input$independent, collapse = "+"), 
                    paste("+", paste(input$added, collapse = "+")), 
                    paste("+", "f(effect, model = ", input$tempeffect, ")"),
                    paste("+", "f(areaID, model = ", input$speffect, ", graph = g)"))
  }else {
    formula = paste("Count ~ 1 + ", paste(input$independent, collapse = "+"), 
                    paste("+", "f(effect, model = ", input$tempeffect, ")"),
                    paste("+", "f(areaID, model = ", input$speffect, ", graph = g)"))
  }
  
  # Fit INLA model
  if(input$factor == "Yes"){
    
    model <- inla(as.formula(formula), data = df, family = distribution(),
                  control.family = list(link = "log"),
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE), verbose = T)
  }else {
    
    model <- inla(as.formula(formula), data = df, family = distribution(),
                  control.family = list(link = "log"),
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), verbose = T)
  }
  return(model)
  
})

# Summary output of the INLA model
fitsum <- eventReactive(input$summary, {
  round(fitsummary()$summary.fixed[,c(1:3,5)], digits = 6)})
output$summary <- renderPrint({return(fitsum())})

output$val <- renderPrint({
  if(is.null(fitsummary())){return(NULL)}
                            
  str1 <- paste("DIC:", round(fitsummary()$dic$dic, digits = 2))
  str2 <- paste("WAIC:", round(fitsummary()$waic$waic, digits = 2))
  str3 <- paste("Marginal log likelihood:", round(fitsummary()$mlik[1], digits = 2))
  HTML(paste(str1, str2, str3, sep = '<br/>'))
})

}

shinyApp(ui, server)
