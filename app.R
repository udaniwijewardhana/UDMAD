## app.R ##
library(shiny)
library(rgdal); library(mapview); library(spdep); library(sf)
library(INLA); library(DT); 
library(plyr); library(dplyr); library(tidyverse)

################################################################################################################
# Shiny App for Seasonal Abundance Models using Areal Data
################################################################################################################

# First the user needs to upload the data csv file into the application and 
# then select whether normalize the numerical predictors or not.
# The data file should include only:
#         1. Trend - Detected Trend
#         2. Count - Species count
#         3. Area - Area Name 
#         4. areaID - Area ID
#         5. Year - Detected Year
#         6. Season - Detected Season
#         7. seasonID - Season ID
#         with or without predictor variables (numeric/factor).
# The above names are case sensitive.
# Factor variables should add as factor/character without numeric/integer.
# A sample format of the data can be found in https://github.com/uwijewardhana/UDMAD.
# Data should be ordered according to factor levels as in sample "Sample Data.csv".

### Shiny User Interface ###

# By default the file size limit is 5MB. Here limit is 70MB.
options(shiny.maxRequestSize = 70*1024^2)

# Increase memory limit
memory.size(max = FALSE)

ui <- fluidPage(
  
  titlePanel(strong("UDMAD - A shiny App for Spatio-temporal modelling of areal data", titleWidth = 450)),
  
  sidebarLayout(
    
    sidebarPanel(div(style='height:950px; overflow: scroll',
                     
                     # Loading the data file
                     fileInput("file", "Choose data CSV File", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                     selectInput("prednorm", "Numeric predictors normalization:", choices=c("rnorm", "stand", "none"), selected = "none"),
                     
                     # Input the 4 files to generate the shapefile (.dbf, .prj, .shp and .shx)
                     fileInput("filemap", "Upload all map files: shp, dbf, shx and prj", accept=c('.shp','.dbf','.sbn','.sbx','.shx','.prj'), multiple=TRUE),
                     
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
    
    mainPanel(fluidRow(column(6, mapviewOutput("mapview")),
                       column(6, plotOutput("g"))),
              fluidRow(column(6, div(style='height:590px; overflow: scroll', DT::dataTableOutput("contents"))),
                       column(6, verbatimTextOutput("summary"))))))

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
    
    output$mapview<-renderMapview({mapview(map)})
  })
  
  g <- reactive({
    shpdf <- input$filemap
    if(is.null(shpdf)){return(NULL)}
    
    previouswd <- getwd()
    uploaddirectory <- dirname(shpdf$datapath[1])
    setwd(uploaddirectory)
    setwd(previouswd)
    
    map <- st_read(paste(uploaddirectory, shpdf$name[grep(pattern="*.shp$", shpdf$name)], sep="/")) # Delete_null_obj=TRUE)
    sf::sf_use_s2(FALSE)
    p <- poly2nb(map) # Create adjacency matrix
    td <- tempdir()
    nb2INLA(paste(td, "map.adj", sep="/"), p)  # Convert to INLA 
    g <- inla.read.graph(filename = paste(td, "map.adj", sep="/"))
    return(g)
  })
  
  output$g <- renderPlot({
              if(is.null(input$filemap)){return(NULL)}
              plot(g())})
  
  # Read Data CSV file
  
  filedata1 <- reactive({
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
    x <- cbind(z,fac, Area = x$Area)
    return(x)
    }else {
    return(x) 
    }
    
  })
  
  # normalize numeric predictor variables (if applicable)
  
  filedata2 <- reactive({
    req(input$file)
    x <- filedata1()
    
    y = dplyr::select_if(x, is.numeric)
    if(ncol(y)>5){
      q = subset(y, select = -c(Count, Trend))
      q <- unique(q)
      p = subset(q, select = -c(Year, seasonID, areaID))
    }else {p = NULL}
    
    if(!is.null(p)){
      for(i in 1:ncol(p)){
        if(input$prednorm == "rnorm"){p[,i] <- round(rnorm(p[,i]), digits = 6)
        } else if(input$prednorm == "stand"){p[,i] <- round(scale(p[,i]), digits = 6)
        } else {p[,i] <- round(p[,i], digits = 6)}}
    } 
    
    p = cbind(p, Year = q$Year, seasonID = q$seasonID, areaID = q$areaID)
    p$ID = paste(p$Year, p$seasonID, p$areaID, sep = "-", collapse = NULL)
    p = subset(p, select = -c(Year, seasonID, areaID))
    
    x$ID =  paste(x$Year, x$seasonID, x$areaID, sep = "-", collapse = NULL)
    z = x %>% select_if(negate(is.numeric))
    j = dplyr::select_if(x, is.numeric)
    j = subset(j, select = c(Year, seasonID, areaID, Count, Trend))
    x = cbind(z,j)
    
    x <- left_join(x, p, by = "ID")
    x = subset(x, select = -c(ID))
    return(x)
  })
  
  # Output of the data table
  
  output$contents <- DT::renderDataTable({
    req(input$file)
    df <- filedata1()
    return(DT::datatable(df, options = list(scrollX = TRUE)))
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
      if(b != "")
        interacts[[b]] <- a
    })})
  
  # Checkbox list of all numeric variables to use 
  independent <- reactive({
    if(!is.null(input$file)){
      inFile <- input$file
      
      x <- as.data.frame(read.csv(inFile$datapath, fileEncoding="UTF-8-BOM"))
      df = x[ , !(names(x) %in% c("Count", "seasonID", "areaID", "Area"))]
      return(names(df))
    }
  })
  
  output$independent <- renderUI({checkboxGroupInput("independent", "Independent (Predictor) Variables:", independent())})
  
  # Variables to Add to the List of Combinations  
  makeInteract <- reactive({
    if(!is.null(input$file)){
      inFile <- input$file
      
      x <- as.data.frame(read.csv(inFile$datapath, fileEncoding="UTF-8-BOM"))
      df = x[ , !(names(x) %in% c("Count", "seasonID", "areaID", "Area"))]    
      return(names(df))
    }
  })
  
  output$makeInteract1 <- renderUI({selectInput("makeInteract1", "Variable1 For Interaction:", makeInteract())})
  output$makeInteract2 <- renderUI({selectInput("makeInteract2", "Variable2 For Interaction:", makeInteract())})
  
  # Create dataframe for regression with predictor variables (numeric/factor)
  Final <- reactive({
    req(input$file)
    x <- filedata1()
    if (is.null(x)){return(NULL)}
    df <- filedata2()
    
    if(input$distribution == "Poisson Hurdle" | input$distribution == "Negative Binomial Hurdle"){
      df$Count[df$Count == 0] <- NA
    }else {
      df$Count = df$Count
    }
    
    y = df[ , (names(df) %in% c(input$independent))]
    
    if(input$factor == "No"){
      
      z = cbind(y, Area = df$areaID, effect = df$Trend, Count = df$Count, seasonID = df$seasonID)
      
      if("Year" %in% colnames(y)){z = z
      }else {z$Year = df$Year}
      
      z$ID = paste(z$Year, z$seasonID, z$Area, sep = "-", collapse = NULL)
      zz = subset(z, select = -c(Count))
      Final <- unique(zz)
      
      p <- aggregate(Count ~ Year + seasonID + Area, z, FUN = sum)
      p$ID = paste(p$Year, p$seasonID, p$Area, sep = "-", collapse = NULL)
      p = subset(p, select = -c(Year, seasonID, Area))
      
      Final <- left_join(Final, p, by = "ID")
      Final = subset(Final, select = -c(ID))
      
    }else {
      
      w = dplyr::select_if(y, is.factor)
      w = cbind(w, Area = df$areaID, seasonID = df$seasonID, Year = df$Year)
      w$ID <- apply(w, 1, function(x) paste(x[!is.na(x)], collapse = "-"))
      Final = unique(w)
      w$Count = df$Count
      ww <- aggregate(Count ~ ID, w, FUN = sum)
      Final <- left_join(Final, ww, by = "ID")
      
      c = dplyr::select_if(y, is.numeric)
      if("Year" %in% colnames(c)){c = subset(c, select = -c(Year))
      }else {c = c}
      c$ID = w$ID; c$effect = df$Trend
      c = unique(c)
      
      Final = merge(x = Final, y = c, by = "ID", all.y = TRUE)
      Final = subset(Final, select = -c(ID))
    }
    return(Final)
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
  
  # formula
  formula <- reactive({
    if(!is.null(input$added)){
      formula = paste("Count ~ 1 +", paste(input$independent, collapse = "+"), 
                      paste("+", paste(input$added, collapse = "+")), 
                      paste("+", "f(effect, model = ", input$tempeffect, ")"),
                      paste("+", "f(Area, model = ", input$speffect, ", graph = ", g(), ")"))
    }else {
      formula = paste("Count ~ 1 + ", paste(input$independent, collapse = "+"), 
                      paste("+", "f(effect, model = ", input$tempeffect, ")"),
                      paste("+", "f(Area, model = ", input$speffect, ", graph = ", g(), ")"))
    }
    return(formula)
  })
  
  # Fit SDM using R-INLA
  fitsummary <- reactive({
    req(input$file)
    if(is.null(filedata2())){return(NULL)}
    df <- Final()
    
    if(input$factor == "Yes"){
      
      model <- inla(as.formula(formula()), data = df, family = distribution(),
                    control.family = list(link = "log"),
                    control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE), verbose = T)
    }else {
      
      model <- inla(as.formula(formula()), data = df, family = distribution(),
                    control.family = list(link = "log"),
                    control.compute = list(dic = TRUE, cpo = TRUE), verbose = T)
    }
    return(summary(model))
  })
  
  # Summary output of SDM
  fitsum <- eventReactive(input$summary, {fitsummary()})
  output$summary <- renderPrint({return(fitsum())})
  
}

shinyApp(ui, server)