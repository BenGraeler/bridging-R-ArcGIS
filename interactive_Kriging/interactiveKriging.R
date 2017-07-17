################################################################################
### Interactive kriging bridge
################################################################################

tool_exec <- function(in_params, out_params){

  arc.progress_label("Loading packages ...")
  arc.progress_pos(10)
  
  library(sp)
  library(raster)
  library(shiny)
  library(rgdal)
  library(gstat)
  
################################################################################
### Define input and output parameters
################################################################################
  
  spdf_data <- in_params[[1]]
  sp_target <- in_params[[2]]
  var <- in_params[[3]]

  out_raster <- out_params[[1]]
  
  arc.progress_label("Loading data ...")
  arc.progress_pos(20)
  
################################################################################
### Load Data 
################################################################################
 
  spdf_esri <- arc.open(spdf_data)
  spdf_esri_sub <- arc.select(spdf_esri, fields = "*")
  
  target_esri <- arc.open(sp_target)
  target_esri_sub <- arc.select(target_esri, fields = "*")
  
  spdf <- arc.data2sp(spdf_esri_sub)
  target <- arc.data2sp(target_esri_sub)
  
################################################################################
### Pre-process Data 
################################################################################
  
  variance <- var(spdf[[var]], na.rm = T)
  diagonal <- spDists(coordinates(t(spdf@bbox)),
                      longlat = !is.projected(spdf))[1, 2]
  
  vgmFit <- vgm(0.9*variance, "Sph", diagonal/2, 0.1*variance)
  empVgm <- variogram(eval(parse(text=paste(var, "~ 1"))), spdf)
  oldABinput <- 0
  
################################################################################
### Load App
################################################################################
  
  arc.progress_label("Starting App ...")
  arc.progress_pos(40)
  
  vgmChoices <- as.list(as.character(vgm()[-1,1]))
  names(vgmChoices) <- as.character(vgm()[-1,2])
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(uiOutput("vgmTypeSel"),
                   uiOutput("psillSli"), 
                   uiOutput("rangeSli"),
                   uiOutput("nuggetSli"),
                   actionButton("autofit", "Autofit", value = FALSE),
                   tags$button(id = 'close',
                               type = "button",
                               class = "btn action-button",
                               onclick = "setTimeout(function(){window.close();},500);",
                               "Close window to proceed.")),
      mainPanel(plotOutput("vgmPlot"),
                textOutput("rmse"))
    )
  )
  
  server <- function(input, output) {
    observe({  
      vgmMod <- reactive({
        if (!any(is.null(c(input$psill, input$vgmType, input$range, input$nugget))))
          curVgm <- vgm(input$psill, input$vgmType, input$range, input$nugget)
        else 
          curVgm <- vgmFit
        
        if (input$autofit > oldABinput) {
          cat("automatic fit \n", input$autofit)
          curVgm <- fit.variogram(empVgm, curVgm)
          print(curVgm)
          oldABinput <<- input$autofit
        } 
        
        curVgm
      })
      
     
      vgmFit <<- vgmMod()
      
      
      output$vgmPlot <- renderPlot(plot(empVgm, vgmFit))
      
      output$rmse <- renderText(paste("Current RMSE:", 
                                      format(sqrt(mean((variogramLine(vgmFit, dist_vector = empVgm$dist) - empVgm$gamma)^2)), 
                                             digits = 6)))
      
      output$vgmTypeSel <- renderUI(selectInput("vgmType", "Select variogram family",
                                                choices =  vgmChoices,
                                                selected = as.character(vgmFit$model[2])))
      
      output$psillSli <- renderUI(sliderInput("psill", "Select the partial sill", 
                                              min = 0, max = as.numeric(format(ceiling(2*variance), digits=4)), 
                                              value = vgmFit$psill[2]))
      
      output$rangeSli <- renderUI(sliderInput("range", "Select the spatial range", 
                                              min = max(0.0001, floor(diagonal)/100), 
                                              max = as.numeric(format(ceiling(2*diagonal), digits = 4)), 
                                              value = vgmFit$range[2], 
                                              step = max(0.0001, floor(diagonal)/1000)))
      
      output$nuggetSli <- renderUI(sliderInput("nugget", "Select the nugget",
                                               min = 0, max = as.numeric(format(ceiling(2*variance), digits = 4)),
                                               value = vgmFit$psill[1]))
      
      if (input$close > 0) {
        stopApp()
      }
    })
  }

  runApp(list(ui = ui, server = server))
  
  arc.progress_label("Performing kriging ...")
  arc.progress_pos(60)
  
################################################################################
### Perform interpolation
################################################################################  
  
  res <- krige(eval(parse(text=paste(var, "~ 1"))), spdf, target, vgmFit)
  
  gridded(res) <- TRUE
  print(str(res))
  
################################################################################
### Write Output
################################################################################
  
  arc.progress_label("Writing Output...")
  arc.progress_pos(80)
  
  if(!is.null(out_raster) && out_raster != "NA")
    writeRaster(raster(res), out_raster, overwrite=TRUE)

  arc.progress_pos(100)
}