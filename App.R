# Preparation ----
library(shiny)
# library(rsconnect)
# UI
library(shinyWidgets)
library(shinyjs)
library(waiter) # wait while running
library(DT)
# library(shinydashboard)
library(config)

# Manipulate data
library(dplyr)
library(data.table)

# Meta analysis
library(meta)
library(metafor)
library(effsize)

# Plot
library(UpSetR)
library(ggpubr)
library(ggrepel)
library(treemapify)
library(gridExtra)
library(grid)
library(patchwork)
library(cowplot)
library(gg.gap)

# Multithreads
library(snowfall)
library(parallel)

# DROMA Packages
library(DROMA.Set)
library(DROMA.R)

## Debug
# library(reactlog)

# Load ----
config_list <- config::get(
  # config = "test"
  # Default is test mode
)

# Initialize database connection
db_path <- config_list$db_path
connectDROMADatabase(db_path)

# Welcome notification
str1 <- "Nice to meet you."
str2 <- "Very welcome to my refactored version (v1.0) —2025/11/17"
str3 <- "You can visit https://github.com/mugpeng/DROMA_DB to reach the toturial."
modal_notification <- modalDialog(
  # p("Nice to meet you. \n, test"),
  HTML(paste(str1, str2, str3, sep = '<br/>')),
  title = "Update Notification",
  footer = tagList(
    actionButton("close_modal", "Close")
  )
)

## Modules----
source("Modules/DrugOmicPair.R")
source("Modules/BatchFeature.R")
source("Modules/DrugFeature.R")
source("Modules/StatAnno.R")
source("Functions/DatabaseHelpers.R")
source("Functions/PlottingHelpers.R")
source("Modules/GlobalSetting.R")

# UI ----
ui <- tagList(
  tags$head(
    tags$title("DROMA(Drug Response Omics association MAp)"),
  ),
  useShinyjs(),  # Enable shinyjs
  # Global Settings Module
  uiGlobalSetting("GlobalSetting"),
  autoWaiter(html = spin_loader(), color = transparent(0.5)),
  navbarPage("DROMA(Drug Response Omics association MAp)",
             ## Drugs-omics pairs analysis ----
             tabPanel("Drugs-omics Pairs Analysis",
                      uiDrugOmicPair("DrugOmicPair")
             ),
             ## Features database significant analysis ----
             tabPanel("Batch Features Associations Analysis",
                      uiBatchFeature("BatchFeature")
             ),
             ## Drug Feature Analysis ----
             tabPanel("Drug Feature Analysis",
                      uiDrugFeature("DrugFeature")
             ),
             ## Statistics and Annotations ----
             tabPanel("Statistics and Annotations",
                      uiStatAnno("StatAnno")
             ),
             ## Contact ----
             tabPanel("Contact",
                      fluidPage(
                        strong("Feel free to talk with me if you find any bugs or have any suggestions. :)"),
                        p(""),
                        p("Email: mugpeng@outlook.com"),
                        p("Email: yc47680@um.edu.mo"),
                        p("github: https://github.com/mugpeng"),
                        p("You can visit https://github.com/mugpeng/DROMA_DB to reach the toturial.")
                      ))
  )
)

# Server ----
server <- function(input, output, session) {
  # Some setup ----
  showModal(modal_notification) # notification
  observeEvent(input$close_modal, {
    removeModal()
  })
  
  # Initialize Global Settings Module
  callModule(serverGlobalSetting, "GlobalSetting")
  
  # stop warn
  storeWarn <- getOption("warn")
  options(warn = -1) 
  # Drugs-omics pairs analysis ----
  callModule(serverDrugOmicPair, "DrugOmicPair")
  # Features database significant analysis ----
  callModule(serverBatchFeature, "BatchFeature")
  # Drug Feature Analysis ----
  callModule(serverDrugFeature, "DrugFeature")
  # Statistics and Annotations ----
  callModule(serverStatAnno, "StatAnno")
}

# Run ----
shinyApp(ui = ui, server = server)

