uiStatAnno <- function(id){
  ns <- NS(id)
  fluidPage(
    column(12,
           navlistPanel(
             tabPanel("Overall Drug Information",
                      tabsetPanel(
                        tabPanel("Drug and Sample Counts",
                                 plotOutput(ns("p_count_drugandsample_sum_with_gap2")),
                                 plotOutput(ns("p_count_drugandsample_facet_with_gap"))
                        ),
                        tabPanel("Molecular Characteristics",
                                 plotOutput(ns("p_mol_character"))
                        ),
                        tabPanel("Drug and Sample Overlap Counts",
                                 plotOutput(ns("p_overlap_drug")),
                                 plotOutput(ns("p_overlap_sample"))
                        ),
                        tabPanel("Drug and Sample Annotation visulization",
                                 plotOutput(ns("p_tumor_bubble")),
                                 plotOutput(ns("p_drug_moa"))
                        )
                      )
             ),
             tabPanel("Annotation",
                      tabsetPanel(
                        tabPanel("Sample",
                                 DT::dataTableOutput(ns("sample_anno")),
                                 br(),
                                 wellPanel(
                                   fluidRow(
                                     column(6,
                                            prettyRadioButtons(
                                              inputId = ns("sample_download_type"),
                                              label = "Select download format:",
                                              choices = c(
                                                "Data (RDS)" = "data_rds",
                                                "Data (CSV)" = "data_csv"
                                              ),
                                              selected = "data_csv",
                                              icon = icon("check"), 
                                              animation = "jelly",
                                              status = "primary"
                                            )
                                     ),
                                     column(6,
                                            downloadBttn(
                                              outputId = ns("download_sample"),
                                              label = "Download",
                                              style = "gradient",
                                              color = "default",
                                              block = TRUE,
                                              size = "sm"
                                            )
                                     )
                                   )
                                 )
                        ),
                        tabPanel("Drug",
                                 DT::dataTableOutput(ns("drug_anno")),
                                 br(),
                                 wellPanel(
                                   fluidRow(
                                     column(6,
                                            prettyRadioButtons(
                                              inputId = ns("drug_download_type"),
                                              label = "Select download format:",
                                              choices = c(
                                                "Data (RDS)" = "data_rds",
                                                "Data (CSV)" = "data_csv"
                                              ),
                                              selected = "data_csv",
                                              icon = icon("check"), 
                                              animation = "jelly",
                                              status = "primary"
                                            )
                                     ),
                                     column(6,
                                            downloadBttn(
                                              outputId = ns("download_drug"),
                                              label = "Download",
                                              style = "gradient",
                                              color = "default",
                                              block = TRUE,
                                              size = "sm"
                                            )
                                     )
                                   )
                                 )
                        ),
                      )),
           )),
  )
}

serverStatAnno <- function(input, output, session){
  ns <- session$ns
  
  # Get available projects
  projects <- listDROMAProjects()
  
  # Generate statistical plots using DROMA_R
  stat_plots <- reactive({
    tryCatch({
      generateStatisticalPlots(
        projects = projects$project_name,
        plot_types = "all",
        use_gap_plots = TRUE
      )
    }, error = function(e) {
      showNotification(paste("Error generating plots:", e$message), type = "error")
      return(list())
    })
  })
  
  # Get sample and drug annotations from database
  sample_annotations <- reactive({
    tryCatch({
      # Query sample metadata from database
      conn <- get("droma_db_connection", envir = .GlobalEnv)
      sample_query <- "SELECT * FROM sample_anno"
      sample_data <- DBI::dbGetQuery(conn, sample_query)
      sample_data
    }, error = function(e) {
      showNotification(paste("Error loading sample annotations:", e$message), type = "error")
      return(data.frame())
    })
  })
  
  drug_annotations <- reactive({
    tryCatch({
      # Query treatment metadata from database
      conn <- get("droma_db_connection", envir = .GlobalEnv)
      drug_query <- "SELECT * FROM drug_anno"
      drug_data <- DBI::dbGetQuery(conn, drug_query)
      drug_data
    }, error = function(e) {
      showNotification(paste("Error loading drug annotations:", e$message), type = "error")
      return(data.frame())
    })
  })
  
  # Plot outputs ----
  output$p_count_drugandsample_facet_with_gap <- renderPlot({
    req(stat_plots())
    if ("counts" %in% names(stat_plots()) && "detailed" %in% names(stat_plots()$counts)) {
      stat_plots()$counts$detailed
    }
  })
  
  output$p_count_drugandsample_sum_with_gap2 <- renderPlot({
    req(stat_plots())
    if ("counts" %in% names(stat_plots()) && "summary" %in% names(stat_plots()$counts)) {
      stat_plots()$counts$summary
    }
  })
  
  output$p_mol_character <- renderPlot({
    req(stat_plots())
    if ("molecular" %in% names(stat_plots())) {
      stat_plots()$molecular
    }
  })
  
  output$p_overlap_drug <- renderPlot({
    req(stat_plots())
    if ("overlaps" %in% names(stat_plots()) && "drugs" %in% names(stat_plots()$overlaps)) {
      stat_plots()$overlaps$drugs
    }
  })
  
  output$p_overlap_sample <- renderPlot({
    req(stat_plots())
    if ("overlaps" %in% names(stat_plots()) && "samples" %in% names(stat_plots()$overlaps)) {
      stat_plots()$overlaps$samples
    }
  })
  
  output$p_tumor_bubble <- renderPlot({
    req(stat_plots())
    if ("tumor_types" %in% names(stat_plots())) {
      stat_plots()$tumor_types
    }
  })
  
  output$p_drug_moa <- renderPlot({
    req(stat_plots())
    if ("drug_moa" %in% names(stat_plots())) {
      stat_plots()$drug_moa
    }
  })
  
  # Table outputs ----
  output$drug_anno <- DT::renderDataTable({ 
    req(drug_annotations())
    DT::datatable(
      drug_annotations(),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: left; color: black; font-size: 14px;',
        htmltools::strong("Drug Annotation Data")
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv')
      ),
      extensions = 'Buttons',
      rownames = FALSE,
      filter = 'top'
    )
  })
  
  output$sample_anno <- DT::renderDataTable({ 
    req(sample_annotations())
    DT::datatable(
      sample_annotations(),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: left; color: black; font-size: 14px;',
        htmltools::strong("Sample Annotation Data")
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv')
      ),
      extensions = 'Buttons',
      rownames = FALSE,
      filter = 'top'
    )
  })
  
  # Download handlers for sample annotations
  output$download_sample <- downloadHandler(
    filename = function() {
      type <- input$sample_download_type
      base_name <- "DROMA_Sample_Annotations"
      
      switch(type,
             "data_rds" = paste0(base_name, ".rds"),
             "data_csv" = paste0(base_name, ".csv")
      )
    },
    content = function(filename) {
      type <- input$sample_download_type
      
      switch(type,
             "data_rds" = saveRDS(sample_annotations(), filename),
             "data_csv" = write.csv(sample_annotations(), filename, row.names = FALSE)
      )
    }
  )
  
  # Download handlers for drug annotations
  output$download_drug <- downloadHandler(
    filename = function() {
      type <- input$drug_download_type
      base_name <- "DROMA_Drug_Annotations"
      
      switch(type,
             "data_rds" = paste0(base_name, ".rds"),
             "data_csv" = paste0(base_name, ".csv")
      )
    },
    content = function(filename) {
      type <- input$drug_download_type
      
      switch(type,
             "data_rds" = saveRDS(drug_annotations(), filename),
             "data_csv" = write.csv(drug_annotations(), filename, row.names = FALSE)
      )
    }
  )
}

