# Modules/DrugFeature.R
# UI Component
uiDrugFeature <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    fluidRow(
      column(4,
             # Drug Selection Panel
             wellPanel(
               # Direct selectizeInput for drug selection
               selectizeInput(
                 ns("select_drug"), "Drug Selection:", choices = NULL,
                 options = list(
                   placeholder = 'Please select a drug',
                   onInitialize = I('function() { this.setValue(""); }')
                 )),
               
               # Data type and tumor type filters
               selectInput(inputId = ns("filter_data_type"), 
                           "Filter by data type:", 
                           choices = c("All" = "all",
                                       "Cell Lines" = "CellLine",
                                       "Patient-Derived Cells" = "PDC",
                                       "Patient-Derived Organoids" = "PDO",
                                       "Patient-Derived Xenografts" = "PDX"
                           ), selected = "all"),
               
               selectInput(inputId = ns("filter_tumor_type"), 
                           "Filter by tumor type:", 
                           choices = c("All" = "all",
                                       "Aerodigestive Tract Cancer" = "aerodigestive tract cancer",
                                       "Bladder Cancer" = "bladder cancer",
                                       "Breast Cancer" = "breast cancer",
                                       "Cervical Cancer" = "cervical cancer",
                                       "Choriocarcinoma" = "choriocarcinoma",
                                       "Endometrial Cancer" = "endometrial cancer",
                                       "Gastrointestinal Cancer" = "gastrointestinal cancer",
                                       "Haematopoietic/Lymphoid Cancer" = "haematopoietic/lymphoid cancer",
                                       "Kidney Cancer" = "kidney cancer",
                                       "Liver Cancer" = "liver cancer",
                                       "Lung Cancer" = "lung cancer",
                                       "Nasopharyngeal Cancer" = "nasopharyngeal cancer",
                                       "Nervous System Cancer" = "nervous system cancer",
                                       "Non-Cancer" = "non-cancer",
                                       "Ovarian Cancer" = "ovarian cancer",
                                       "Pancreatic Cancer" = "pancreatic cancer",
                                       "Prostate Cancer" = "prostate cancer",
                                       "Retinoblastoma" = "retinoblastoma",
                                       "Sarcoma" = "sarcoma",
                                       "Skin Cancer" = "skin cancer",
                                       "Stomach Cancer" = "stomach cancer",
                                       "Testicular Cancer" = "testicular cancer",
                                       "Thyroid Cancer" = "thyroid cancer",
                                       "Uterine Cancer" = "uterine cancer",
                                       "Vulvar Cancer" = "vulvar cancer"
                           ), selected = "all"),
               
               # Add "Compare by which" selection
               hr(),
               h4("Comparison Options"),
               selectInput(inputId = ns("compare_by"), 
                           "Compare by which:", 
                           choices = c("TumorType", "Gender", "FullEthnicity", "Age"),
                           selected = "TumorType"),
               
               # Download section
               hr(), # Add a horizontal rule for separation
               prettyRadioButtons(
                 inputId = ns("download_type"),
                 label = "Select download format:",
                 choices = c(
                   "Data (RDS)" = "data_rds",
                   "Data (CSV)" = "data_csv"
                 ),
                 selected = "data_csv",
                 icon = icon("check"), 
                 animation = "jelly",
                 status = "primary"
               ), 
               downloadBttn(
                 outputId = ns("download_content"),
                 label = "Download",
                 style = "gradient",
                 color = "default",
                 block = TRUE,
                 size = "sm"
               )
             )
      ),
      
      column(8,
             # Results Panels with tabs for different visualizations
             tabsetPanel(
               tabPanel("Drug Data", 
                        h4("Drug Sensitivity Data"),
                        DT::dataTableOutput(ns("drug_table"))),
               
               tabPanel("Annotated Data", 
                        h4("Drug Sensitivity with Sample Annotations"),
                        DT::dataTableOutput(ns("annotated_drug_table"))),
               
               tabPanel("Comparison", 
                        h4("Drug Sensitivity by Selected Attribute"),
                        uiOutput(ns("comparison_ui")),
                        plotOutput(ns("comparison_plot"), height = "500px")
               )
             )
      )
    )
  )
}

# Server Component
serverDrugFeature <- function(input, output, session) {
  ns <- session$ns
  
  # Get available projects from database
  projects <- listDROMAProjects()
  
  # Create MultiDromaSet (cached)
  multi_dromaset <- reactive({
    createMultiDromaSetFromDatabase(project_names = projects)
  })
  
  # Update drug selection choices
  observe({
    tryCatch({
      drugs_list <- listDROMATreatments(projects = projects)
      drugs_choices <- unique(drugs_list$TreatmentName)
      
      updateSelectizeInput(session = session, inputId = 'select_drug',
                           choices = drugs_choices, server = TRUE,
                           options = list(placeholder = 'Please select a drug'))
    }, error = function(e) {
      showNotification(paste("Error loading drugs:", e$message), type = "error")
    })
  })
  
  # Get drug sensitivity data with annotations
  drug_sensitivity_data <- reactive({
    shiny::validate(
      shiny::need(input$select_drug != "", "Please select a drug.")
    )
    
    tryCatch({
      # Use getDrugSensitivityData from DROMA_R which includes annotations
      getDrugSensitivityData(
        dromaset_object = multi_dromaset(),
        select_drugs = input$select_drug,
        data_type = input$filter_data_type,
        tumor_type = input$filter_tumor_type,
        overlap_only = FALSE,
        include_annotations = TRUE
      )
    }, error = function(e) {
      showNotification(paste("Error loading drug data:", e$message), type = "error")
      return(NULL)
    })
  })
  
  # Drug sensitivity table (basic data)
  output$drug_table <- DT::renderDataTable({
    req(drug_sensitivity_data())
    data <- drug_sensitivity_data()
    
    # Select basic columns for display
    basic_cols <- c("ProjectID", "SampleID", "TreatmentName", "raw_value", "zscore_value")
    display_data <- data[, intersect(basic_cols, names(data)), drop = FALSE]
    
    DT::datatable(
      display_data,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        dom = 'Bfrtip'
      ),
      caption = "Drug sensitivity data - showing both raw and Z-score normalized values",
      rownames = FALSE
    ) %>%
      DT::formatRound(columns = c("raw_value", "zscore_value"), digits = 3)
  })
  
  # Annotated drug table (with sample metadata)
  output$annotated_drug_table <- DT::renderDataTable({
    req(drug_sensitivity_data())
    data <- drug_sensitivity_data()
    
    DT::datatable(
      data,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        dom = 'Bfrtip'
      ),
      caption = "Drug sensitivity with sample annotations - showing both raw and Z-score normalized values",
      rownames = FALSE
    ) %>%
      DT::formatRound(columns = c("raw_value", "zscore_value"), digits = 3)
  })
  
  # Dynamic UI for comparison settings
  output$comparison_ui <- renderUI({
    req(input$compare_by, drug_sensitivity_data())
    
    # Check if the comparison variable is continuous (numeric)
    if (input$compare_by %in% names(drug_sensitivity_data()) && 
        is.numeric(drug_sensitivity_data()[[input$compare_by]])) {
      # For continuous variables
      fluidRow(
        column(6,
               sliderInput(ns("group_bins"), 
                           "Number of groups:",
                           min = 2, max = 10, value = 4, step = 1)
        ),
        column(6,
               checkboxInput(ns("show_jitter"), "Show Groups Boxplot", value = TRUE)
        )
      )
    } else {
      # For categorical variables - no additional controls needed
      NULL
    }
  })
  
  # Comparison visualization
  output$comparison_plot <- renderPlot({
    req(drug_sensitivity_data(), input$compare_by)
    data <- drug_sensitivity_data()
    
    # Check if comparison variable exists
    if (!input$compare_by %in% names(data)) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = paste("Variable", input$compare_by, "not found in data")) + 
               theme_void())
    }
    
    # Always use z-score value for plots
    plot_data <- data
    plot_data$value <- plot_data$zscore_value
    value_type_label <- "Z-score normalized drug sensitivity values"
    
    # Handle missing values in the comparison variable
    plot_data <- plot_data[!is.na(plot_data[[input$compare_by]]), ]
    
    if (nrow(plot_data) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No data available for this comparison") + 
               theme_void())
    }
    
    # Source helper function for plotting
    if (!exists("createDrugComparisonPlot")) {
      source("Functions/DrugFeatureHelpers.R", local = TRUE)
    }
    
    # Different visualization based on the type of the comparison variable
    if (is.numeric(data[[input$compare_by]])) {
      return(createDrugComparisonPlot(
        data = plot_data, 
        comparison_var = input$compare_by,
        value_column = "value",
        value_label = value_type_label,
        num_bins = input$group_bins,
        show_groups_boxplot = input$show_jitter
      ))
    } else {
      # For categorical variables
      return(createDrugComparisonPlot(
        data = plot_data,
        comparison_var = input$compare_by,
        value_column = "value",
        value_label = value_type_label
      ))
    }
  })
  
  # Download handler
  output$download_content <- downloadHandler(
    filename = function() {
      type <- input$download_type
      drug_name_safe <- gsub("[^A-Za-z0-9_]", "_", input$select_drug)
      base_name <- paste0("DrugFeature_", drug_name_safe)
      
      switch(type,
             "data_rds" = paste0(base_name, "_data.rds"),
             "data_csv" = paste0(base_name, "_data.csv")
      )
    },
    content = function(filename) {
      type <- input$download_type
      data_to_save <- drug_sensitivity_data()
      
      switch(type,
             "data_rds" = saveRDS(data_to_save, filename),
             "data_csv" = write.csv(data_to_save, filename, row.names = FALSE)
      )
    }
  )
}

