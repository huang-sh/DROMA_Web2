uiBatchFeature <- function(id){
  ns <- NS(id)
  fluidPage(
    useWaiter(), # for overall loading
    useSweetAlert(), # for notifications
    
    # Add progress elements
    tags$head(
      tags$style(HTML("
        #progressBox {
          background-color: #f8f9fa;
          padding: 15px;
          border-radius: 5px;
          margin-bottom: 20px;
        }
        .progress {
          margin-bottom: 10px;
        }
      "))
    ),
    fluidRow(
      # Select Features 1 ----
      column(4,
             selectInput(inputId = ns("select_features1"), 
                         "Please select the feature type:", 
                         choices = c("Copy Number Data" = "cnv",
                                     "DNA Methylation" = "meth",
                                     "Gene Fusion" = "fusion",
                                     "Gene Mutation" = "mutation_gene",
                                     "Gene Site Mutation" = "mutation_site",
                                     "mRNA Expression" = "mRNA",
                                     "Protein RPPA Expression" = "proteinrppa",
                                     "Protein MS Expression" = "proteinms",
                                     "Drug Sensitivity" = "drug"
                         ), selected = "proteinms"
             )),
      # Select specific feature ----
      column(4,
             selectizeInput(
               ns("select_specific_feature"), "Features Selection:", choices = NULL,
               options = list(
                 placeholder = 'Please select a feature',
                 onInitialize = I('function() { this.setValue(""); }'), selected = ""
               ))),
      # Select Features 2 ----
      column(4,
             selectInput(inputId = ns("select_features2"), 
                         "Please select the second feature:", 
                         choices = c("Copy Number Data" = "cnv",
                                     "DNA Methylation" = "meth",
                                     "Gene Fusion" = "fusion",
                                     "Gene Mutation" = "mutation_gene",
                                     "Gene Site Mutation" = "mutation_site",
                                     "mRNA Expression" = "mRNA",
                                     "Protein RPPA Expression" = "proteinrppa",
                                     "Protein MS Expression" = "proteinms",
                                     "Drug Sensitivity" = "drug"
                         ), selected = "drug"
             )),
    ),
    # Add data_type and tumor_type filters
    fluidRow(
      column(6,
             selectInput(inputId = ns("data_type"), 
                         "Filter by data type:", 
                         choices = c("All" = "all",
                                     "Cell Lines" = "CellLine",
                                     "Patient-Derived Cells" = "PDC",
                                     "Patient-Derived Organoids" = "PDO",
                                     "Patient-Derived Xenografts" = "PDX"
                         ), selected = "all"
             )),
      column(6,
             selectInput(inputId = ns("tumor_type"), 
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
                         ), selected = "all"
             ))
    ),
    fluidRow(
      column(12,
             hidden(
               div(id = ns("progressBox"),
                   h4("Analysis Progress"),
                   progressBar(ns("progressBar"), value = 0, 
                               title = "Calculating...",
                               display_pct = TRUE) 
               )
             )
      )
    ),
    # Output results ----
    wellPanel(
      column(12,
             tabsetPanel(
               tabPanel("Volcano Plot",
                        plotOutput(ns("volcano_plot"), height = "600px"))
             )),
      h5(".")
    ),
    # Download section
    fluidRow(
      column(6,
             prettyRadioButtons(
               inputId = ns("download_type"),
               label = "Select download format:",
               choices = c(
                 "PDF" = "pdf",
                 "CSV Results" = "csv",
                 "Plot Data" = "plot_data"
               ),
               selected = "pdf",
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
      ),
      column(6,
             h4(strong("NOTEs:")),
             h5("1. The analysis compares one feature against all features in the selected dataset type."),
             h5("2. The xaxis and yaxis are effect size and p value calculated from Meta analysis results."),
             h5("3. Multiple cores can speed up computation but uses more memory.")
      )
    )
  )
}

serverBatchFeature <- function(input, output, session){
  ns <- session$ns
  
  # Get available projects from database
  projects <- listDROMAProjects(exclude_clinical = T)
  
  # Filter projects to only those with drug data (drug list doesn't change, so compute once)
  projects_with_drug <- listDROMAProjects(feature_type = "drug", show_names_only = TRUE, exclude_clinical = T)
  filtered_projects_drug <- projects[projects$project_name %in% projects_with_drug, ]
  
  # Create MultiDromaSet (cached)
  multi_dromaset <- reactive({
    db_path <- config::get()$db_path
    createMultiDromaSetFromDatabase(project_names = projects$project_name, db_path = db_path)
  })
  
  # Track z-score changes
  zscore_tracker <- reactiveVal(if(exists("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)) 
                               base::get("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)$timestamp 
                               else Sys.time())
  
  # Update tracker when global state changes
  observe({
    if(exists("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)) {
      zscore_tracker(base::get("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)$timestamp)
    }
    invalidateLater(1000) # Check every second
  })
  
  # Add reactive values for progress tracking
  progress_vals <- reactiveValues(
    start_time = NULL,
    features_done = 0,
    total_features = 0,
    estimated_time = "Calculating..."
  )
  
  # Select Features reactive ----
  features_search_sel <- reactiveValues()
  
  observeEvent(input$select_features1, {
    tryCatch({
      # Get features list from database
      if (input$select_features1 == "drug") {
        # Handle multiple projects: loop through and combine results
        all_drugs <- character(0)
        for (proj in filtered_projects_drug$project_name) {
          drugs_vec <- listDROMAFeatures(projects = proj, feature_type = "drug")
          all_drugs <- c(all_drugs, drugs_vec)
        }
        features_search_sel$features <- unique(all_drugs)
      } else {
        # Filter projects to only those with this feature type
        projects_with_feature <- listDROMAProjects(feature_type = input$select_features1, show_names_only = TRUE, exclude_clinical = T)
        filtered_projects <- projects[projects$project_name %in% projects_with_feature, ]
        # Handle multiple projects: loop through and combine results
        all_features <- character(0)
        for (proj in filtered_projects$project_name) {
          features_vec <- listDROMAFeatures(
            projects = proj,
            feature_type = input$select_features1
          )
          all_features <- c(all_features, features_vec)
        }
        features_search_sel$features <- unique(all_features)
      }
      
      updateSelectizeInput(session = session, 
                           inputId = 'select_specific_feature',
                           label = 'Features Selection:', 
                           choices = features_search_sel$features, 
                           server = TRUE,
                           selected = "")
    }, error = function(e) {
      showNotification(paste("Error loading features:", e$message), type = "error")
    })
  })
  
  # Calculate results ----
  results <- reactive({
    # React to z-score changes
    zscore_tracker()
    
    req(input$select_features1, input$select_features2, input$select_specific_feature)
    
    # Show progress box
    shinyjs::show("progressBox")
    progress_vals$start_time <- Sys.time()
    progress_vals$features_done <- 0
    
    # Simple progress indication
    updateProgressBar(session = session,
                      id = ns("progressBar"),
                      value = 50,
                      title = "Processing...")
    
    # Check if z-score is enabled
    zscore_enabled <- TRUE
    if(exists("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)) {
      zscore_enabled <- isTRUE(base::get("GLOBAL_ZSCORE_STATE", envir = .GlobalEnv)$enabled)
    }
    
    # Use optimal number of cores
    used_core <- ifelse(parallel::detectCores()/2 > 8, 8, parallel::detectCores()/2)
    
    # Explicitly get feature list for feature2
    feature2_list <- NULL
    tryCatch({
      if (input$select_features2 == "drug") {
        # Handle multiple projects: loop through and combine results
        all_drugs <- character(0)
        for (proj in filtered_projects_drug$project_name) {
          drugs_vec <- listDROMAFeatures(projects = proj, feature_type = "drug")
          all_drugs <- c(all_drugs, drugs_vec)
        }
        feature2_list <- unique(all_drugs)
      } else {
        # Filter projects to only those with this feature type
        projects_with_feature <- listDROMAProjects(feature_type = input$select_features2, show_names_only = TRUE, exclude_clinical = T)
        filtered_projects <- projects[projects$project_name %in% projects_with_feature, ]
        # Handle multiple projects: loop through and combine results
        all_features <- character(0)
        for (proj in filtered_projects$project_name) {
          features_vec <- listDROMAFeatures(
            projects = proj,
            feature_type = input$select_features2
          )
          all_features <- c(all_features, features_vec)
        }
        feature2_list <- unique(all_features)
      }
    }, error = function(e) {
      showNotification(paste("Error retrieving feature list:", e$message), type = "error")
      return(NULL)
    })
    
    req(feature2_list)
    
    # Use batchFindSignificantFeatures from DROMA_R
    tryCatch({
      results <- batchFindSignificantFeatures(
        dromaset_object = multi_dromaset(),
        feature1_type = input$select_features1,
        feature1_name = input$select_specific_feature,
        feature2_type = input$select_features2,
        feature2_name = feature2_list,  # Pass explicit list
        data_type = input$data_type,
        tumor_type = input$tumor_type,
        overlap_only = FALSE,
        cores = used_core,
        show_progress = TRUE,
        test_top_n = NULL  # NULL means all features
      )
      
      # Update progress to complete
      updateProgressBar(session = session,
                        id = ns("progressBar"),
                        value = 100,
                        title = "Complete!")
      
      # Hide progress box when done
      shinyjs::hide("progressBox")
      
      # Show completion message
      if (!is.null(results) && nrow(results) > 0) {
        sendSweetAlert(
          session = session,
          title = "Analysis Complete",
          text = sprintf("Processed %d features in %s",
                         nrow(results),
                         format(difftime(Sys.time(), progress_vals$start_time), digits = 2)),
          type = "success"
        )
      }
      
      results
    }, error = function(e) {
      # Hide progress box on error
      shinyjs::hide("progressBox")
      showNotification(paste("Error in batch analysis:", e$message), type = "error")
      return(NULL)
    })
  })
  
  # Render volcano plot ----
  p_volcano <- reactive({
    req(results())
    plotMetaVolcano(results(),
                    es_t = 0.2,
                    P_t = 0.001,
                    label = TRUE,
                    top_label_each = 5,
                    title = paste(input$select_features1, input$select_specific_feature,
                                  "vs", input$select_features2))
  })
  
  output$volcano_plot <- renderPlot({
    p_volcano()
  })
  
  # Download handler ----
  output$download_content <- downloadHandler(
    filename = function() {
      type <- input$download_type
      base_name <- paste0(input$select_features1, "_", 
                          input$select_specific_feature, "_vs_",
                          input$select_features2)
      
      switch(type,
             "pdf" = paste0(base_name, "_volcano.pdf"),
             "csv" = paste0(base_name, "_results.csv"),
             "plot_data" = paste0(base_name, "_data.rds"))
    },
    content = function(file) {
      type <- input$download_type
      
      switch(type,
             "pdf" = {
               ggsave(file, 
                      plot = p_volcano(),
                      width = 10, height = 8)
             },
             "csv" = {
               write.csv(results(), file, row.names = FALSE)
             },
             "plot_data" = {
               saveRDS(results(), file)
             })
    }
  )
}

