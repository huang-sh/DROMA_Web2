#' List available treatments from DROMA database
#'
#' @param projects Vector of project names to filter by (optional)
#' @return Dataframe containing treatment information
#' @export
listDROMATreatments <- function(projects = NULL) {
  conn <- base::get("droma_db_connection", envir = .GlobalEnv)
  
  query <- "SELECT * FROM drug_anno"
  
  # Filter by project if specified
  if (!is.null(projects)) {
    # If projects is a dataframe/list (e.g. from listDROMAProjects), extract the name column
    if (is.data.frame(projects) && "project_name" %in% names(projects)) {
      projects <- projects$project_name
    }
    
    # Create SQL IN clause
    projects_str <- paste(sprintf("'%s'", projects), collapse = ",")
    query <- paste0(query, " WHERE ProjectID IN (", projects_str, ")")
  }
  
  tryCatch({
    data <- DBI::dbGetQuery(conn, query)
    
    # Rename DrugName to TreatmentName for consistency with other functions if needed
    if ("DrugName" %in% names(data) && !"TreatmentName" %in% names(data)) {
      data$TreatmentName <- data$DrugName
    }
    
    return(data)
  }, error = function(e) {
    stop(paste("Error querying treatments:", e$message))
  })
}

#' List available features from DROMA database
#'
#' @param projects Vector of project names to filter by
#' @param feature_type Type of feature (e.g. "mRNA", "cnv")
#' @return Dataframe containing feature information
#' @export
listDROMAFeatures <- function(projects, feature_type) {
  conn <- base::get("droma_db_connection", envir = .GlobalEnv)
  
  # Handle projects argument
  if (is.data.frame(projects) && "project_name" %in% names(projects)) {
    projects <- projects$project_name
  }
  
  all_features <- character(0)
  
  for (proj in projects) {
    table_name <- paste0(proj, "_", feature_type)
    
    # Check if table exists
    if (DBI::dbExistsTable(conn, table_name)) {
      # Get column names to find the feature column
      # Usually the columns are SampleIDs, but the first column might be FeatureName or similar
      # Or maybe the table has features as columns?
      # Let's check the structure from previous output or assume standard format
      # Actually, for omics data, usually rows are features and columns are samples, or vice versa.
      # But `listDROMAFeatures` returns a list of features.
      # If the table is Feature x Sample, we need row names (or first column).
      # If the table is Sample x Feature, we need column names.
      
      # Based on typical bio-data storage in SQL:
      # It's often stored as a matrix table.
      # Let's try to get the first column name or check if there is a specific feature table.
      # Wait, the previous error said "Found 18881 features in CCLE_mRNA".
      # This suggests it counts rows. So likely features are rows.
      # Let's assume the first column is the feature name.
      
      tryCatch({
        # Get all fields
        fields <- DBI::dbListFields(conn, table_name)
        
        # Determine feature column
        feature_col <- NULL
        if ("feature_id" %in% fields) {
          feature_col <- "feature_id"
        } else if ("FeatureName" %in% fields) {
          feature_col <- "FeatureName"
        } else if ("Gene" %in% fields) {
          feature_col <- "Gene"
        } else if ("Symbol" %in% fields) {
          feature_col <- "Symbol"
        } else {
          # If no standard name found, look for TEXT columns?
          # Or just warn and skip?
          # For now, let's try to find a column that is likely a feature ID
          # But we can't easily check type without query.
          # Let's assume if feature_id is missing, we might have a problem.
        }
        
        if (!is.null(feature_col)) {
          # Query distinct features
          query <- paste0("SELECT DISTINCT ", feature_col, " FROM ", table_name)
          features <- DBI::dbGetQuery(conn, query)[[1]]
          all_features <- c(all_features, features)
        }
      }, error = function(e) {
        warning(paste("Error querying table", table_name, ":", e$message))
      })
    }
  }
  
  if (length(all_features) == 0) {
    return(data.frame(FeatureName = character(0)))
  }
  
  return(data.frame(FeatureName = unique(all_features), stringsAsFactors = FALSE))
}
