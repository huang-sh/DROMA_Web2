# DROMA Web Application (v1.0) - Refactored

This is a refactored version of the DROMA Shiny web application that uses the **DROMA.R** and **DROMA.Set** R packages for backend functionality.

## Overview

DROMA (Drug Response Omics association MAp) is an interactive web application for analyzing drug-omics associations across multiple cancer datasets. This refactored version replaces custom backend logic with standardized R package functions while maintaining the original UI and module structure.

## Key Changes from Previous Version

### Architecture Improvements

1. **Database-Driven Data Loading**: All data is now loaded from SQLite database using DROMA.Set functions
2. **Standardized Analysis**: Uses DROMA.R package functions for all statistical analyses
3. **Z-score Normalization**: Handled through function parameters (`zscore=TRUE/FALSE`) instead of global data manipulation
4. **MultiDromaSet Support**: Uses MultiDromaSet objects for multi-project analysis
5. **Simplified Code**: Removed redundant custom functions in favor of package functions

### Backend Function Mapping

| Old Custom Function | New Package Function | Package |
|---------------------|---------------------|---------|
| `selFeatures()` | `loadMultiProjectMolecularProfiles()` | DROMA.Set |
| `pairDrugOmic()` | `analyzeDrugOmicPair()` | DROMA.R |
| `BatchFindSigFeaturesPlus()` | `batchFindSignificantFeatures()` | DROMA.R |
| `process_drug_data()` | `getDrugSensitivityData()` | DROMA.R |
| Custom meta-analysis | Built into `analyzeDrugOmicPair()` | DROMA.R |
| Custom plot functions | `createForestPlot()`, `plotMetaVolcano()` | DROMA.R |
| Data loading scripts | `createMultiDromaSetFromDatabase()` | DROMA.Set |

## Project Structure

```
DROMA_Web2/
├── App.R                  # Main Shiny application
├── config.yml             # Configuration file (database paths)
├── README.md              # This file
├── Modules/               # Shiny modules (UI + Server)
│   ├── DrugOmicPair.R    # Drug-Omic pair analysis module
│   ├── DrugFeature.R     # Drug feature analysis module
│   ├── BatchFeature.R    # Batch feature analysis module
│   ├── StatAnno.R        # Statistics and annotations module
│   └── GlobalSetting.R   # Global settings module
└── Functions/             # Helper functions
    └── DrugFeatureHelpers.R  # Drug comparison plot helpers
```

## Prerequisites

### Required R Packages

1. **DROMA Packages** (must be installed):
   - `DROMA.Set` - Data management and database access
   - `DROMA.R` - Analysis and visualization functions

2. **Shiny and UI packages**:
   - shiny
   - shinyWidgets
   - shinyjs
   - waiter
   - DT

3. **Data manipulation**:
   - dplyr
   - data.table

4. **Plotting**:
   - ggplot2
   - ggpubr
   - ggrepel
   - treemapify
   - gridExtra
   - grid
   - patchwork
   - cowplot
   - gg.gap
   - UpSetR

5. **Meta-analysis**:
   - meta
   - metafor
   - effsize

6. **Other**:
   - parallel
   - snowfall

### Database Setup

The application requires a DROMA SQLite database file. Update the `config.yml` file with the correct path.

## Installation

1. **Install DROMA packages** (if not already installed):

```r
# Install DROMA.Set
devtools::install_github("your-repo/DROMA_Set")

# Install DROMA.R
devtools::install_github("your-repo/DROMA_R")
```

2. **Install other required packages**:

```r
# Install CRAN packages
install.packages(c("shiny", "shinyWidgets", "shinyjs", "waiter", "DT",
                   "dplyr", "data.table", "ggplot2", "ggpubr", "ggrepel",
                   "treemapify", "gridExtra", "grid", "patchwork", "cowplot",
                   "gg.gap", "UpSetR", "meta", "metafor", "effsize",
                   "parallel", "snowfall", "config"))
```

3. **Configure database path**:

Edit `config.yml` to point to your DROMA SQLite database.

4. **Run the application**:

```r
shiny::runApp("DROMA_Web2")
```

## Module Descriptions

### 1. Drug-Omics Pairs Analysis (`DrugOmicPair.R`)

Analyzes associations between individual drugs and molecular features across multiple datasets.

**Key Features**:
- Select molecular type (mRNA, CNV, methylation, protein, mutations, fusions)
- Select specific feature and drug
- Filter by data type and tumor type
- Automatic meta-analysis across projects
- Forest plot and correlation visualizations
- Download results in multiple formats

**Backend**: Uses `analyzeDrugOmicPair()` from DROMA.R with MultiDromaSet objects.

### 2. Drug Feature Analysis (`DrugFeature.R`)

Explores drug sensitivity data with sample annotations and comparisons.

**Key Features**:
- View drug sensitivity data with annotations
- Compare drug response by tumor type, gender, ethnicity, age
- Interactive data tables
- Download annotated data

**Backend**: Uses `getDrugSensitivityData()` from DROMA.R.

### 3. Batch Features Associations Analysis (`BatchFeature.R`)

Performs batch analysis comparing one feature against all features of another type.

**Key Features**:
- Compare one feature vs. all features of selected type
- Parallel processing support
- Volcano plot visualization
- Meta-analysis across datasets
- Progress tracking

**Backend**: Uses `batchFindSignificantFeatures()` from DROMA.R.

### 4. Statistics and Annotations (`StatAnno.R`)

Displays summary statistics and dataset annotations.

**Key Features**:
- Drug and sample count visualizations
- Molecular characteristics overview
- Sample and drug overlap analysis
- Tumor type and drug MOA visualizations
- Downloadable annotation tables

**Backend**: Uses `generateStatisticalPlots()`, `generateOverlapPlots()`, `generateCountPlots()` from DROMA.R.

### 5. Global Settings (`GlobalSetting.R`)

Controls application-wide settings.

**Key Features**:
- Toggle z-score normalization on/off
- Applies settings across all modules
- Floating settings button

## Usage

### Basic Workflow

1. **Start the application**: Run `shiny::runApp("DROMA_Web2")`
2. **Configure settings**: Click "Global Settings" to enable/disable z-score normalization
3. **Select module**: Navigate to the analysis module of interest
4. **Configure analysis**: Select features, drugs, and filters
5. **View results**: Plots and tables are generated automatically
6. **Download**: Use download buttons to save results

### Z-score Normalization

The application supports z-score normalization controlled via Global Settings:

- **Enabled** (default): Data is normalized using z-scores for better cross-dataset comparison
- **Disabled**: Raw data values are used

This setting affects all data loading operations through the `zscore` parameter in DROMA.Set and DROMA.R functions.

## Important Implementation Notes

### Parameter Names

The refactored version uses standardized parameter names:
- `feature_type` (not `molecular_type`)
- `select_features` (not `features`)
- `zscore = TRUE/FALSE` (parameter-based normalization)

### Data Loading

All data loading goes through:
1. `createMultiDromaSetFromDatabase()` - Creates MultiDromaSet objects
2. `loadMultiProjectMolecularProfiles()` - Loads molecular data
3. `loadMultiProjectTreatmentResponse()` - Loads drug response data

### Meta-analysis

Meta-analysis is automatically performed when using MultiDromaSet objects with appropriate DROMA.R functions (e.g., `analyzeDrugOmicPair()`).

## Troubleshooting

### Database Connection Issues

If you see "No database connection found" errors:
1. Check that `config.yml` has the correct database path
2. Ensure the database file exists and is readable
3. Verify DROMA.Set package is installed

### Missing Package Errors

Install all required packages as listed in Prerequisites section.

### Performance Issues

For large batch analyses:
- The application automatically uses optimal number of CPU cores
- Progress bars show analysis status
- Consider using filters to reduce dataset size

## Development Notes

### Adding New Modules

1. Create module file in `Modules/` directory
2. Define `ui<ModuleName>()` and `server<ModuleName>()` functions
3. Source the module in `App.R`
4. Add to UI navigation in `App.R`

### Adding Helper Functions

Place shared utility functions in `Functions/` directory and source as needed in modules.

## Version History

- **v1.0** (2025-11-17): Refactored to use DROMA.R and DROMA.Set packages
- **v0.3** (2025-05-13): Original version with custom backend functions

## License

[Add license information]

## Contact

For bug reports or feature requests:
- Email: mugpeng@outlook.com
- Email: yc47680@um.edu.mo
- GitHub: https://github.com/mugpeng/DROMA_DB

## Acknowledgments

This application is part of the DROMA project. Visit https://droma01.github.io/DROMA/ for more information.

