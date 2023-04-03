
# The All-Seeing Application

if (!require('rsconnect')) { 
  install.packages(c('rsconnect', 'pacman')) }

# Packing the essentials
pacman::p_load(
  
  # Wrangling
  stringr,
  readr,
  magrittr,
  dplyr,
  plotly,
  data.table,
  tidyr,
  
  # Plotting
  ggplot2,
  
  # Shiny
  shiny,
  shinydashboard,
  shinyWidgets,
  shinyjs,
)

# Max file size is 1 GB:
options(shiny.maxRequestSize = 1000 * 1024 ^ 2)

# Dashboard
ui <- dashboardPage(
  dashboardHeader(title = "Murphy Lab All-See"),
  
  # Sidebar
  dashboardSidebar(
    tags$head(tags$style(
      HTML(
        '.main-header > .navbar {        margin-left: 300px;      }      .main-header .logo {         width: 300px;   font-family: Source Sans Pro,sans-serif;font-size: 18px;   }      .content-wrapper, .main-footer, .right-side {          margin-left: 300px;      }    '
      )
    )),
    width = 300,
    sidebarMenu(
      menuItem(
        "Upload peptides",
        tabName = "Uploads",
        icon = icon("arrow-up", lib = "glyphicon"),
        badgeLabel = "Step 1a",
        badgeColor = "green"
      ),
      menuItem(
        "Peptide insights",
        tabName = "QCDash",
        icon = icon("search"),
        badgeLabel = "Step 1b",
        badgeColor = "blue"
      ),
      menuItem(
        "Upload proteins",
        tabName = "Pep-See",
        icon = icon("arrow-up", lib = "glyphicon"),
        badgeLabel = "Step 2a",
        badgeColor = "purple"
      ),
      menuItem(
        "Protein insights",
        icon = icon("search"),
        tabName = "widgets",
        badgeLabel = "Step 2b",
        badgeColor = "red"
      )
    )
  ),
  
  # Each tab has a corresponding tabItem
  dashboardBody(useShinyjs(),
                tabItems(
                  tabItem(tabName = "Uploads",
                          # To place multiple boxes in the same row, nest boxes within the same fluidRow
                          fluidRow(
                            # Tabs are comprised of boxes like this one
                            box(
                              title = "Upload peptide-spectrum match",
                              fileInput(
                                "psm",
                                "Pick your file",
                                # Permissible file types
                                accept = c(
                                  "text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv",
                                  "text/tab-separated-values",
                                  "application/vnd.ms-excel"
                                ),
                                multiple = TRUE
                              ),
                              checkboxInput("header", "Header?", TRUE))
                          ),
                          # Loading bar pb4
                          fluidRow(
                            box(
                              title = "Precisely what percentage is loaded:",
                              status = "primary",
                              progressBar(
                                id = "pb4",
                                value = 0,
                                display_pct = T
                              )
                            )
                          )), 
                  tabItem(tabName = "QCDash",
                          fluidRow(box(
                            title = "Figures", status = "primary", tableOutput("contents")
                          )),
                          fluidRow(
                            box(
                              # Summed intensity 
                              title = "Summed Intensities",
                              status = "primary",
                              plotlyOutput("intensiTMT")
                            ),
                            box(
                              # Delta mass
                              title = "Delta Masses",
                              status = "primary",
                              plotlyOutput("massHisteria")
                            ),
                            box(
                              # Relative intensity
                              title = "Relative Intensities", 
                              status = "primary", 
                              plotlyOutput("seeBars")
                              ),
                            box(
                              # Individual intensity
                              title = "Individual Intensities",
                              status = "primary",
                              plotlyOutput("seeRainbow16")
                            )
                            
                          )),
                  tabItem(tabName = "Pep-See",
                          class = "visible",
                          fluidRow(
                            box(
                              title = "Upload protein results",
                              fileInput(
                                "protein",
                                "Pick your file",
                                accept = c(
                                  "text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv",
                                  "text/tab-separated-values",
                                  "application/vnd.ms-excel"
                                ),
                                multiple = TRUE
                              ),
                              checkboxInput("header", "Header?", TRUE)
                            )
                          )),
                  # Protein results tab is updated after file is uploaded
                  tabItem(tabName = "widgets",
                          uiOutput("widgets_tab_content"))
                )
  )
)


server <- function(input, output, session) {
  
  # Create peptides data frame
  peptide_data <- reactive({
    req(input$psm, input$header, file.exists(input$psm$datapath))
    
    # Get the file extension of the uploaded file
    ext <- tools::file_ext(input$psm$name)
    
    # Check if the file extension is csv
    if(ext == "csv") {
      # If the file is a csv, read it into R using read_csv function
      localPep <- read_csv(input$psm$datapath, col_names = input$header)
    } 
    # If the file extension is tsv
    else if(ext == "tsv") {
      # If the file is a tsv, read it into R using read_tsv function
      localPep <- read_tsv(input$psm$datapath, col_names = input$header)
    } 
    # If the file extension is neither csv nor tsv
    else {
      # Stop the execution and display an error message if the file type is unsupported.
      stop("Unsupported file type. Please upload a CSV or TSV file.")
    }
    
    # Update the progress bar with an ID of "pb4" and a value of 15.
    updateProgressBar(id = "pb4", value = 15)
    
    localPep
  })
  
  # Correction factors for Paula's data
  correction_factors <- reactive({
    sums <- peptide_data() %>% 
      select(matches("(MDA|MCF|sample|TMT|tmt|Sample)")) %>%
      colSums()
    # correction factors from the column sums
    max(sums)/sums
  })
  
  # Unnormalized peptide data
  peptide_data_unnorm <- reactive({
    req(peptide_data())
    peptide_data() %>% 
      select(matches("(MDA|MCF|sample|TMT|tmt|Sample)"))
  })
  
  # Peptide data normalized using correction factors
  peptide_data_norm <- reactive({
    req(peptide_data_unnorm(), correction_factors())
    sweep(peptide_data_unnorm(), 2, correction_factors(), "*")
  })
  
  
  ## (A Paula exclusive) Select columns containing MCF and MDA values into n_peptide_data
  n_peptide_data <- reactive({
    peptide_data() %>%
      select(-c(MCF_ctrl_1:MDA_KD_4)) %>%
      cbind(peptide_data_norm())
  })
  
  ## (Another Paula exclusive) Create new columns that store the mean values
  n_peptide_data <- reactive({
    n_peptide_data() %>%
      mutate(MCF_log_KD_ctrl = log2((MCF_KD_1 + MCF_KD_2 + MCF_KD_3 + MCF_KD_4) /
                                      (MCF_ctrl_1 + MCF_ctrl_2 + MCF_ctrl_3 + MCF_ctrl_4)
      )) %>%
      mutate(MDA_log_KD_ctrl = log2((MDA_KD_1 + MDA_KD_2 + MDA_KD_3 + MDA_KD_4) /
                                      (MDA_ctrl_1 + MDA_ctrl_2 + MDA_ctrl_3 + MDA_ctrl_4)
      ))
  })
  
  
  # Creates data frame for protein results
  protein_data <- reactive({
    req(input$protein, input$header, file.exists(input$protein$datapath))
    
    ext <- tools::file_ext(input$protein$name)
    if(ext == "csv") {
      localPep <- read_csv(input$protein$datapath, col_names = input$header)
    } else if(ext == "tsv") {
      localPep <- read_tsv(input$protein$datapath, col_names = input$header)
    } else {
      stop("Unsupported file type. Please upload a CSV or TSV file.")
    }
    
    localPep
  })
  
  # Protein Correction Factors
  protein_correction_factors <- reactive({
    sums <- protein_data() %>% 
      select(matches("(MDA|MCF|sample|TMT|tmt|Sample)")) %>%
      colSums()
    ##correction factors from the column sums
    max(sums)/sums
  })
  
  # Unnormalized and k means clustered data
  
  # peptide_data_unnorm <- reactive({
  #   req(peptide_data())
  #   peptide_data() %>% 
  #     select(MCF_ctrl_1:MDA_KD_4)
  # })
  
  protein_data_unnorm <- reactive({
    req(protein_data())

    Proteo <- protein_data()

    Proteo$Gene <- make.unique(Proteo$Gene)

    Proteo$Gene[is.na(Proteo$Gene)] <- "Unknown"

    print(colnames(Proteo))

    # # Select columns and set row names from the "Gene" column
    Proteo <- Proteo %>%
      select(matches("(Gene|MDA|MCF|sample|TMT|tmt|Sample)")) %>%
      column_to_rownames("Gene")

    Proteo %>%
      select(matches("(MDA|MCF|sample|TMT|tmt|Sample)"))
    # %>%
    #   column_to_rownames("Gene")

    # Select columns for k-means clustering
    clustering_data <- Proteo %>%
      select(matches("(MDA|MCF|sample|TMT|tmt|Sample)"))

    # Perform k-means clustering
    set.seed(123)
    k <- 10 # Set the number of clusters
    kmeans_result <- kmeans(clustering_data, centers = k)

    # Add cluster assignments to the data frame
    Proteo$cluster <- kmeans_result$cluster

    Proteo
  })
  
  # Unnormalized and hierarchically clustered data
  
  # protein_data_unnorm <- reactive({
  #   # Load and preprocess protein_data
  #   req(protein_data())
  #   Proteo <- protein_data()
  #   
  #   # Process Gene column
  #   Proteo$Gene <- make.unique(Proteo$Gene)
  #   Proteo$Gene[is.na(Proteo$Gene)] <- "Unknown"
  #   
  #   # Print column names
  #   print(colnames(Proteo))
  #   
  #   # Select and process columns
  #   Proteo <- Proteo %>%
  #     select(matches("(Gene|MDA|MCF|sample|TMT|tmt|Sample)")) %>%
  #     column_to_rownames("Gene")
  #   
  #   # Final selection
  #   Proteo <- Proteo %>%
  #     select(matches("(MDA|MCF|sample|TMT|tmt|Sample)"))
  #   
  #   # Apply hierarchical clustering
  #   distance_matrix <- dist(Proteo, method = "euclidean")
  #   clustering_result <- hclust(distance_matrix, method = "average")
  #   
  #   # Reorder rows based on clustering results
  #   clustered_data <- Proteo[clustering_result$order, ]
  #   
  #   # Return clustered data
  #   clustered_data
  # })
  # 
  
  # Normalize protein results
  protein_data_norm <- reactive({
    req(protein_data_unnorm(), protein_correction_factors())
    sweep(protein_data_unnorm(), 2, protein_correction_factors(), "*")
  })
  
  n_protein_data <- reactive({
    protein_data() %>% 
      select(-c(MCF_ctrl_1:MDA_KD_4)) %>% 
      cbind(protein_data_norm())
  })
  
  ## make some new columns that store the mean values 
  n_protein_data <- reactive({
    n_protein_data() %>%
      mutate(MCF_log_KD_ctrl = log2((MCF_KD_1 + MCF_KD_2 + MCF_KD_3 + MCF_KD_4) /
                                      (MCF_ctrl_1 + MCF_ctrl_2 + MCF_ctrl_3 + MCF_ctrl_4)
      )) %>%
      mutate(MDA_log_KD_ctrl = log2((MDA_KD_1 + MDA_KD_2 + MDA_KD_3 + MDA_KD_4) /
                                      (MDA_ctrl_1 + MDA_ctrl_2 + MDA_ctrl_3 + MDA_ctrl_4)
      ))
  })
  
  
  # Define reactive expression for flipped dataframe
  protein_data_flipped <- reactive({
    temp_prot <- as.data.frame(t(protein_data()))
    names(temp_prot) <- temp_prot["Gene", ] # make row Gene into column names
    picked_prot <-
      temp_prot[str_which(rownames(temp_prot), "MDA|MCF|sample|TMT|tmt|Sample"), ]
    
    
    proats <-
      protein_data() %>% select(matches("(Gene|MDA|MCF|sample|TMT|tmt|Sample)"))
    protein_data_unnorm_t <-
      t(proats) # flip the data frame on its side
    
    colnames(protein_data_unnorm_t) <-
      protein_data_unnorm_t["Gene", ] # make row Gene into column names
    
    protein_data_unnormz <-
      protein_data_unnorm_t[-which(names(protein_data_unnorm_t) == "Gene"), ]
    
    picked_prot
  })
  
  
  ## Heat map data wrangling
  
  scaled_pdn_df <- reactive({
    pdn_df <- protein_data_norm()
    
    # Make intensities relative
    normalized_pdn_df <- t(apply(pdn_df, 1, function(row) row / max(row)))
    normalized_pdn_df[is.na(normalized_pdn_df) | is.nan(normalized_pdn_df) | is.infinite(normalized_pdn_df)] <- 0
    
    # Scale the rows
    scaled_pdn_df <- t(scale(t(normalized_pdn_df)))
    scaled_pdn_df[is.na(scaled_pdn_df) | is.nan(scaled_pdn_df) | is.infinite(scaled_pdn_df)] <- 0
    
    return(scaled_pdn_df)
  })
  
  
  
  uploaded_file <- reactive({
    !is.null(input$protein)
  })
  
  # If no gene is searched for in the proteins search on the Protein Insights page, set default as "peace" -> will default to first protein in dataframe
  # If a gene is searched for, set that as variable "findation", to be used later by output$foundation
  findation <- reactive({
    if (input$query == "") {
    placehold <- "peace"
    placehold
    } else {
      placehold <- input$query
      print(placehold)
      placehold
    }
  })
  
  # Render the "Widgets" tab content only if a file has been uploaded
  output$widgets_tab_content <- renderUI({
    if (uploaded_file()) {
      tabItem(
        tabName = "Search",
        fluidRow(
          box(
            title = "Find your protein",
            background = "light-blue",
            "Just search for it."
          )
        ),
        fluidRow(
          box(
            status = "primary",
            selectInput(
              'query',
              'Search proteome by gene symbol',
              choices = protein_data_flipped() %>% colnames()
            )
          )
        ),
        fluidRow(
          box(
            status = "primary",
            plotlyOutput("foundation")
          )
        ),
        fluidRow(
          box(
            status = "primary",
            plotlyOutput("protein_heatmap")
          ),
          box(
            status = "primary",
            plotlyOutput("mda_volcano")
          )
        ),
      )
    }
  })
  
  
  
  
  ## create the factor vector
  f1 <- reactive({
    factor(c(0, 0, 0, 0, 1, 1, 1, 1))
  })
  
  ## select the MCF columns and perform t-tests
  Stats_MCF <- reactive({
    p_MCF <- apply(n_peptide_data()[, c("MCF_ctrl_1", "MCF_ctrl_2", "MCF_KD_1", "MCF_KD_2", "MCF_ctrl_3", "MCF_KD_3", "MCF_ctrl_4", "MCF_KD_4")], 1, function(x) t.test(x[1:4], x[5:8])$p.value)
    p.adjust(p_MCF, method = "BH")
  })
  
  ## select the MDA columns and perform t-tests
  Stats_MDA <- reactive({
    p_MDA <- apply(n_peptide_data()[, c("MDA_ctrl_1", "MDA_ctrl_2", "MDA_KD_1", "MDA_KD_2", "MDA_ctrl_3", "MDA_KD_3", "MDA_ctrl_4", "MDA_KD_4")], 1, function(x) t.test(x[1:4], x[5:8])$p.value)
    p.adjust(p_MDA, method = "BH")
  })
  
  ## create data frames from the results
  Stats_MCF <- reactive({
    data.frame(p.value = Stats_MCF())
  })
  
  Stats_MDA <- reactive({
    data.frame(p.value = Stats_MDA())
  })
  
  ## select the p.values column
  Stats_MCF <- reactive({
    select(Stats_MCF(), p.value)
  })
  
  Stats_MDA <- reactive({
    select(Stats_MDA(), p.value)
  })
  
  #### COMBINE STATS WITH DATA
  low_cutoff <- reactiveVal(-0.5)
  
  ## MCF-7
  n_peptide_data_statistics <- reactive({
    n_peptide_data() %>%
      slice(1:500) %>% ## only take the first 500 rows
      cbind(Stats_MCF()) %>% ## takes the pvalues and binds them to the dataset a new row
      rename(MCF_pval=p.value) %>%
      mutate(MCF_adj_p=p.adjust(MCF_pval,method="BH")) %>% # adjust the p values for multiple-hypothesis testing
      mutate(MCF_log_adj_p=-log10(MCF_adj_p)) %>%
      mutate(MCF_log_p=-log10(MCF_pval)) %>%
      mutate(MCF_significant=case_when(MCF_log_p>1.3 & (MCF_log_KD_ctrl<low_cutoff() | MCF_log_KD_ctrl>0.5) ~ "significant", TRUE ~ "not significant"))
  })
  
  output$mcf_significant_table <- renderTable({
    n_peptide_data_statistics() %>% group_by(MCF_significant) %>%
      summarise(n=n())
  })
  
  ## MDA 231
  n_peptide_data_stats <- reactive({
    n_peptide_data_statistics() %>%
      cbind(Stats_MDA()) %>% ## takes the pvalues and binds them to the dataset a new row
      rename(MDA_pval=p.value) %>%
      mutate(MDA_adj_p=p.adjust(MDA_pval,method="BH")) %>% # adjust the p values for multiple-hypothesis testing
      mutate(MDA_log_adj_p=-log10(MDA_adj_p)) %>%
      mutate(MDA_log_p=-log10(MDA_pval)) %>%
      mutate(MDA_significant=case_when(MDA_log_p>1.3 & (MDA_log_KD_ctrl<low_cutoff() | MDA_log_KD_ctrl>0.5) ~ "significant", TRUE ~ "not significant")) %>% group_by(MDA_significant) %>%
      summarise(n=n())
  })
  
  
  output$mda_significant_table <- renderTable({
    n_peptide_data_stats() %>% group_by(MDA_significant) %>%
      summarise(n=n())
  })
  

  
  
  # modpep <- reactive ({peptide_data() %>% pull("Modified.Peptide")})
  # 
  # # N-TERMINI LABELLED
  # NLabel <- reactive({
  #   updateProgressBar(id = "pb4", value = 20)
  #   hasTag <- str_detect(modpep(), "305")
  #   labelledNumber <- as.numeric(hasTag)
  #   labelledNumber[is.na(labelledNumber)] <- 0
  #   updateProgressBar(id = "pb4", value = 35)
  #   sum(labelledNumber) / length(hasTag) %>% round(digits = 2)
  # })
  
  Nmod <- reactive({
    if ("Modified Peptide" %in% colnames(peptide_data())) {
      peptide_data() %>% pull("Modified Peptide")
    } else {
      peptide_data() %>% pull("Razor Assigned Modifications")
    }
  })
  
  # N-termini labelling efficiency
  NLabel <- reactive ({
    updateProgressBar(id = "pb4", value = 20)
    hasTag <- str_detect(Nmod(), "N-term\\(304\\.2071\\)|305")
    labelledNumber <- as.numeric(hasTag)
    labelledNumber[is.na(labelledNumber)] <- 0
    updateProgressBar(id = "pb4", value = 35)
    sum(labelledNumber) / length(hasTag) %>% round(digits = 2)
  })
  
  # # LYSINES LABELLED
  # KLabel <- reactive({
  #   hasTag <- str_detect(modpep(), "432")
  #   labelledNumber <- as.numeric(hasTag)
  #   labelledNumber[is.na(labelledNumber)] <- 0
  #   updateProgressBar(id = "pb4", value = 40)
  #   sum(labelledNumber) / sum(str_count(modpep(), "K")) %>% round(digits = 2)
  # })
  # 
  # # LABELLING EFFICIENCY
  # labellingEfficiency <- reactive({
  #   hasTag <- str_detect(modpep(), "432|305")
  #   labelledNumber <- as.numeric(hasTag)
  #   labelledNumber[is.na(labelledNumber)] <- 0
  #   updateProgressBar(id = "pb4", value = 55)
  #   sum(labelledNumber) / length(modpep()) %>% round(digits = 2)
  # })
  # 
  # # DIGEST EFFICIENCY
  # digestEfficiency <- reactive({
  #   digpep <- peptide_data() %>% pull("Peptide")
  #   
  #   hasntCut <- str_detect(digpep, "R\\w|K\\w")
  #   cutNumber <- as.numeric(hasntCut)
  #   cutNumber[is.na(cutNumber)] <- 0
  #   updateProgressBar(id = "pb4", value = 60)
  #   digestEfficiency <- 1 - (sum(cutNumber) / length(hasntCut))
  # })
  
  modpep <- reactive({
    if ("Modified Peptide" %in% names(peptide_data())) {
      peptide_data() %>% pull("Modified Peptide")
    } else {
      character()
    }
  })
  
  KLabel <- reactive({
    if (length(modpep()) > 0) {
      hasTag <- str_detect(modpep(), "432")
      labelledNumber <- as.numeric(hasTag)
      labelledNumber <- replace(labelledNumber, is.na(labelledNumber) | is.null(labelledNumber), 0)
      # labelledNumber[is.na(labelledNumber)] <- 0
      # updateProgressBar(id = "pb4", value = 40)
      totalK <- str_count(modpep(), "K")

      sum(labelledNumber) / sum(replace(totalK, is.na(totalK) | is.null(totalK), 0)) %>% round(digits = 2)
      # sum(labelledNumber) / sum(str_count(modpep(), "K")) %>% round(digits = 2)
    } else {
      0
    }
  })
  
  labellingEfficiency <- reactive({
    if (length(modpep()) > 0) {
      hasTag <- str_detect(modpep(), "432|305")
      labelledNumber <- as.numeric(hasTag)
      labelledNumber[is.na(labelledNumber)] <- 0
      updateProgressBar(id = "pb4", value = 55)
      sum(labelledNumber) / length(modpep()) %>% round(digits = 2)
    } else {
      0
    }
  })
  
  digestEfficiency <- reactive({
    if (!("Peptide" %in% names(peptide_data()))) {
      return(0)
    }
    
    digpep <- peptide_data() %>% pull("Peptide")
    if (length(digpep) > 0) {
      hasntCut <- str_detect(digpep, "R\\w|K\\w")
      cutNumber <- as.numeric(hasntCut)
      cutNumber[is.na(cutNumber)] <- 0
      updateProgressBar(id = "pb4", value = 60)
      1 - (sum(cutNumber) / length(hasntCut))
    } else {
      0
    }
  })
  
  
  # EFFICIENCIES LIST - contains digest and labelling efficiencies
  
  # efficiencies <- reactive({
  #   effiList <- data.frame(
  #     Digest = c(round(digestEfficiency() * 100, digits = 2)),
  #     Labelling = c(labellingEfficiency() * 100),
  #     KLabelling = c(KLabel() * 100),
  #     NLabelling = c(NLabel() * 100)
  #   )
  #   effiList
  # })
  
  efficiencies <- reactive({
    effiList <- data.frame(
      Digest = c(round(digestEfficiency() * 100, digits = 2)),
      Labelling = c(labellingEfficiency() * 100),
      KLabelling = c(KLabel() * 100),
      NLabelling = c(NLabel() * 100)
    )
    print("Hey")
    effiList
  })
  
  # sampling <- 
  #   # reactive({
  #   str_which(names(peptide_data_norm), "MDA|MCF|sample|TMT|tmt|Sample")
  # # })
  
  sampling <- reactive({
    print(peptide_data_unnorm())
    str_which(names(peptide_data_norm()), "MDA|MCF|sample|TMT|tmt|Sample") })
  
  # samplingL <- reactive({
  #   length(sampling())
  # })
  
  # Sum TMT intensities
  tmtConcatAll <- reactive({
    rowSums(peptide_data_unnorm()[, c(as.numeric(ncol(peptide_data_unnorm()))-(length(sampling()))):as.numeric(ncol(peptide_data_unnorm()))], na.rm = FALSE)
  })
  ##Got rid of sampling()-1
  
  # Delta masses
  # deltaMass <- reactive({
  #   req(peptide_data())
  #   if ("Delta Mass" %in% colnames(peptide_data())) {
  #     peptide_data()$Delta.Mass
  #   } else {
  #     NULL
  #   }
  # })
  
  deltaMass <- reactive({
    req(peptide_data())
    if ("Delta Mass" %in% colnames(peptide_data())) {
      peptide_data() %>% 
        select("Delta Mass") %>% 
        pull()
    } else {
      NULL
    }
  })
  
  # TMTensity <- reactive({
  #   req(peptide_data())
  #   
  #   preptides <-
  #     as.data.frame(peptide_data() %>% select(sampling()))
  #   
  #   for (k in 1:ncol(preptides)) {
  #     # Get 90th percentile intensity score/filtering for top 10 percent
  #     thresholdIntensity <- quantile(preptides[, k], probs = 0.9)
  #     sorted <- sort(preptides[, k], decreasing = TRUE)
  #     cutAndSorted <- sorted[sorted > thresholdIntensity && sorted > 0]
  #     
  #     names <- paste("cutAndSorted", k, sep = "")
  #     assign(names, cutAndSorted)
  #   }
  #   
  #   seeIt<- list()
  #   
  #   seeIt <- dplyr::bind_cols(mget(paste0('cutAndSorted', 1:samplingL())))
  #   
  #   return(seeIt)
  # })
  
  TMTensity <- reactive({
    req(peptide_data_unnorm())
    
    peptides <-
      as.data.frame(peptide_data_unnorm() %>% select(sampling()))
    
    for (k in 1:ncol(peptides)) {
      # Get 90th percentile intensity score/filtering for top 10 percent
      thresholdIntensity <- quantile(peptides[, k], probs = 0.9)
      sorted <- sort(peptides[, k], decreasing = TRUE)
      cutAndSorted <- sorted[sorted > thresholdIntensity && sorted > 0]
      
      names <- paste("cutAndSorted", k, sep = "")
      assign(names, cutAndSorted)
    }
    
    seeIt <- list()
    
    seeIt <- dplyr::bind_cols(mget(paste0('cutAndSorted', 1:length(sampling()))))
    
    return(seeIt)
  })
  
  # Bars <- reactive({
  #   req(peptide_data())
  #   preptides <- 
  #     as.data.frame(peptide_data() %>% select(sampling()))
  #   tags <- as.vector(colnames(preptides))
  #   updateProgressBar(id = "pb4", value = 65)
  #   howIntense <- vector()
  #   
  #   for (k in 1:ncol(preptides)) {
  #     howIntense[k] <-
  #       sum(eval(parse(text = paste(
  #         "preptides$", tags[k], sep = ""
  #       ))))
  #   }
  #   
  #   updateProgressBar(id = "pb4", value = 80)
  #   
  #   joe <-
  #     as.data.frame(dplyr::bind_cols(colnames(preptides), howIntense))
  #   colnames(joe) <- c("TMT", "Intensity")
  #   
  #   
  #   updateProgressBar(id = "pb4", value = 100)
  #   
  #   return(joe)
  # })
  
  Bars <- reactive({
    req(peptide_data_unnorm())
    preptides <- 
      as.data.frame(peptide_data_unnorm() %>% select(sampling()))
    tags <- as.vector(colnames(preptides))
    # updateProgressBar(id = "pb4", value = 65)
    howIntense <- vector()
    
    for (k in 1:ncol(preptides)) {
      howIntense[k] <-
        sum(eval(parse(text = paste(
          "preptides$", tags[k], sep = ""
        ))))
    }
    
    howIntense <- colSums(preptides)
    
    # updateProgressBar(id = "pb4", value = 80)
    
    joe <-
      as.data.frame(dplyr::bind_cols(colnames(preptides), howIntense))
    colnames(joe) <- c("TMT", "Intensity")
    
    # updateProgressBar(id = "pb4", value = 100)
    print(joe)
    return(joe)
  })
  
  # loopy <- reactive({
  #   req(peptide_data())
  #   for (k in 1:samplingL()) {
  #   # Get 90th percentile intensity score/filtering for top 10 percent
  #   thresholdIntensity <- quantile(preptides[, k], probs = 0.9)
  #   sorted <- sort(preptides[, k], decreasing = TRUE)
  #   cutAndSorted <- sorted[sorted > thresholdIntensity]
  #   
  #   temporary <- paste("add_histogram(x = TMTensity()$cutAndSorted1,  name = \"TMT126\") %>%", , sep = "")
  #   }
  # })
  
  output$contents <- renderTable(digits = 3, colnames = TRUE, {
    req(digestEfficiency())
    req(labellingEfficiency())
    
    x <- c("Variable", "Value")
    unique_peptides <- if ("Peptide" %in% names(peptide_data())) {
      unique(peptide_data()$Peptide) %>% length()
    } else {
      0
    }
    
    unique_protein_ids <- if ("Protein ID" %in% names(peptide_data())) {
      unique(peptide_data() %>% 
               select("Protein ID") %>% 
               pull()) %>% length()
    } else {
      0
    }
    
    results <-
      data.frame(
        c("Digest efficiency", "Labelling efficiency", "Lysines labelled", "N-termini labelled", "Unique peptide sequences", "Unique UniProt protein IDs"),
        c(paste(efficiencies()$Digest, "%", sep = ""),
          paste(efficiencies()$Labelling %>% round(digits = 2), "%", sep = ""),
          paste(efficiencies()$KLabelling %>% round(digits = 2), "%", sep = ""),
          paste(efficiencies()$NLabelling %>% round(digits = 2), "%", sep = ""),
          unique_peptides,
          unique_protein_ids
        )
      )
    updateProgressBar(id = "pb4", value = 100)
    colnames(results) <- x
    results
  })
  
  
  
  
  # output$intensiTMT <- renderPlotly({
  #   plot_ly(alpha = 0.62,
  #           xbins = list(
  #             start = -10000,
  #             size = 9000,
  #             end = 1800000
  #           )) %>% add_histogram(x = tmtConcatAll(), name = "") %>% layout(barmode = "stack",
  #                                                                          yaxis = list(
  #                                                                            type = "linear", 
  #                                                                            zeroline = TRUE,
  #                                                                            title = "# PSM's"
  #                                                                          ),
  #                                                                          xaxis = list(
  #                                                                            zeroline = TRUE, 
  #                                                                            rangeslider = list(type = "date"),
  #                                                                            title = "Intensity"
  #                                                                          ))
  # })
  
  output$intensiTMT <- renderPlotly({
    plot_ly(alpha = 0.62,
            xbins = list(
              start = -10000,
              size = 30000,
              end = 10000000
            )) %>% add_histogram(x = tmtConcatAll(), name = "") %>% layout(barmode = "stack",
                                                                           yaxis = list(
                                                                             type = "linear", 
                                                                             zeroline = TRUE,
                                                                             title = "# PSM's"
                                                                           ),
                                                                           xaxis = list(
                                                                             zeroline = TRUE, 
                                                                             rangeslider = list(type = "date"),
                                                                             title = "Intensity"
                                                                           ))
  })
  
  # output$massHisteria <- renderPlotly({
  #   plot_ly(alpha = 0.62,
  #           xbins = list(
  #             start = -0.8,
  #             size = 0.001,
  #             end = 0.8
  #           )) %>% add_histogram(x = deltaMass(),  name = "") %>% layout(barmode = "stack",
  #                                                                        yaxis = list(
  #                                                                          type = "log",
  #                                                                          zeroline = FALSE,
  #                                                                          nticks = 5,
  #                                                                          title = "# PSM's"
  #                                                                        ),
  #                                                                        xaxis = list(
  #                                                                          title = "Delta mass (domain restriction ~0)"
  #                                                                        )
  #           )
  # })
  
  massHisteriaPlot <- reactive({
    if (!is.null(deltaMass())) {
      plot_ly(alpha = 0.62,
              xbins = list(
                start = -0.8,
                size = 0.001,
                end = 0.8
              )) %>% add_histogram(x = deltaMass(),  name = "") %>% layout(barmode = "stack",
                                                                           yaxis = list(
                                                                             type = "log",
                                                                             zeroline = FALSE,
                                                                             nticks = 5,
                                                                             title = "# PSM's"
                                                                           ),
                                                                           xaxis = list(
                                                                             title = "Delta mass (domain restriction ~0)"
                                                                           )
              )
    }
  })
  
  output$massHisteria <- renderPlotly({
    massHisteriaPlot()
  })
  
  
  # output$seeBars <- renderPlotly({
  #   plot_ly(
  #     x = Bars()$TMT,
  #     y = Bars()$Intensity / max(Bars()$Intensity),
  #     type = 'bar'
  #   ) %>%
  #     layout(
  #       title = "TMT Intensities (Percentile Rank = 90)",
  #       xaxis = list(title = ""),
  #       yaxis = list(title = "Relative intensity")
  #     )
  # })
  
  # output$seeRainbow16 <- renderPlotly({
  #   if (sum(str_detect(colnames(peptide_data()), "sample"))  == 10) {
  #     plot_ly(alpha = 0.66,
  #             xbins = list(
  #               start = 0,
  #               size = 100000,
  #               end = 3900000
  #             )) %>%
  #       add_histogram(x = TMTensity()$cutAndSorted1,  name = "TMT126") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted2,  name = "TMT127N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted3,  name = "TMT127C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted4,  name = "TMT128N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted5,  name = "TMT128C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted6,  name = "TMT129N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted7,  name = "TMT129C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted8,  name = "TMT130N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted9,  name = "TMT130C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted10,  name = "TMT131N") %>%
  #       # add_histogram(x = TMTensity()$cutAndSorted11,  name = "TMT131C") %>%
  #       # add_histogram(x = TMTensity()$cutAndSorted12,  name = "TMT132N") %>%
  #       # add_histogram(x = TMTensity()$cutAndSorted13,  name = "TMT132C") %>%
  #       # add_histogram(x = TMTensity()$cutAndSorted14,  name = "TMT133N") %>%
  #       # add_histogram(x = TMTensity()$cutAndSorted15,  name = "TMT133C") %>%
  #       # add_histogram(x = TMTensity()$cutAndSorted16,  name = "TMT134") %>%
  #       # for (k in 1:samplingL()) {
  #       #   add_histogram(x = TMTensity()$cutAndSorted1)
  #       # } %>%
  #       layout(
  #         barmode = "stack",
  #         xaxis = list(
  #           zeroline = TRUE, 
  #           title = "Summed intensities"
  #         ),
  #         yaxis = list(
  #           type = "log",
  #           zeroline = TRUE,
  #           nticks = 5,
  #           title = "# PSM's"
  #         )
  #       )
  #   } else if (sum(str_detect(colnames(peptide_data()), "sample"))  == 16) {
  #     plot_ly(alpha = 0.66,
  #             xbins = list(
  #               start = 0,
  #               size = 100000,
  #               end = 3900000
  #             )) %>%
  #       add_histogram(x = TMTensity()$cutAndSorted1,  name = "TMT126") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted2,  name = "TMT127N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted3,  name = "TMT127C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted4,  name = "TMT128N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted5,  name = "TMT128C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted6,  name = "TMT129N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted7,  name = "TMT129C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted8,  name = "TMT130N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted9,  name = "TMT130C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted10,  name = "TMT131N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted11,  name = "TMT131C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted12,  name = "TMT132N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted13,  name = "TMT132C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted14,  name = "TMT133N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted15,  name = "TMT133C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted16,  name = "TMT134") %>%
  #       # for (k in 1:samplingL()) {
  #       #   add_histogram(x = TMTensity()$cutAndSorted1)
  #       # } %>%
  #       layout(
  #         barmode = "stack",
  #         xaxis = list(
  #           zeroline = TRUE, 
  #           title = "Summed intensities"
  #         ),
  #         yaxis = list(
  #           type = "log",
  #           zeroline = TRUE,
  #           nticks = 5,
  #           title = "# PSM's"
  #         )
  #       )
  #   } else {
  #     plot_ly(alpha = 0.66,
  #             xbins = list(
  #               start = 0,
  #               size = 100000,
  #               end = 3900000
  #             )) %>%
  #       add_histogram(x = TMTensity()$cutAndSorted1,  name = "TMT126") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted2,  name = "TMT127N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted3,  name = "TMT127C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted4,  name = "TMT128N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted5,  name = "TMT128C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted6,  name = "TMT129N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted7,  name = "TMT129C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted8,  name = "TMT130N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted9,  name = "TMT130C") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted10,  name = "TMT131N") %>%
  #       add_histogram(x = TMTensity()$cutAndSorted11,  name = "TMT131C") %>%
  #       # for (k in 1:samplingL()) {
  #       #   add_histogram(x = TMTensity()$cutAndSorted1)
  #       # } %>%
  #       layout(
  #         barmode = "stack",
  #         xaxis = list(
  #           zeroline = TRUE, 
  #           title = "Summed intensities"
  #         ),
  #         yaxis = list(
  #           type = "log",
  #           zeroline = TRUE,
  #           nticks = 5,
  #           title = "# PSM's"
  #         )
  #       )
  #   }
  #   
  # } 
  # )
  
  output$seeBars <- renderPlotly({
    plot_ly(
      x = Bars()$TMT,
      y = Bars()$Intensity / max(Bars()$Intensity),
      type = 'bar'
    ) %>%
      layout(
        title = "TMT Intensities (Percentile Rank = 90)",
        xaxis = list(title = ""),
        yaxis = list(title = "Relative intensity")
      )
  })
  
  output$seeRainbow16 <- renderPlotly({
    if (sum(str_detect(colnames(peptide_data_unnorm()), "sample"))  == 10) {
      plot_ly(alpha = 0.66,
              xbins = list(
                start = 0,
                size = 100000,
                end = 3900000
              )) %>%
        add_histogram(x = TMTensity()$cutAndSorted1,  name = "TMT126") %>%
        add_histogram(x = TMTensity()$cutAndSorted2,  name = "TMT127N") %>%
        add_histogram(x = TMTensity()$cutAndSorted3,  name = "TMT127C") %>%
        add_histogram(x = TMTensity()$cutAndSorted4,  name = "TMT128N") %>%
        add_histogram(x = TMTensity()$cutAndSorted5,  name = "TMT128C") %>%
        add_histogram(x = TMTensity()$cutAndSorted6,  name = "TMT129N") %>%
        add_histogram(x = TMTensity()$cutAndSorted7,  name = "TMT129C") %>%
        add_histogram(x = TMTensity()$cutAndSorted8,  name = "TMT130N") %>%
        add_histogram(x = TMTensity()$cutAndSorted9,  name = "TMT130C") %>%
        add_histogram(x = TMTensity()$cutAndSorted10,  name = "TMT131N") %>%
        # add_histogram(x = TMTensity()$cutAndSorted11,  name = "TMT131C") %>%
        # add_histogram(x = TMTensity()$cutAndSorted12,  name = "TMT132N") %>%
        # add_histogram(x = TMTensity()$cutAndSorted13,  name = "TMT132C") %>%
        # add_histogram(x = TMTensity()$cutAndSorted14,  name = "TMT133N") %>%
        # add_histogram(x = TMTensity()$cutAndSorted15,  name = "TMT133C") %>%
        # add_histogram(x = TMTensity()$cutAndSorted16,  name = "TMT134") %>%
        # for (k in 1:samplingL()) {
        #   add_histogram(x = TMTensity()$cutAndSorted1)
        # } %>%
        layout(
          barmode = "stack",
          xaxis = list(
            zeroline = TRUE, 
            title = "Summed intensities"
          ),
          yaxis = list(
            type = "log",
            zeroline = TRUE,
            nticks = 5,
            title = "# PSM's"
          )
        )
    } else if (sum(str_detect(colnames(peptide_data_unnorm()), "sample"))  == 16) {
      plot_ly(alpha = 0.66,
              xbins = list(
                start = 0,
                size = 100000,
                end = 3900000
              )) %>%
        add_histogram(x = TMTensity()$cutAndSorted1,  name = "TMT126") %>%
        add_histogram(x = TMTensity()$cutAndSorted2,  name = "TMT127N") %>%
        add_histogram(x = TMTensity()$cutAndSorted3,  name = "TMT127C") %>%
        add_histogram(x = TMTensity()$cutAndSorted4,  name = "TMT128N") %>%
        add_histogram(x = TMTensity()$cutAndSorted5,  name = "TMT128C") %>%
        add_histogram(x = TMTensity()$cutAndSorted6,  name = "TMT129N") %>%
        add_histogram(x = TMTensity()$cutAndSorted7,  name = "TMT129C") %>%
        add_histogram(x = TMTensity()$cutAndSorted8,  name = "TMT130N") %>%
        add_histogram(x = TMTensity()$cutAndSorted9,  name = "TMT130C") %>%
        add_histogram(x = TMTensity()$cutAndSorted10,  name = "TMT131N") %>%
        add_histogram(x = TMTensity()$cutAndSorted11,  name = "TMT131C") %>%
        add_histogram(x = TMTensity()$cutAndSorted12,  name = "TMT132N") %>%
        add_histogram(x = TMTensity()$cutAndSorted13,  name = "TMT132C") %>%
        add_histogram(x = TMTensity()$cutAndSorted14,  name = "TMT133N") %>%
        add_histogram(x = TMTensity()$cutAndSorted15,  name = "TMT133C") %>%
        add_histogram(x = TMTensity()$cutAndSorted16,  name = "TMT134") %>%
        # for (k in 1:samplingL()) {
        #   add_histogram(x = TMTensity()$cutAndSorted1)
        # } %>%
        layout(
          barmode = "stack",
          xaxis = list(
            zeroline = TRUE, 
            title = "Summed intensities"
          ),
          yaxis = list(
            type = "log",
            zeroline = TRUE,
            nticks = 5,
            title = "# PSM's"
          )
        )
    } else {
      plot_ly(alpha = 0.66,
              xbins = list(
                start = 0,
                size = 100000,
                end = 3900000
              )) %>%
        add_histogram(x = TMTensity()$cutAndSorted1,  name = "TMT126") %>%
        add_histogram(x = TMTensity()$cutAndSorted2,  name = "TMT127N") %>%
        add_histogram(x = TMTensity()$cutAndSorted3,  name = "TMT127C") %>%
        add_histogram(x = TMTensity()$cutAndSorted4,  name = "TMT128N") %>%
        add_histogram(x = TMTensity()$cutAndSorted5,  name = "TMT128C") %>%
        add_histogram(x = TMTensity()$cutAndSorted6,  name = "TMT129N") %>%
        add_histogram(x = TMTensity()$cutAndSorted7,  name = "TMT129C") %>%
        add_histogram(x = TMTensity()$cutAndSorted8,  name = "TMT130N") %>%
        add_histogram(x = TMTensity()$cutAndSorted9,  name = "TMT130C") %>%
        add_histogram(x = TMTensity()$cutAndSorted10,  name = "TMT131N") %>%
        add_histogram(x = TMTensity()$cutAndSorted11,  name = "TMT131C") %>%
        # for (k in 1:samplingL()) {
        #   add_histogram(x = TMTensity()$cutAndSorted1)
        # } %>%
        layout(
          barmode = "stack",
          xaxis = list(
            zeroline = TRUE, 
            title = "Summed intensities"
          ),
          yaxis = list(
            type = "log",
            zeroline = TRUE,
            nticks = 5,
            title = "# PSM's"
          )
        )
    }
  }
  )
  
  ## make a volcano MCF-7
  output$mcf_volcano <- renderPlotly({
    n_peptide_data_statistics() %>%
      arrange(desc(MCF_significant)) %>%
      ggplot(aes(x=MCF_log_KD_ctrl,
                 y=MCF_log_p,
                 description=`Gene`,
                 color=MCF_significant)) +
      geom_point(alpha=0.7,size=1) +
      # theme_light() +
      # scale_colour_manual(values=plain_cols1) +
      theme(axis.text=element_text(size=12)) %>%
      ggplotly()
  })
  
  output$mda_volcano <- renderPlotly({
    n_peptide_data_stats() %>%
      arrange(desc(MDA_significant)) %>% 
      ggplot(aes(x=MDA_log_KD_ctrl,
                 y=MDA_log_p,
                 description=`Gene`,
                 color=MDA_significant))+
      geom_point(alpha=0.7)+
      # theme_light()+
      # scale_colour_manual(values=plain_cols1)+
      theme(axis.text=element_text(size=12)) %>%
      ggplotly()
  })
  
  findation <- reactive({
    if (input$query == '') {
      searching <- dplyr::select(protein_data_flipped(), "LDHA")
      searching <- t(searching)
      searching[,1:4]
      print(input$query)
    } else {
      searching <- dplyr::select(as.data.frame(protein_data_flipped()), as.character(input$query))
      searching <- t(searching)
      print("hi")
      print(input$query)
    }
    searching
  })
  
  
  # findation <- reactive({
  #   searching <- dplyr::select(TMTs, as.character(input$query))
  #   searching <- t(searching)
  #   searching
  # })
  
  output$foundation <- renderPlotly({
    plot_ly(
      x = colnames(findation()),
      y = as.numeric(findation()),
      type = 'bar',
      color = str_extract(colnames(findation()), "[^_]+")
      # marker = list(
      #   color = c(
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#78cdff',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#3c8dbc',
      #     '#828282',
      #     '#828282',
      #     '#828282',
      #     '#828282',
      #     '#333',
      #     '#333',
      #     '#333',
      #     '#333'
      #   )
      # )
    ) %>%
      layout(
        title = input$query,
        xaxis = list(title = ""),
        yaxis = list(title = "Intensity")
      )
  })
  
  output$protein_heatmap <- renderPlotly({
    plot_ly(
      x = colnames(scaled_pdn_df()),
      y = rownames(scaled_pdn_df()),
      z = scaled_pdn_df(),
      type = "heatmap",
      # hover
      colorscale = "Blues"
    ) %>%
      layout(
        title = "Your file, heat mapped",
        xaxis = list(title = "", side = "bottom"),
        yaxis = list(showticklabels = FALSE, ticks = "none")
      )
  })
  
}

shinyApp(ui, server)