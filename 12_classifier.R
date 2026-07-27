library(sesame)
library(sesameData)
library(caret)
library(dplyr)
library(here)

###clone figshare repository including the given folder structure to your local directory
###populate raw data subfolders with data from either GEO or figshare itself
###refer to the README.txt file at https://github.com/MerkLab/Meningioma_DNA_methylation_project for further information

####the first part of the code serves to generate a random forest model to predict METHlow and METHhigh meningiomas based on significantly differentially methylated probes
###this model is based on samples from the discovery cohort
###data-driven thresholds for prediction confidence are calculated
###random forest model is saved as an .RDS file
###beginning with line 122, all code and corresponding R packages used to generate the webapp classifier are shown

###final rf model on differentially methylated probes EPIC compatible
#######repeat running the rf model, check accuracy and  check probability distributions for model that works only on differentially methylated probes
#get all significantly dif methylated probes between clusters, but only EPIC probes
setwd(here("data", "12_Classifier"))
probes_hyper = read.csv(file = "probes_hyper.csv", header = T)
probes_hypo = read.csv(file = "probes_hypo.csv", header = T)
probes_merge = rbind(probes_hyper, probes_hypo)
sig_probes=probes_merge$Probe_ID
sig_EPIC=mLiftOver(sig_probes, "EPIC")

#subset to 5432 probes differentially methylated and present on EPIC, and prepare for model
setwd(here("data", "general_datasets"))
betas.disc=readRDS(file = "Betas_discovery_preprocessed_corrected_EPIC_match.RDS")
disc.sig = as.data.frame(betas.disc)
disc.sig = disc.sig[rownames(disc.sig) %in% sig_EPIC,]

targets_disc= read.csv(file="meta_discovery_current.csv", header = T)
rownames(targets_disc) = targets_disc$X
all(colnames(disc.sig) == rownames(targets_disc))

dataset <- t(disc.sig)
dataset <- as.data.frame(dataset)
all(rownames(dataset) == rownames(targets_disc))

condition=factor(targets_disc$cluster_1000_new_2)
condition <- factor(condition, levels = c("2", "1"), labels = c("METH_high", "METH_low"))
dataset <- cbind(dataset,condition)


#train model with sig probes
ctrl <- trainControl(
  method = "cv",                 # cross-validation
  number = 10,                    # number of folds
  classProbs = TRUE,             # required for ROC
  summaryFunction = twoClassSummary,
  savePredictions = "final" 
)

#retrain the model
rf_model_sig <- train(condition~., data=dataset, method="rf", metric="ROC", trControl=ctrl)

rf_model_sig$finalModel  # contains the randomForest object
rf_model_sig$levels 

oo_pred_sig <- rf_model_sig$pred %>%
  # keep only predictions corresponding to best mtry
  filter(mtry == rf_model_sig$bestTune$mtry) %>%
  arrange(rowIndex)

# 1️⃣ Inspect distribution
summary(oo_pred_sig$METH_high)
ggplot(oo_pred_sig, aes(x = METH_high, fill = obs)) +
  geom_density(alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(title = "OOF Probability Distribution", x = "Probability METH_high", y = "Density")

# 2️⃣ Define data-driven confidence thresholds
# For example, take 10th percentile of METH_high probabilities for true positives
# as lower bound of confident METH_high
high_conf_thresh_sig <- quantile(oo_pred_sig$METH_high[oo_pred_sig$obs=="METH_high"], 0.1)
# Take 90th percentile of METH_high probabilities for true negatives as upper bound of confident METH_low
low_conf_thresh_sig <- quantile(oo_pred_sig$METH_high[oo_pred_sig$obs=="METH_low"], 0.9)

# Intermediate zone is between these thresholds
low_conf_thresh_sig
high_conf_thresh_sig

# 3️⃣ Assign final calls based on thresholds
oo_pred_sig <- oo_pred_sig %>%
  mutate(final_call = case_when(
    METH_high >= high_conf_thresh_sig        ~ "Confident METH_high",
    METH_high <= low_conf_thresh_sig         ~ "Confident METH_low",
    TRUE                                 ~ "Intermediate"
  ))


# 4️⃣ Inspect results
table(oo_pred_sig$final_call, oo_pred_sig$obs)

# 5️⃣ Optional: visualize thresholds on probability distribution
ggplot(oo_pred_sig, aes(x = METH_high, fill = obs)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = c(low_conf_thresh_sig, high_conf_thresh_sig), linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(title = "OOF Probability Distribution with Confidence Thresholds",
       x = "Probability METH_high", y = "Density")


class(oo_pred_sig$obs)
str(oo_pred_sig$obs)
str(oo_pred_sig$pred)

obs_sig = as.numeric(oo_pred_sig$obs)
pred_sig = as.numeric(oo_pred_sig$pred)

xtab = table(pred_sig, obs_sig)
cm_sig = caret::confusionMatrix(xtab)
print(cm_sig)

setwd(here("results"))
saveRDS(rf_model_sig, file = "rf_classifier_meningioma_methylation_cluster_final.rds")
#use this one for app



#####the following is the code used for the web app classifier
###this code is hosted on the Hertie Institute server

library(shiny)
library(caret)
library(sesame)
library(sesameData)
library(conumee)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(png)
library(shinycssloaders)
library(randomForest)

options(shiny.maxRequestSize = 100 * 1024^2)

# =====================================================================
# Configuration
# =====================================================================
MIN_IDAT_FILE_SIZE <- 1 * 1024^2    # 1 MB
MAX_IDAT_FILE_SIZE <- 50 * 1024^2   # 50 MB per channel
MIN_FEATURE_MATCH_WARN <- 0.90      # Warn below this; still run (matches original imputation behavior)
MIN_FEATURE_MATCH_HARD <- 0.50      # Hard-fail only when most features are missing

CONF_THRESH_HIGH <- 0.548           # Confident
CONF_THRESH_LOW  <- 0.2062            # Intermediate confidence (below = Unclassified)

# =====================================================================
# Upload validation
# =====================================================================
validate_idat_uploads <- function(red_file, grn_file) {
  if (is.null(red_file) || is.null(grn_file)) {
    stop("Please upload both Red and Green IDAT files.")
  }
  
  if (!grepl("_Red\\.idat$", red_file$name)) {
    stop("The Red channel filename must end exactly with '_Red.idat'.")
  }
  
  if (!grepl("_Grn\\.idat$", grn_file$name)) {
    stop("The Green channel filename must end exactly with '_Grn.idat'.")
  }
  
  red_sample <- sub("_Red\\.idat$", "", basename(red_file$name))
  grn_sample <- sub("_Grn\\.idat$", "", basename(grn_file$name))
  
  if (!identical(red_sample, grn_sample)) {
    stop("The Red and Green IDAT files do not belong to the same sample.")
  }
  
  file_sizes <- c(red_file$size, grn_file$size)
  
  if (any(is.na(file_sizes)) || any(file_sizes < MIN_IDAT_FILE_SIZE)) {
    stop("One or both uploaded files are too small to be valid IDAT files.")
  }
  
  if (any(file_sizes > MAX_IDAT_FILE_SIZE)) {
    stop("One or both uploaded files exceed the allowed size of 50 MB per file.")
  }
  
  invisible(red_sample)
}

validate_idat_content <- function(idat_prefix) {
  sset <- tryCatch(
    openSesame(idat_prefix, func = NULL),
    error = function(e) {
      stop(
        "The uploaded files could not be read as a valid IDAT pair.",
        call. = FALSE
      )
    }
  )
  
  intensities <- tryCatch(
    as.data.frame(totalIntensities(sset)),
    error = function(e) {
      stop(
        "The uploaded files do not contain readable IDAT intensity data.",
        call. = FALSE
      )
    }
  )
  
  probe_ids <- rownames(intensities)
  
  if (is.null(probe_ids) || length(probe_ids) == 0L) {
    stop("No probe identifiers were found in the uploaded files.", call. = FALSE)
  }
  
  platform <- tryCatch(
    inferPlatformFromProbeIDs(probe_ids),
    error = function(e) NA_character_
  )
  
  if (length(platform) != 1L || is.na(platform) ||
      !platform %in% c("EPIC", "EPICv2")) {
    stop(
      "Only valid Illumina EPIC or EPICv2 IDAT files are supported.",
      call. = FALSE
    )
  }
  
  invisible(platform)
}

# ---- Load caret model ----
model <- readRDS("rf_classifier_meningioma_methylation_cluster_final.rds")

# ---- Extract feature names used in training ----
model_features <- model$finalModel$xNames

data.c <- readRDS("data.c.rds")
anno <- readRDS("anno.rds")

# ---- Label mapping ----
group_labels_html <- c(
  "METH_high" = "METH<sup>high</sup>",
  "METH_low"  = "METH<sup>low</sup>"
)
group_labels_text <- c(
  "METH_high" = "METH^high",
  "METH_low"  = "METH^low"
)

# ---- Logging helper ----
log_run <- function(type = "uploaded", sample_id = NA, extra = NA) {
  app_dir <- normalizePath(".", winslash = "/", mustWork = FALSE)
  log_file <- file.path(app_dir, "log.txt")
  
  if (!file.exists(log_file)) {
    tryCatch({
      file.create(log_file)
    }, error = function(e) {
      return(NULL)
    })
  }
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  session_id <- shiny::getDefaultReactiveDomain()$token
  
  log_line <- paste(
    timestamp, type, sample_id, 1, session_id, extra,
    sep = "\t"
  )
  
  tryCatch({
    cat(log_line, "\n", file = log_file, append = TRUE)
  }, error = function(e) {
    message("Logging failed: ", e$message)
  })
}

# ---- Shared helpers ----
collapse_duplicate_probe_values <- function(x) {
  if (!any(duplicated(names(x)))) {
    return(x)
  }
  tapply(x, names(x), mean)
}

get_confidence_display <- function(prob) {
  if (prob >= CONF_THRESH_HIGH) {
    list(label = "Confident", color = "#006400")
  } else if (prob >= CONF_THRESH_LOW) {
    list(label = "Intermediate confidence", color = "#DAA520")
  } else {
    list(label = "Unclassified", color = "#B22222")
  }
}

deduplicate_intensity_rows <- function(intensities) {
  if (!any(duplicated(rownames(intensities)))) {
    return(intensities)
  }
  
  probe_ids <- rownames(intensities)
  numeric_cols <- vapply(intensities, is.numeric, logical(1L))
  mat <- as.matrix(intensities[, numeric_cols, drop = FALSE])
  
  aggregated <- rowsum(mat, probe_ids, reorder = FALSE)
  aggregated <- aggregated / as.vector(table(probe_ids)[rownames(aggregated)])
  
  result <- as.data.frame(aggregated)
  rownames(result) <- rownames(aggregated)
  result
}

harmonize_betas_to_EPIC <- function(beta) {
  platform <- tryCatch(
    inferPlatformFromProbeIDs(names(beta)),
    error = function(e) NA_character_
  )
  
  if (is.na(platform) || length(platform) != 1L) {
    stop("Could not determine methylation array platform from probe IDs.", call. = FALSE)
  }
  
  if (platform == "EPIC") {
    return(collapse_duplicate_probe_values(beta))
  }
  
  if (platform == "EPICv2") {
    lift <- sesameDataGet("liftOver.EPICv2ToEPIC")
    lift <- lift[
      !(duplicated(lift$ID_target) | duplicated(lift$ID_target, fromLast = TRUE)),
    ]
    
    mapping <- lift
    beta <- beta[names(beta) %in% mapping$ID_source]
    names(beta) <- mapping$ID_target[match(names(beta), mapping$ID_source)]
    return(collapse_duplicate_probe_values(beta))
  }
  
  stop("Unsupported platform: ", platform)
}

process_IDAT_common_EPIC_EPICv2_totalIntensities <- function(idat_prefix) {
  epic <- sesameDataGet("EPIC.probeInfo")
  lift <- sesameDataGet("liftOver.EPICv2ToEPIC")
  lift <- lift[
    !(duplicated(lift$ID_target) | duplicated(lift$ID_target, fromLast = TRUE)),
  ]
  epic_probes <- epic$mapped.probes.hg19@ranges@NAMES
  common_epic_ids <- intersect(epic_probes, lift$ID_target)
  
  sset <- openSesame(idat_prefix, func = NULL)
  intensities <- as.data.frame(totalIntensities(sset))
  platform <- inferPlatformFromProbeIDs(rownames(intensities))
  
  if (platform == "EPICv2") {
    mapping <- lift %>% filter(ID_target %in% common_epic_ids)
    intensities <- subset(intensities, rownames(intensities) %in% mapping$ID_source)
    rownames(intensities) <- mapping$ID_target[
      match(rownames(intensities), mapping$ID_source)
    ]
    intensities <- deduplicate_intensity_rows(intensities)
  } else if (platform == "EPIC") {
    intensities <- subset(intensities, rownames(intensities) %in% common_epic_ids)
  } else {
    stop("Unsupported platform: ", platform)
  }
  
  intensities
}

align_betas_for_prediction <- function(beta) {
  beta <- harmonize_betas_to_EPIC(beta)
  
  matched <- model_features %in% names(beta)
  match_frac <- mean(matched)
  
  if (match_frac < MIN_FEATURE_MATCH_HARD) {
    stop(
      sprintf(
        "Only %.1f%% of model features matched after platform harmonization (minimum %.0f%% required).",
        100 * match_frac,
        100 * MIN_FEATURE_MATCH_HARD
      ),
      call. = FALSE
    )
  }
  
  beta_aligned <- beta[model_features]
  missing_features <- model_features[!matched]
  
  if (length(missing_features) > 0L) {
    beta_aligned[missing_features] <- NA
  }
  
  beta_aligned[is.na(beta_aligned)] <- median(beta_aligned, na.rm = TRUE)
  
  list(
    newdata = as.data.frame(t(beta_aligned)),
    feature_match_pct = 100 * match_frac,
    low_feature_match = match_frac < MIN_FEATURE_MATCH_WARN
  )
}

run_classification <- function(prefix) {
  beta <- openSesame(prefix, BPPARAM = BiocParallel::SerialParam())
  beta <- betasCollapseToPfx(beta)
  beta <- collapse_duplicate_probe_values(beta)
  beta_for_plot <- beta
  
  aligned <- align_betas_for_prediction(beta)
  
  pred <- predict(model, aligned$newdata)
  pred_prob <- predict(model, aligned$newdata, type = "prob")
  pred_class <- as.character(pred)
  
  list(
    beta_for_plot = beta_for_plot,
    pred_label = unname(group_labels_html[pred_class]),
    pred_label_text = unname(group_labels_text[pred_class]),
    pred_conf = pred_prob[1, pred_class, drop = TRUE],
    feature_match_pct = aligned$feature_match_pct,
    low_feature_match = aligned$low_feature_match
  )
}

compute_cnv_segments <- function(prefix) {
  query <- process_IDAT_common_EPIC_EPICv2_totalIntensities(prefix)
  data.q <- CNV.load(query)
  x <- CNV.fit(data.q, data.c, anno)
  x <- CNV.bin(x)
  x <- CNV.detail(x)
  CNV.segment(x)
}

draw_cnv_plot <- function(x, main_title) {
  CNV.genomeplot(x, main = main_title)
  grid.text(
    "log2 ratios of total intensities",
    x = unit(0.005, "npc"),
    y = unit(0.5, "npc"),
    rot = 90,
    gp = gpar(fontsize = 12, fontface = "bold")
  )
}

clear_results <- function(
    beta_for_plot,
    pred_label,
    pred_conf,
    pred_label_text,
    feature_match_pct,
    idat_base,
    last_run_type,
    current_prefix,
    output
) {
  beta_for_plot(NULL)
  pred_label(NULL)
  pred_conf(NULL)
  pred_label_text(NULL)
  feature_match_pct(NULL)
  idat_base(NULL)
  last_run_type(NULL)
  current_prefix(NULL)
  
  output$pred_result <- renderUI(NULL)
  output$beta_plot <- renderPlot(NULL)
  output$cnv_plot <- renderPlot(NULL)
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Meningioma Methylation Cluster Classifier"),
  p(
    "Welcome to our Meningioma Methylation Cluster Classifier compatible with EPIC and EPICv2 DNA methylation arrays.",
    style = "color:grey; margin-top:-8px;"
  ),
  p(
    "This web application will process your idat files and predict methylation cluster assignment based on our model (Accuracy: 98.7%, 95% CI: 96.25-99.73; Sensitivity: 94.12%; Specificity: 100%).",
    style = "color:grey; margin-top:-8px;"
  ),
  p(
    "β-value distribution of the sample is provided for QC. Genome-wide copy number events provide further prognostic information.",
    style = "color:grey; margin-top:-8px;"
  ),
  p(
    "This classifier is solely intended to be used in combination with meningioma tumor samples. We recommend the usage of the Heidelberg CNS tumor methylation classifier prior to our tool to ensure correct diagnosis.",
    style = "color:grey; margin-top:-8px;"
  ),
  p(
    "Please refer to www.hih-tuebingen.de/datenschutzerklaerung for our privacy policy.",
    style = "color:grey; margin-top:-8px;"
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("red_idat", "Upload Red Channel IDAT (.idat)", accept = ".idat"),
      fileInput("grn_idat", "Upload Green Channel IDAT (.idat)", accept = ".idat"),
      actionButton("run", "Run Analysis", class = "btn-primary"),
      actionButton("run_example", "Run Example Data", class = "btn-success"),
      downloadButton("download_pdf", "Download Results as PDF"),
      helpText(
        "Upload both Red and Green IDAT files from the same sample before Run Analysis, or pick Run Example Data. Both analysis as well as pdf report might take up to 1 minute, only press once."
      )
    ),
    mainPanel(
      uiOutput("pred_result") %>% withSpinner(color = "#007BFF"),
      h4("β-value Density Plot"),
      plotOutput("beta_plot", height = "150px") %>% withSpinner(color = "#28a745"),
      h4("CNV Genome Plot"),
      plotOutput("cnv_plot", height = "375px") %>% withSpinner(color = "#FF5733")
    )
  )
)

# ---- SERVER ----
server <- function(input, output, session) {
  
  beta_for_plot <- reactiveVal(NULL)
  pred_label <- reactiveVal(NULL)
  pred_label_text <- reactiveVal(NULL)
  pred_conf <- reactiveVal(NULL)
  feature_match_pct <- reactiveVal(NULL)
  idat_base <- reactiveVal(NULL)
  last_run_type <- reactiveVal(NULL)
  current_prefix <- reactiveVal(NULL)
  
  render_results_ui <- function(
    idat_name,
    pred_label_value,
    pred_conf_value,
    feature_match_value
  ) {
    conf <- get_confidence_display(pred_conf_value)
    cluster_color <- if (grepl("high", pred_label_value)) "red" else "blue"
    
    output$pred_result <- renderUI({
      HTML(paste0(
        "<div style='font-size:16px; font-weight:bold; margin-bottom:10px;'>",
        "Results report for ", idat_name, "</div>",
        "<div style='font-size:22px; font-weight:bold; color:", cluster_color, ";'>",
        "Predicted Meningioma Methylation Cluster: ", pred_label_value,
        "<br>",
        "<span style='font-size:20px; color:", conf$color, "; font-weight:bold;'>",
        "Probability: ",
        sprintf("%.1f%% (%s)", pred_conf_value * 100, conf$label),
        "</span>",
        "<br>",
        "<span style='font-size:16px; color:#444; font-weight:normal;'>",
        sprintf("Model features matched: %.1f%%", feature_match_value),
        "</span>",
        "</div>"
      ))
    })
  }
  
  make_beta_plot <- function(beta_values) {
    df <- data.frame(Beta = as.numeric(beta_values))
    ggplot(df, aes(x = Beta)) +
      geom_density(fill = "skyblue", alpha = 0.5) +
      theme_minimal(base_size = 14) +
      labs(x = expression(beta * "-value"), y = "Density") +
      xlim(0, 1)
  }
  
  run_analysis_pipeline <- function(prefix, sample_name, run_type) {
    withProgress(message = paste0("Processing ", run_type, " data..."), value = 0, {
      incProgress(0.15, detail = "Checking IDAT content and platform...")
      validate_idat_content(prefix)
      
      idat_base(sample_name)
      last_run_type(run_type)
      current_prefix(prefix)
      
      if (run_type == "uploaded") {
        log_run(type = "uploaded", sample_id = sample_name)
      } else {
        log_run(type = "example", sample_id = sample_name)
      }
      
      incProgress(0.3, detail = "Loading IDATs and computing β-values...")
      cls <- run_classification(prefix)
      
      beta_for_plot(cls$beta_for_plot)
      pred_label(cls$pred_label)
      pred_label_text(cls$pred_label_text)
      pred_conf(cls$pred_conf)
      feature_match_pct(cls$feature_match_pct)
      
      if (isTRUE(cls$low_feature_match)) {
        showNotification(
          sprintf(
            "Warning: only %.1f%% of model features matched. Results may be unreliable.",
            cls$feature_match_pct
          ),
          type = "warning",
          duration = 12
        )
      }
      
      render_results_ui(
        idat_base(),
        pred_label(),
        pred_conf(),
        feature_match_pct()
      )
      output$beta_plot <- renderPlot(make_beta_plot(beta_for_plot()))
      
      incProgress(0.9, detail = "Computing CNV...")
      cnv_x <- compute_cnv_segments(prefix)
      output$cnv_plot <- renderPlot({
        draw_cnv_plot(cnv_x, basename(prefix))
      })
    })
  }
  
  observeEvent(input$run, {
    req(input$red_idat, input$grn_idat)
    
    tryCatch({
      withProgress(message = "Processing uploaded files...", value = 0, {
        incProgress(0.05, detail = "Validating filenames and file sizes...")
        
        sample_name <- validate_idat_uploads(
          input$red_idat,
          input$grn_idat
        )
        
        safe_session_id <- gsub("[^A-Za-z0-9_-]", "_", session$token)
        tmp_dir <- file.path(tempdir(), "uploaded_IDATs", safe_session_id)
        dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
        
        red_name <- basename(input$red_idat$name)
        grn_name <- basename(input$grn_idat$name)
        
        red_path <- file.path(tmp_dir, red_name)
        grn_path <- file.path(tmp_dir, grn_name)
        
        red_copied <- file.copy(
          input$red_idat$datapath,
          red_path,
          overwrite = TRUE
        )
        grn_copied <- file.copy(
          input$grn_idat$datapath,
          grn_path,
          overwrite = TRUE
        )
        
        if (!isTRUE(red_copied) || !isTRUE(grn_copied)) {
          stop("The uploaded files could not be copied for processing.")
        }
        
        prefix <- file.path(tmp_dir, sample_name)
      })
      
      run_analysis_pipeline(prefix, sample_name, "uploaded")
    }, error = function(e) {
      log_run(type = "error", sample_id = NA, extra = conditionMessage(e))
      
      clear_results(
        beta_for_plot,
        pred_label,
        pred_conf,
        pred_label_text,
        feature_match_pct,
        idat_base,
        last_run_type,
        current_prefix,
        output
      )
      
      showNotification(
        paste("Analysis failed:", conditionMessage(e)),
        type = "error",
        duration = 12
      )
    })
  })
  
  observeEvent(input$run_example, {
    tryCatch({
      example_prefix <- "example_data"
      app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
      red_path <- file.path(app_dir, paste0(example_prefix, "_Red.idat"))
      grn_path <- file.path(app_dir, paste0(example_prefix, "_Grn.idat"))
      
      if (!file.exists(red_path) || !file.exists(grn_path)) {
        stop(
          "Example IDAT files not found in the app directory: ",
          app_dir,
          call. = FALSE
        )
      }
      
      safe_session_id <- gsub("[^A-Za-z0-9_-]", "_", session$token)
      tmp_dir <- file.path(tempdir(), "example_IDATs", safe_session_id)
      dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
      
      red_copied <- file.copy(
        red_path,
        file.path(tmp_dir, basename(red_path)),
        overwrite = TRUE
      )
      grn_copied <- file.copy(
        grn_path,
        file.path(tmp_dir, basename(grn_path)),
        overwrite = TRUE
      )
      
      if (!isTRUE(red_copied) || !isTRUE(grn_copied)) {
        stop("Example IDAT files could not be copied for processing.")
      }
      
      prefix <- file.path(tmp_dir, example_prefix)
      run_analysis_pipeline(prefix, example_prefix, "example")
    }, error = function(e) {
      log_run(type = "error", sample_id = "example", extra = conditionMessage(e))
      
      clear_results(
        beta_for_plot,
        pred_label,
        pred_conf,
        pred_label_text,
        feature_match_pct,
        idat_base,
        last_run_type,
        current_prefix,
        output
      )
      
      showNotification(
        paste("Example analysis failed:", conditionMessage(e)),
        type = "error",
        duration = 12
      )
    })
  })
  
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste0("Meningioma_Results_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(
        beta_for_plot(),
        pred_label(),
        pred_label_text(),
        pred_conf(),
        feature_match_pct(),
        idat_base(),
        last_run_type(),
        current_prefix()
      )
      
      df <- data.frame(Beta = as.numeric(beta_for_plot()))
      p_beta <- ggplot(df, aes(x = Beta)) +
        geom_density(fill = "skyblue", alpha = 0.5) +
        theme_minimal(base_size = 13) +
        labs(x = expression(beta * "-value"), y = "Density") +
        scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)), limits = c(0, 1))
      
      conf <- get_confidence_display(pred_conf())
      cluster_color <- if (pred_label_text() == "METH^high") "red" else "blue"
      
      cluster_grob <- textGrob(
        paste0(
          "Predicted Meningioma Methylation Cluster: ",
          pred_label_text()
        ),
        x = unit(0.07, "npc"),
        just = "left",
        gp = gpar(fontsize = 16, col = cluster_color, fontface = "bold")
      )
      conf_grob <- textGrob(
        paste0(
          "Probability: ",
          round(pred_conf() * 100, 1),
          "% (",
          conf$label,
          ")"
        ),
        x = unit(0.07, "npc"),
        just = "left",
        gp = gpar(fontsize = 16, col = conf$color, fontface = "bold")
      )
      match_grob <- textGrob(
        sprintf("Model features matched: %.1f%%", feature_match_pct()),
        x = unit(0.07, "npc"),
        just = "left",
        gp = gpar(fontsize = 14, col = "#444444")
      )
      pred_grob <- arrangeGrob(cluster_grob, conf_grob, match_grob, ncol = 1)
      
      prefix_path <- current_prefix()
      cnv_x <- compute_cnv_segments(prefix_path)
      
      cnv_plot_file <- tempfile(fileext = ".png")
      png(
        cnv_plot_file,
        width = 7.2 * 1.2,
        height = 3.9 * 1.2,
        units = "in",
        res = 300
      )
      draw_cnv_plot(cnv_x, idat_base())
      dev.off()
      
      cnv_image <- rasterGrob(readPNG(cnv_plot_file), interpolate = TRUE)
      cnv_grob <- editGrob(
        cnv_image,
        vp = viewport(x = 0.5, y = 0.5, width = 0.94, height = 1.0)
      )
      
      top_title_grob <- textGrob(
        "Meningioma Methylation Cluster Classifier",
        x = 0.5,
        y = 1,
        just = "center",
        gp = gpar(fontsize = 20, fontface = "bold")
      )
      
      title_grob <- textGrob(
        paste0("Results report for ", idat_base()),
        x = unit(0.05, "npc"),
        just = "left",
        gp = gpar(fontsize = 18, fontface = "bold")
      )
      timestamp_grob <- textGrob(
        paste("Generated on:", Sys.time()),
        x = unit(0.95, "npc"),
        just = "right",
        gp = gpar(fontsize = 10, col = "grey40", fontface = "italic")
      )
      
      pdf(file, width = 8.5, height = 11)
      grid.arrange(
        top_title_grob,
        title_grob,
        pred_grob,
        p_beta,
        cnv_grob,
        timestamp_grob,
        ncol = 1,
        heights = unit(c(0.5, 0.7, 1.2, 2.0, 5.0, 1.0), "in")
      )
      dev.off()
    }
  )
  
  session$onSessionEnded(function() {
    safe_session_id <- gsub("[^A-Za-z0-9_-]", "_", session$token)
    session_upload_dir <- file.path(
      tempdir(),
      "uploaded_IDATs",
      safe_session_id
    )
    
    if (dir.exists(session_upload_dir)) {
      unlink(session_upload_dir, recursive = TRUE, force = TRUE)
    }
    
    example_upload_dir <- file.path(tempdir(), "example_IDATs", safe_session_id)
    if (dir.exists(example_upload_dir)) {
      unlink(example_upload_dir, recursive = TRUE, force = TRUE)
    }
  })
}

# ---- Run App ----
shinyApp(ui, server)


