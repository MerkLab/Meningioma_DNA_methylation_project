library(sesame)
library(sesameData)
library(caret)
library(dplyr)

###final rf model on differentially methylated probes EPIC compatible
#######repeat running the rf model, check accuracy and  check probability distributions for model that works only on differentially methylated probes
#get all significantly dif methylated probes between clusters, but only EPIC probes
probes_hyper = read.csv(file = "probes_hyper.csv", header = T)
probes_hypo = read.csv(file = "probes_hypo.csv", header = T)
probes_merge = rbind(probes_hyper, probes_hypo)
sig_probes=probes_merge$Probe_ID
sig_EPIC=mLiftOver(sig_probes, "EPIC")

#subset to 5432 probes differentially methylated and present on EPIC, and prepare for model
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

saveRDS(rf_model_sig, file = "rf_classifier_meningioma_methylation_cluster_final.rds")
#use this one for app



#####the following is the code used for the web app classifier

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

# ---- Logging helper (log.txt already exists) ----
log_run <- function(type = "uploaded", sample_id = NA) {
  # Determine app directory robustly
  app_dir <- normalizePath(".", winslash = "/", mustWork = FALSE)
  log_file <- file.path(app_dir, "log.txt")
  
  # If the file doesn't exist, create it (only if possible)
  if (!file.exists(log_file)) {
    tryCatch({
      file.create(log_file)
    }, error = function(e) {
      return(NULL)
    })
  }
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  session_id <- shiny::getDefaultReactiveDomain()$token
  
  log_line <- paste(timestamp, type, sample_id, 1, session_id, sep = "\t")
  
  tryCatch({
    cat(log_line, "\n", file = log_file, append = TRUE)
  }, error = function(e) {
    message("Logging failed: ", e$message)
  })
}

# ---- Helper function ----
process_IDAT_common_EPIC_EPICv2_totalIntensities <- function(idat_prefix) {
  epic <- sesameDataGet("EPIC.probeInfo")
  lift <- sesameDataGet("liftOver.EPICv2ToEPIC")
  lift <- lift[!(duplicated(lift$ID_target) | duplicated(lift$ID_target, fromLast = TRUE)), ]
  epic_probes <- epic$mapped.probes.hg19@ranges@NAMES
  common_epic_ids <- intersect(epic_probes, lift$ID_target)
  sset <- openSesame(idat_prefix, func = NULL)
  intensities <- as.data.frame(totalIntensities(sset))
  platform <- inferPlatformFromProbeIDs(rownames(intensities))
  if (platform == "EPICv2") {
    mapping <- lift %>% filter(ID_target %in% common_epic_ids)
    intensities <- subset(intensities, rownames(intensities) %in% mapping$ID_source)
    rownames(intensities) <- mapping$ID_target[match(rownames(intensities), mapping$ID_source)]
  } else if (platform == "EPIC") {
    intensities <- subset(intensities, rownames(intensities) %in% common_epic_ids)
  } else {
    stop("Unsupported platform: ", platform)
  }
  intensities
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Meningioma Methylation Cluster Classifier"),
  p("Welcome to our Meningioma Methylation Cluster Classifier compatible with EPIC and EPICv2 DNA methylation arrays.", style = "color:grey; margin-top:-8px;"),
  p("This web application will process your idat files and predict methylation cluster assignment based on our model (Accuracy: 98.7%, 95% CI: 96.25-99.73; Sensitivity: 94.12%; Specificity: 100%).", style = "color:grey; margin-top:-8px;"),
  p("β-value distribution of the sample is provided for QC. Genome-wide copy number events provide further prognostic information.", style = "color:grey; margin-top:-8px;"),
  p("This classifier is solely intended to be used in combination with meningioma tumor samples. We recommend the usage of the Heidelberg CNS tumor methylation classifier prior to our tool to ensure correct diagnosis.", style = "color:grey; margin-top:-8px;"),
  p("Please refer to www.hih-tuebingen.de/datenschutzerklearung for our privacy policy.", style = "color:grey; margin-top:-8px;"),
  sidebarLayout(
    sidebarPanel(
      fileInput("red_idat", "Upload Red Channel IDAT (.idat)", accept = ".idat"),
      fileInput("grn_idat", "Upload Green Channel IDAT (.idat)", accept = ".idat"),
      actionButton("run", "Run Analysis", class = "btn-primary"),
      actionButton("run_example", "Run Example Data", class = "btn-success"),
      downloadButton("download_pdf", "Download Results as PDF"),
      helpText("Upload both Red and Green IDAT files from the same sample before Run Analysis, or pick Run Example Data. Both analysis as well as pdf report might take up to 1 minute, only press once.")
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
  pred_conf <- reactiveVal(NULL)
  idat_base <- reactiveVal(NULL)
  last_run_type <- reactiveVal(NULL)
  
  # ---- Uploaded data ----
  observeEvent(input$run, {
    req(input$red_idat, input$grn_idat)
    withProgress(message = "Processing uploaded files...", value = 0, {
      
      tmp_dir <- file.path(tempdir(), "uploaded_IDATs")
      dir.create(tmp_dir, showWarnings = FALSE)
      
      red_path <- file.path(tmp_dir, input$red_idat$name)
      grn_path <- file.path(tmp_dir, input$grn_idat$name)
      
      file.copy(input$red_idat$datapath, red_path, overwrite = TRUE)
      file.copy(input$grn_idat$datapath, grn_path, overwrite = TRUE)
      
      prefix <- sub("_Red.idat$|_Grn.idat$", "", red_path)
      idat_base(basename(prefix))
      last_run_type("uploaded")
      
      # ---- LOG RUN ----
      log_run(type = "uploaded", sample_id = idat_base())
      
      incProgress(0.3, detail = "Loading IDATs and computing β-values...")
      beta <- openSesame(prefix, BPPARAM = BiocParallel::SerialParam())
      beta <- betasCollapseToPfx(beta)
      if (any(duplicated(names(beta)))) beta <- tapply(beta, names(beta), mean)
      beta_for_plot(beta)
      
      beta_aligned <- beta[model_features]
      missing_features <- setdiff(model_features, names(beta_aligned))
      if (length(missing_features) > 0) beta_aligned[missing_features] <- NA
      beta_aligned[is.na(beta_aligned)] <- median(beta_aligned, na.rm = TRUE)
      newdata <- as.data.frame(t(beta_aligned))
      
      incProgress(0.7, detail = "Predicting cluster...")
      pred <- predict(model, newdata)
      pred_prob <- predict(model, newdata, type = "prob")
      
      pred_label(group_labels_html[as.character(pred)])
      pred_conf(pred_prob[1, as.character(pred), drop = TRUE])
      
      render_results_ui(idat_base(), pred_label(), pred_conf())
      output$beta_plot <- renderPlot(make_beta_plot(beta_for_plot()))
      
      incProgress(0.9, detail = "Computing CNV...")
      output$cnv_plot <- renderPlot(make_cnv_plot(prefix))
    })
  })
  
  # ---- Example data ----
  observeEvent(input$run_example, {
    withProgress(message = "Running Example Data...", value = 0, {
      
      example_prefix <- "example_data"
      red_path <- paste0(example_prefix, "_Red.idat")
      grn_path <- paste0(example_prefix, "_Grn.idat")
      
      if (!file.exists(red_path) || !file.exists(grn_path)) {
        showNotification("Example IDAT files not found!", type = "error")
        return()
      }
      
      tmp_dir <- file.path(tempdir(), "example_IDATs")
      dir.create(tmp_dir, showWarnings = FALSE)
      
      file.copy(red_path, file.path(tmp_dir, basename(red_path)), overwrite = TRUE)
      file.copy(grn_path, file.path(tmp_dir, basename(grn_path)), overwrite = TRUE)
      
      prefix <- sub("_Red.idat$|_Grn.idat$", "", file.path(tmp_dir, basename(red_path)))
      idat_base(basename(prefix))
      last_run_type("example")
      
      beta <- openSesame(prefix, BPPARAM = BiocParallel::SerialParam())
      beta <- betasCollapseToPfx(beta)
      beta_for_plot(beta)
      
      beta_aligned <- beta[model_features]
      beta_aligned[is.na(beta_aligned)] <- median(beta_aligned, na.rm = TRUE)
      newdata <- as.data.frame(t(beta_aligned))
      
      pred <- predict(model, newdata)
      pred_prob <- predict(model, newdata, type = "prob")
      
      pred_label(group_labels_html[as.character(pred)])
      pred_conf(pred_prob[1, as.character(pred), drop = TRUE])
      
      render_results_ui(idat_base(), pred_label(), pred_conf())
      output$beta_plot <- renderPlot(make_beta_plot(beta_for_plot()))
      output$cnv_plot <- renderPlot(make_cnv_plot(prefix))
    })
  })
  
  # ---- Helper functions ----
  render_results_ui <- function(idat_name, pred_label_value, pred_conf_value) {
    output$pred_result <- renderUI({
      
      cluster_html <- pred_label_value
      cluster_color <- if (grepl("high", cluster_html)) "red" else "blue"
      
      conf_color <- if (pred_conf_value >= 0.548) {
        "#006400"
      } else if (pred_conf_value >= 0.2062) {
        "#DAA520"
      } else {
        "#B22222"
      }
      
      conf_label <- if (pred_conf_value >= 0.548) {
        "Confident"
      } else if (pred_conf_value >= 0.2062) {
        "Intermediate confidence"
      } else {
        "Unclassified"
      }
      
      HTML(paste0(
        "<div style='font-size:16px; font-weight:bold; margin-bottom:10px;'>",
        "Results report for ", idat_name, "</div>",
        
        "<div style='font-size:22px; font-weight:bold; color:", cluster_color, ";'>",
        "Predicted Meningioma Methylation Cluster: ", cluster_html,
        "<br>",
        "<span style='font-size:20px; color:", conf_color, "; font-weight:bold;'>",
        "Probability: ",
        sprintf('%.1f%% (%s)', pred_conf_value * 100, conf_label),
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
  
  make_cnv_plot <- function(prefix) {
    query <- process_IDAT_common_EPIC_EPICv2_totalIntensities(prefix)
    data.q <- CNV.load(query)
    x <- CNV.fit(data.q, data.c, anno)
    x <- CNV.bin(x)
    x <- CNV.detail(x)
    x <- CNV.segment(x)
    CNV.genomeplot(x, main = basename(prefix))
    
    grid.text(
      "log2 ratios of total intensities",
      x = unit(0.005, "npc"),
      y = unit(0.5, "npc"),
      rot = 90,
      gp = gpar(fontsize = 12, fontface = "bold")
    )
  }
  
  # ---- PDF Export (unchanged) ----
  output$download_pdf <- downloadHandler(
    filename = function() { paste0("Meningioma_Results_", Sys.Date(), ".pdf") },
    content = function(file) {
      req(beta_for_plot(), pred_label(), pred_conf(), idat_base(), last_run_type())
      
      df <- data.frame(Beta = as.numeric(beta_for_plot()))
      p_beta <- ggplot(df, aes(x = Beta)) +
        geom_density(fill = "skyblue", alpha = 0.5) +
        theme_minimal(base_size = 13) +
        labs(x = expression(beta * "-value"), y = "Density") +
        scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)), limits = c(0, 1))
      
      cluster_html <- pred_label()
      cluster_text <- gsub("<sup>", "^", gsub("</sup>", "", cluster_html))
      cluster_color <- if (cluster_text == "METH^high") "red" else "blue"
      conf_label <- if (pred_conf() >= 0.512) "Confident" else if (pred_conf() >= 0.2348) "Intermediate confidence" else "Unclassified"
      conf_color <- if (pred_conf() >= 0.512) "#006400" else if (pred_conf() >= 0.2348) "#DAA520" else "#B22222"
      
      cluster_grob <- textGrob(
        paste0("Predicted Meningioma Methylation Cluster: ", cluster_text),
        x = unit(0.07, "npc"), just = "left",
        gp = gpar(fontsize = 16, col = cluster_color, fontface = "bold")
      )
      conf_grob <- textGrob(
        paste0("Probability: ", round(pred_conf() * 100, 1), "% (", conf_label, ")"),
        x = unit(0.07, "npc"), just = "left",
        gp = gpar(fontsize = 16, col = conf_color, fontface = "bold")
      )
      pred_grob <- arrangeGrob(cluster_grob, conf_grob, ncol = 1)
      
      base_dir <- if (last_run_type() == "uploaded") "uploaded_IDATs" else "example_IDATs"
      prefix_path <- file.path(tempdir(), base_dir, idat_base())
      
      cnv_plot_file <- tempfile(fileext = ".png")
      query <- process_IDAT_common_EPIC_EPICv2_totalIntensities(prefix_path)
      data.q <- CNV.load(query)
      x <- CNV.fit(data.q, data.c, anno)
      x <- CNV.bin(x)
      x <- CNV.detail(x)
      x <- CNV.segment(x)
      png(cnv_plot_file, width = 7.2 * 1.2, height = 3.9 * 1.2, units = "in", res = 300)
      CNV.genomeplot(x, main = idat_base())
      dev.off()
      cnv_image <- rasterGrob(readPNG(cnv_plot_file), interpolate = TRUE)
      
      cnv_grob <- grobTree(
        textGrob("log2 ratios of total intensities", x = unit(0.02, "npc"), y = unit(0.5, "npc"),
                 rot = 90, gp = gpar(fontsize = 11, fontface = "bold", col = "black")),
        editGrob(cnv_image, vp = viewport(x = 0.5, y = 0.5, width = 0.94, height = 1.2))
      )
      top_title_grob <- textGrob(
        "Meningioma Methylation Cluster Classifier",
        x = 0.5, y = 1, just = "center",
        gp = gpar(fontsize = 20, fontface = "bold")
      )
      
      title_grob <- textGrob(paste0("Results report for ", idat_base()), x = unit(0.05, "npc"), just = "left",
                             gp = gpar(fontsize = 18, fontface = "bold"))
      timestamp_grob <- textGrob(paste("Generated on:", Sys.time()), x = unit(0.95, "npc"), just = "right",
                                 gp = gpar(fontsize = 10, col = "grey40", fontface = "italic"))
      
      pdf(file, width = 8.5, height = 11)
      grid.arrange(
        top_title_grob,title_grob, pred_grob, p_beta, cnv_grob, timestamp_grob,
        ncol = 1,
        heights = unit(c(0.5,0.7, 1.0, 2.0, 5.0, 1.0), "in")
      )
      dev.off()
    }
  )
}

# ---- Run App ----
shinyApp(ui, server)


