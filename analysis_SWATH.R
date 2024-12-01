# Load necessary libraries
library(gridExtra)
library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(limma)

# Read data
rawswath <- read.xlsx("path/to/your/SWATH-3groups-data.xlsx")

# Data preprocessing
swath_dt <- rawswath %>%
  dplyr::mutate(tosplit = Replicate, todes = Protein.Description) %>%
  tidyr::separate(tosplit, sep = "_", into = c("trtsplit", "sw", "ba", "01", "no" )) %>%
  tidyr::separate(todes, sep = " OS=", into = c("Protein.Name", "toos")) %>%
  tidyr::separate(toos, sep = " OX=", into = c("Organism", "toox")) %>%
  tidyr::separate(toox, sep = " GN=", into = c("OX", "togn")) %>%
  tidyr::separate(togn, sep = " PE=", into = c("GN", "tope")) %>%
  tidyr::separate(tope, sep = " SV=", into = c("PE", "SV")) %>%
  dplyr::mutate(Protein = gsub("_.*", "", Protein.Name)) %>%
  dplyr::mutate(Protein.IDs = gsub(".*\\|", "", Protein)) %>%
  tidyr::separate(trtsplit, sep = "-", into = c("treatments", "P", "subrep"))

# Data summarization
swath_dt_summary <- swath_dt %>%
  dplyr::group_by(treatments, P, Protein.IDs) %>%
  dplyr::summarise(total_area_sums = sum(Total.Area, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(treatments, P) %>%
  dplyr::mutate(sumeachtrt = sum(total_area_sums, na.rm = T),
                total_area_sums_divided = total_area_sums/sum(total_area_sums, na.rm = T)) %>%
  dplyr::ungroup() %>%
  mutate(trt = gsub("\\.N.*", "", treatments),
         total_area_sums_log2 = log2(total_area_sums)) %>%
  dplyr::mutate(trt_p = paste(trt, "_", P, sep = ""))

# Normalization
intensity_data <- swath_dt_summary %>%
  dplyr::select(trt_p, Protein.IDs, total_area_sums_log2) %>%
  tidyr::pivot_wider(names_from = trt_p, values_from = total_area_sums_log2)

intensity_data <- as.data.frame(intensity_data)
rownames(intensity_data) <- intensity_data$Protein.IDs
intensity_data <- intensity_data[,-1]

loess_normalized_data <- normalizeCyclicLoess(as.matrix(intensity_data))
loess_normalized_df <- as.data.frame(loess_normalized_data)
loess_normalized_df$Protein.IDs <- rownames(intensity_data)

loess_normalized_df_long <- loess_normalized_df %>%
  pivot_longer(Control_P18:Dimethoate_P20, names_to = "treatments", values_to = "Normalized.Log2.Total.Area")

swath_dt_summary_norm <- merge(swath_dt_summary, loess_normalized_df_long, 
                               by.x = c("Protein.IDs", "trt_p"), 
                               by.y = c("Protein.IDs", "treatments"))

# Prepare data for analysis
swath_dt <- swath_dt_summary_norm %>%
  mutate(Protein.IDs = factor(Protein.IDs),
         trt = factor(treatments, levels = c("Control", "Dimethoate")))

control_group <- "Control"
trt_groups <- unique(swath_dt$trt)
trt_groups <- trt_groups[trt_groups != control_group]

# Proteins of interest
prot_oi <- c("PRKAB1", "RPTOR", "MTOR", "SIRT1", "BECN1", "RHEB")

# Analysis for Dimethoate
trt_id <- "Dimethoate"
data_subset <- subset(swath_dt, trt %in% c(control_group, trt_id))

# Linear model
lm_results <- lm(Normalized.Log2.Total.Area ~ Protein.IDs*trt, data = data_subset)
summary_results <- summary(lm_results)

# ANOVA
aov_lm <- aov(Normalized.Log2.Total.Area ~ Protein.IDs * trt, data = data_subset)
emmeans(aov_lm, "trt", by = "Protein.IDs", pairwise = TRUE, adjust = "tukey")

# Plotting
p_list <- list()

for (protein in prot_oi) {
  subset_data <- data_subset[(data_subset$GN %in% protein),]
  
  model <- lm(Normalized.Log2.Total.Area ~ trt, data = subset_data)
  summary_model <- summary(model)
  coef_summary <- summary_model$coefficients
  
  p_value <- coef_summary[grepl("trtDimethoate", rownames(coef_summary)),]$`Pr(>|t|)`
  
  significance_level <- 0.05
  letters <- if (p_value < significance_level) c("a", "b") else c("a", "a")
  names(letters) <- unique(subset_data$trt)
  
  p <- ggplot(subset_data, aes(x = trt, y = Normalized.Log2.Total.Area)) +
    geom_boxplot(aes(color = trt), lwd = 1, size = 10) +
    geom_jitter(aes(color = trt), width = 0.2, size = 3) +
    stat_summary(aes(color = trt), fun = mean, geom = "point", shape = 18, size = 5) +
    annotate("text", x = names(letters), y = max(subset_data$Normalized.Log2.Total.Area) + 0.05, 
             label = letters, vjust = -0.5, color = "black") +
    labs(title = paste(protein), x = "Treatments", y = "Normalized Log2 Total Area Sum", color = "Treatment") +
    scale_color_manual(values = c("Control" = "#619CFF", "Dimethoate" = "#F8766D")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.position = "none")
  
  p_list[[protein]] <- p
}

# Combine plots
combined_plot <- grid.arrange(grobs = p_list, ncol = 3, nrow = 2)

# Save combined plot
ggsave(plot = combined_plot, "path/to/your/output/six_prot_dimethoate.png", 
       height = 6, width = 12, dpi = 400, unit = "in")
