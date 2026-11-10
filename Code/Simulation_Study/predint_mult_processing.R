library(tidyverse)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(here)
library(ggsci)
library(wesanderson)

# set working directory to Code and Data
# setwd("/Code_and_Data")


# Explanation of variables ------------------------------------------------

# This script processes simulation results from .rds files. It calculates coverage
# probabilities, equi-tailed probabilities, and confidence interval widths for
# various statistical methods. The final output is a set of long-format
# ("tidy") data frames suitable for analysis and plotting.

# cerrnas: count error: "NAs found in the working weights variable 'wz'"
# cerrfin: count error: "Some elements in the working weights variable 'wz' are not finite"
# cwarwz: count warning: "diagonal elements of the working weights variable 'wz' have been replaced by 1.819e-12"
# cwarcon: count warning: "convergence not obtained in 30 IRLS iterations"
# cwarany: count warning: any other warning
# cwarzc: count warning: check if there are any deleted columns due to zero counts  "Deleted 2 columns of the response matrix due to zero counts"

# Define constants to avoid repetition ------------------------------------

# Grouping variables used throughout the analysis
GROUPING_VARS <- c(
  "probs", "phi", "k", "n", "m", "quant_min", "quant_max", "alpha",
  "dispersion", "alternative", "tol", "n_bisec", "C", "minpn",
  "minpnk", "minp", "maxp"
)

COVERAGE_METHODS <- c(
  "coverage_pointwise", "coverage_bonferroni", "coverage_bisection",
  "coverage_bisection_asym", "coverage_bonf_bisec", "coverage_mvn",
  "coverage_percentile_bt1", "coverage_percentile_bt_gemini",
  "coverage_bayesian_mcmc_cauchy_mean", "coverage_bayesian_mcmc_cauchy_marg", "coverage_bayesian_mcmc_cauchy_scs",
  "coverage_bayesian_mcmc_beta_mean", "coverage_bayesian_mcmc_beta_marg","coverage_bayesian_mcmc_beta_scs"
)



result_dirs <- c(".\\Code_and_Data\\Results\\intermediate_results\\threecats\\",
                 ".\\Code_and_Data\\Results\\intermediate_results\\fivecats\\",
                 ".\\Code_and_Data\\Results\\intermediate_results\\tencats\\"
)

# Load and Summarize Data in a Single Step -----------------------------

df_tot_covprob_long <- NULL
df_tot_eqt_long <- NULL
df_tot_widths_long <- NULL
df_tot_asymmetry <- NULL

for (i in 1:length(result_dirs)) {
  
df_summary <- list.files(
  path = result_dirs[i],
  pattern = "\\.rds$",
  full.names = TRUE
) %>%
  map_dfr(readRDS) %>%
  group_by(across(all_of(GROUPING_VARS))) %>%
  summarise(
    total_sim = n(),
    
    # Calculate coverage probabilities
    across(all_of(COVERAGE_METHODS),
           ~ mean(as.integer(.), na.rm = TRUE),
           .names = "covprob_{.col}"),
    # Calculate equi-tailed probabilities
    across(matches("_V\\d+$"), 
           ~ (n() - sum(. == 1, na.rm = TRUE)) / n(),
           .names = "{.col}_eqt_u"),
    across(matches("_V\\d+$"), 
           ~ (n() - sum(. == -1, na.rm = TRUE)) / n(),
           .names = "{.col}_eqt_l"),
    
    # calculations for the error probabilities
    across(matches("_V\\d+$"), 
           ~ sum(. == 1, na.rm = TRUE) / n(),
           .names = "{.col}_error_u"),
    across(matches("_V\\d+$"), 
           ~ sum(. == -1, na.rm = TRUE)  / n(),
           .names = "{.col}_error_l"),
    
    # Calculate mean of width statistics
    across(ends_with(c("_avg_width", "_median_width", "_sd_width", "_q25_width", "_q75_width")),
           ~ mean(., na.rm = TRUE)),
    
    # Summarise error/warning counts
    across(c(czercounts, cerrnas, cerrfin, cwarwz, cwarcon, cwarany, cwarzc),
           ~ sum(., na.rm = TRUE)),
    
    .groups = "drop" # Drop grouping for subsequent operations
  )

# Reshape Data to Tidy (Long) Format ------------------------------------

# Pivot coverage probabilities to long format
df_covprob_long <- df_summary %>%
  dplyr::select(all_of(GROUPING_VARS), total_sim, starts_with("covprob_")) %>%
  pivot_longer(
    cols = starts_with("covprob_"),
    names_to = "covprob_type",
    names_prefix = "covprob_",
    values_to = "covprob"
  )

df_eqt_long <- df_summary %>%
  dplyr::select(all_of(GROUPING_VARS), total_sim, ends_with(c("_eqt_l", "_eqt_u"))) %>%
  pivot_longer(
    cols = ends_with(c("_eqt_l", "_eqt_u")),
    names_to = c("method", ".value"),
    names_pattern = "(.*)_(eqt_[ul])"
  ) %>%
  extract(
    method,
    into = c("base_method", "CategoryNumber"),
    regex = "(.+)_V(\\d+)",
    remove = FALSE, 
    convert = TRUE  
  ) %>%
  rowwise() %>%
  mutate(
    prob = {
      all_probs <- str_extract_all(probs, "\\d+\\.\\d+")[[1]]
      as.numeric(all_probs[CategoryNumber])
    }
  ) %>%
  ungroup()

  
df_asymmetry <- df_summary %>%
  dplyr::select(all_of(GROUPING_VARS), total_sim, ends_with(c("_error_l", "_error_u"))) %>%
  
  pivot_longer(
    cols = ends_with(c("_error_l", "_error_u")),
    names_to = c("method", ".value"),
    names_pattern = "(.*)_(error_[ul])" 
  ) %>%
  
  mutate(
    asymmetry_ratio = ifelse(error_l > 0 & error_u > 0, error_u / error_l, NA_real_)
  ) %>% 
  extract(
    method,
    into = c("base_method", "CategoryNumber"),
    regex = "(.+)_V(\\d+)",
    remove = FALSE, 
    convert = TRUE  
  ) %>%
  rowwise() %>%
  mutate(
    prob = {
      all_probs <- str_extract_all(probs, "\\d+\\.\\d+")[[1]]
      as.numeric(all_probs[CategoryNumber])
    }
  ) %>%
  ungroup()

# Pivot all width statistics into a single tidy data frame
df_widths_long <- df_summary %>%
  dplyr::select(all_of(GROUPING_VARS), total_sim, ends_with("_width")) %>%
  pivot_longer(
    cols = ends_with("_width"),
    names_to = c("method", "statistic"),
    names_pattern = "(.*)_(avg|median|sd|q25|q75)_width",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = statistic,
    values_from = value,
    names_glue = "{statistic}_width" 
  )

df_tot_covprob_long <- rbind(df_tot_covprob_long, df_covprob_long)
df_tot_eqt_long <- rbind(df_tot_eqt_long, df_eqt_long)
df_tot_widths_long <- rbind(df_tot_widths_long, df_widths_long)
df_tot_asymmetry <- rbind(df_tot_asymmetry, df_asymmetry)
}

# Some checks
unique(df_tot_covprob_long$total_sim)
unique(df_tot_covprob_long$probs)
unique(df_tot_covprob_long$k)
unique(df_tot_covprob_long$n)


# Plots -------------------------------------------------------------------
# Define Standard Plotting Parameters
my_base_size <- 10 # font size

# Colors

# Okabe-Ito Palette (Colorblind Safe)
# Order: Orange, Sky Blue, Bluish Green, Yellow, Blue, Vermillion, Reddish Purple, Grey
oi_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Map Categories (C) to specific colors (Plot 1)
cols_categories <- c("3" = oi_palette[2],  # Sky Blue
                     "5" = oi_palette[1],  # Orange
                     "10" = oi_palette[6]) # Vermillion

## Coverage Probability ---------------------------------------------------

new_covprob_labels <- c(
  "coverage_pointwise" = "Pointwise",
  "coverage_bonferroni" = "Bonferroni",
  "coverage_bisection" = "Sym. Calibration",
  "coverage_bisection_asym" = "Asym. Calibration",
  "coverage_bonf_bisec" = "Marg. Calibration",
  "coverage_mvn" = "MVN",
  "coverage_percentile_bt1" = "Rank-Based SCS",
  "coverage_percentile_bt_gemini" = "Max. Abs. Res.",
  "coverage_bayesian_mcmc_cauchy_mean" = "B: Cauchy (mean)", 
  "coverage_bayesian_mcmc_cauchy_marg" = "B: Cauchy (marg)", 
  "coverage_bayesian_mcmc_cauchy_scs" = "B: Cauchy (SCS)",
  "coverage_bayesian_mcmc_beta_mean" = "B: Beta (mean)", 
  "coverage_bayesian_mcmc_beta_marg" = "B: Beta (marg)",
  "coverage_bayesian_mcmc_beta_scs" = "B: Beta (SCS)"
)

df_tot_covprob_long <- df_tot_covprob_long %>% 
  mutate(covprob_type = fct_relevel(covprob_type, 
                                    "coverage_pointwise", 
                                    "coverage_mvn", 
                                    "coverage_bonferroni",
                                    "coverage_bisection",
                                    "coverage_bisection_asym",
                                    "coverage_bonf_bisec",
                                    "coverage_percentile_bt_gemini",
                                    "coverage_percentile_bt1",
                                    "coverage_bayesian_mcmc_cauchy_mean", 
                                    "coverage_bayesian_mcmc_cauchy_marg", 
                                    "coverage_bayesian_mcmc_cauchy_scs",
                                    "coverage_bayesian_mcmc_beta_mean", 
                                    "coverage_bayesian_mcmc_beta_marg",
                                    "coverage_bayesian_mcmc_beta_scs"))

### complete Plot -------------------------------------------------------
df_tot_covprob_long <- df_tot_covprob_long %>%
  mutate(rho = (phi - 1) / (n - 1),
    sum_rho = ifelse(rho <= 0.1, 0.1, rho))
unique(df_tot_covprob_long$sum_rho)

plot_covprob_all <- df_tot_covprob_long %>% 
  #filter(total_sim > 300) %>%
  filter(covprob_type != "coverage_pointwise" & 
           covprob_type != "coverage_mvn" & 
           covprob_type != "coverage_bonferroni") %>%
  filter(C < 10) %>%
  #filter(covprob_type != "coverage_bayesian_mcmc_gamma_marg" & covprob_type != "coverage_bayesian_mcmc_cauchy_marg") %>%
  #filter(covprob_type == "coverage_bayesian_mcmc_cauchy_scs") %>%
  ggplot( aes(x = log10(minpnk), y = covprob))+ 
  geom_point(aes(color = sum_rho, shape = factor(n)), size = 2)+
  #facet_grid(phi ~.)+ 
  facet_grid(phi ~ covprob_type, labeller = labeller(covprob_type = new_covprob_labels))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  geom_hline(yintercept = 0.95) +
  geom_hline(yintercept = 0.95+sqrt(0.95*(1-0.95)/500),linetype="dashed") +
  geom_hline(yintercept = 0.95-sqrt(0.95*(1-0.95)/500),linetype="dashed") +
  ylab("Coverage Probability") +
  labs(color ="Probability Vector", shape = "Cluster Size (n)") +
  theme_bw(base_size = 18)+
  ylim(0.75, 1)

plot_covprob_all

### Splitted in 2 parts ---------------------------------------------------
# Used in Manuscript Figure 1

#write.csv(df_tot_covprob_long, "predint_covprob.csv", row.names = FALSE)

# Get the names of your facet levels in their current order
facet_levels <- levels(factor(df_tot_covprob_long$covprob_type))

plot_covprob_top <- df_tot_covprob_long %>%
  filter(covprob_type %in% facet_levels[1:5]) %>%
  arrange(desc(C)) %>%
  ggplot(aes(x = minpn, y = covprob)) +
  geom_point(aes(color = factor(C), shape = factor(k)), size = 1.25, alpha = 1) +
  facet_grid(
    phi ~ covprob_type,
    labeller = labeller(
      phi = function(value) parse(text = paste0("\u03A6 = ", value)),
      covprob_type = new_covprob_labels
    )
  ) +
  geom_hline(yintercept = 0.95) +
  geom_hline(yintercept = 0.95 + sqrt(0.95 * (1 - 0.95) / 1000), linetype = "dashed") +
  geom_hline(yintercept = 0.95 - sqrt(0.95 * (1 - 0.95) / 1000), linetype = "dashed") +
  scale_x_log10(
    name = expression(paste(min((pi[c] * n)))),
    labels = scales::label_number(drop0trailing = TRUE)
  ) +
  coord_cartesian(ylim = c(0.75, 1)) +
  ylab("Simultaneous Coverage Probability") +
  labs(colour = "No. of categories (C)", shape = "No. of clusters (k)") +
  theme_bw(base_size = my_base_size) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_color_manual(values = cols_categories, name = "No. of categories (C)")

plot_covprob_bottom <- df_tot_covprob_long %>%
  filter(covprob_type %in% facet_levels[c(6, 7, 8, 11, 14)]) %>%
  arrange(desc(C)) %>%
  ggplot(aes(x = minpn, y = covprob)) +
  geom_point(aes(color = factor(C), shape = factor(k)), size = 1.25, alpha = 1) +
  facet_grid(
    phi ~ covprob_type,
    labeller = labeller(
      phi = function(value) parse(text = paste0("\u03A6 = ", value)),
      covprob_type = new_covprob_labels
    )
  ) +
  geom_hline(yintercept = 0.95) +
  geom_hline(yintercept = 0.95 + sqrt(0.95 * (1 - 0.95) / 1000), linetype = "dashed") +
  geom_hline(yintercept = 0.95 - sqrt(0.95 * (1 - 0.95) / 1000), linetype = "dashed") +
  scale_x_log10(
    name = expression(paste(min((pi[c] * n)))),
    labels = scales::label_number(drop0trailing = TRUE)
  ) +
  coord_cartesian(ylim = c(0.75, 1)) +
  ylab("Simultaneous Coverage Probability") +
  labs(colour = "No. of categories (C)", shape = "No. of clusters (k)") +
  theme_bw(base_size = my_base_size) +
  scale_color_manual(values = cols_categories, name = "No. of categories (C)")


plot_covprob_all_split <- plot_covprob_top / plot_covprob_bottom + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

plot_covprob_all_split

# --- Save --- 
# Figure 1 in Manuscript
ggsave("./Code_and_Data/Figures/Budig_predint_multnom_Figure_1.pdf",
       plot = plot_covprob_all_split,
       width = 16,       # Full text width
       height = 18,      # Max text height
       units = "cm",
       device = cairo_pdf
      )         


# Only three categories
plot_covprob_top <- df_tot_covprob_long %>%
  filter(covprob_type %in% facet_levels[1:5]) %>%
  filter(C == 3) %>%
  arrange(desc(C)) %>%
  ggplot(aes(x = minpn, y = covprob)) +
  geom_point(aes(color = factor(k)), size = 1.25, alpha = 1) +
  facet_grid(
    phi ~ covprob_type,
    labeller = labeller(
      phi = function(value) parse(text = paste0("\u03A6 = ", value)),
      covprob_type = new_covprob_labels
    )
  ) +
  geom_hline(yintercept = 0.95) +
  geom_hline(yintercept = 0.95 + sqrt(0.95 * (1 - 0.95) / 1000), linetype = "dashed") +
  geom_hline(yintercept = 0.95 - sqrt(0.95 * (1 - 0.95) / 1000), linetype = "dashed") +
  scale_x_log10(
    name = expression(paste(min((pi[c] * n))))
  ) +
  coord_cartesian(ylim = c(0.75, 1)) +
  ylab("Simultaneous Coverage Probability") +
  labs(colour = "No. of clusters (k)") +
  theme_bw(base_size = my_base_size) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) 

plot_covprob_bottom <- df_tot_covprob_long %>%
  filter(covprob_type %in% facet_levels[c(6, 7, 8, 11, 14)]) %>%
  filter(C == 3) %>%
  arrange(desc(C)) %>%
  ggplot(aes(x = minpn, y = covprob)) +
  geom_point(aes(color = factor(k)), size = 1.25, alpha = 1) +
  facet_grid(
    phi ~ covprob_type,
    labeller = labeller(
      phi = function(value) parse(text = paste0("\u03A6 = ", value)),
      covprob_type = new_covprob_labels
    )
  ) +
  geom_hline(yintercept = 0.95) +
  geom_hline(yintercept = 0.95 + sqrt(0.95 * (1 - 0.95) / 1000), linetype = "dashed") +
  geom_hline(yintercept = 0.95 - sqrt(0.95 * (1 - 0.95) / 1000), linetype = "dashed") +
  scale_x_log10(
    name = expression(paste(min((pi[c] * n))))
  ) +
  coord_cartesian(ylim = c(0.75, 1)) +
  ylab("Coverage Probability") +
  labs(colour = "No. of clusters (k)") +
  theme_bw(base_size = my_base_size)


plot_covprob_all_split <- plot_covprob_top / plot_covprob_bottom + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

plot_covprob_all_split


## Equal Tail Probability -------------------------------------------------
pattern <- "bisection|bisec|percentile|bayesian_mcmc_cauchy_scs"

new_eqt_labels <- c(
  "pointwise" = "Pointwise",
  "bonferroni" = "Bonferroni",
  "bisection" = "Sym. Calibration",
  "bisection_asym" = "Asym. Calibration",
  "bonf_bisec" = "Marg. Calibration",
  "mvn" = "MVN",
  "percentile_bt_gemini" = "Max. Abs. Res.",
  "percentile_bt1" = "Rank-Based SCS",
  "bayesian_mcmc_cauchy_mean" = "B: Cauchy (mean)", 
  "bayesian_mcmc_cauchy_marg" = "B: Cauchy (marg)", 
  "bayesian_mcmc_cauchy_scs" = "B: Cauchy (SCS)",
  "bayesian_mcmc_beta_mean" = "B: Beta (mean)", 
  "bayesian_mcmc_beta_marg" = "B: Beta (marg)",
  "bayesian_mcmc_beta_scs" = "B: Beta (SCS)"
)

df_tot_eqt_long <- df_tot_eqt_long %>% 
  mutate(base_method = fct_relevel(base_method, 
                                    "pointwise", 
                                    "mvn", 
                                    "bonferroni",
                                    "bisection",
                                    "bisection_asym",
                                    "bonf_bisec",
                                    "percentile_bt1",
                                    "percentile_bt_gemini",
                                    "bayesian_mcmc_cauchy_mean", 
                                    "bayesian_mcmc_cauchy_marg", 
                                    "bayesian_mcmc_cauchy_scs",
#                                    "bayesian_mcmc_gamma_mean", 
#                                    "bayesian_mcmc_gamma_marg",
#                                    "bayesian_mcmc_gamma_scs",
                                    "bayesian_mcmc_beta_mean", 
                                    "bayesian_mcmc_beta_marg",
                                    "bayesian_mcmc_beta_scs"))

base_brbg_palette <- brewer.pal(n = 11, name = "BrBG")
base_brbg_palette <- brewer.pal(n = 11, name = "RdBu")
base_brbg_palette <- brewer.pal(n = 11, name = "Spectral")
base_brbg_palette <- brewer.pal(n = 11, name = "RdYlGn")
base_brbg_palette <- brewer.pal(n = 11, name = "RdYlBu")

color_ramp_brbg <- colorRampPalette(base_brbg_palette)

precise_palette <- color_ramp_brbg(100)
pal <- wes_palette("Zissou1", 100, type = "continuous")
## Equal Tails Plot -------------------------------------------------------

eqt_lu_facet_levels <- levels(factor(df_tot_eqt_long$base_method))

hline_data_l <- df_tot_eqt_long %>% 
  filter(grepl(pattern, base_method)) %>%
  #filter(base_method %in% eqt_lu_facet_levels[1:6]) %>%
  distinct(phi, base_method) %>%
  mutate(y_intercept = 1 - 0.05 / 3 / 2) 

hline_data_u <- df_tot_eqt_long %>%
  filter(grepl(pattern, base_method)) %>%
  #filter(base_method %in% eqt_lu_facet_levels[7:12]) %>%
  distinct(phi, base_method) %>% 
  mutate(y_intercept = 1 - 0.05 / 3 / 2) 

plot_eqt_l <- df_tot_eqt_long %>% 
  filter(grepl(pattern, base_method)) %>% 
  filter(C == 3) %>%
  ggplot(aes(x = minpn, y = eqt_l)) +
  geom_point(aes(color = prob, shape = factor(k)), size = 1) +
  ggh4x::facet_nested(
  phi ~ "Lower Border" + base_method,
  labeller = labeller(
    phi = function(value) parse(text = paste0("\u03A6 = ", value)),
    base_method = new_eqt_labels
  ),
  strip = ggh4x::strip_nested(
    text_x = list(
      element_text(size = 7.2)    # Method names
    ),
    by_layer_x = TRUE
  )
) +
  geom_hline(data = hline_data_l, aes(yintercept = y_intercept)) +
 scale_x_log10(
    labels = scales::label_number(drop0trailing = TRUE)
  ) +
  scale_color_gradientn(colours = pal) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  coord_cartesian(ylim = c(0.9, NA)) +
  ylab("Marginal Coverage Probability") +
  labs(colour = expression(paste("Category Probabilities ", pi[c])), shape = "No. of clusters (k)")

plot_eqt_u <- df_tot_eqt_long %>% 
  filter(grepl(pattern, base_method)) %>% 
  filter(C == 3) %>%
  ggplot(aes(x = minpn, y = eqt_u)) +
  geom_point(aes(color = prob, shape = factor(k)), size = 1) +
  ggh4x::facet_nested(
  phi ~ "Upper Border" + base_method,
  labeller = labeller(
    phi = function(value) parse(text = paste0("\u03A6 = ", value)),
    base_method = new_eqt_labels
  ),
  strip = ggh4x::strip_nested(
    text_x = list(
      element_text(size = 7.2)    # Method names
    ),
    by_layer_x = TRUE
  )
) +
  geom_hline(data = hline_data_u, aes(yintercept = y_intercept)) +
  scale_x_log10(
    name = expression(paste(min((pi[c] * n)))),
    labels = scales::label_number(drop0trailing = TRUE)
  ) +
  scale_color_gradientn(colours = pal) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim = c(0.9, NA)) +
  ylab("Marginal Coverage Probability") +
  labs(colour = expression(paste("Category Probabilities ", pi[c])), shape = "No. of clusters (k)")


plot_eqt_lu_split <- plot_eqt_l / plot_eqt_u + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
plot_eqt_lu_split

# --- Save ---
# Figure 2 in Manuscript
# ggsave("./Code_and_Data/Figures/Budig_predint_multnom_Figure_2.pdf",
#       plot = plot_eqt_lu_split,
#       width = 16,       # Full text width
#       height = 18,      # Max text height
#       units = "cm",
#       device = cairo_pdf )      
     
# Alternative (used in manuscript)

# Combine eqt_l and eqt_u into a long 'border' column
df_eqt_combined <- df_tot_eqt_long %>% 
  filter(grepl(pattern, base_method), C == 3) %>%
  pivot_longer(
    cols = c(eqt_l, eqt_u),
    names_to = "border",
    values_to = "coverage_prob"
  ) %>%
  mutate(
    border = factor(border, 
                    levels = c("eqt_l", "eqt_u"), 
                    labels = c("Lower Border", "Upper Border"))
  )

# Single reference line dataset
hline_data <- df_eqt_combined %>% 
  distinct(phi, base_method, border) %>% 
  mutate(y_intercept = 1 - 0.05 / 3 / 2)


plot_eqt_lu_split2 <- ggplot(df_eqt_combined, aes(x = minpn, y = coverage_prob)) +
  geom_point(aes(color = prob, shape = factor(k)), size = 1) +
  geom_hline(data = hline_data, aes(yintercept = y_intercept)) +
  ggh4x::facet_nested(
    border + phi ~ base_method,  # Moves border to the row strips
    labeller = labeller(
      phi = function(value) parse(text = paste0("\u03A6 = ", value)),
      base_method = new_eqt_labels
    ),
    strip = ggh4x::strip_nested(
      text_x = list(element_text(size = 7.2)),
      by_layer_y = TRUE
    )
  ) +
  scale_x_log10(
    name = expression(paste(min((pi[c] * n)))),
    labels = scales::label_number(drop0trailing = TRUE)
  ) +
  scale_color_gradientn(colours = pal) +
  coord_cartesian(ylim = c(0.9, NA)) +
  ylab("Marginal Coverage Probability") +
  labs(
    colour = expression(paste("Category Probabilities ", pi[c])),
    shape = "No. of clusters (k)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave("./Code_and_Data/Figures/Budig_predint_multnom_Figure_2.pdf",
       plot = plot_eqt_lu_split2,
       width = 16,       # Full text width
       height = 18,      # Max text height
       units = "cm",
       device = cairo_pdf )      
      

## Asymmetry Plot ---------------------------------------------------------
# Not used in Manuscript
# Create the plot
asymmetry_plot <- df_tot_asymmetry %>% 
  filter(grepl(pattern, base_method)) %>%
  #filter(C == 5 | C == 3) %>%
  #filter( C == 3) %>%
  #filter( C == 5) %>%
  #filter( C == 10) %>%
  ggplot(aes(x = log10(minpnk), y = asymmetry_ratio)) +
  geom_point(aes(color = prob, shape = factor(C)), alpha = 0.7) +
  
  # Add a horizontal line at y=1 to indicate perfect symmetry
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
  
  # Use a log scale for the y-axis
  scale_y_log10(breaks = c(0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32)) +
  
  # Facet by overdispersion (phi) and the method
  facet_grid(phi ~ base_method, labeller = labeller(base_method = new_eqt_labels,
                                                    phi =  function(value) parse(text = paste0("\u03A6 = ", value)))) +
  # scale_color_viridis_d() +
  # Use the high-resolution palette here
  scale_color_gradientn(colors = precise_palette)+
  labs(
    #title = "Asymmetry Ratio of Prediction Interval Tail Probabilities",
    #subtitle = "Ratio = P(y > Upper) / P(y < Lower)",
    x = expression(paste(log10(min( (pi[c]*n*K) )))),
    y = "Asymmetry Ratio",
    color = expression(paste("Category Probability (", pi[c], ")")),
    shape = "No. of categories (C)"
  ) +
  theme_bw() +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# Display the plot
print(asymmetry_plot)

ggsave("C:\\Users\\Budig\\Google_Drive\\Uni\\Phd\\03_Predint_multinomial\\Code\\Code_and_Data\\Figures\\plot_asymmetry_ratio.png",
       plot = asymmetry_plot,
       width = 10, height = 8,
       dpi = 900)

### Tryouts ---------------------------------------------------------------
# Not used in Manuscripts
# First, create a summary data frame
df_tot_covprob_summary <- df_tot_covprob_long %>%
  # Group by the variables you want on your axes and facets
  group_by(covprob_type, phi, minpnk, C) %>%
  # Calculate the mean coverage probability for each group
  summarise(mean_covprob = mean(covprob, na.rm = TRUE), .groups = "drop")

# Now, create the plot
df_tot_covprob_summary %>%
  filter(covprob_type != "coverage_pointwise" & 
           covprob_type != "coverage_mvn" & 
           covprob_type != "coverage_bonferroni" & 
           covprob_type != "coverage_nonparametric_bt" &
           covprob_type != "coverage_bayesian_mcmc_gamma") %>%
  ggplot( aes(x = log10(minpnk), y = mean_covprob, color = covprob_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.95, linetype = "solid") +
  # Use facets for the dispersion parameter
  facet_grid(C ~ phi, labeller = labeller(phi = label_parsed)) +
  labs(
    y = "Mean Coverage Probability",
    x = expression(paste(log10(min(m[g]*pi[gc])))),
    color = "Method"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")


ggplot(df_tot_covprob_long, aes(x = covprob_type, y = covprob)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  # Facet by the most important conditions
  facet_grid(phi ~ n, labeller = label_both) +
  labs(
    title = "Distribution of Coverage Probabilities by Method",
    x = "Method",
    y = "Coverage Probability"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Angle labels if they overlap
    legend.position = "none" # Color fill is for clarity, legend isn't essential
  )


# First, create a summary data frame
df_heatmap <- df_tot_covprob_long %>%
  group_by(covprob_type, phi, n) %>%
  summarise(mean_covprob = mean(covprob, na.rm = TRUE), .groups = "drop")

# Create the heatmap
ggplot(df_heatmap, aes(x = covprob_type, y = factor(phi), fill = mean_covprob)) +
  geom_tile(color = "white") + # Add white lines between tiles
  geom_text(aes(label = round(mean_covprob, 3)), size = 3, color = "white") + # Show the value
  facet_wrap(~ n, labeller = label_bquote(n == .(n))) + # Facet by cluster size
  scale_fill_viridis(
    name = "Mean\nCoverage",
    limits = c(0.8, 1.0) # Set limits to standardize color scale
  ) +
  labs(
    title = "Heatmap of Mean Coverage Probability",
    x = "Method", y = expression(paste("Dispersion ", Phi))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df_tot_covprob_long %>%
  mutate(coverage_error = covprob - 0.95) %>%
  ggplot(aes(x = covprob_type, y = coverage_error, fill = covprob_type)) +
  geom_boxplot() +
  # The new target is 0
  geom_hline(yintercept = 0, linetype = "solid", color = "red") +
  facet_grid(factor(C)~ phi) +
  labs(
    title = "Error in Coverage Probability (Deviation from 0.95)",
    x = "Method",
    y = "Coverage Error (Coverage - 0.95)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggplot(df_widths_long, aes(x = method, y = avg_width, fill = method)) +
  geom_violin(trim = FALSE) +
  # A log scale is often necessary if widths vary widely
  scale_y_log10() +
  facet_wrap(~ phi, labeller = label_bquote(phi == .(phi))) +
  labs(
    title = "Distribution of Average Interval Widths by Method",
    x = "Method",
    y = "Average Interval Width (Log Scale)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# First, calculate ranks for each scenario
df_ranks <- df_tot_covprob_long %>%
  # Group by each unique simulation setting
  group_by(phi, n, probs, minpnk) %>%
  # Rank methods by how close their coverage is to 0.95
  mutate(rank = rank(abs(covprob - 0.95), ties.method = "min")) %>%
  ungroup()

# Plot the distribution of ranks for each method
ggplot(df_ranks, aes(x = rank, fill = covprob_type)) +
  # A frequency polygon is great for this
  geom_freqpoly(aes(color = covprob_type), binwidth = 1, linewidth = 1.2) +
  labs(
    title = "Frequency of Performance Ranks by Method",
    x = "Performance Rank (1 is best)",
    y = "Number of Scenarios",
    color = "Method"
  ) +
  theme_bw()

df_tot_eqt_long <- df_tot_eqt_long %>%
  mutate(asymmetry_ratio = ifelse(eqt_l > 0, eqt_u / eqt_l, NA_real_))

# Create the plot
asymmetry_plot <- ggplot(df_tot_eqt_long, aes(x = log10(minpnk), y = ifelse(eqt_l > 0, eqt_u / eqt_l, NA_real_)
)) +
  geom_point(aes(color = prob, shape = factor(C)), alpha = 0.7) +
  
  # Add a horizontal line at y=1 to indicate perfect symmetry
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
  
  # Use a log scale for the y-axis
  #scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  
  # Facet by overdispersion (phi) and the method
  facet_grid(phi ~ base_method) +
  
  # Apply themes and labels
  #scale_color_viridis_d() +
  scale_color_gradientn(colors = precise_palette)+
  labs(
    #title = "Asymmetry Ratio of Prediction Interval Tail Probabilities",
    #subtitle = "Ratio = P(y > Upper) / P(y < Lower)",
    x = expression(paste(log10(min( (m[g]*pi[gc]) )))),
    y = "Asymmetry Ratio",
    color = expression(paste("Category Probability (", pi[c], ")"))
  ) +
  theme_bw() +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# Display the plot
print(asymmetry_plot)

