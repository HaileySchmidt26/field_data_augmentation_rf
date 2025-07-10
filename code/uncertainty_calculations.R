# A Novel Approach to Field Data Augmentation with Remote Sensing and Machine Learning in Rangelands
# R Code by Javier Osorio Leyton & Hailey E. Schmidt
# Part III: Estimating Uncertainty
# -------------------------------------------------------
# load libraries
library(ggplot2)
library(dplyr)

# assuming your classification probabilities are stored as dataframes from Part I
# rf_probs = full vis, rf_probs_subvis = subset vis, rf_probs_gt = ground truth, rf_probs_aug = augmented
rf_probs$source <- "Full VIS"
rf_probs_subvis$source <- "Subset VIS"
rf_probs_gt$source <- "Ground Truth"
rf_probs_aug$source <- "Augmented"

# add max probability as a column
rf_probs$max_prob <- apply(rf_probs[, 1:(ncol(rf_probs)-2)], 1, max)
rf_probs_subvis$max_prob <- apply(rf_probs_subvis[, 1:(ncol(rf_probs_subvis)-2)], 1, max)
rf_probs_gt$max_prob <- apply(rf_probs_gt[, 1:(ncol(rf_probs_gt)-2)], 1, max)
rf_probs_aug$max_prob <- apply(rf_probs_aug[, 1:(ncol(rf_probs_aug)-2)], 1, max)

# function to compute entropy for a single vector of probabilities
compute_entropy <- function(prob_vec) {
  -sum(prob_vec * log2(prob_vec + 1e-10))  # small constant to avoid log(0)
}

# apply to each row of the class probability columns (excluding 'source' and any non-prob columns)
rf_probs$entropy <- apply(rf_probs[, 1:(ncol(rf_probs)-2)], 1, compute_entropy)
rf_probs_subvis$entropy <- apply(rf_probs_subvis[, 1:(ncol(rf_probs_subvis)-2)], 1, compute_entropy)
rf_probs_gt$entropy <- apply(rf_probs_gt[, 1:(ncol(rf_probs_gt)-2)], 1, compute_entropy)
rf_probs_aug$entropy <- apply(rf_probs_aug[, 1:(ncol(rf_probs_aug)-2)], 1, compute_entropy)

# combine rows and change to long format for plotting
all_probs <- bind_rows(rf_probs, rf_probs_subvis, rf_probs_gt, rf_probs_aug)
#write.csv(all_probs, "all_probs.csv", row.names = FALSE) 

plot_df <- all_probs %>%
  select(source, max_prob, entropy) %>%
  pivot_longer(cols = c(max_prob, entropy),
               names_to = "metric",
               values_to = "value")

# rename for clarity in the graph
plot_df$metric <- factor(plot_df$metric,
                         levels = c("entropy", "max_prob"),
                         labels = c("Prediction Uncertainty (Entropy)",
                                    "Prediction Confidence (Max Probability)"))

# plot results as a two panel graph 
# reorder and rename
plot_df$source <- factor(plot_df$source,
                         levels = c("Full VIS", "Subset VIS", "Ground Truth", "Augmented"),
                         labels = c("1",
                                    "2",
                                    "3",
                                    "4"))

ggplot(plot_df, aes(x = source, y = value, fill = source)) +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.4) +
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Model Confidence and Uncertainty Across Treatments",
       x = "Data Treatment",
       y = "Value") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold")) +
  scale_fill_brewer(palette = "Accent")
