##------------------------------------------------------------
# This snippet creates bar plots to visualize mean activity across hippocampal subfields for BOLD and VASO data
# under memory vs math contrast 
#the required input is *_activity.txt that is generated via the bash script in folder voxelwise_glm
##------------------------------------------------------------

setwd('Desktop/VASO_HC/Figs_script/HC_subfield_activity/')
library(ggplot2)
library(dplyr)

## function to calculate means and SEs
calculate_mean_se <- function(file_paths, n_sub) {
  data_list <- lapply(file_paths, function(file) read.table(file, header = FALSE, sep = "\t")$V1)
  
  # Calculate mean and SE for each region
  means <- sapply(data_list, mean)
  ses <- sapply(data_list, function(data) sd(data) / sqrt(n_sub))
  
  return(list(means = means, ses = ses))
}

N_sub = 6

## Define file paths for BOLD and VASO data
bold_files <- c("bold_Sub_activity.txt", "bold_CA1_activity.txt", "bold_CA2_activity.txt", 
                "bold_CA3_activity.txt", "bold_DG_activity.txt", "bold_SRLM_activity.txt")

vaso_files <- c("vaso_Sub_activity.txt", "vaso_CA1_activity.txt", "vaso_CA2_activity.txt", 
                "vaso_CA3_activity.txt", "vaso_DG_activity.txt", "vaso_SRLM_activity.txt")

# Calculate mean and SE for BOLD and VASO
bold_stats <- calculate_mean_se(bold_files, N_sub)
vaso_stats <- calculate_mean_se(vaso_files, N_sub)

# Create data frames for plotting
bold_df <- data.frame(
  Region = c("Sub", "CA1", "CA2", "CA3", "CA4/DG", "SRLM"),
  Mean = bold_stats$means,
  SE = bold_stats$ses
)

vaso_df <- data.frame(
  Region = c("Sub", "CA1", "CA2", "CA3", "CA4/DG", "SRLM"),
  Mean = vaso_stats$means,
  SE = vaso_stats$ses
)


## plot BOLD
p <- ggplot(bold_df, aes(x = Region, y = Mean)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#2ca02c", color = "white", size = 1.2) +  # All bars green
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, color = "white", size = 1.2) +
  labs(
    x = "HC subregion",
    y = "Mean BOLD Signal change"
  ) +
  ylim(0, 2.5) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "white", size = 16),
    axis.title = element_text(color = "white", size = 18, face = "bold"),
    axis.line = element_line(color = "white", size = 1.2),
    axis.ticks = element_line(color = "white", size = 1.2),
    plot.title = element_text(color = "white", size = 20, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave("HC_subfield_activity_bold.pdf", plot = p, width = 15, height = 16, units = "in", dpi = 300)

###### Look for between group statistical difference 

bold_data <- bind_rows(
  lapply(1:6, function(i) {
    data.frame(
      Signal = read.table(bold_files[i], header = FALSE, sep = "\t")$V1,
      Region = c("Sub", "CA1", "CA2", "CA3", "CA4/DG", "SRLM")[i]
    )
  })
)

anova_result_bold <- aov(Signal ~ Region, data = bold_data)
summary(anova_result_bold)  # Non-significant




## plot VASO
p1 <- ggplot(vaso_df, aes(x = Region, y = Mean)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#2ca02c", color = "white", size = 1.2) +  # All bars green
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, color = "white", size = 1.2) +
  labs(
    x = "HC subregion",
    y = "Mean VASO Signal change"
  ) +
  ylim(0, 2.5) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "white", size = 16),
    axis.title = element_text(color = "white", size = 18, face = "bold"),
    axis.line = element_line(color = "white", size = 1.2),
    axis.ticks = element_line(color = "white", size = 1.2),
    plot.title = element_text(color = "white", size = 20, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave("HC_subfield_activity_vaso.pdf", plot = p1, width = 15, height = 16, units = "in", dpi = 300)

###### Look for between group statistical difference 

vaso_data <- bind_rows(
  lapply(1:6, function(i) {
    data.frame(
      Signal = read.table(vaso_files[i], header = FALSE, sep = "\t")$V1,
      Region = c("Sub", "CA1", "CA2", "CA3", "CA4/DG", "SRLM")[i]
    )
  })
)

# Perform one-way ANOVA for VASO
anova_result_vaso <- aov(Signal ~ Region, data = vaso_data)
summary(anova_result_vaso)  # Non-significant
