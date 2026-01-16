##------------------------------------------------------------
# This snippet creates bar plots to visualize mean activity across hippocampal subfields (for conjuction ROIs representing 
# significant HC voxels for BOLD + VASO under memory vs math contrast at p < 0.05 ) and across entire anatomically defined ROIs regardless of statistical thresholding.
# the required input is *_activity.txt that is generated via the bash script in folder voxelwise_glm
##------------------------------------------------------------

setwd('HC_subfield_activity_conjunctionROIs/')
library(ggplot2)
library(dplyr)
library(gdata)

## function to calculate means and SEs

calculate_mean_se_correct_sem <- function(files, n_sub) {
  data_list <- lapply(files, function(f) read.table(f, header = FALSE))
  
  means_list <- lapply(data_list, function(x) x$V1)  # subject means
  
  means <- sapply(means_list, mean)                  # group mean
  
  # SEM calculated from between-subject SD of means
  sems <- sapply(means_list, function(x) sd(x) / sqrt(n_sub))
  
  return(list(means = means, ses = sems))
}

N_sub = 6


## Define file names for BOLD and VASO data
bold_files <- c("bold_Sub_activity.txt", "bold_CA1_activity.txt", "bold_CA2_activity.txt", 
                "bold_CA3_activity.txt", "bold_DG_activity.txt", "bold_SRLM_activity.txt")

vaso_files <- c("vaso_Sub_activity.txt", "vaso_CA1_activity.txt", "vaso_CA2_activity.txt", 
                "vaso_CA3_activity.txt", "vaso_DG_activity.txt", "vaso_SRLM_activity.txt")

regions <- c("Sub", "CA1", "CA2", "CA3", "CA4/DG", "SRLM")

# Calculate mean and SE for BOLD and VASO
bold_stats <- calculate_mean_se_correct_sem(bold_files, N_sub)
vaso_stats <- calculate_mean_se_correct_sem(vaso_files, N_sub)

# Create data frames for plotting
bold_df <- data.frame(
  Region = regions,
  Mean = bold_stats$means,
  SE = bold_stats$ses
)


vaso_df <- data.frame(
  Region = regions,
  Mean = vaso_stats$means,
  SE = vaso_stats$ses
)

## plot BOLD
p <- ggplot(bold_df, aes(x = Region, y = Mean)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#2ca02c", color = "white", size = 3) +  # All bars green
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, color = "white", size = 3) +
  labs(
    x = "HC subregion",
    y = "Mean BOLD Signal change"
  ) +
  ylim(0, 5.5) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "white", size = 40),
    axis.title = element_text(color = "white", size = 40, face = "bold"),
    axis.line = element_line(color = "white", size = 3),
    axis.ticks = element_line(color = "white", size = 3),
    plot.title = element_text(color = "white", size = 30, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave("HC_subfield_activity_bold.pdf", plot = p, width = 15, height = 16, units = "in", dpi = 300)

###### Look for between subfield statistical difference 

# Read means per subject
bold_subject_means <- lapply(bold_files, function(f) read.table(f, header = FALSE)[,1])

# Combine into one data frame
bold_df_subject <- as.data.frame(bold_subject_means)
colnames(bold_df_subject) <- regions

# Add subject ID
bold_df_subject$Subject <- 1:nrow(bold_df_subject)



friedman_result_bold <- friedman.test(as.matrix(bold_df_subject[, regions]))
friedman_result_bold








## plot VASO
p1 <- ggplot(vaso_df, aes(x = Region, y = Mean)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#2ca02c", color = "white", size = 3) +  # All bars green
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, color = "white", size = 3) +
  labs(
    x = "HC subregion",
    y = "Mean VASO Signal change"
  ) +
  ylim(0, 5.5) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "white", size = 40),
    axis.title = element_text(color = "white", size = 40, face = "bold"),
    axis.line = element_line(color = "white", size = 3),
    axis.ticks = element_line(color = "white", size = 3),
    plot.title = element_text(color = "white", size = 30, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave("HC_subfield_activity_vaso.pdf", plot = p1, width = 15, height = 16, units = "in", dpi = 300)


###### Look for between subfield statistical difference 

vaso_subject_means <- lapply(vaso_files, function(f) read.table(f, header = FALSE)[,1])

# Combine into one data frame
vaso_df_subject <- as.data.frame(vaso_subject_means)
colnames(vaso_df_subject) <- regions

# Add subject ID
vaso_df_subject$Subject <- 1:nrow(vaso_df_subject)

friedman_result_vaso <- friedman.test(as.matrix(vaso_df_subject[, regions]))
friedman_result_vaso


########################################### Do the same but this time for entire anatomical ROIs
setwd('../HC_whole_subfield_activity/')
keep(calculate_mean_se_correct_sem, N_sub,sure = T)


bold_files <- c("bold_Sub_activity.txt", "bold_CA1_activity.txt", "bold_CA2_activity.txt", 
                "bold_CA3_activity.txt", "bold_DG_activity.txt", "bold_SRLM_activity.txt")

vaso_files <- c("vaso_Sub_activity.txt", "vaso_CA1_activity.txt", "vaso_CA2_activity.txt", 
                "vaso_CA3_activity.txt", "vaso_DG_activity.txt", "vaso_SRLM_activity.txt")

regions <- c("Sub", "CA1", "CA2", "CA3", "CA4/DG", "SRLM")

# Calculate mean and SE for BOLD and VASO
bold_stats <- calculate_mean_se_correct_sem(bold_files, N_sub)
vaso_stats <- calculate_mean_se_correct_sem(vaso_files, N_sub)

# Create data frames for plotting
bold_df <- data.frame(
  Region = regions,
  Mean = bold_stats$means,
  SE = bold_stats$ses
)


vaso_df <- data.frame(
  Region = regions,
  Mean = vaso_stats$means,
  SE = vaso_stats$ses
)

## plot BOLD
p <- ggplot(bold_df, aes(x = Region, y = Mean)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#2ca02c", color = "white", size = 3) +  # All bars green
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, color = "white", size = 3) +
  labs(
    x = "HC subregion",
    y = "Mean BOLD Signal change"
  ) +
  ylim(-0.4, 0.8) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "white", size = 40),
    axis.title = element_text(color = "white", size = 40, face = "bold"),
    axis.line = element_line(color = "white", size = 3),
    axis.ticks = element_line(color = "white", size = 3),
    plot.title = element_text(color = "white", size = 30, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave("HC_subfield_activity_bold1.pdf", plot = p, width = 15, height = 16, units = "in", dpi = 300)

###### Look for between subfield statistical difference 

# Read means per subject
bold_subject_means <- lapply(bold_files, function(f) read.table(f, header = FALSE)[,1])

# Combine into one data frame
bold_df_subject <- as.data.frame(bold_subject_means)
colnames(bold_df_subject) <- regions

# Add subject ID
bold_df_subject$Subject <- 1:nrow(bold_df_subject)



friedman_result_bold <- friedman.test(as.matrix(bold_df_subject[, regions]))
friedman_result_bold








## plot VASO
p1 <- ggplot(vaso_df, aes(x = Region, y = Mean)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#2ca02c", color = "white", size = 3) +  # All bars green
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, color = "white", size = 3) +
  labs(
    x = "HC subregion",
    y = "Mean VASO Signal change"
  ) +
  ylim(-0.4, 0.8) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "white", size = 40),
    axis.title = element_text(color = "white", size = 40, face = "bold"),
    axis.line = element_line(color = "white", size = 3),
    axis.ticks = element_line(color = "white", size = 3),
    plot.title = element_text(color = "white", size = 30, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave("HC_subfield_activity_vaso.pdf", plot = p1, width = 15, height = 16, units = "in", dpi = 300)


###### Look for between subfield statistical difference 

vaso_subject_means <- lapply(vaso_files, function(f) read.table(f, header = FALSE)[,1])

# Combine into one data frame
vaso_df_subject <- as.data.frame(vaso_subject_means)
colnames(vaso_df_subject) <- regions

# Add subject ID
vaso_df_subject$Subject <- 1:nrow(vaso_df_subject)

friedman_result_vaso <- friedman.test(as.matrix(vaso_df_subject[, regions]))
friedman_result_vaso

