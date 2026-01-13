# Load libraries
library(ggplot2)

# Read the data files
bold <- scan("BOLD_HC.txt")
vaso <- scan("VASO_HC.txt")

# Make a data frame in long format
df <- data.frame(
  CNR = c(bold, vaso),
  Contrast = factor(rep(c("BOLD","VASO"), each = length(bold)))
)

# Violin plot with points and median
p <- ggplot(df, aes(x = Contrast, y = CNR, fill = Contrast)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.4) +
  geom_point(size = 4, alpha = 0.6, colour = "white") +
  stat_summary(fun = median, geom = "point", size = 5, shape = 23, alpha = 0.6,fill = "black") +
  theme_classic(base_size = 14) +
  labs(
    #title = "Voxel-wise hippocampal CNR (memory > math)",
    y = "CNR",
    x = ""
  ) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black", size = 20),
        axis.title = element_text(color = "black", size = 20),
        #plot.title = element_text(color = "black", size = 20, face = "bold", hjust = 0.5)
        )

ggsave("CNR_comp.pdf", plot = p, width = 6, height = 6, units = "in", dpi = 300)


## check for between group significant differences

wilcox.test(bold, vaso, paired = TRUE)
