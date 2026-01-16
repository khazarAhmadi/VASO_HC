##------------------------------------------------------------
# This snippet creates box plots to visualize mean tSNR values in the hippocampus for BOLD and VASO data, 
# acquired across three scanners
##------------------------------------------------------------


library(ggplot2)
library(readxl)


## starting with bold

BOLD <- read_excel('Scanner_tSNR_comp_bold.xlsx')
View(BOLD)


p <- ggplot(BOLD, aes(x = Scanner, y = HC_tSNR, fill = Scanner)) +
  geom_boxplot(size = 0.7, color = "white") +
  geom_point(color = 'white', size = 8, shape = 21, stroke = 1.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 12,
               fill = "red", color = "red", stroke = 2) +
  ylim(0, 36) +
  scale_fill_brewer(palette = "Set2") +  # <- sharper, more vivid colors
  theme(
    text = element_text(size = 20, color = "white"),
    panel.border = element_rect(fill = NA, color = "white", size = 2),
    axis.text.x = element_text(size = 30, color = "white"),
    axis.text.y = element_text(size = 30, color = "white"),
    panel.background = element_rect(fill = 'black'),
    plot.background = element_rect(fill = 'black', colour = 'black'),
    legend.background = element_rect(fill = "black"),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

ggsave("HC_tSNR_bold_scanners.pdf", plot = p, width = 15, height = 16, units = "in", dpi = 300)


## apply the same for vaso 

VASO <- read_excel('Scanner_tSNR_comp_vaso.xlsx')
View(VASO)

p1 <- ggplot(VASO, aes(x = Scanner, y = HC_tSNR, fill = Scanner)) +
  geom_boxplot(size = 0.7, color = "white") +
  geom_point(color = 'white', size = 8, shape = 21, stroke = 1.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 12,
               fill = "red", color = "red", stroke = 2) +
  ylim(0, 36) +
  scale_fill_brewer(palette = "Set2") +  # <- sharper, more vivid colors
  theme(
    text = element_text(size = 20, color = "white"),
    panel.border = element_rect(fill = NA, color = "white", size = 2),
    axis.text.x = element_text(size = 30, color = "white"),
    axis.text.y = element_text(size = 30, color = "white"),
    panel.background = element_rect(fill = 'black'),
    plot.background = element_rect(fill = 'black', colour = 'black'),
    legend.background = element_rect(fill = "black"),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

ggsave("HC_tSNR_vaso_scanners.pdf", plot = p1, width = 15, height = 16, units = "in", dpi = 300)
