#################################################################### 
# Snippet to plot subfield-specific tSNR values for BOLD and VASO and check for statistical difference 
####################################################################

library(readxl)
library(nlme)
library(emmeans)
library(ggplot2)

df <- read_excel('subfields_tSNR.xlsx')

df_bold <- subset(df, contrast == 'BOLD')
df_vaso <- subset(df,contrast =='VASO')

subfields <- c("Sub", "CA1", "CA2", "CA3", "CA4DG")

df_long_b <- cbind(
  df_bold[c("subject", "Run", "contrast")],
  stack(df_bold[subfields])
)

names(df_long_b)[4:5] <- c("tSNR", "Subfield")


df_long_v <- cbind(
  df_vaso[c("subject", "Run", "contrast")],
  stack(df_vaso[subfields])
)

names(df_long_v)[4:5] <- c("tSNR", "Subfield")


### First plot tSNR across subfields for BOLD
p <- ggplot(df_long_b, aes(x = Subfield, y = tSNR, fill = Subfield)) +
  
  geom_boxplot(size = 0.7, color = "white") +
  geom_point(color = 'white', size = 8, shape = 21, stroke = 1.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 12,
               fill = "red", color = "red", stroke = 2) +
  ylim(0, 45) +
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
ggsave("Subfields_tSNR_bold.pdf", plot = p, width = 15, height = 16, units = "in", dpi = 300)


#### Now plot the same for VASO 

pv <- ggplot(df_long_v, aes(x = Subfield, y = tSNR, fill = Subfield)) +
  
  geom_boxplot(size = 0.7, color = "white") +
  geom_point(color = 'white', size = 8, shape = 21, stroke = 1.2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 12,
               fill = "red", color = "red", stroke = 2) +
  ylim(0, 45) +
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
ggsave("Subfields_tSNR_vaso.pdf", plot = pv, width = 15, height = 16, units = "in", dpi = 300)



#################################### check for statistical differences using mixed-effect models 
model <- lme(
  tSNR ~ Subfield,
  random = ~ 1 | subject,
  data = df_long_b,
  method = "REML"
)

summary(model)
emmeans(model, pairwise ~ Subfield, adjust = "tukey")


modelv <- lme(
  tSNR ~ Subfield,
  random = ~ 1 | subject,
  data = df_long_v,
  method = "REML"
)

summary(modelv)
emmeans(modelv, pairwise ~ Subfield, adjust = "tukey")
