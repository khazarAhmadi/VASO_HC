library(R.matlab)
library(nlme)

#### read z-transformed signal change values per bin from *.mat files

sub <- readMat("All_cont_sub_Z.mat")
bold_sub <- sub$ALL.cont.sub.Z
sub_v <- readMat("All_cont_sub_Z_v.mat")
vaso_sub <- sub_v$ALL.cont.sub.Z.v

ca1 <- readMat('All_cont_ca1_Z.mat')
bold_ca1 <- ca1$ALL.cont.ca1.Z

ca1_v <- readMat('All_cont_ca1_Z_v.mat')
vaso_ca1 <- ca1_v$ALL.cont.ca1.Z.v

ca2 <- readMat('All_cont_ca2_Z.mat')
bold_ca2 <- ca2$ALL.cont.ca2.Z

ca2_v <- readMat("All_cont_ca2_Z_v.mat")
vaso_ca2 <- ca2_v$ALL.cont.ca2.Z.v

depth <- 1:ncol(bold_sub)
n_sub <- nrow(bold_sub)


#### build long data frames and merge them

df_bold_sub <- data.frame(
  Subject = rep(1:n_sub, each = length(depth)),
  Depth = rep(depth, times = n_sub),
  Contrast = "BOLD",
  Signal = as.vector(t(bold_sub))
)


df_bold_ca1 <- data.frame(
  Subject = rep(1:n_sub, each = length(depth)),
  Depth = rep(depth, times = n_sub),
  Contrast = "BOLD",
  Signal = as.vector(t(bold_ca1))
)



df_bold_ca2 <- data.frame(
  Subject = rep(1:n_sub, each = length(depth)),
  Depth = rep(depth, times = n_sub),
  Contrast = "BOLD",
  Signal = as.vector(t(bold_ca2))
)


df_vaso_sub <- data.frame(
  Subject = rep(1:n_sub, each = length(depth)),
  Depth = rep(depth, times = n_sub),
  Contrast = "VASO",
  Signal = as.vector(t(vaso_sub))
)


df_vaso_ca1 <- data.frame(
  Subject = rep(1:n_sub, each = length(depth)),
  Depth = rep(depth, times = n_sub),
  Contrast = "VASO",
  Signal = as.vector(t(vaso_ca1))
)



df_vaso_ca2 <- data.frame(
  Subject = rep(1:n_sub, each = length(depth)),
  Depth = rep(depth, times = n_sub),
  Contrast = "VASO",
  Signal = as.vector(t(vaso_ca2))
)




df_sub <- rbind(df_bold_sub, df_vaso_sub)

df_sub$DepthCat <- factor(
  cut(df_sub$Depth,
      breaks = c(0,10,20,30), 
      labels = c("SRLM","Inner","Outer")) ## srlm the first 10 bins, inner bins 11:20 and outer bins 21:30
)




df_ca1 <- rbind(df_bold_ca1, df_vaso_ca1)

df_ca1$DepthCat <- factor(
  cut(df_ca1$Depth,
      breaks = c(0,10,20,30), 
      labels = c("SRLM","Inner","Outer")) 
)





df_ca2 <- rbind(df_bold_ca2, df_vaso_ca2)

df_ca2$DepthCat <- factor(
  cut(df_ca2$Depth,
      breaks = c(0,10,20,30), 
      labels = c("SRLM","Inner","Outer")) 
)



#### fit mixed-effect model with interaction effect of contrast and depth, including random effect of subjects and ...
# aoutocorrelations across depths 


model_sub_cont <- lme(
  Signal ~ Contrast * Depth,                       # interaction between contrast and continuous depth
  random = ~1 | Subject,                           # subject random intercept
  correlation = corAR1(form = ~ Depth | Subject/Contrast),  # AR1 within each subject and contrast profile
  data = df_sub,
  method = "REML"
)

summary(model_sub_cont)



model_ca1_cont <- lme(
  Signal ~ Contrast * Depth,                       
  random = ~1 | Subject,                          
  correlation = corAR1(form = ~ Depth | Subject/Contrast), 
  data = df_ca1,
  method = "REML"
)

summary(model_ca1_cont)



model_ca2_cont <- lme(
  Signal ~ Contrast * Depth,                      
  random = ~1 | Subject,                          
  correlation = corAR1(form = ~ Depth | Subject/Contrast),  
  data = df_ca2,
  method = "REML"
)

summary(model_ca2_cont)
