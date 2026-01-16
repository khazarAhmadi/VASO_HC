## ---------------------------------------------------------
## Sensitivity analysis to obtain minimum detectable effect size 
## ---------------------------------------------------------
# obtain standardized MDES (one-sample Cohen’s d) for 80% power at α = 0.05

n <- 6
power <- 0.8
sig_level <- 0.05
pwr.t.test(n = n,
           d = NULL,
           sig.level = sig_level,
           power = power,
           type = "one.sample",
           alternative = "two.sided")


# this returns 1.43 as the standardized MDES


d_min <- 1.43 

## ---------------------------------------------------------
## function to calculate group mean, SD and the MDES in raw contrast units, calculated as d_min x SD 
## ---------------------------------------------------------

compute_stats <- function(filename, d_min) {
  x <- read.table(filename, header = FALSE)$V1
  
  mean_x <- mean(x)
  sd_x   <- sd(x)
  mdes_raw <- d_min * sd_x
  
  data.frame(
    Mean = mean_x,
    SD = sd_x,
    MDES_raw = mdes_raw
  )
}

## ---------------------------------------------------------
## Load data and compute statistics
## ---------------------------------------------------------

bold_gm  <- compute_stats("BOLD_GM_sensitivity_mean.txt", d_min)
vaso_gm  <- compute_stats("VASO_GM_sensitivity_mean.txt", d_min)
bold_hc  <- compute_stats("BOLD_HC_sensitivity_mean.txt", d_min)
vaso_hc  <- compute_stats("VASO_HC_sensitivity_mean.txt", d_min)

## ---------------------------------------------------------
## Combine into one table
## ---------------------------------------------------------

results <- rbind(
  data.frame(Imaging = "BOLD", ROI = "GM", bold_gm),
  data.frame(Imaging = "VASO", ROI = "GM", vaso_gm),
  data.frame(Imaging = "BOLD", ROI = "HC", bold_hc),
  data.frame(Imaging = "VASO", ROI = "HC", vaso_hc)
)

print(results)
