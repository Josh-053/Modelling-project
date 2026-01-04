data <- read.table("dataset.txt", 
                  sep = " ",        # space-separated
                  na.strings = ".", # treat "." as NA
                  header = TRUE)

library(ggplot2)
library(tidyverse)

data <- data %>%
  mutate(DOSE = as.integer(DOSE),
         EVID = case_when(
           is.na(CONC) & DOSE>0 ~ 1,
           CONC>=0 ~ 0,
           is.na(CONC) & is.na(DOSE) ~ 2))

data <- data %>% group_by(ID) %>%  
  fill(c(BW),.direction = "downup") %>% 
  ungroup()

data$AMT = data$DOSE * data$BW

write.table(data, file = "dataset_clean.txt", quote = FALSE, sep = " ", na = ".")
write.csv(data, file = "dataset_clean.csv", quote = FALSE, row.names = FALSE, na = ".")






library(ggplot2)
setwd('/home/ali/Desktop/Modeling/project')
dat <- read.table (file='run01_sdtab', skip=1, header=T)
ggplot (dat %>% filter(CONC!=0), aes(x = IPRED,
                 y = CONC)) +
  geom_point (color='blue') +
  geom_abline(aes(slope=1, intercept=0), linetype='dashed') +
  theme_light()

## check %BLQ
conc <- na.omit(data$CONC) # omit missing values from the concentrations
LLOQ <- 1.0 # given: LLOQ = 1.0 mg/L from ELISA
perc_BLQ <- mean(conc < LLOQ) * 100
perc_BLQ

## recommended BLQ handling method: M7+ (little bias, even at higher %BLQ)
## also recommended: M3 method after building final model

===vergeet even de lijnen hieronder, moet ze nog nakijken!!===

## filter conc < 1.0 mg/L
data_M7_plus <- data %>% 
  mutate(CONC = ifelse(CONC < 1, 0, CONC))

## WARNING: filtering the dataset will remove data points 

TALD <- data_filtered %>%
  mutate(
    number_of_doses = pmin(floor(TIME / 14) + 1, 4), # Calculate number of doses administered (max 4)
    tald = TIME - (number_of_doses - 1) * 14
  )

# ggplot(TALD,
#        aes(x = tald,
#            y = CONC))+
#   geom_point(size = 0.5) +
#   scale_y_log10()+
#   stat_summary_bin(bins = 9, # binning measured data point
#                    fun.y = mean, # arithmetic mean
#                    geom = "line")+
#   labs(x = "Time After Last Dose (days)",
#        y = "Plasma Concentration of Free PMX001 (")

## determined pattern: one (Monolix)

TVCL = 0.157 L/d /24 h/d

bioavailability = 58% (https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1016/j.clpt.2004.12.212)
CL = CL/F * F = 0.0187 * 0.58 = 0.010846