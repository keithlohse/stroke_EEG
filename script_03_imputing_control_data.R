# Title: Brain Age - Imputing Control Data
# Author: Keith Lohse, PhD, PStat
# Date: 2023-10-26

library(tidyverse); library(patchwork)
library(psych); library(car); library(mice)
library(factoextra); library(corrplot)

getwd()
setwd("C:/Users/lohse/Box/Lohse Brain Age EEG/data/")
list.files()

#Read Data
DATA <- read.csv("data_CLEAN_COMBINED.csv", stringsAsFactors = TRUE)
summary(DATA$subID)
summary(DATA$Channels)

DATA <- DATA %>% select(-X) 

DATA %>% group_by(dataset, subID) %>%
  summarize(count = n(),
         age_ave = mean(age),
         sex=sex[1]) %>%
  group_by(dataset) %>%
  summarize(count = n(),
            age_mean = mean(age_ave),
            age_sd = sd(age_ave),
         female=sum(sex=="F"),
         male = sum(sex=="M"))
  


# Add rows for subjects with missing channels
all_subs <- levels(DATA$subID)
all_channels <- levels(DATA$Channels)
all_dat <- expand.grid(all_subs, all_channels)
all_dat <- all_dat %>% rename(subID = "Var1", Channels="Var2") %>%
  arrange(subID, Channels)
head(all_dat)

DATA <- merge(x=all_dat, y=DATA, by=c("subID", "Channels"), all=TRUE)


DATA <- DATA %>% group_by(subID) %>%
  fill(age, dataset, sex, .direction="downup") 

# Imputing missing Exponents ----
# Multivariate imputation through chained equations ----
# Conditional Multiple Imputation: Conditional MI, as indicated in its name, 
# follows an iterative procedure, modeling the conditional distribution of a 
# certain variable given the other variables. This technique allows users to be 
# more flexible as a distribution is assumed for each variable rather than the 
# whole dataset.
head(DATA)
EXP_WIDE <- DATA %>% select(subID, Channels, Exponent) %>%
  pivot_wider(names_from = Channels, values_from = Exponent)

head(EXP_WIDE)

md.pattern(EXP_WIDE[,2:29])
md.pairs(EXP_WIDE[,2:29])

set.seed(21)
imp1 <- mice(EXP_WIDE[,2:29], 
             m = 5,
             method="pmm")
imp1
imp1$imp

imp_tot <- complete(imp1, "broad")
imp_tot

# Combining Imputed Slopes with other variables
head(DATA)
head(imp_tot)
IMP_LONG <- imp_tot %>% mutate(subID = EXP_WIDE$subID) %>%
  pivot_longer(cols=C3.1:T8.5 , 
               names_to=c("Channels", "iter"),
               names_pattern = "(.*)\\.(.*)", # regular expressions are the worst!
               values_to = "Imp_Exp") %>%
  arrange(subID, Channels, iter)

IMP_LONG <- IMP_LONG %>% group_by(subID, Channels) %>%
  summarize(Imp_Exp = mean(Imp_Exp))

DATA <- merge(x=IMP_LONG, y=DATA, by=c("subID", "Channels"), all=TRUE)



# Imputing missing Offsets ----
# Multivariate imputation through chained equations
# Conditional Multiple Imputation: Conditional MI, as indicated in its name, 
# follows an iterative procedure, modeling the conditional distribution of a 
# certain variable given the other variables. This technique allows users to be 
# more flexible as a distribution is assumed for each variable rather than the 
# whole dataset.
head(DATA)
OFF_WIDE <- DATA %>% select(subID, Channels, Offset ) %>%
  pivot_wider(names_from = Channels, values_from = Offset)

head(OFF_WIDE)

md.pattern(OFF_WIDE[,2:29])
md.pairs(OFF_WIDE[,2:29])

set.seed(21)
imp1 <- mice(OFF_WIDE[,2:29], 
             m = 5,
             method="pmm")
imp1
imp1$imp

imp_tot <- complete(imp1, "broad")
imp_tot

# Combining Imputed Slopes with other variables
head(DATA)
head(imp_tot)
IMP_LONG <- imp_tot %>% mutate(subID = OFF_WIDE$subID) %>%
  pivot_longer(cols=C3.1:T8.5 , 
               names_to=c("Channels", "iter"),
               names_pattern = "(.*)\\.(.*)", # regular expressions are the worst!
               values_to = "Imp_Off") %>%
  arrange(subID, Channels, iter)

IMP_LONG <- IMP_LONG %>% group_by(subID, Channels) %>%
  summarize(Imp_Off = mean(Imp_Off))

DATA <- merge(x=IMP_LONG, y=DATA, by=c("subID", "Channels"), all=TRUE)
head(DATA)


# Saving Data ----
write.csv(DATA, "./data_IMPUTED_CONTROLS.csv")
rm(list = ls())
