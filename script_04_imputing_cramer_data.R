# Title: Brain Age - Adults with Stroke (Cramer Data)
# Author: Keith Lohse, PhD, PStat
# Date: 2023-10-26

# 1.0 Import data ----
library(tidyverse); library(patchwork)
library(psych); library(car); library(mice)
library(factoextra); library(corrplot)

getwd()
setwd("C:/Users/kelop/Box/Lohse Brain Age EEG/data/")

list.files()

# Read in the merged data files for the Cramer dataset
CRAMER <- read.csv("cramer_MERGED.csv", stringsAsFactors = TRUE)
head(CRAMER)

hist(CRAMER[CRAMER$Channels=="E126",]$CF1, breaks=25)


# Read in imputed control data for reference
CTL_DATA <- read.csv("data_IMPUTED_CONTROLS.csv", stringsAsFactors = TRUE)
head(CTL_DATA)



# Why do we have 203 Cz channels? 
summary(CRAMER$Channels)
CRAMER <- CRAMER %>% group_by(subID, Channels) %>% dplyr::slice(1)
summary(CRAMER$Channels)
hist(CRAMER[CRAMER$Channels=="E126",]$CF1, breaks=25)



# 2.0 Imputed Data and Relabel High Density EEG ----
# Impute missing values in Cramer data ----
# Convert 256 Hydrogel electrode labels to 10-20 labels, 
# drop the unused electrodes from the Cramer data:
head(CRAMER)
summary(CRAMER$Channels)
summary(CTL_DATA$Channels)

# Imputation using MICE
head(CRAMER)
all_subs <- levels(CRAMER$subID)
all_channels <- levels(CRAMER$Channels)
all_dat <- expand.grid(all_subs, all_channels)
all_dat <- all_dat %>% rename(subID = "Var1", Channels="Var2") %>%
  arrange(subID, Channels)
head(all_dat)

DATA <- merge(x=all_dat, y=CRAMER, by=c("subID", "Channels"), all=TRUE)
summary(DATA$Channel)

# MICE Imputation of Exponents ----
head(DATA)
EXP_WIDE <- DATA %>% dplyr::select(subID, Channels, Exponent) %>%
  pivot_wider(names_from = Channels, values_from = Exponent)

head(EXP_WIDE)

md.pattern(EXP_WIDE[,2:258])
md.pairs(EXP_WIDE[,2:258])

set.seed(21)
imp1 <- mice(EXP_WIDE[,2:258], 
             m = 5,
             method="pmm")
imp1
imp1$imp

imp_tot <- complete(imp1, "broad")
imp_tot

# Combining Imputed Slopes with other variables
head(DATA)
head(imp_tot)
EXP_LONG <- imp_tot %>% mutate(subID = EXP_WIDE$subID) %>%
  pivot_longer(cols=Cz.1:E99.5, 
               names_to=c("Channels", "iter"),
               names_pattern = "(.*)\\.(.*)", # regular expressions are the worst!
               values_to = "Imp_Exp") %>%
  arrange(subID, Channels, iter)

EXP_LONG <- EXP_LONG %>% group_by(subID, Channels) %>%
  summarize(Imp_Exp = mean(Imp_Exp))
head(EXP_LONG)
head(DATA)

DATA <- merge(x=EXP_LONG, y=DATA, by=c("subID", "Channels"), all=TRUE)


# MICE Imputation of Offsets ----
head(DATA)
OFF_WIDE <- DATA %>% dplyr::select(subID, Channels, Offset) %>%
  pivot_wider(names_from = Channels, values_from = Offset)

head(OFF_WIDE)

md.pattern(OFF_WIDE[,2:258])
md.pairs(OFF_WIDE[,2:258])

set.seed(21)
imp1 <- mice(OFF_WIDE[,2:258], 
             m = 5,
             method="pmm")
imp1
imp1$imp

imp_tot <- complete(imp1, "broad")
imp_tot

# Combining Imputed Slopes with other variables
head(DATA)
head(imp_tot)
OFF_LONG <- imp_tot %>% mutate(subID = EXP_WIDE$subID) %>%
  pivot_longer(cols=Cz.1:E99.5, 
               names_to=c("Channels", "iter"),
               names_pattern = "(.*)\\.(.*)", # regular expressions are the worst!
               values_to = "Imp_Off") %>%
  arrange(subID, Channels, iter)

OFF_LONG <- OFF_LONG %>% group_by(subID, Channels) %>%
  summarize(Imp_Off = mean(Imp_Off))
head(OFF_LONG)
head(DATA)

DATA <- merge(x=OFF_LONG, y=DATA, by=c("subID", "Channels"), all=TRUE)

head(DATA)


####

write.csv(DATA, "./data_CRAMER_IMPUTED.csv")

# Downsample to select channels and rename ----
# See: https://www.egi.com/images/HydroCelGSN_10-10.pdf

DATA <- read.csv("./data_CRAMER_IMPUTED.csv", header=TRUE,
                 stringsAsFactors = TRUE)

levels(DATA$Channels)
DATA %>% group_by(Channels) %>% summarize(Count=n())

list.files()
CONTROLS <- read.csv("./data_IMPUTED_CONTROLS.csv", header=TRUE,
                 stringsAsFactors = TRUE)
summary(CONTROLS$Channels)
unique(CONTROLS$Channels)
# List of common electrodes from Pathania, Euler, and MPI
keep_list <- c(unique(CONTROLS$Channels)) 

# rename the 256 channels to closest 10-20 electrode
DATA <- DATA %>% 
  #select(-X) %>%
  mutate(Channels = fct_recode(Channels,
                           Cz = "Cz", C3 = "E59", C4 = "E183", F3 = "E36",
                           F4 = "E224", F7 = "E47", F8 = "E2", Fp1 = "E37",
                           Fp2 = "E18", Fpz = "E26", Fz = "E21", O1 = "E116",
                           O2 = "E150", P3 = "E87", P4 = "E153", T7 = "E69",
                           T8 = "E202", P7 = "E96", P8 = "E170", Pz = "E101",
                           Poz = "E119", F2 = "E5", FC5 = "E49", Ft10 = "E219",
                           C6 = "E194", Ft9 = "E67", F6 = "E222", FT8 = "E211",
                           AF8 = "E10", CpZ = "E81", CP6 = "E172", C5 = "E64",
                           CP4 = "E164", P10 = "E169", F9 = "E252", P1 = "E88",
                           P5 = "E86", AF3 = "E34", C1 = "E44", PO8 = "E161",
                           AF4 = "E12", Afz = "E20", TP8 = "E179", FC3 = "E42",
                           CP3 = "E66", P6 = "E162", PO3 = "E109", C2 = "E185",
                           FC1 = "E24", PO4 = "E140", Oz = "E126", CP2 = "E143",
                           FC2 = "E207", CP1 = "E79", TP9 = "E94", F1 = "E29",
                           Fcz = "E15", TP10 = "E190", F10 = "E226", P2 = "E142",
                           F5 = "E48", P9 = "E106", FC4 = "E206", CP5 = "E76",
                           FC6 = "E213", Afz = "E27", Po7 = "E97", AF7 = "E46",
                           Afz = "E26", Tp7 = "E84", FT7 = "E62", T9 = "E68",
                           T10 = "E210")) %>%
  filter(Channels %in% keep_list) %>%
  mutate(Channels = factor(Channels))

DATA %>% group_by(subID) %>% summarise(count=n())

head(DATA)

keep_list[keep_list %in% levels(DATA$Channels)==TRUE]
keep_list[keep_list %in% levels(DATA$Channels)==FALSE]

write.csv(DATA, "./data_CRAMER_IMPUTED_SUBSET.csv")







