# Title: Brain Age - Merging Datasets
# Author: Keith Lohse, PhD, PStat
# Date: 2023-10-26
library(tidyverse)

getwd()
setwd("C:/Users/lohse/Box/Lohse Brain Age EEG/data/")

list.files()

# Read in the merged data files for each data set
PATHANIA <- read.csv("pathania_MERGED.csv", stringsAsFactors = TRUE)
head(PATHANIA)

EULER <- read.csv("euler_MERGED.csv", stringsAsFactors = TRUE)
head(EULER)

MPI <- read.csv("mpi_MERGED.csv", stringsAsFactors = TRUE)
head(MPI)
# MPI <- MPI %>% mutate(CF1 = ifelse(is.na(CF1)==TRUE, 0, CF1),
#                PW1 = ifelse(is.na(PW1)==TRUE, 0, PW1),
#                BW1 = ifelse(is.na(BW1)==TRUE, 0, BW1))

# Common channels across all datasets
pathania_channels <- c(levels(PATHANIA$Channel))
euler_channels <- c(levels(EULER$Channel))
mpi_channels <- c(levels(MPI$Channel))

# Keep the pathania channels that are also in the euler channels
keep <- pathania_channels[pathania_channels %in% euler_channels]

# Keep the shared channels that are also in the mpi channels
keep <- keep[keep %in% mpi_channels]

# We will filter out the unwanted channels later. 



# Combining the cleaned data into one ----
colnames(PATHANIA)
colnames(EULER)
colnames(MPI)


CLEAN <- rbind(PATHANIA %>% select(subID, file_id, dataset, Channels, r.2, Error,
                                Offset, Exponent, CF1, PW1, BW1, CF2, PW2, 
                                BW2, CF3, PW3, BW3, CF4, PW4, BW4, CF5, PW5,
                                BW5, CF6, PW6, BW6, age, sex), 
               EULER%>% select(subID, file_id, dataset, Channel, r.2, Error,
                               Offset, Exponent, CF1, PW1, BW1, CF2, PW2, 
                               BW2, CF3, PW3, BW3, CF4, PW4, BW4, CF5, PW5,
                               BW5, CF6, PW6, BW6, age, sex) %>%
                 rename(Channels=Channel), 
               MPI%>% select(subID, file_id, dataset, Channels, r.2, Error,
                             Offset, Exponent, CF1, PW1, BW1, CF2, PW2, 
                             BW2, CF3, PW3, BW3, CF4, PW4, BW4, CF5, PW5,
                             BW5, CF6, PW6, BW6, age, sex))
head(CLEAN)
summary(CLEAN$Channels)
summary(CLEAN$dataset)
CLEAN <- CLEAN %>% filter(Channels %in% keep) %>% mutate(Channels = factor(Channels))
summary(CLEAN$Channels)

summary(factor(CLEAN$Channel))

getwd()
write.csv(CLEAN, "./data_CLEAN_COMBINED.csv")
