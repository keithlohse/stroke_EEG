# Title: Brain Age - Merging Individual Participant Data
# Author: Keith Lohse, PhD, PStat
# Date: 2023-10-268
library(tidyverse)

getwd()
setwd("C:/Users/lohse/Box/Lohse Brain Age EEG/")

list.files()
list.files("./data/processed EEG data/")
list.files("./data/processed EEG data/2. Asher Processed Pathania FOOOF DATA/")

# 1.0 Merging Pathania Data ----
setwd("./data/processed EEG data/2. Asher Processed Pathania FOOOF DATA/")
file_names <- list.files()
file_names

file_names<-file_names[grep("fooof_results.csv", file_names)]



# filling in column names for missing variables
var_list <- c("Channels",	"r.2",	"Error", "Offset" , "Exponent",
              "CF1",	"PW1",	"BW1",	"CF2",	"PW2",	"BW2",
              "CF3",	"PW3",	"BW3",	"CF4",	"PW4",	"BW4",
              "CF5",	"PW5",	"BW5",  "CF6",	"PW6",	"BW6")

subject <- read.csv("./oa22_ec_fooof_results.csv",
                    header=TRUE, 
                    stringsAsFactors = TRUE)
head(subject)
var_list %in% names(subject)
needed <- var_list[var_list %in% names(subject) == FALSE]
needed

for (i in needed) {
  print(i)
  subject$New_Col <- c(NA) 
  names(subject)[names(subject) == 'New_Col'] <- eval(i)
}

subject



# Putting an if else statement inside of our for-loop
for(name in file_names) {
  print(name)
  subject <- read.csv(name,
                      header=TRUE, 
                      stringsAsFactors = TRUE)
  
  needed <- var_list[var_list %in% names(subject) == FALSE]
  
  for (i in needed) {
    print(i)
    subject$New_Col <- c(NA) 
    names(subject)[names(subject) == 'New_Col'] <- eval(i)
  }
  
  # if the MASTER data set doesn't exist, create it
  if (!exists("MASTER")){
    MASTER <- data.frame(subject)
    MASTER$file_id <- name
    
    
  } else {
    # Create the temporary data file:
    temp_dataset <- data.frame(subject)
    temp_dataset$file_id <-  name
    
    MASTER<-rbind(MASTER, temp_dataset)
    
    # Remove or "empty" the temporary data set
    rm(temp_dataset)
  }
}

head(MASTER)

# move the file ID and Hz columns to the front of the dataset
MASTER <- MASTER %>% relocate(file_id)
head(MASTER)

# Break the file id at the underscore
str_split(MASTER$file_id, "_")[[1]]
str_sub(map_chr(str_split(MASTER$file_id, "_"), 1), start=-4, end=-1)
map_chr(str_split(MASTER$file_id, "_"), 2)

MASTER$subID <- str_sub(map_chr(str_split(MASTER$file_id, "_"), 1), start=-4, end=-1)
MASTER$dataset <- c("pathania")
MASTER$condition <- map_chr(str_split(MASTER$file_id, "_"), 2)

MASTER <- MASTER %>% relocate(file_id, subID, dataset, condition) %>%
  filter(condition == "eo") %>%
  select(-condition)
head(MASTER)

# Merging in the demographic data ----
setwd("C:/Users/lohse/Box/Lohse Brain Age EEG/data/")
list.files()

DEMO <- read.csv("./RBANS_aging_study_08262021.csv",
                    header=TRUE, 
                    stringsAsFactors = TRUE)
head(DEMO)
DEMO <- DEMO %>% select(subID, age, sex)


summary(DEMO$subID)
summary(MASTER$subID)
MASTER$subID <- factor(MASTER$subID)
summary(MASTER$subID)

MASTER <- merge(x=MASTER, y=DEMO, by="subID")
head(MASTER)

# Recode the sex variable to match other datasets
MASTER$sex <- ifelse(MASTER$sex=="f", "F", "M")

getwd()
write.csv(MASTER, "pathania_MERGED.csv")

rm(list = ls())






# 2.0 Merging Euler Data ----
setwd("C:/Users/lohse/Box/Lohse Brain Age EEG/")
setwd("./data/processed EEG data/euler/")
file_names <- list.files()
file_names

var_list <- c("Channels",	"r.2",	"Error", "Offset" , "Exponent",
              "CF1",	"PW1",	"BW1",	"CF2",	"PW2",	"BW2",
              "CF3",	"PW3",	"BW3",	"CF4",	"PW4",	"BW4",
              "CF5",	"PW5",	"BW5",  "CF6",	"PW6",	"BW6")


# Putting an if else statement inside of our for-lopp
for(name in file_names) {
  print(name)
  subject <- read.csv(name,
                      header=TRUE, 
                      stringsAsFactors = TRUE)
  
  
  needed <- var_list[var_list %in% names(subject) == FALSE]
  
  for (i in needed) {
    print(i)
    subject$New_Col <- c(NA) 
    names(subject)[names(subject) == 'New_Col'] <- eval(i)
  }
  
  # if the MASTER data set doesn't exist, create it
  if (!exists("MASTER")){
    MASTER <- data.frame(subject)
    MASTER$file_id <- name
    
    
  } else {
    # Create the temporary data file:
    temp_dataset <- data.frame(subject)
    temp_dataset$file_id <-  name
    
    MASTER<-rbind(MASTER, temp_dataset)
    
    # Remove or "empty" the temporary data set
    rm(temp_dataset)
  }
}

head(MASTER)

# move the file ID and Hz columns to the front of the dataset
MASTER <- MASTER %>% relocate(file_id)
head(MASTER)

# Break the file id at the underscore
str_split(MASTER$file_id, "_")[[1]]
map_chr(str_split(map_chr(str_split(MASTER$file_id, "_"), 1), "channel"), 2)

MASTER$subID <- map_chr(str_split(map_chr(str_split(MASTER$file_id, "_"), 1), "channel"), 2)
MASTER$dataset <- c("euler")

MASTER <- MASTER %>% relocate(file_id, subID, dataset)
head(MASTER)

# Merging in the demographic data ----
setwd("C:/Users/lohse/Box/Lohse Brain Age EEG/data/")
list.files()

DEMO <- read.csv("./euler_demographic_data.csv",
                 header=TRUE, 
                 stringsAsFactors = TRUE)
head(DEMO)
DEMO <- DEMO %>% select(Participant, dmage, dmsex) %>%
  rename(subID = Participant) %>%
  rename(age = dmage) %>%
  rename(sex = dmsex)


summary(DEMO$subID)
summary(MASTER$subID)
MASTER$subID <- factor(MASTER$subID)
summary(MASTER$subID)

MASTER <- merge(x=MASTER, y=DEMO, by="subID")
head(MASTER)

write.csv(MASTER, "euler_MERGED.csv")
rm(list = ls())






# 3.0 Merging MPI Data ----
setwd("C:/Users/lohse/Box/Lohse Brain Age EEG/")
list.files("./data/processed EEG data/")
setwd("./data/processed EEG data/1. Asher Processed MPI FOOOF DATA/")
file_names <- list.files()
file_names

# Extract only those filenames that match the pattern:
file_names<-file_names[grep("fooof_results.csv", file_names)]
file_names <- file_names[-1]

var_list <- c("Channels",	"r.2",	"Error", "Offset" , "Exponent",
              "CF1",	"PW1",	"BW1",	"CF2",	"PW2",	"BW2",
              "CF3",	"PW3",	"BW3",	"CF4",	"PW4",	"BW4",
              "CF5",	"PW5",	"BW5",  "CF6",	"PW6",	"BW6")


# Putting an if else statement inside of our for-lopp
for(name in file_names) {
  print(name)
  subject <- read.csv(name,
                      header=TRUE, 
                      stringsAsFactors = TRUE)
  
  needed <- var_list[var_list %in% names(subject) == FALSE]
  
  for (i in needed) {
    print(i)
    subject$New_Col <- c(NA) 
    names(subject)[names(subject) == 'New_Col'] <- eval(i)
  }
  
  # if the MASTER data set doesn't exist, create it
  if (!exists("MASTER")){
    MASTER <- data.frame(subject)
    MASTER$file_id <- name
    
    
  } else {
    # Create the temporary data file:
    temp_dataset <- data.frame(subject)
    temp_dataset$file_id <-  name
    
    MASTER<-rbind(MASTER, temp_dataset)
    
    # Remove or "empty" the temporary data set
    rm(temp_dataset)
  }
}

head(MASTER)

# move the file ID and Hz columns to the front of the dataset
MASTER <- MASTER %>% relocate(file_id)
head(MASTER)

# Break the file id at the underscore
str_split(MASTER$file_id, "_")[[1]]
map_chr(str_split(MASTER$file_id, "_"), 1)

MASTER$subID <- map_chr(str_split(MASTER$file_id, "_"), 1)

MASTER$dataset <- c("mpi")

MASTER <- MASTER %>% relocate(file_id, subID, dataset)
head(MASTER)

# Merging in the demographic data ----
setwd("C:/Users/lohse/Box/Lohse Brain Age EEG/data/")
list.files()

DEMO <- read.csv("./mpi_demographic.csv",
                 header=TRUE, 
                 stringsAsFactors = TRUE)
head(DEMO)
DEMO <- DEMO %>% select(subID, gender, Age) %>%
  rename(subID = subID) %>%
  rename(age = Age) %>%
  rename(sex = gender)


summary(DEMO$subID)
summary(MASTER$subID)
MASTER$subID <- factor(MASTER$subID)
summary(MASTER$subID)

MASTER <- merge(x=MASTER, y=DEMO, by="subID")
head(MASTER)

write.csv(MASTER, "mpi_MERGED.csv")

rm(list = ls())










# 4.0 Merging Cramer Data ----
setwd("./processed EEG data/3. Asher Processed CRAMER FOOOF Data/")
file_names <- list.files()
file_names

file_names<-file_names[grep("fooof_results.csv", file_names)]
file_names<-file_names[-1]


# filling in column names for missing variables
var_list <- c("Channels",	"r.2",	"Error", "Offset" , "Exponent",
              "CF1",	"PW1",	"BW1",	"CF2",	"PW2",	"BW2",
              "CF3",	"PW3",	"BW3",	"CF4",	"PW4",	"BW4",
              "CF5",	"PW5",	"BW5",  "CF6",	"PW6",	"BW6")


# Putting an if else statement inside of our for-loop
for(name in file_names) {
  print(name)
  subject <- read.csv(name,
                      header=TRUE, 
                      stringsAsFactors = TRUE)
  
  needed <- var_list[var_list %in% names(subject) == FALSE]
  
  for (i in needed) {
    print(i)
    subject$New_Col <- c(NA) 
    names(subject)[names(subject) == 'New_Col'] <- eval(i)
  }
  
  # if the MASTER data set doesn't exist, create it
  if (!exists("MASTER")){
    MASTER <- data.frame(subject)
    MASTER$file_id <- name
    
    
  } else {
    # Create the temporary data file:
    temp_dataset <- data.frame(subject)
    temp_dataset$file_id <-  name
    
    MASTER<-rbind(MASTER, temp_dataset)
    
    # Remove or "empty" the temporary data set
    rm(temp_dataset)
  }
}

head(MASTER)

# move the file ID and Hz columns to the front of the dataset
MASTER <- MASTER %>% relocate(file_id)
head(MASTER)

# Break the file id at the underscore
str_split(MASTER$file_id, "_")[[1]]
str_sub(map_chr(str_split(MASTER$file_id, "_"), 1), start=-4, end=-1)
map_chr(str_split(MASTER$file_id, "_"), 2)

MASTER$subID <- map_chr(str_split(MASTER$file_id, "_"), 1)
MASTER$dataset <- c("cramer")

MASTER <- MASTER %>% relocate(file_id, subID, dataset)
head(MASTER)

# # Merging in the demographic data ----
setwd("C:/Users/lohse/Box/Lohse Brain Age EEG/data/")
list.files()
# 
DEMO <- read.csv("./cramer_demographic_data.csv",
                 header=TRUE,
                 stringsAsFactors = TRUE) %>% 
  filter(is.na(EEG_Directory_Name)==FALSE)

head(DEMO)
summary(DEMO$subID)
summary(DEMO$gender)
summary(DEMO$age)


summary(MASTER$subID)
MASTER$subID <- factor(MASTER$subID)
summary(MASTER$subID)

MASTER <- merge(x=MASTER, y=DEMO, by="subID")
head(MASTER)

getwd()
write.csv(MASTER, "cramer_MERGED.csv")
 
rm(list = ls())
