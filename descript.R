setwd(getwd())

#load packages
library(readxl)
library(tidyverse)

#load data
data <- read_excel("data/data.xlsx")
ma_id <- read_excel("data/info_ma.xlsx")

#modify Probably Low (PL) and Probably High (PH) RoB to either low (L) or high risk of bias (H)
data$BLINDING_PATIENTS_PERSONNEL <- ifelse(data$BLINDING_PATIENTS_PERSONNEL == 'PH',"H", data$BLINDING_PATIENTS_PERSONNEL)
data$BLINDING_PATIENTS_PERSONNEL <- ifelse(data$BLINDING_PATIENTS_PERSONNEL == 'PL',"L", data$BLINDING_PATIENTS_PERSONNEL)

data$RANDOM_SEQUENCE_GENERATION <- ifelse(data$RANDOM_SEQUENCE_GENERATION == 'PH',"H", data$RANDOM_SEQUENCE_GENERATION)
data$RANDOM_SEQUENCE_GENERATION <- ifelse(data$RANDOM_SEQUENCE_GENERATION == 'PL',"L", data$RANDOM_SEQUENCE_GENERATION)

data$ALLOCATION_CONCEALMENT <- ifelse(data$ALLOCATION_CONCEALMENT == 'PH',"H", data$ALLOCATION_CONCEALMENT)
data$ALLOCATION_CONCEALMENT <- ifelse(data$ALLOCATION_CONCEALMENT == 'PL',"L", data$ALLOCATION_CONCEALMENT)

data$INCOMPLETE_DATA <- ifelse(data$INCOMPLETE_DATA == 'PH',"H", data$INCOMPLETE_DATA)
data$INCOMPLETE_DATA <- ifelse(data$INCOMPLETE_DATA == 'PL',"L", data$INCOMPLETE_DATA)

#modify unavailable or unclear blinding manual assessment to high RoB (H)
data$BLINDING_MANUAL <- ifelse(data$BLINDING_MANUAL == 'UNA',"U", data$BLINDING_MANUAL)
data$BLINDING_MANUAL <- ifelse(data$BLINDING_MANUAL == 'U',"H", data$BLINDING_MANUAL)

# modify Unclear RoB (U) to high (H)
data$BLINDING_PATIENTS_PERSONNEL <- ifelse(data$BLINDING_PATIENTS_PERSONNEL == 'U',"H", data$BLINDING_PATIENTS_PERSONNEL)
data$RANDOM_SEQUENCE_GENERATION <- ifelse(data$RANDOM_SEQUENCE_GENERATION == 'U',"H", data$RANDOM_SEQUENCE_GENERATION)
data$ALLOCATION_CONCEALMENT <- ifelse(data$ALLOCATION_CONCEALMENT == 'U',"H", data$ALLOCATION_CONCEALMENT)
data$INCOMPLETE_DATA <- ifelse(data$INCOMPLETE_DATA == 'U',"H", data$INCOMPLETE_DATA)

# sum total patients
data %>% 
  mutate(sum = SUBJECTS_EXP + SUBJECTS_CONTROL) %>%
  select(sum) %>% 
  sum()

# summarize patient counts per RCTs according to blinding
count_L <- data %>% 
  mutate(sum = SUBJECTS_EXP + SUBJECTS_CONTROL) %>% 
  select(BLINDING_MANUAL, sum) %>% 
  filter(BLINDING_MANUAL == "L") 
summary(count_L$sum)

count_H <- data %>% 
  mutate(sum = SUBJECTS_EXP + SUBJECTS_CONTROL) %>% 
  select(BLINDING_MANUAL, sum) %>% 
  filter(BLINDING_MANUAL == "H") 
summary(count_H$sum)

# summarize patient counts per meta-analysis
count_MA <- data %>% 
  group_by(ID_RS) %>% 
  mutate(sum = sum(SUBJECTS_EXP + SUBJECTS_CONTROL)) %>% 
  group_by(ID_RS) %>% 
  select(ID_RS, sum) %>% 
  distinct() %>% 
  filter(!is.na(sum))

min <- data %>% 
  group_by(ID_RS) %>% 
  mutate(x_min = min(SUBJECTS_EXP + SUBJECTS_CONTROL)) %>% 
  select(ID_RS, x_min) %>% 
  distinct() %>% 
  filter(!is.na(x_min))

med <- data %>% 
  group_by(ID_RS) %>% 
  mutate(x_med = median(SUBJECTS_EXP + SUBJECTS_CONTROL)) %>% 
  select(ID_RS, x_med) %>% 
  distinct() %>% 
  filter(!is.na(x_med))

max <- data %>% 
  group_by(ID_RS) %>% 
  mutate(x_max = max(SUBJECTS_EXP + SUBJECTS_CONTROL)) %>% 
  select(ID_RS, x_max) %>% 
  distinct() %>% 
  filter(!is.na(x_max))

table <- 
  left_join(count_MA, max) %>% 
  left_join(min) %>% 
  left_join(med) %>% 
  mutate(num = str_c(sum, x_max, sep = " ("),
         num = str_c(num, x_med, sep = ", "), 
         num = str_c(num, x_min, sep = ", "),
         num = str_c(num, "", sep = ")"))

# extract meta-analysis Author + reference in the paper
table$studlab <- as.character(ma_id$Author)

# count RCTs per meta-analysis
count_RCT <- data %>% group_by(ID_RS) %>% count()
summary(count_RCT)
# number of blinded and non-blinded studies
data %>% group_by(BLINDING_MANUAL) %>% 
  count()

# counts per publication date
data %>% group_by(BLINDING_MANUAL) %>% 
  count(YEAR <= 2019 & YEAR >= 2010) %>% 
  mutate(n_percent = ifelse(BLINDING_MANUAL == "H", (n/200)*100, (n/267)*100)) %>% 
  filter(`YEAR <= 2019 & YEAR >= 2010` == TRUE)

data %>% group_by(BLINDING_MANUAL) %>% 
  count(YEAR <= 2009 & YEAR >= 2000) %>% 
  mutate(n_percent = ifelse(BLINDING_MANUAL == "H", (n/200)*100, (n/267)*100)) %>% 
  filter(`YEAR <= 2009 & YEAR >= 2000` == TRUE)

data %>% group_by(BLINDING_MANUAL) %>% 
  count(YEAR <= 1999 & YEAR >= 1990) %>% 
  mutate(n_percent = ifelse(BLINDING_MANUAL == "H", (n/200)*100, (n/267)*100)) %>% 
  filter(`YEAR <= 1999 & YEAR >= 1990` == TRUE)

data %>% group_by(BLINDING_MANUAL) %>% 
  count(YEAR <= 1989) %>% 
  mutate(n_percent = ifelse(BLINDING_MANUAL == "H", (n/200)*100, (n/267)*100)) %>% 
  filter(`YEAR <= 1989` == TRUE)

#median of publication per group
data %>% group_by(BLINDING_MANUAL) %>% 
  summarise(Mean=mean(YEAR), Max=max(YEAR), Min=min(YEAR), Median=median(YEAR), Std=sd(YEAR))

# options(pillar.sigfig = 6)
            
# counts per sample size / blinded
data %>% 
  mutate(sample_size = SUBJECTS_EXP + SUBJECTS_CONTROL) %>% 
  filter(BLINDING_MANUAL == "L") %>% 
  count(sample_size < 100) %>% 
  mutate(n_percent = (n/267)*100)

data %>% 
  mutate(sample_size = SUBJECTS_EXP + SUBJECTS_CONTROL) %>% 
  filter(BLINDING_MANUAL == "L") %>% 
  count(sample_size < 200 & sample_size >= 100) %>% 
  mutate(n_percent = (n/267)*100)

data %>% 
  mutate(sample_size = SUBJECTS_EXP + SUBJECTS_CONTROL) %>% 
  filter(BLINDING_MANUAL == "L") %>% 
  count(sample_size >= 200) %>% 
  mutate(n_percent = (n/267)*100)

# counts per sample size / non-blinded
data %>% 
  mutate(sample_size = SUBJECTS_EXP + SUBJECTS_CONTROL) %>% 
  filter(BLINDING_MANUAL == "H") %>% 
  count(sample_size < 100) %>% 
  mutate(n_percent = (n/200)*100)

data %>% 
  mutate(sample_size = SUBJECTS_EXP + SUBJECTS_CONTROL) %>% 
  filter(BLINDING_MANUAL == "H") %>% 
  count(sample_size < 200 & sample_size >= 100) %>% 
  mutate(n_percent = (n/200)*100)

data %>% 
  mutate(sample_size = SUBJECTS_EXP + SUBJECTS_CONTROL) %>% 
  filter(BLINDING_MANUAL == "H") %>% 
  count(sample_size >= 200) %>% 
  mutate(n_percent = (n/200)*100)


# counts per bias
data %>% group_by(BLINDING_MANUAL) %>% 
  count(RANDOM_SEQUENCE_GENERATION) %>% 
  mutate(n_percent = ifelse(BLINDING_MANUAL == "H", (n/200)*100, (n/267)*100)) %>% 
  filter(!is.na(BLINDING_MANUAL))

data %>% group_by(BLINDING_MANUAL) %>% 
  count(ALLOCATION_CONCEALMENT) %>% 
  mutate(n_percent = ifelse(BLINDING_MANUAL == "H", (n/200)*100, (n/267)*100)) %>% 
  filter(!is.na(BLINDING_MANUAL))

data %>% group_by(BLINDING_MANUAL) %>% 
  count(INCOMPLETE_DATA) %>% 
  mutate(n_percent = ifelse(BLINDING_MANUAL == "H", (n/200)*100, (n/267)*100)) %>% 
  filter(!is.na(BLINDING_MANUAL))

data %>% group_by(BLINDING_MANUAL) %>% 
  count(BLINDING_PATIENTS_PERSONNEL) %>% 
  mutate(n_percent = ifelse(BLINDING_MANUAL == "H", (n/200)*100, (n/267)*100)) %>% 
  filter(!is.na(BLINDING_MANUAL))

