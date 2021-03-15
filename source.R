## Tested with R 3.6.0 or MRAN (R Open) 4.0.2

setwd(getwd()) #your working directory with all R files

#Load Packages Snapshots from the Checkpoint Server for Reproducibility
library(checkpoint) 
checkpoint("2020-10-11") #CRAN Snapshot date used

library(meta) #General Package for Meta-Analysis 
library(readxl) #Read Excel Files
library(writexl) #Write Excel Files
library(tidyverse) #General Data Manipulation/Analysis
 