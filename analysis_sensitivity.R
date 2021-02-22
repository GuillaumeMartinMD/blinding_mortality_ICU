setwd(getwd())

library(meta)
library(readxl)
library(tidyverse)

## Read data
data <- read_excel("data/data.xlsx")

## Convert Probably Low to Low Risk (L) and Probably High and Unclear to High Risk (H)
# and convert  Low to letter "A" for ranking reasons in R
data$BLINDING_PATIENTS_PERSONNEL <- ifelse(data$BLINDING_PATIENTS_PERSONNEL == 'PH',"H", data$BLINDING_PATIENTS_PERSONNEL)
data$BLINDING_PATIENTS_PERSONNEL <- ifelse(data$BLINDING_PATIENTS_PERSONNEL == 'PL',"L", data$BLINDING_PATIENTS_PERSONNEL)
data$BLINDING_PATIENTS_PERSONNEL <- ifelse(data$BLINDING_PATIENTS_PERSONNEL == 'U',"H", data$BLINDING_PATIENTS_PERSONNEL)
data$BLINDING_PATIENTS_PERSONNEL <- ifelse(data$BLINDING_PATIENTS_PERSONNEL == 'L',"A", data$BLINDING_PATIENTS_PERSONNEL)

data$RANDOM_SEQUENCE_GENERATION <- ifelse(data$RANDOM_SEQUENCE_GENERATION == 'PH',"H", data$RANDOM_SEQUENCE_GENERATION)
data$RANDOM_SEQUENCE_GENERATION <- ifelse(data$RANDOM_SEQUENCE_GENERATION == 'PL',"L", data$RANDOM_SEQUENCE_GENERATION)
data$RANDOM_SEQUENCE_GENERATION <- ifelse(data$RANDOM_SEQUENCE_GENERATION == 'U',"H", data$RANDOM_SEQUENCE_GENERATION)
data$RANDOM_SEQUENCE_GENERATION <- ifelse(data$RANDOM_SEQUENCE_GENERATION == 'L',"A", data$RANDOM_SEQUENCE_GENERATION)

data$ALLOCATION_CONCEALMENT <- ifelse(data$ALLOCATION_CONCEALMENT == 'PH',"H", data$ALLOCATION_CONCEALMENT)
data$ALLOCATION_CONCEALMENT <- ifelse(data$ALLOCATION_CONCEALMENT == 'PL',"L", data$ALLOCATION_CONCEALMENT)
data$ALLOCATION_CONCEALMENT <- ifelse(data$ALLOCATION_CONCEALMENT == 'U',"H", data$ALLOCATION_CONCEALMENT)
data$ALLOCATION_CONCEALMENT <- ifelse(data$ALLOCATION_CONCEALMENT == 'L',"A", data$ALLOCATION_CONCEALMENT)

data$INCOMPLETE_DATA <- ifelse(data$INCOMPLETE_DATA == 'PH',"H", data$INCOMPLETE_DATA)
data$INCOMPLETE_DATA <- ifelse(data$INCOMPLETE_DATA == 'PL',"L", data$INCOMPLETE_DATA)
data$INCOMPLETE_DATA <- ifelse(data$INCOMPLETE_DATA == 'U',"H", data$INCOMPLETE_DATA)
data$INCOMPLETE_DATA <- ifelse(data$INCOMPLETE_DATA == 'L',"A", data$INCOMPLETE_DATA)

## Convert "NA" to "Unclear", then Unclear (U) to High (H), and Low to letter "A" for ranking reasons in R 
data$BLINDING_MANUAL <- ifelse(data$BLINDING_MANUAL == 'UNA',"U", data$BLINDING_MANUAL)
data$BLINDING_MANUAL <- ifelse(data$BLINDING_MANUAL == 'U',"H", data$BLINDING_MANUAL)
data$BLINDING_MANUAL <- ifelse(data$BLINDING_MANUAL == 'L',"A", data$BLINDING_MANUAL)

## Convert data as numeric
data <- as.data.frame(data)
data$ID_RS <- as.numeric(data$ID_RS)
data$OUTCOME_EXP <- as.numeric(data$OUTCOME_EXP)
data$SUBJECTS_EXP <- as.numeric(data$SUBJECTS_EXP)
data$OUTCOME_CTRL <- as.numeric(data$OUTCOME_CTRL)
data$SUBJECTS_CONTROL <- as.numeric(data$SUBJECTS_CONTROL)

## Create dataframe with needed data for the meta-analysis, by group
data_MA <- data %>% group_by(ID_RS) %>% 
    mutate(event.e = OUTCOME_EXP, n.e = SUBJECTS_EXP, event.c = OUTCOME_CTRL, n.c = SUBJECTS_CONTROL, 
           blind = BLINDING_MANUAL, sequence = RANDOM_SEQUENCE_GENERATION, 
           allocation = ALLOCATION_CONCEALMENT, incomplete = INCOMPLETE_DATA, 
           sample = SUBJECTS_EXP + SUBJECTS_CONTROL) %>% 
    select(event.e, n.e, event.c, n.c, ID_RS, blind, sequence, allocation, incomplete, sample)

data_MA <- data_MA %>% group_by(ID_RS) %>% 
  filter (n() > 3) ## remove meta-analysis with < 4 RCT (2 and 24)


#### metabin loop to recompute meta_analysis OR and adjust on RoB + sample ####

## Create vectors for loop
n <- 34 #number of meta_analysis

logROR.reg.sequence <- vector('numeric', n)
logSE.reg.sequence <- vector('numeric', n)

logROR.reg.allocation <- vector('numeric', n)
logSE.reg.allocation <- vector('numeric', n)

logROR.reg.incomplete <- vector('numeric', n)
logSE.reg.incomplete <- vector('numeric', n)

logROR.reg.sample <- vector('numeric', n)
logSE.reg.sample <- vector('numeric', n)

##metabin loop
for (i in 1:n){
  assign(paste0("df", i), as.data.frame(split(data_MA, data_MA$ID_RS)[[i]])) ##split dataframe in n df according to the n groups, with df[i] as name
  event.e <- get(paste0("df", i))[, "event.e"] ##extract data and assign to variable to be used in metabin
  n.e <- get(paste0("df", i))[, "n.e"]
  event.c <- get(paste0("df", i))[, "event.c"]
  n.c <- get(paste0("df", i))[, "n.c"]
  sequence <- get(paste0("df", i))[, "sequence"]
  allocation <- get(paste0("df", i))[, "allocation"]
  incomplete <- get(paste0("df", i))[, "incomplete"]
  sample <- get(paste0("df", i))[, "sample"]
  
  meta <- metabin(event.e, n.e, event.c, n.c, data = get(paste0("df", i)), ##metabin function
            method="MH",
            sm="OR",
            incr="TACC", allincr=F,
            addincr=F, allstudies=F,
            MH.exact=F, RR.cochrane=gs("RR.cochrane"),
            level=gs("level"), level.comb=gs("level.comb"),
            comb.fixed=F, comb.random=T,
            hakn=F,
            method.tau="DL",
            tau.preset=NULL, TE.tau=NULL,
            method.bias="peters",
            backtransf=T, pscale = 1,
            print.CMH=T,
            warn=gs("warn"))
  
## adjustement on sequence generation
    if (length(unique((get(paste0("df", i))$sequence))) == 2 & length(get(paste0("df", i))$sequence) > 3)
  {
    logROR.reg.sequence[i] <- metareg(x = meta, blind + sequence, method.tau = "REML")$b[2]
    logSE.reg.sequence[i] <- metareg(x = meta, blind + sequence, method.tau = "REML")$se[2]
  }
  reslt.reg_sequence <- cbind(logROR.reg.sequence, logSE.reg.sequence)
  

## adjustement on allocation concealment
      if (length(unique((get(paste0("df", i))$allocation))) == 2 & length(get(paste0("df", i))$allocation) > 3)
  {
    logROR.reg.allocation[i] <- metareg(x = meta, blind + allocation, method.tau = "REML")$b[2]
    logSE.reg.allocation[i] <- metareg(x = meta, blind + allocation, method.tau = "REML")$se[2]
  }
  reslt.reg_allocation <- cbind(logROR.reg.allocation, logSE.reg.allocation)
  
## adjustement on incomplete data  
      if (length(unique((get(paste0("df", i))$incomplete))) == 2 & length(get(paste0("df", i))$incomplete) > 3)
  {
    logROR.reg.incomplete[i] <- metareg(x = meta, blind + incomplete, method.tau = "REML")$b[2]
    logSE.reg.incomplete[i] <- metareg(x = meta, blind + incomplete, method.tau = "REML")$se[2]
  }
  reslt.reg_incomplete <- cbind(logROR.reg.incomplete, logSE.reg.incomplete)

    ## adjustement on sample size
  logROR.reg.sample[i] <- metareg(x = meta, blind + sample, method.tau = "REML")$b[2]
  logSE.reg.sample[i] <- metareg(x = meta, blind + sample, method.tau = "REML")$se[2]
  reslt.reg_sample <- cbind(logROR.reg.sample, logSE.reg.sample)
  
}

##convert to data.frame
reslt.reg_sequence <- as.data.frame(reslt.reg_sequence)
reslt.reg_allocation <- as.data.frame(reslt.reg_allocation)
reslt.reg_incomplete <- as.data.frame(reslt.reg_incomplete)
reslt.reg_sample <- as.data.frame(reslt.reg_sample)

#### meta-epidemiological analysis ####

## adjustment on sequence generation
metaepi_sequence <- metagen(reslt.reg_sequence$logROR.reg.sequence, 
                            reslt.reg_sequence$logSE.reg.sequence, 
                            data = reslt.reg_sequence, 
                            method.tau="DL", sm = "OR", 
                            backtransf = T, hakn = F, comb.fixed=F)

## adjustment on allocation concealment
metaepi_allocation <- metagen(reslt.reg_allocation$logROR.reg.allocation,
                              reslt.reg_allocation$logSE.reg.allocation, 
                              data = reslt.reg_allocation, 
                              method.tau="DL", sm = "OR", 
                              backtransf = T, hakn = F, comb.fixed=F)

## adjustment on incomplete data
metaepi_incomplete <- metagen(reslt.reg_incomplete$logROR.reg.incomplete, 
                              reslt.reg_incomplete$logSE.reg.incomplete, 
                              data = reslt.reg_incomplete, 
                              method.tau="DL", sm = "OR", 
                              backtransf = T, hakn = F, comb.fixed=F)

## adjustment on sample size
metaepi_sample <- metagen(reslt.reg_sample$logROR.reg.sample, 
                          reslt.reg_sample$logSE.reg.sample, 
                          data = reslt.reg_sample, 
                          method.tau="DL", sm = "OR", 
                          backtransf = T, hakn = F, comb.fixed=F)

metaepi_sample$k <- paste0(n) # trick to remove .00 after number of studies

metaepi_sensitivity <-  metabind(metaepi_sequence, #bind multiple analysis
                                 metaepi_allocation, 
                                 metaepi_incomplete, 
                                 metaepi_sample,
                                 name = c("Random sequence generation",
                                          "Allocation concealment", 
                                          "Incomplete outcome data", 
                                          "Sample size"))

metaepi_sensitivity$studlab <- c("Random sequence generation",
                                 "Allocation concealment", 
                                 "Incomplete outcome data", 
                                 "Sample size")
## forest plot
forest(metaepi_sensitivity,
       studlab = T,
       smlab = "",
       digits = 2,
       rightcols = c("effect", "ci"), rightlabs = c("ROR", "CI (95%)"),
       leftlabs = c("Adjustment on :", "Number of contributive\n meta-analyses"), 
       leftcols = c("byvar", "k"),
       # just = "left",
       label.right="Non-blinded less beneficial", col.label.right="black",
       label.left="Non-blinded more beneficial", col.label.left="black",
       squaresize = 0.8, 
       colgap = "1mm",
       col.square = "cornflowerblue", col.inside = "black",  
       col.by = "white",
       xlab.pos = -0.1, 
       fs.study=10, fs.axis = 10,
       fs.study.label=10, fs.heading=10, ff.heading="bold", plotwidth = "8cm",
       ff.axis="plain", ff.random = "bold",  ff.xlab = "plain",
       ff.smlab="plain",  ff.study="plain", ff.study.label="plain",
       ff.fixed="plain", ff.hetstat="plain"
       )

dev.print(device = pdf, 
          file = "output/sensitivity_forest.pdf", 
          width = 10, height = 5)

dev.print(device = tiff, res=300,
          file = "output/sensitivity_forest.tiff", 
          width = 2700, height = 1400)

