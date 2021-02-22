setwd(getwd())

library(meta)
library(readxl)
library(writexl)
library(tidyverse)

## Read data
data <- read_excel("data/data.xlsx")
ma_id <- read_excel("data/info_MA.xlsx")

## Convert Probably Low to Low Risk (L) and Probably High and Unclear to High Risk (H)
data$BLINDING_PATIENTS_PERSONNEL <- ifelse(data$BLINDING_PATIENTS_PERSONNEL == 'PH',"H", data$BLINDING_PATIENTS_PERSONNEL)
data$BLINDING_PATIENTS_PERSONNEL <- ifelse(data$BLINDING_PATIENTS_PERSONNEL == 'PL',"L", data$BLINDING_PATIENTS_PERSONNEL)
data$BLINDING_PATIENTS_PERSONNEL <- ifelse(data$BLINDING_PATIENTS_PERSONNEL == 'U',"H", data$BLINDING_PATIENTS_PERSONNEL)

data$RANDOM_SEQUENCE_GENERATION <- ifelse(data$RANDOM_SEQUENCE_GENERATION == 'PH',"H", data$RANDOM_SEQUENCE_GENERATION)
data$RANDOM_SEQUENCE_GENERATION <- ifelse(data$RANDOM_SEQUENCE_GENERATION == 'PL',"L", data$RANDOM_SEQUENCE_GENERATION)
data$RANDOM_SEQUENCE_GENERATION <- ifelse(data$RANDOM_SEQUENCE_GENERATION == 'U',"H", data$RANDOM_SEQUENCE_GENERATION)

data$ALLOCATION_CONCEALMENT <- ifelse(data$ALLOCATION_CONCEALMENT == 'PH',"H", data$ALLOCATION_CONCEALMENT)
data$ALLOCATION_CONCEALMENT <- ifelse(data$ALLOCATION_CONCEALMENT == 'PL',"L", data$ALLOCATION_CONCEALMENT)
data$ALLOCATION_CONCEALMENT <- ifelse(data$ALLOCATION_CONCEALMENT == 'U',"H", data$ALLOCATION_CONCEALMENT)

data$INCOMPLETE_DATA <- ifelse(data$INCOMPLETE_DATA == 'PH',"H", data$INCOMPLETE_DATA)
data$INCOMPLETE_DATA <- ifelse(data$INCOMPLETE_DATA == 'PL',"L", data$INCOMPLETE_DATA)
data$INCOMPLETE_DATA <- ifelse(data$INCOMPLETE_DATA == 'U',"H", data$INCOMPLETE_DATA)

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
           allocation = ALLOCATION_CONCEALMENT, 
           incomplete = INCOMPLETE_DATA) %>% 
    select(event.e, n.e, event.c, n.c, ID_RS, blind, sequence, allocation,
           incomplete, SHORT_TERM)

#### metabin loop to recompute meta_analysis OR ####

## Create vectors for loop
n <- 36 #number of meta_analysis
OR <- vector("numeric", n)
zval <- vector("numeric", n)
pval <- vector("numeric", n)
CI_low <- vector("numeric", n)
CI_high <- vector("numeric", n)
I2 <- vector('numeric', n)
tau <- vector('numeric', n)
pvalQ <- vector('numeric', n)
logROR.reg <- vector('numeric', n)
logSE.reg <- vector('numeric', n)

for (i in 1:n){
  assign(paste0("df", i), as.data.frame(split(data_MA, data_MA$ID_RS)[[i]])) ##split dataframe in n df according to the n groups, with df[i] as name
  event.e <- get(paste0("df", i))[, "event.e"] ##extract data and assign to variable to be used in metabin
  n.e <- get(paste0("df", i))[, "n.e"]
  event.c <- get(paste0("df", i))[, "event.c"]
  n.c <- get(paste0("df", i))[, "n.c"]
  blind <- get(paste0("df", i))[, "blind"]
  sequence <- get(paste0("df", i))[, "sequence"]
  allocation <- get(paste0("df", i))[, "allocation"]
  incomplete <- get(paste0("df", i))[, "incomplete"]

  meta <- metabin(event.e, n.e, event.c, n.c, ##metabin function
            method="MH",
            sm="OR",
            incr="TACC", allincr=F,
            addincr=F, allstudies=F,
            MH.exact=F, 
            comb.fixed=F, comb.random=T,
            hakn=F,
            method.tau="DL",
            tau.preset=NULL, TE.tau=NULL,
            method.bias="peters",
            backtransf=T, pscale = 1,
            print.CMH=T,
            keepdata = T)
  
  ## vectors storing data
  OR[i] <- exp(meta$TE.random)
  CI_low[i] <- exp(meta$lower.random)
  CI_high[i] <- exp(meta$upper.random)
  CI <- c(rbind(CI_low, CI_high))
  CI <- matrix(unlist(CI), ncol = 2, byrow = TRUE)
  OR_CI <- round(cbind(OR, CI), 2)
  colnames(OR_CI) <- c("OR", "lower", "upper")
  I2[i] <- round((meta$I2 * 100), 1)
  pvalQ[i] <- round(meta$pval.Q, 3)
  tau[i] <- meta$tau
  results <- cbind (OR_CI, I2, pvalQ, tau) ##final dataframe with all needed data
  
  #metaregression
  logROR.reg[i] <- metareg(x = meta, blind, method.tau = "REML")$b[2]
  logSE.reg[i] <- metareg(x = meta, blind, method.tau = "REML")$se[2]
  reslt.reg <- cbind(logROR.reg, logSE.reg)

}

#save to excel to check consistency with published study results 
results <- data.frame(results) 
results$studlab <- ma_id$Author 
results_tbl <- results %>% as_tibble()  
write_xlsx(results_tbl, path = "output/meta_analyses_OR.xlsx", 
           col_names = TRUE, format_headers = TRUE)

reslt.reg <- data.frame(reslt.reg) #convert to data.frame 

## Count number of studies blinded or not-blinded within meta-analysis
counts_blind <- data_MA %>% 
  group_by(ID_RS) %>% 
  filter(blind == "A") %>% 
  count()

counts_unblind <- data_MA %>% 
  group_by(ID_RS) %>% 
  filter(blind == "H") %>% 
  count() 

counts <- tibble(ID_RS = counts_blind$ID_RS, 
                 counts_blind = as.character(counts_blind$n), 
                 counts_unblind = as.character(counts_unblind$n)) %>% 
  filter(!is.na(ID_RS)) %>% 
  mutate(name = ma_id$Author)

## extract number of RCTs within meta-analysis
reslt.reg$L <- counts$counts_blind
reslt.reg$H <- counts$counts_unblind

## Extract study time points and add to regression results
term_MA <- ma_id$SHORT_TERM
reslt.reg$term <- term_MA

## Extract study names and add to regression results
studlab <- as.character(ma_id$Author) 
reslt.reg$studlab <- studlab 

## Change label for study time-points
reslt.reg$term <- replace(reslt.reg$term, reslt.reg$term==0, 
                          "Long-term mortality assessment time-point")
reslt.reg$term <- replace(reslt.reg$term, reslt.reg$term==1, 
                          "Short-term mortality assessment time-point")

#### meta-epidemiological analysis ####

metaepi <- metagen(reslt.reg$logROR.reg, 
                   reslt.reg$logSE.reg, 
                   data = reslt.reg, 
                   method.tau="DL", 
                   sm = "OR", 
                   byvar = reslt.reg$term,
                   backtransf = T, 
                   hakn = F, comb.fixed=F, 
                   studlab = reslt.reg$studlab)

interact <- metareg(x = metaepi, term, method.tau = "REML") #interaction test
interact$pval[2]

## forest plot
forest(metaepi, 
       sortvar = -TE, 
       studlab = reslt.reg$studlab, 
       print.byvar = F, 
       digits = 2,
       smlab = " ", 
       xlab.pos = -0.1,
       xlim=c(0.1, 10), 
       rightcols = c("effect", "ci", "w.random"), 
       rightlabs = c("ROR", "CI (95%)", "Weight"), colgap = "1mm",
       leftlabs = c("Study", "Blinded \nRCTs", "Non-blinded \nRCTs"),
       leftcols = c("studlab", "L", "H"), 
       label.right="Non-blinded less beneficial", col.label.right="black",
       label.left="Non-blinded more beneficial", col.label.left="black",
       plotwidth = "8cm", squaresize = 1.1,
       lty.random=2, 
       col.square = "lightblue", col.inside = "black",  col.diamond = "cornflowerblue", 
       col.random = "cornflowerblue", col.by = "gray20",  
       fs.study=10, fs.axis = 10, fs.heading=10, fs.random = 11, fs.study.label=10, 
       ff.study="plain", ff.study.label="plain", ff.heading="bold", 
       ff.axis="plain", ff.random = "bold", 
       ff.xlab = "plain", ff.smlab="plain", ff.fixed="plain", ff.hetstat="plain"
       )

dev.print(device = pdf, 
          file = "output/main_forest.pdf", 
          width = 12, height = 12)

dev.print(device = tiff, res=300,
          file = "output/main_forest.tiff", 
          width = 2900, height = 3200)

