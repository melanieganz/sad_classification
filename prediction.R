rm(list = ls())
require('randomForest')
require('mets')
require('parallel')
require('ROCR')
require('e1071')

out.folder <- '/data1/patrick/fmri/hvi_trio/sad_classification/'

####
## LOAD AND ORGANIZE DATA
####

###
## Load fMRI data
###

dd0 <- read.csv('/data1/patrick/fmri/hvi_trio/sad_classification/sad_fmri.csv')

# Convert variables to character
dd0$group <- as.character(dd0$group)
dd0$neutral_files <- as.character(dd0$neutral_files)
dd0$angry_files <- as.character(dd0$angry_files)
dd0$fear_files <- as.character(dd0$fear_files)
dd0$mr.id <- as.character(dd0$mr.id)

# Define first scan based on mr.id
dd0$scan_order <- NA
for (i in unique(dd0$cimbi.id)){
    dd0[which(dd0$cimbi.id == i), 'scan_order'] <- order(dd0[dd0$cimbi.id == i,'mr.id'])
}

# Data included in Camilla manuscript
# dd0 <- dd0[dd0$camilla_dataset == 1,]

###
## Load SB data
###
 
dd_sb <- read.csv('/data1/patrick/fmri/hvi_trio/sad_classification/DBproject_SumWin_SB_MelPat.csv')

# Rename columnes
colnames(dd_sb)[colnames(dd_sb) == 'CIMBI.ID'] <- 'cimbi.id'
colnames(dd_sb)[colnames(dd_sb) == 'Person.status'] <- 'group'
colnames(dd_sb)[colnames(dd_sb) == 'Gender'] <- 'sex'
colnames(dd_sb)[colnames(dd_sb) == "Age.at.SB.scan"] <- 'age.sb'
colnames(dd_sb)[colnames(dd_sb) == "SB.scan.date"] <- 'sb.scan.date'
colnames(dd_sb)[colnames(dd_sb) == "SB.ID"] <- 'sb.id'
colnames(dd_sb)[colnames(dd_sb) == "HighBinding_SB_BPnd_NonPV_GM"] <- 'sb.hb'
colnames(dd_sb)[colnames(dd_sb) == "Neocortex_SB_BPnd_NonPV_GM"] <- 'sb.neo'

# Determine genotype status
dd_sb$httlpr2 <- factor(dd_sb$SLC6A4.5HTTLPR == 'll', labels = c('sx', 'll'))
dd_sb$sert2 <- factor(dd_sb$SLC6A4.5HTTLPR == 'll' & dd_sb$"SLC6A4.5HTTLPR.A.G..l.allele." == 'AA', labels = c('sx', 'lala'))

# Set variable type
dd_sb$sb.id <- as.character(dd_sb$sb.id)
dd_sb$group <- as.character(dd_sb$group)

# Determine season for each scan based on sb.scan.date
dd_sb$season <- NA
for (i in seq(nrow(dd_sb))){
    if (as.numeric(format(as.Date(dd_sb[i,'sb.scan.date'], '%d/%m/%Y'), '%m')) %in% seq(4,8)){
        dd_sb[i,'season'] <- 'S'
    } else if((as.numeric(format(as.Date(dd_sb[i,'sb.scan.date'], '%d/%m/%Y'), '%m')) %in% c(seq(10,12), seq(2)))){
        dd_sb[i,'season'] <- 'W'
    }
}
dd_sb$season <- factor(dd_sb$season, levels = c('S', 'W'))

# Define first scan based on sb.id
dd_sb$scan_order <- NA
for (i in unique(dd_sb$cimbi.id)){
    dd_sb[which(dd_sb$cimbi.id == i), 'scan_order'] <- order(dd_sb[dd_sb$cimbi.id == i,'sb.id'])
}

dd_sb <- cbind(dd_sb[,c('cimbi.id', 'group', 'season', 'sb.id', 'sex', 'age.sb', 'sb.hb', 'sb.neo', 'httlpr2', 'sert2', 'scan_order')])

###
## Load DASB data
###

dd_dasb <- read.csv('/data1/patrick/fmri/hvi_trio/sad_classification/DBproject_SumWin_DASB_MelPat.csv')

# Rename columnes
colnames(dd_dasb)[colnames(dd_dasb) == 'CIMBI.ID'] <- 'cimbi.id'
colnames(dd_dasb)[colnames(dd_dasb) == 'Person.status'] <- 'group'
colnames(dd_dasb)[colnames(dd_dasb) == 'Gender'] <- 'sex'
colnames(dd_dasb)[colnames(dd_dasb) == "Age.at.DASB.scan"] <- 'age.dasb'
colnames(dd_dasb)[colnames(dd_dasb) == "DASB.scan.date"] <- 'dasb.scan.date'
colnames(dd_dasb)[colnames(dd_dasb) == "DASB.ID"] <- 'dasb.id'
colnames(dd_dasb)[colnames(dd_dasb) == "HighBinding_DASB_BPnd_NonPV_GM"] <- 'dasb.hb'
colnames(dd_dasb)[colnames(dd_dasb) == "GlobNeoCort_DASB_BPnd_NonPV_GM"] <- 'dasb.neo'

# Determine genotype status
dd_dasb$httlpr2 <- factor(dd_dasb$SLC6A4.5HTTLPR == 'll', labels = c('sx', 'll'))
dd_dasb$sert2 <- factor(dd_dasb$SLC6A4.5HTTLPR == 'll' & dd_dasb$"SLC6A4.5HTTLPR.A.G..l.allele." == 'AA', labels = c('sx', 'lala'))

# Set variable type
dd_dasb$dasb.id <- as.character(dd_dasb$dasb.id)
dd_dasb$group <- as.character(dd_dasb$group)

# Determine season for each scan based on sb.scan.date
dd_dasb$season <- NA
for (i in seq(nrow(dd_dasb))){
    if (as.numeric(format(as.Date(dd_dasb[i,'dasb.scan.date'], '%d/%m/%Y'), '%m')) %in% seq(4,8)){
        dd_dasb[i,'season'] <- 'S'
    } else if((as.numeric(format(as.Date(dd_dasb[i,'dasb.scan.date'], '%d/%m/%Y'), '%m')) %in% c(seq(10,12), seq(2)))){
        dd_dasb[i,'season'] <- 'W'
    }
}
dd_dasb$season <- factor(dd_dasb$season, levels = c('S', 'W'))

# Define first scan based on sb.id
dd_dasb$scan_order <- NA
for (i in unique(dd_dasb$cimbi.id)){
    dd_dasb[which(dd_dasb$cimbi.id == i), 'scan_order'] <- order(dd_dasb[dd_dasb$cimbi.id == i,'dasb.id'])
}

dd_dasb <- cbind(dd_dasb[,c('cimbi.id', 'group', 'season', 'dasb.id', 'sex', 'age.dasb', 'dasb.hb', 'dasb.neo', 'httlpr2', 'sert2', 'scan_order')])

###
## Load Neuropsych data
####

dd_neuropsych <- read.csv('/data1/Ganz/Project14/DBproject_SumWin_Neuropsy_MelPat.csv',sep = ";")

# Rename columnes
colnames(dd_neuropsych)[colnames(dd_neuropsych) == 'CIMBI.ID'] <- 'cimbi.id'
colnames(dd_neuropsych)[colnames(dd_neuropsych) == 'Person.status'] <- 'group'
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "Date.of.neuropsychological.examination"] <- 'date.neuropsych'
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "Date.of.NEO.P.IR.examination"] <- 'date.neopir'

#Covariates
colnames(dd_neuropsych)[colnames(dd_neuropsych) == 'Gender'] <- 'sex'
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "Age.at.neuropsych"] <- 'age.neuropsych'
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "Rist.-.Index"] <- 'IQ'
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "MDI"] <- 'MDI'

# Livs paper
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "SDMT.'total'.score"] <- 'neuropsych.SDMT' #cognitive processing speed
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "Letter-number.sequencing.total.score"] <- 'neuropsych.LNS' 
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "SRT.-.Total.mean.reaction.latency"] <- 'neuropsych.SRT' # motor speed
                        
# Deas paper                        
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "N:.Neuroticism"] <- 'neopir.neuroticism'
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "E: Extraversion"] <- 'neopir.extraversion'
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "O: Openness"] <- 'neopir.openness'
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "A: Agreeableness"] <- 'neopir.agreeableness'
colnames(dd_neuropsych)[colnames(dd_neuropsych) == "C: Conscientiousness"] <- 'neopir.conscientiousness'

# Determine season for testing based on testing dates
dd_neuropsych$season.neuropsych <- NA
for (i in seq(nrow(dd_neuropsych))){
  if (as.numeric(format(as.Date(dd_neuropsych[i,'date.neuropsych'], '%d-%m-%Y'), '%m')) %in% seq(4,8)){
    dd_neuropsych[i,'season.neuropsych'] <- 'S'
  } 
  else if((as.numeric(format(as.Date(dd_neuropsych[i,'date.neuropsych'], '%d-%m-%Y'), '%m')) %in% c(seq(10,12), seq(2)))){
    dd_neuropsych[i,'season.neuropsych'] <- 'W'
  }
}
dd_neuropsych$season.neuropsych <- factor(dd_neuropsych$season.neuropsych, levels = c('S', 'W'))

dd_neuropsych$season.neopir <- NA
for (i in seq(nrow(dd_neuropsych))){
  if (as.numeric(format(as.Date(dd_neuropsych[i,'date.neopir'], '%d-%m-%Y'), '%m')) %in% seq(4,8)){
    dd_neuropsych[i,'season.neopir'] <- 'S'
  } 
  else if((as.numeric(format(as.Date(dd_neuropsych[i,'date.neopir'], '%d-%m-%Y'), '%m')) %in% c(seq(10,12), seq(2)))){
    dd_neuropsych[i,'season.neopir'] <- 'W'
  }
}
dd_neuropsych$season.neopir <- factor(dd_neuropsych$season.neopir, levels = c('S', 'W'))


####
## DEFINE FUNCTIONS
####

### Function name: fx_scramble
### Function returns:
###     (1) permuted dataframe
### Function description:
###     Permutes a "dataframe" based on some column name ("to.scramble") but ensures that shuffled column name is consistent for each person ("unique.id")
fx_scramble <- function(dataframe, unique.id = 'cimbi.id', to.scramble = 'group'){
    
    # Obtain cimbi.id
    ids <- unique(dataframe[,unique.id])
    # Set of group assignments
    group_list <- c()
    for (id in ids){
        group_list <- c(group_list, unique(dataframe[dataframe[,unique.id] == id, to.scramble]))
    }
    
    # Permuted group assignment
    new_group_list <- sample(group_list, replace = F)
    
    # Update group column to reflect permuted group assignment
    for (id in ids){
        new_id <- new_group_list[which(id == ids)]
        dataframe[dataframe$cimbi.id == id, to.scramble] <- new_id
    }
    return(dataframe)
}

### Function name: fx_sample
### Function returns:
###     (1) train ("train"), dataframe
###     (2) test ("test"), dataframe
###     (3) HC train ids ("hc_train_list"), dataframe
###     (4) HC test ids ("hc_test_list"), dataframe
###     (5) SAD train ids ("sad_train_list"), dataframe
###     (6) SAD test ids ("sad_test_list"), dataframe
### Function description:
###     Removes rows from an input "dataframe" based on "criteria" and columns except "outvar" and "predvar". Train and Test dataframes are created based on random sampling with "n" datasets  from each group (Healthy Control & Case) allocated to Train dataframe. That is, Train data.frame is balanced across groups.  Test dataframe consists of remaining datasets (not necessarily balanced).
fx_sample <- function(dataframe, n, criteria, outvar='group', predvar=NULL){
    
    # Check that predvar is specified
    if (is.null(predvar)){
        stop(paste0('ERROR: Predictor variables not specified!'))
    }
    
    # Check that criteria specified is supported
    criteria_set <- c('summer', 'winter', 'first', 'random')
    if (!tolower(criteria) %in% criteria_set){
        stop(paste0('Incorrect criteria. Choose from: ', paste0(criteria_set, collapse = ', ')))
    }
    
    # Subset data.frame based on criteria
    if (criteria == 'summer'){
        df_sub <- dataframe[dataframe[,'season'] == 'S',]
    } else if (criteria == 'winter'){
        df_sub <- dataframe[dataframe[,'season'] == 'W',]
    } else if (criteria == 'first'){
        df_sub <- dataframe[dataframe[,'scan_order'] == 1,]
    } else if (criteria == 'random'){
        df_sub <- data.frame()
        for (i in unique(dataframe$cimbi.id)){
            df_sub <- rbind(df_sub,
                            dataframe[sample(which(dataframe$cimbi.id == i), 1),])
        }
    }
    
    # Define train and test datasets
    hc_list <- which(df_sub$group == 'Healthy Control')
    hc_train <- sample(hc_list, n)
    hc_test <- hc_list[!hc_list %in% hc_train]
    
    sad_list <- which(df_sub$group == 'Case')
    sad_train <- sample(sad_list, n)
    sad_test <- sad_list[!sad_list %in% sad_train]
    
    # Define train and test data.frames
    df_train <- data.frame(rbind(cbind(group = df_sub[hc_train,'group'], df_sub[hc_train, predvar]),
                                 cbind(group = df_sub[sad_train,'group'], df_sub[sad_train, predvar])
                                 ))
    df_test <- data.frame(rbind(cbind(group = df_sub[hc_test,'group'], df_sub[hc_test, predvar]),
                                 cbind(group = df_sub[sad_test,'group'], df_sub[sad_test, predvar])
    ))
    
    # List of column names to pull
    col.list <- c('cimbi.id', 'mr.id', 'dasb.id', 'sb,id', 'group')
    col.get <- col.list[col.list %in% colnames(dataframe)]
    if (!length(col.get)){
        stop(paste0('Could not find any subject identifiable columns. Choose from: ', paste0(col.list, collapse = ', ')))
    }
    
    return(list(train = df_train, 
                test = df_test, 
                hc_train_list = df_sub[hc_train, col.get],
                hc_test_list = df_sub[hc_test, col.get],
                sad_train_list = df_sub[sad_train, col.get],
                sad_test_list = df_sub[sad_test, col.get]
                ))
}

### Function name: fx_model
### Function returns:
###     (1) predicted class ("pred.class")
###     (2) predicted probability ("pred.prob")
###     (3) observed ("actual") and
###     (4) model ("model")
###     (5) type of model ("type")
### Function description:
###     Takes output from fx_sample as input ("df_list"). Runs selected model where "outvar" is predicted based on "predvar" (predictor variables).

fx_model <- function(df_list, outvar='group', predvar=NULL, model.type=NULL){
    
    # Check that predvar and model.type are specified
    if(is.null(predvar) | is.null(model.type)){
        stop(paste0('ERROR: Predictor variables or model.type not specified'))
    }
    
    # Check that model specified is supported
    model.set <- c('logistic', 'rf', 'svm')
    if(!tolower(model.type) %in% model.set){
        stop(paste0('Specify appropriate model type. Choose from: ', paste0(model.set, collapse = ', ')))
    } else {model.type <- tolower(model.type)}
    
    # model
    model.formula <- as.formula(paste(outvar, '~ .'))
    
    # Apply model based on model.type
    if (model.type == 'logistic'){
        model <- glm(model.formula, data = df_list$train, family = 'binomial')
    } else if (model.type == 'rf'){
        model <- randomForest(model.formula, data = df_list$train)
    } else if (model.type == 'svm'){
        model <- svm(model.formula, data = df_list$train)
    } else {stop(paste0('Cannot apply model.type: ', model.type, '! How did you get this far?!'))}
    
    # Predictors for unseen (test) datasets
    x.test <- df_list[['test']][,predvar]
    
    # Use model to predict test datasets
    if (model.type == 'logistic'){
        pred.prob <- predict(model, newdata = x.test, type = 'resp')
        lev <- c('Healthy Control', 'Case')
        pred.class <- factor(lev[as.numeric(pred.prob > 0.5) + 1], levels = c('Healthy Control', 'Case'))
    } else if (model.type == 'rf'){
        pred.class <- predict(model, newdata = x.test, type = 'class')
        pred.prob <- predict(model, newdata = x.test, type = 'prob')[,'Case']
    } else if (model.type == 'svm'){
        pred.class <- predict(model, x.test)
        pred.prob <- NULL
    } else {stop(paste0('Cannot evaluate performance of model.type: ', model.type, '! How did you get this far?!'))}
    
    return(list(pred.class = pred.class, pred.prob = pred.prob, actual = df_list[['test']][,outvar], model = model, type = model.type))
    
}

### Function name: fx_cTable
### Function returns:
###     (1) contingency table
### Function description:
###     Computes and returns contingency table. Function takes either (1) a list containing many models evaluated ("many" = T) or (2) a single model object ("many" = F).
fx_cTable <- function(modelObj, many = T, groups = c('Healthy Control', 'Case')){
    
    if (!many){
        cTable <- table(modelObj$pred.class, modelObj$actual)
    } else {
        nObj <- length(modelObj)
        cTable <- matrix(0, nrow = 2, ncol = 2)
        rownames(cTable) <- groups
        colnames(cTable) <- groups
        for (i in seq(nObj)){
            cTable <- table(modelObj[[i]]$pred.class, modelObj[[i]]$actual) + cTable
        }
    }
    return (cTable)
}

### Function name: fx_modelPerf
### Function returns:
###     (1) Various performance metrics
###     (2) Contingency table from which (1) is computed
### Function description:
###     Computes and returns performance measures associated with a contingency table. Function takes either (1) a list containing ("pred") and ("actual") or (2) an already computed contingency table (in which case, make.c_table = F).  Metrics computed and returned along with contingency table.

fx_modelPerf <- function(modelObj, make.c_table = T){
    
    if (make.c_table){
        # Contingency table
        c_table <- table(modelObj$pred.class, modelObj$actual)
    } else {c_table <- modelObj}
    
    # performance measures
    perf <- list()
    
    # reference
    perf$positive <- 'Case'
    perf$negative <- 'Healthy Control'
    
    # true positive
    perf$TP <- c_table[perf$positive, perf$positive]
    # false positive
    perf$FP <- c_table[perf$positive, perf$negative]
    # true negative
    perf$TN <- c_table[perf$negative, perf$negative]
    # false negative
    perf$FN <- c_table[perf$negative, perf$positive]
    # sensitivity
    perf$sensitivity <- perf$TP/(perf$TP + perf$FN)
    # specificity
    perf$specificity <- perf$TN/(perf$TN + perf$FP)
    # F1
    perf$F1 <- (2*perf$TP)/((2*perf$TP) + perf$FP + perf$FN)
    # accuracy
    perf$accuracy <- (perf$TP + perf$TN)/sum(c_table)
    # contingency table
    perf$c_table <- c_table
    
    return(perf)
}

### Function name: fx_internalPerm
### Function returns:
###     (1) Object of model performance metrics
### Function description:
###     In light of a discussion with Melanie, this function performs a shuffling of group labels in "dd.frame" and then a random train/test split "rsplit" times.  This aligns how the null distribution is derived and how the observed performance measures are assessed.
fx_internalPerm <- function(dd.frame, rsplit, measure = 'accuracy', model.type){

    # Check that model specified is supported
    model.set <- c('logistic', 'rf', 'svm')
    if(!tolower(model.type) %in% model.set){
        stop(paste0('Specify appropriate model type. Choose from: ', paste0(model.set, collapse = ', ')))
    } else {model.type <- tolower(model.type)}
    
    dd.scramble <- fx_scramble(dd.frame)
    modelOutput <- mclapply(seq(rsplit), function(i) {
        set.seed(i)
        fx_model(fx_sample(dd.scramble, nTrainSize, splitType, predvar = predvar), predvar = predvar, model.type = model.type)},
        mc.cores = 20)
    modelTable <- fx_cTable(modelOutput)
    modelTablePerf <- fx_modelPerf(modelTable, make.c_table = F)
    return(modelTablePerf)
}

### Function name: fx_nullComparison
### Function returns:
###     (1) Observed null distribution values
###     (2) Observed value
###     (3) Null p-value
###     (A) Histogram plot of 1-3
### Function description:
###     Returns a histogram of observed values (null distribution) and vertical line of single value (observed value) and null-derived p-value. Function takes (1) either list or array of values describing null distribution ("permObject"), (2) either list or value describing observed measure ("obsObject") and (3) name of measure ("measure") used to read if inputs are list.  Plot object (i.e., histogram of null distribution), observed value and permutation-derived p-value is returned
fx_nullComparison <- function(permObject, obsObject, measure = 'accuracy'){
    
    # Check that measure is among those in list
    measure_list <- c('sensitivity', 'specificity', 'f1', 'accuracy')
    if (!tolower(measure) %in% measure_list){
        stop(paste('Please choose measure from:', paste0(measure_list, collapse = ', ')))
    }
    
    # Number of permutations
    n <- length(permObject)
    
    # Determine whether permObject is list or (assumed) array. If list, return array of specified measure
    if (typeof(permObject) == 'list'){
        permDistribution <- unlist(sapply(seq(n), function(i) permObject[[i]][measure]))
    } else {permDistribution <- permObject}
    
    # Determine whether obsObject is list or (assumed) array. If list, return specified measure
    if(typeof(obsObject) == 'list'){
        obsMeasure <- unlist(obsObject[measure])
    } else {obsMeasure <- obsObject}
    
    # Permutation derived p-value
    permP <- sum(permDistribution > obsMeasure)/length(permDistribution)
    
    # Generate histogram
    hist(permDistribution, las = 1, main = paste0('Null distribution: ', measure))
    mtext(paste0('obs. perf: ', signif(obsMeasure, 3), '; perm. p-value = ', signif(permP,3)))
    title(sub = paste0('N = ', n, ' permutations'))
    abline(v = obsMeasure, lty = 2, col = 'red')
    
    out <- list()
    out$permValues <- permDistribution
    out$obsValue <- obsMeasure
    out$p.value <- permP
    
    return(out)
}

### Function name: fx_popROC
### Function returns:
###     (1) prediction object (ROCR related)
###     (2) performance object (ROCR related)
### Function description:
###     Something of a wrapper on the ROCR package.  Derives measures that can be plotted as ROC curve based on set of resamples (or permutations).  Where permutation test is performed, the true "actual.class" is not known from modelObj.  Therefore, "perm.test" must be specified with a reference model object that is assumed to have class specified correctly.  In principle, any permutation test should be a reference for an observed model, meaning an observed model with correct labels should be available.
fx_popROC <- function(modelObj, measure = 'tpr', x.measure = 'fpr', perm.test = F){
    
    label.ordering <- c('Healthy Control', 'Case')
    
    # Number of runs in modelObj
    nruns <- length(modelObj)
    
    # Extract all estimated probabilities from modelObj
    est.prob <- unlist(lapply(seq(nruns), function(i) modelObj[[i]]$pred.prob))
    
    # Determine subject-specific names and compute subject-specific mean probability
    set.names <- unique(names(est.prob))
    mean.prob <- sapply(set.names, function(i) mean(est.prob[names(est.prob) == i]))
    
    # If ROC for permutation is being derived, model object with class correctly specified is used to derive "actual.class"
    if(!identical(perm.test,F)){
        all.class <- names(perm.test$predObj@predictions[[1]])
        actual.class <- unlist(lapply(names(mean.prob), function(i) perm.test$predObj@labels[[1]][names(mean.prob[i]) == all.class]))    
    } else {
        all.class <- unlist(lapply(seq(nruns), function(i) modelObj[[i]]$actual))
        actual.class <- unlist(lapply(set.names, function(i) unique(all.class[names(est.prob) == i])))
    }
    
    # Compute prediction and performance objects
    predObj <- prediction(mean.prob, actual.class, label.ordering = label.ordering)
    perfObj <- performance(predObj, measure, x.measure)
    
    return(list(predObj = predObj, perfObj = perfObj))
}

####
## RUN MODELS
####

###
## Define predictor variables
###

# fMRI
predvar <- c("angry_lamy", "angry_ramy", "fear_lamy", "fear_ramy", "neutral_lamy", "neutral_ramy")
dd <- dd0
dd[,predvar] <- scale(dd0[,predvar])

# # # SB
# predvar <- c('sb.hb', 'sb.neo')
# dd <- dd_sb
# dd[,predvar] <- scale(dd[,predvar])
# 
# # 
# # # DASB
# predvar <- c('dasb.hb', 'dasb.neo')
# dd <- dd_sb
# dd[,predvar] <- scale(dd[,predvar])

###
## Example code - single-run (based on fMRI data)
###

a.scramble <- fx_scramble(dd)
a <- fx_sample(dd, 12, 'winter', predvar=predvar)
b.logit <- fx_model(a,predvar=predvar, model.type = 'logistic')
b.rF <- fx_model(a,predvar=predvar, model.type = 'rf')
c <- fx_modelPerf(b.logit)

###
## Example code - evaluate performance (based on fMRI data)
###

# Train/Test split
nTrainSize <- 12
# Data type to use
splitType <- 'winter'

## Repeated random splits to derive mean performance
# n random splits
rsplit <- 10

# Evaluate performance for each split
rF_rsplit <- mclapply(seq(rsplit), function(i) fx_model(fx_sample(dd, nTrainSize, splitType, predvar=predvar), predvar=predvar, model.type = 'rf'), mc.cores = 20)

# Contingency table across splits
rF_rsplitTable <- fx_cTable(rF_rsplit)

# Performance measures across splits
rsplitPerf <- fx_modelPerf(rF_rsplitTable, make.c_table = F)

# Generate ROC
rF.popROC <- fx_popROC(rF_rsplit)
plot(rF.popROC$perfObj)

## Permute labels to derive null distribution

###
## Random split resampling within each permutation
###

# n permutations
perm <- 100

permOutput <- lapply(seq(perm), function(i) fx_internalPerm(dd))
fx_nullComparison(permOutput, rsplitPerf)

###
## Old permutation method
##

### Description: For each permutation ("perm") group labels are shuffled and then ids are divided into train/test groups ("fx_sample").  A model is fit ("fx_model") with "nTrainSize" individuals from each group in the train dataset based on "splitType".  Model is evaluated with unseen observations to derive performance measure.  Performing this "perm" times derives a null distribution against which true observations can be compared.

### Why might this be a problem? Through discussions between Patrick and Melanie, Melanie suggested that this is an imperfect null distribution because each permutation does not assess accuracy for each dataset, only those that happened to be assigned to test group.  As an alternative, she suggested shuffling group labels and then evaluating model "rsplit" times to derive an accuracy measure for a given permutation.  This process would be performed "perm" times to derive a null distribution.

### PMF 20170519

# n permutations
perm <- 100

# Evaluate performance for each permutation
rF_perm <- mclapply(seq(perm),function(i) fx_model(fx_sample(dd, nTrainSize, splitType, predvar=predvar, perm = T), predvar=predvar, model.type = 'rf'),mc.cores = 20)

# Permutation contingency table
rF_permTable <- fx_cTable(rF_perm)

# Permutation performance measures
rF_permPerf <- lapply(seq(perm),
                  function(i) fx_modelPerf(rF_perm[[i]]))
rF_perm.popROC <- fx_popROC(rF_perm, perm.test = rF.popROC)

plot(rF.popROC$perfObj, col = 'red')
plot(rF_perm.popROC$perfObj, add = T, col = 'black')

# Performance measures vs. null distribution
perfSummary <- fx_nullComparison(rF_permPerf, rsplitPerf, measure = 'accuracy')

