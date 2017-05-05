rm(list = ls())
library('randomForest')
library('mets')

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
dd0 <- dd0[dd0$camilla_dataset == 1,]

# List of predictor variable names
predvar <- c("angry_lamy", "angry_ramy", "fear_lamy", "fear_ramy", "neutral_lamy", "neutral_ramy")


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

predvar <- c('sb.hb', 'sb.neo')

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

predvar <- c('dasb.hb', 'dasb.neo')

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

#### Function returns list with 6 objects: (1) train and (2) test dataframes along with HC and SAD train and test group (3-6). Removes rows from an input "dataframe" based on "criteria" and columns except "outvar" and "predvar". Train and Test data.frames are created based on random sampling with "n" datasets  from each group (Healthy Control & Case) allocated to Train data.frame. That is, Train data.frame is balanced across groups.  Test data.frame consists of remaining datasets (not necessarily balanced).

fx_sample <- function(dataframe,n,criteria,outvar='group',predvar){

    # Check that criteria specified is supported
    criteria_set <- c('summer', 'winter', 'first', 'random')
    if (!criteria %in% criteria_set){
        stop(paste('Incorrect criteria. Choose from: ', paste0(criteria_set, collapse = ', ')))
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
    
    return(list(train = df_train, 
                test = df_test, 
                hc_train_list = df_sub[hc_train, c('cimbi.id','mr.id', 'group')],
                hc_test_list = df_sub[hc_test, c('cimbi.id','mr.id', 'group')],
                sad_train_list = df_sub[sad_train, c('cimbi.id','mr.id', 'group')],
                sad_test_list = df_sub[sad_test, c('cimbi.id','mr.id', 'group')]
                ))
}

### Function returns list of 3 objects: (1) observed and (2) predicted group status for datasets in Test data.frame and (3) related rF model. Function takes output from fx_sample as input ("df_list"). Runs randomForest model where "outvar" is predicted based on "predvar" (predictor variables).

fx_rF <- function(df_list,outvar='group',predvar){
    
    # model
    rF_formula <- as.formula(paste(outvar, '~ .'))
    rF <- randomForest(rF_formula, proximity = T, importance = T, data = df_list$train)
    
    # Predictors for unseen (test) datasets
    x.test <- df_list[['test']][,predvar]
    
    #Use training model to predict unseen (test) datasets
    pred <- predict(rF, newdata = x.test, type = 'class')
    
    return(list(pred = pred, actual = df_list[['test']][,outvar], rF = rF))
}

### Function returns list of 3 objects: (1) observed and (2) predicted group status for datasets in Test data.frame and (3) related logistic model. Function takes output from fx_sample as input ("df_list"). Runs logistic model where "outvar" is predicted based on "predvar" (predictor variables).

fx_logit <- function(df_list,outvar='group',predvar){
    
    # model
    logit_formula <- as.formula(paste(outvar, '~ .'))
    l <- glm(logit_formula, data = df_list$train, family = 'binomial')
    
    # Predictors for unseen (test) datasets
    x.test <- df_list[['test']][,predvar]
    
    #Use training model to predict unseen (test) datasets
    pred <- predict(l, newdata = x.test, type = 'response')
    lev <- levels(df_list[['train']][,outvar])
    pred.fact <- factor(lev[round(pred)+1], levels = c('Healthy Control', 'Case'))
    
    return (list(pred = pred.fact, actual = df_list[['test']][,outvar], logit = l))
}

####
## RUN MODELS
####

## Example code
# MR
a <- fx_sample(dd0, 12, 'winter', predvar=predvar)
b <- fx_rF(a,predvar=predvar)
c <- fx_logit(a,predvar=predvar)

a <- fx_sample(dd_sb, 5, 'winter', predvar=predvar)
a <- fx_sample(dd_dasb, 15, 'winter', predvar=predvar)

# Evaluate performance
perm <- 1000
rF_out <- lapply(seq(perm), 
             function(i) fx_rF(fx_sample(dd0, 12, 'winter', predvar=predvar),predvar=predvar))

rF_table <- matrix(0, nrow = 2, ncol = 2)
rownames(rF_table) <- c('Healthy Control', 'Case')
colnames(rF_table) <- c('Healthy Control', 'Case')
for (i in seq(length(rF_out))){
    rF_table <- table(rF_out[[i]]$pred, rF_out[[i]]$actual) + rF_table
}
print(rF_table)

logit_out <- lapply(seq(perm), 
                 function(i) fx_logit(fx_sample(dd0, 12, 'winter', predvar=predvar),predvar=predvar))
logit_table <- matrix(0, nrow = 2, ncol = 2)
rownames(logit_table) <- c('Healthy Control', 'Case')
colnames(logit_table) <- c('Healthy Control', 'Case')
for (i in seq(length(tmp))){
    logit_table <- table(logit_out[[i]]$pred, logit_out[[i]]$actual) + logit_table
}
print(logit_table)


####
## OLD STUFF (MAYBE USEFUL ONE DAY)
####

### Function returns performance measures associated with prediction models
### Takes as input number of permutations ("perm")

fx_modelRun <- function(perm){
    
    # Performance tables
    rF_perf <- logit_perf <- 
        matrix(0, ncol = 2, nrow = 2, 
               dimnames = list(c('Case', 'HC'), c('Case', 'HC')))
    
    # Iteration-by-iteration percent correct.
    pcorr <- data.frame(rF_perf = rep(0,perm), logit_perf = rep(0,perm))
    
    for (i in seq(perm)){
        
        # Create Train and Test data.frames
        df_list <- fx_sample(dd0, 12, 'first', predvar=predvar)
        
        # Perform rF and logistic regression
        rF_out <- fx_rF(df_list, predvar=predvar)
        logit_out <- fx_logit(df_list, predvar=predvar)
        
        # Update rF performance table
        rF_perf['Case','Case'] <- 
            rF_perf['Case','Case'] + sum(rF_out$pred == 'Case' & rF_out$y == 'Case')
        
        rF_perf['Case','HC'] <- 
            rF_perf['Case','HC'] + sum(rF_out$pred == 'Case' & rF_out$y == 'Healthy Control')
        
        rF_perf['HC','Case'] <- 
            rF_perf['HC','Case'] + sum(rF_out$pred == 'Healthy Control' & rF_out$y == 'Case')
        
        rF_perf['HC','HC'] <- 
            rF_perf['HC','HC'] + sum(rF_out$pred == 'Healthy Control' & rF_out$y == 'Healthy Control')
        
        # Update logit performance table
        logit_perf['Case','Case'] <- 
            logit_perf['Case','Case'] + sum(logit_out$pred == 'Case' & logit_out$y == 'Case')
        
        logit_perf['Case','HC'] <- 
            logit_perf['Case','HC'] + sum(logit_out$pred == 'Case' & logit_out$y == 'Healthy Control')
        
        logit_perf['HC','Case'] <- 
            logit_perf['HC','Case'] + sum(logit_out$pred == 'Healthy Control' & logit_out$y == 'Case')
        
        logit_perf['HC','HC'] <- 
            logit_perf['HC','HC'] + sum(logit_out$pred == 'Healthy Control' & logit_out$y == 'Healthy Control')
        # Update performance tables
        rF_corr <- sum(rF_out$pred == 'Case' & rF_out$y == 'Case') + sum(rF_out$pred == 'Healthy Control' & rF_out$y == 'Healthy Control')
        rF_tot <- rF_corr + sum(rF_out$pred == 'Case' & rF_out$y == 'Healthy Control') + sum(rF_out$pred == 'Case' & rF_out$y == 'Healthy Control')
        logit_corr <- sum(logit_out$pred == 'Case' & logit_out$y == 'Case') + sum(logit_out$pred == 'Healthy Control' & logit_out$y == 'Healthy Control')
        logit_tot <- logit_corr + sum(logit_out$pred == 'Case' & logit_out$y == 'Healthy Control') + sum(logit_out$pred == 'Healthy Control' & logit_out$y == 'Case')
        
        # Update iteration-by-iteration performance arrays
        pcorr[i,'rF_perf'] <- signif(rF_corr/rF_tot, 4)*100
        pcorr[i,'logit_perf'] <- signif(logit_corr/logit_tot, 4)*100
        
    }
    
    # Compute overall performance
    rF_pcorr <- signif((sum(diag(rF_perf))/sum(rF_perf))*100, 5)
    logit_pcorr <- signif((sum(diag(logit_perf))/sum(logit_perf))*100, 5)
    
    return(list(pcorr = pcorr, 
                rF_pcorr = rF_pcorr, logit_pcorr = logit_pcorr, 
                rF_perf = rF_perf, logit_perf = logit_perf))
}
