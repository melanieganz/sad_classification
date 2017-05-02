rm(list = ls())
library('randomForest')
library('mets')

####
## LOAD AND ORGANIZE DATA
####

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

dd0 <- dd0[dd0$camilla_dataset == 1,]

# List of predictor variable names
predvar <- c("angry_lamy", "angry_ramy", "fear_lamy", "fear_ramy", "neutral_lamy", "neutral_ramy")


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
a <- fx_sample(dd0, 12, 'winter', predvar=predvar)
b <- fx_rF(a,predvar=predvar)
c <- fx_logit(a,predvar=predvar)

# Evaluate performance
perm <- 1000
rF_out <- lapply(seq(perm), 
             function(i) fx_rF(fx_sample(dd0, 12, 'winter', predvar=predvar),predvar=predvar))

rF_table <- matrix(0, nrow = 2, ncol = 2)
rownames(rF_table) <- c('Healthy Control', 'Case')
colnames(rF_table) <- c('Healthy Control', 'Case')
for (i in seq(length(tmp))){
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