# Define data
param <- list()
param$fMRI <- list()
param$fMRI$df_name <- dd0
param$fMRI$nTrainSize <- 16
param$fMRI$splitType_set <- c('winter', 'summer', 'first', 'random')
param$fMRI$predvar <- c("angry_lamy", "angry_ramy", "fear_lamy", "fear_ramy", "neutral_lamy", "neutral_ramy")
param$fMRI$df_name[,param$fMRI$predvar] <- scale(param$fMRI$df_name[,param$fMRI$predvar])

param$SB <- list()
param$SB$df_name <- dd_sb
param$SB$nTrainSize <- 5
param$SB$splitType_set <- c('winter', 'summer', 'first', 'random')
param$SB$predvar <- c('sb.hb', 'sb.neo')
param$SB$df_name[,param$SB$predvar] <- scale(param$SB$df_name[,param$SB$predvar])

param$DASB <- list()
param$DASB$df_name <- dd_dasb
param$DASB$nTrainSize <- 15
param$DASB$splitType_set <- c('winter', 'summer', 'first', 'random')
param$DASB$predvar <- c('dasb.hb', 'dasb.neo')
param$DASB$df_name[,param$DASB$predvar] <- scale(param$DASB$df_name[,param$DASB$predvar])

# Output directory
top <- '/data1/patrick/fmri/hvi_trio/sad_classification/'

# n random splits
rsplit <- 100

# n permutations
perm <- 10000

# Performance measures
measure_set <- c('specificity', 'sensitivity', 'F1', 'accuracy')
for (name in names(param)){
    
    nTrainSize <- param[[name]][['nTrainSize']]
    dd <- param[[name]][['df_name']]
    predvar <- param[[name]][['predvar']]
    splitType <- 'winter'
    for (splitType in param[[name]][['splitType_set']]){
        
        print(paste0('Working on: ', name, ', ', splitType))
        
        ## Derive observe accuracy
        
        # Evaluate performance for each split
        rF_rsplit <- mclapply(seq(rsplit), function(i) {
            set.seed(i)
            fx_model(fx_sample(dd, nTrainSize, splitType, predvar=predvar), predvar=predvar, model.type = 'rF')},
            mc.cores = 20)
        # Contingency table across splits
        rF_rsplitTable <- fx_cTable(rF_rsplit)
        # Performance measures across splits
        rF_rsplitPerf <- fx_modelPerf(rF_rsplitTable, make.c_table = F)
        
        
        logistic_rsplit <- mclapply(seq(rsplit), function(i) {
            set.seed(i)
            fx_model(fx_sample(dd, nTrainSize, splitType, predvar=predvar), predvar=predvar, model.type = 'logistic')}, 
            mc.cores = 20)
        # Contingency table across splits
        logistic_rsplitTable <- fx_cTable(logistic_rsplit)
        # Performance measures across splits
        logistic_rsplitPerf <- fx_modelPerf(logistic_rsplitTable, make.c_table = F)
        
        ## Derive null distribution
        
        rF_permOutput <- sapply(seq(perm), function(i) fx_internalPerm(dd, rsplit = rsplit, model.type = 'rF')[measure_set])
        pdf(paste0(top, name, '_', splitType, '_', nTrainSize, '_rF.pdf'))
        for (measure in measure_set){
            fx_nullComparison(unlist(rF_permOutput[measure,]), rF_rsplitPerf, measure = measure)
        }
        dev.off()
        
        logistic_permOutput <- sapply(seq(perm), function(i) fx_internalPerm(dd, rsplit = rsplit, model.type = 'logistic')[measure_set])
        pdf(paste0(top, name, '_', splitType, '_', nTrainSize, '_logistic.pdf'))
        for (measure in measure_set){
            fx_nullComparison(unlist(logistic_permOutput[measure,]), logistic_rsplitPerf, measure = measure)
        }
        dev.off()
    }
}

