####
## DEFINE DATA
####

param <- list()

# param$fMRI <- list()
# param$fMRI$df_name <- dd0
# param$fMRI$nTrainSize <- 16
# param$fMRI$splitType_set <- c('winter', 'summer', 'first', 'random')
# param$fMRI$predvar <- c("angry_lamy", "angry_ramy", "fear_lamy", "fear_ramy", "neutral_lamy", "neutral_ramy")
# param$fMRI$df_name[,param$fMRI$predvar] <- scale(param$fMRI$df_name[,param$fMRI$predvar])
# 
# param$SB <- list()
# param$SB$df_name <- dd_sb
# param$SB$nTrainSize <- 5
# param$SB$splitType_set <- c('winter', 'summer', 'first', 'random')
# param$SB$predvar <- c('sb.hb', 'sb.neo')
# param$SB$df_name[,param$SB$predvar] <- scale(param$SB$df_name[,param$SB$predvar])
#  
# param$DASB <- list()
# param$DASB$df_name <- dd_dasb
# param$DASB$nTrainSize <- 15
# param$DASB$splitType_set <- c('winter', 'summer', 'first', 'random')
# param$DASB$predvar <- c('dasb.hb', 'dasb.neo')
# param$DASB$df_name[,param$DASB$predvar] <- scale(param$DASB$df_name[,param$DASB$predvar])
# 
# param$srt <- list()
# param$srt$df_name <- dd_np.srt
# param$srt$nTrainSize <- 20
# param$srt$splitType_set <- c('winter', 'summer', 'first', 'random')
# param$srt$predvar <- c("srt")
# param$srt$df_name[,param$srt$predvar] <- scale(param$srt$df_name[,param$srt$predvar])
# 
# param$sdmt_lns <- list()
# param$sdmt_lns$df_name <- dd_np.sdmt_lns
# param$sdmt_lns$nTrainSize <- 20
# param$sdmt_lns$splitType_set <- c('winter', 'summer', 'first', 'random')
# param$sdmt_lns$predvar <- c('sdmt.correct','lns')
# param$sdmt_lns$df_name[,param$sdmt_lns$predvar] <- scale(param$sdmt_lns$df_name[,param$sdmt_lns$predvar])

param$mixed1 <- list()
param$mixed1$df_name <- dd_np.all
param$mixed1$nTrainSize <- 18
param$mixed1$splitType_set <- c('winter', 'summer', 'random')
param$mixed1$predvar <- c("neopir.neuroticism","neopir.extraversion","neopir.openness","neopir.agreeableness","neopir.conscientiousness", 'sdmt.correct', 'lns', 'srt')
param$mixed1$df_name[,param$mixed1$predvar] <- scale(param$mixed1$df_name[,param$mixed1$predvar])

param$mixed2 <- list()
param$mixed2$df_name <- dd_fmri_np
param$mixed2$nTrainSize <- 8
param$mixed2$splitType_set <- c('winter', 'summer', 'random')
param$mixed2$predvar <- c("neopir.neuroticism","neopir.extraversion","neopir.openness","neopir.agreeableness","neopir.conscientiousness", 'sdmt.correct', 'lns', 'srt', "angry_lamy", "angry_ramy", "fear_lamy", "fear_ramy", "neutral_lamy", "neutral_ramy")
param$mixed2$df_name[,param$mixed2$predvar] <- scale(param$mixed2$df_name[,param$mixed2$predvar])

####
## ANALYSIS
####

# Output directory
if (Sys.getenv("LOGNAME") == 'mganz'){
    top <- '/data1/Ganz/Project14/Rresults'
} else {
    top <-  '/data1/patrick/fmri/hvi_trio/sad_classification/'
}

# n random splits
rsplit <- 100

# n permutations
perm <- 1000

startTime <- Sys.time()
print(startTime)
# Performance measures
measure_set <- c('specificity', 'sensitivity', 'accuracy')
for (name in names(param)){
    
    nTrainSize <- param[[name]]$nTrainSize
    dd <- param[[name]][['df_name']]
    predvar <- param[[name]][['predvar']]

    for (splitType in param[[name]][['splitType_set']]){
        
        print(paste0('Working on: ', name, ', ', splitType))
        
        ### Derive observed accuracy
        
        ## rF
        rF_rsplit <- mclapply(seq(rsplit), function(i) {
#             set.seed(i)
            fx_model(fx_sample(dd, nTrainSize, splitType, predvar=predvar), predvar=predvar, model.type = 'rF')},
            mc.cores = 20)
        # Contingency table across splits
        rF_rsplitTable <- fx_cTable(rF_rsplit)
        # Performance measures across splits
        rF_rsplitPerf <- fx_modelPerf(rF_rsplitTable, make.c_table = F)
        
        ## logistic
        logistic_rsplit <- mclapply(seq(rsplit), function(i) {
            set.seed(i)
            fx_model(fx_sample(dd, nTrainSize, splitType, predvar=predvar), predvar=predvar, model.type = 'logistic')}, 
            mc.cores = 20)
        # Contingency table across splits
        logistic_rsplitTable <- fx_cTable(logistic_rsplit)
        # Performance measures across splits
        logistic_rsplitPerf <- fx_modelPerf(logistic_rsplitTable, make.c_table = F)
        
        ## svm
        svm_rsplit <- mclapply(seq(rsplit), function(i) {
            set.seed(i)
            fx_model(fx_sample(dd, nTrainSize, splitType, predvar=predvar), predvar=predvar, model.type = 'svm')},
            mc.cores = 20)
        # Contingency table across splits
        svm_rsplitTable <- fx_cTable(svm_rsplit)
        # Performance measures across splits
        svm_rsplitPerf <- fx_modelPerf(svm_rsplitTable, make.c_table = F)
        
        ### Derive null distribution
        
        permOutput <- lapply(seq(perm), function(i){

            if (i %% 100 == 0){
                print(paste0('Reached perm: ', i, ' at ', Sys.time()))
            }
            
            # Set seed so that null distribution is reproducible and consistent across models
            set.seed(i)
            dd.perm <- fx_scramble(param[[name]]$df_name)
            
            tmp <- mclapply(seq(rsplit), function(j) {
                
                # Set seed so that rsplit of scrambled data is reproducible and consistent across models
                set.seed(j)
                dd.sample <- fx_sample(dd.perm, nTrainSize, splitType, predvar=predvar)
                
                # Model performance
                logisticObj <- fx_model(dd.sample, predvar=predvar, model.type = 'logistic')
                rFObj <- fx_model(dd.sample, predvar=predvar, model.type = 'rf')
                svmObj <- fx_model(dd.sample, predvar=predvar, model.type = 'svm')
                
                # Return objects for each model type, for each rsplit
                return(list(logisticObj = logisticObj, rFObj = rFObj, svmObj = svmObj))},
                mc.cores = 20)
            
            # Compute average model performance measures across each rsplit
            logisticPerf <- fx_modelPerf(fx_cTable(lapply(seq(rsplit), function(i) tmp[[i]]$logisticObj)), make.c_table = F)
            rfPerf <- fx_modelPerf(fx_cTable(lapply(seq(rsplit), function(i) tmp[[i]]$rFObj)), make.c_table = F)
            svmPerf <- fx_modelPerf(fx_cTable(lapply(seq(rsplit), function(i) tmp[[i]]$svmObj)), make.c_table = F)
            # Return model performance measures for each model type
            return(list(logisticPerf = logisticPerf, rfPerf = rfPerf, svmPerf = svmPerf))
        })
        
        # Visualize observed model performance against model-specific null distributions
        modelPerm.types <- c('logisticPerf', 'rfPerf', 'svmPerf')
        modelObs.types <- c('logistic_rsplitPerf', 'rF_rsplitPerf', 'svm_rsplitPerf')
        modelTypes <- c('logistic', 'rF', 'svm')
        nModelTypes <- length(modelTypes)
        pdf(paste0(top, name, '_', splitType, '_', nTrainSize, '_allModels.pdf'))
        for (j in seq(nModelTypes)){
            for(measure in measure_set){
                nullDistribution <- unlist(sapply(seq(perm), function(i) permOutput[[i]][[modelPerm.types[j]]][measure]))
                fx_nullComparison(nullDistribution, get(modelObs.types[j]), measure = measure, model.type = modelTypes[j])
                }
            }
        dev.off()
    }
}
endTime <- Sys.time()
