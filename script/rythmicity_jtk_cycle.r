
##############################################################################################################
#################### Redo analysis based focused on jtk_cycle algorithm and the mgp imputaion method.
# link to github is found here: https://github.com/mfcovington/jtk-cycle/blob/develop/JTK_CYCLE.R
set.seed(42)
library("dplyr")
library(tidyverse)
library(tidyr)
#source(path) # get imputation module
path <- "script/JTK_CYCLE.R"
source(path)
library(pbapply)

######## Modify and aplly algo to nuclear proteomics dataset ###############

prots <- read.csv("jtk/proteins.txt", sep = "\t")
prots <- as.matrix(prots)
prots <- t(scale(t(prots)))


### The term periods represent the time range of rythmicity. 
### takes the form value_a*sampling_time : value_b*samples_time
### In this case, sampling time is not more than 24 hours according to collected timepoints.
### so we are interested in circadian time between 20 to 24 so 5*4:6*4
## jtk params
n_timepoints <- 6
n_reps <- 3
sampling_time <- 4
periods <- 5:6 # 5*4 : 6*4 = 20 to 24 hour rythms search.
length(periods)

# define the function to run the jtk algorithm
run_jtk <- function(mat) {
    cat("\nRunning JTK for:", "Nuclear Proteins", "\n-------------------------------------")

    # init jtk
    jtkdist(n_timepoints, n_reps)
    jtk.init(periods, sampling_time)
    
    st <- system.time({
        res <- apply(mat,1,function(z) {
            jtkx(z)
            c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
        })
        res <- as.data.frame(t(res))
        bhq <- p.adjust(unlist(res[,1]),"BH")
        res <- cbind(bhq,res)
        colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
        res <- res[order(res$ADJ.P,-res$AMP),]
        # results <- cbind(annot,res,mat)
        # results <- results[order(res$ADJ.P,-res$AMP),]
        })

    savefile <- paste0("jtk/jtk_", "nucProts", ".txt")
    write.table(res, file = savefile, quote = FALSE, sep = "\t", row.names = TRUE)
}

run_jtk(prots)
cat("JTK Cycle complete...")
