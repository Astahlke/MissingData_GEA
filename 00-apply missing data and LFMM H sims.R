# LFMM code with mean imputation of NAs - H1/H5/H9 sims

# BE SURE TO CHECK:
# Code path names throughout

# Line 56: Check sample file: 30/100/500 individuals
# Line 73: sample size
# Line 76: [1]/[2]/[3]/[4]/[5]/[6]/[7] for percentage missing data
# NOTE that this indexing (1-7) requires that all 7 levels of missing data
# are in each missing_inputs folder, or the wrong file may be indexed

# Change output file names (x4) at end of code
# NOTE: for no missing runs, comment out lines 66-74 and 92-94

rm(list=ls())
library(raster)
library(RMTstat)
library(LEA)

sim <- dir(path="/path/to/simulations/", pattern="^H5")  #H1/H5/H9 ##################
sim
repl <- c(0:9)

setwd("/path/to/lfmm/temp/")  # temp directory to hold intermediate lfmm files
# SEARCH FOR THIS AND REPLACE WITH YOUR TEMP DIRECTORY IN LOOP CODE (APPEARS THREE TIMES)

##########################
### Loop starts here 
##########################

  for (s in 1:length(sim)) {
  for (r in 1: length(repl)) {
    
###read in the habitat files
###R1 goes with sim rep0
###Change H9, H5, H1
    fname <- paste("/path/to/Surfaces_Sample/L10H5",r,"_aa.asc", sep="")  #H1/H5/H9 ##################
    qrule <- raster(fname)                                

###read in the snp dataset
    fname1 <- paste("/path/to/simulations/",sim[s],"/R",repl[r],".csv", sep="")
    data <- read.csv(fname1)

###extract coordinates and the habitat value at each coordinate    
    coord <- data[,2:3]
    habitat <- extract(qrule,coord)

###convert UTM - subtract corner coordinate and /1000 to convert from m to km
    xvar <- (data[,2] - 477895.2)/1000         
    yvar <- (data[,3] - 4905702)/1000 

###create the full data set
    fulldata <- cbind(xvar, yvar, habitat, data[,8:207])
    
###subset individuals from the sample100 or sample500.csv files
    samp <- read.csv("/path/to/sample500.csv", header=F)
    samp <- samp[,1]; sampdata <- fulldata[samp,]

###extract the environmental variables alone  
    env <- sampdata[,1:3] 
    env <- scale(env, center=T, scale=T)
    env <- data.matrix(env)
    colnames(env)<-c("xvar","yvar","habitat")
    
###extract the snps alone
    snps <- sampdata[,4:203]
    snps <- snps[,seq(1, ncol(snps), 2)]
    colnames(snps) <- seq(0, 99, by=1)
    snps <- data.matrix(snps)
    
#enter NAs
##in glob2rx(), change "_30 "_100" or "_500"   
    missing_files <- list.files("/path/to/Input Files/", pattern = glob2rx("random_*_500.csv"), full.names = F, recursive = F) 
    
###change missing_files[x] to [1] to [6] for each percentage of missing data
    missing <- read.csv(paste0("/path/to/Input Files/", missing_files[4]),header=FALSE)  
    missing <- as.numeric(missing[,1])  
    snps[missing] <- NA
     
###remove NA rows (individuals) from both snp and env datasets; there are a few replicates where NA individuals so there will be 499 instead of 500 individuals
    subsetNA <- apply(snps,1,function(x) length(unique(x))>1)
    snps <- snps[subsetNA,]                                      
    env <- env[subsetNA,] 
    
###remove monomorphic snps (I don't think we have any, but just in case)
    subsetMM <- apply(snps,2,function(x) length(unique(x))>1)
    snps <- snps[,subsetMM] 

###Impute missing data with colMeans
    
    colMeans(snps, na.rm = TRUE)
    mean(snps[,1], na.rm=TRUE)
    
    cM <- colMeans(snps, na.rm =TRUE)
    indx <- which(is.na(snps), arr.ind=TRUE)
    snps[indx] <- cM[indx[,2]]
    snps[,1]
    
## You shouldn't have to change anything inside the loop from here; but you will have to change the output file names outside of the loop 
    
    
    ##################
    ## CALCULATE K ##
    ##################
    
    # CROSS-ENTROPY K:
    write.geno(snps, "genotypes.geno")
    project = NULL
    project = snmf("genotypes.geno", K=1:20, entropy=T, rep=2, project = "new")
    sumpr <- summary(project); sumxe <- sumpr$crossEntropy; xemean <- sumxe[2,]
    xeK <- which.min(xemean); xeK <- as.vector(xeK)
    
    unlink("/path/to/lfmm/temp/*.*", recursive=T)  ### REMOVE snmf TEMP FILES
    
      #############
      ## LEA - K ##
      #############
      
      write.table(snps, file="snps.lfmm", row.names=F, col.names=F)
      write.table(env, file="env.env", row.names=F, col.names=F)
      
      proj.K = NULL
      proj.K = lfmm("snps.lfmm", "env.env", K=xeK, repetitions=5, project="new", missing.data=F, iterations=10000 , burnin=5000)
      
      # Zscores and p-values #
      
      vars <- c("xvar", "yvar", "hab")
      
      for (v in 1:3) {
        zs <- z.scores(proj.K, K=xeK, d=v); zs.med <- apply(zs, MARGIN=1, median); lambda = median(zs.med^2)/0.456
        pv <- p.values(proj.K, K=xeK, d=v); pv.med <- apply(pv, MARGIN=1, median)
        loc_names <- c(0:99); env_names <- rep(vars[v], 100)
        zpOut <- data.frame(loc_names, env_names, zs.med, pv.med)
        bonf <- ifelse(zpOut[,4] < (.001/(3*100)), TRUE, FALSE) # alpha = .001; 3 vars, 100 markers
        zpOut_sig <- zpOut[bonf==TRUE,]; assign(vars[v], zpOut_sig)
        assign(paste("GIF_",vars[v],"_xeK", sep=""), lambda)
      }
      
      unlink("/path/to/lfmm/temp/*.*", recursive=T)  ### REMOVE lfmm TEMP FILES to avoid confusion
      
      xp.out <- rbind(xvar, yvar, hab)
      x.out <- xp.out
      if (nrow(xp.out) > 0) {
        x.out[5] <- "xe"
        # p.out[5] <- "patt"}
      outs <- x.out
      label1 <- rep(sim[s], nrow(xp.out)); label2 <- rep(repl[r], nrow(xp.out))
      outs <- cbind(label1,label2,outs); outs <- outs[complete.cases(outs),]
      
     if (nrow(outs) > 0) {   
        colnames(outs) <- c("sim","rep","locus","env", "zscore","bonf_pval","K")
        fname <- paste("out_", sim[s], "_R", repl[r], sep=""); assign(fname,outs)  
      }
      
    } else {
      
      ###############
      ## LEA - xeK ##
      ###############
      
      write.table(snps, file="snps.lfmm", row.names=F, col.names=F)
      write.table(env, file="env.env", row.names=F, col.names=F)
      
      proj.K = NULL
      proj.K = lfmm("snps.lfmm", "env.env", K=xeK, repetitions=5, project="new", missing.data=T, iterations=10000 , burnin=5000)
      
      # Zscores and p-values #
      
      vars <- c("xvar", "yvar", "hab")
      for (v in 1:3) {
        zs <- z.scores(proj.K, K=xeK, d=v); zs.med <- apply(zs, MARGIN=1, median); lambda = median(zs.med^2)/0.456
        pv <- p.values(proj.K, K=xeK, d=v); pv.med <- apply(pv, MARGIN=1, median)
        loc_names <- c(0:99); env_names <- rep(vars[v], 100)
        zpOut <- data.frame(loc_names, env_names, zs.med, pv.med)
        bonf <- ifelse(zpOut[,4] < (.001/(3*100)), TRUE, FALSE) # alpha = .001; 3 vars, 100 markers
        zpOut_sig <- zpOut[bonf==TRUE,]; assign(vars[v], zpOut_sig)
        assign(paste("GIF_",vars[v],"_xeK", sep=""), lambda)
      }
      
      x.out <- rbind(xvar, yvar, hab)
      if (nrow(x.out) > 0) {
        x.out[5] <- "xe"}
      
      unlink("/path/to/lfmm/temp/*.*", recursive=T)  ### REMOVE lfmm TEMP FILES to avoid confusion
      
      ## Pull results
      outs <- x.out
      
      label1 <- rep(sim[s], nrow(outs)); label2 <- rep(repl[r], nrow(outs))
      outs <- cbind(label1,label2,outs); outs <- outs[complete.cases(outs),]
      
      if (nrow(outs) > 0) {   
        colnames(outs) <- c("sim","rep","locus","env", "zscore","bonf_pval","K")
        fname <- paste("out_", sim[s], "_R", repl[r], sep=""); assign(fname,outs)  
      }
    }
    
    ## Output the Ks and GIFs too:
    outKs <- cbind(xeK, GIF_xvar_xeK, GIF_yvar_xeK, GIF_hab_xeK)
    outKs <- as.data.frame(outKs)
    outKs <- cbind(sim[s], repl[r], outKs)
    
    colnames(outKs) <- c("sim","rep","xeK","GIFx_xK","GIFy_xK","GIFh_xK")
    fname <- paste("outKs_", sim[s], repl[r], sep="")
    assign(fname, outKs)
    
  }
 }


########################
##### Loop ends here
########################

## save output; there will be 4 output files in total
# 1. Ks and GIFs
# 2. Raw results
# 3. True Positives
# 4. Summary file 

## Follow this output file naming : "_Xpercent_Xind_RawData.csv" with the conditions that you ran

## (1) Ks and GIFs:
save <- ls()[grep("outKs_", ls())]
bar = NA
for (l in 1:length(save)) {
  foo <- get(save[l]); bar <- rbind(foo, bar)
}
bar <- bar[-nrow(bar),]

### CHANGE OUTPUT FILE NAME: 
fname <- paste("/path/to/lfmm/output/", sim[s], "_50percent_500ind_K_GIF.csv", sep="")
write.csv(bar, file=fname, row.names=F)


## (2) Raw results:
save <- ls()[grep("out_", ls())]
bar = NA
for (l in 1:length(save)) {
  foo <- get(save[l]); bar <- rbind(foo, bar)
}
bar <- bar[-nrow(bar),]

### CHANGE OUTPUT FILE NAME: 
fname <- paste("/path/to/lfmm/output/", sim[s], "_50percent_500ind_RawData.csv", sep="")
write.csv(bar, file=fname,  row.names=F)


## (3) True positives:
tp <- bar[bar$locus == 0,]

### CHANGE OUTPUT FILE NAME: 
fname <- paste("/path/to/lfmm/output/", sim[s], "_50percent_500ind_TruePos.csv", sep="")
write.csv(tp, file=fname)


## (4) Summary of true and false positives
fp <- bar[bar$locus != 0,]

K <- "xe"

summary <- matrix(NA, nrow=length(sim), ncol=5)
colnames(summary) <- c("sim", "K", "tp", "fp", "fp.sd")

summary[,1] <- as.vector(sapply(sim, function(x) rep(x,2)))
summary[,2] <- as.vector(rep(K,length(sim)))


for (i in 1:length(sim)) {
  foo <- tp[tp$sim == sim[i],]
  baz <- fp[fp$sim == sim[i],]
  
  for (j in 1:length(K)) {
    
    # true positives    
    bar <- foo[foo$K == K[j],]
    bar1 <- bar[!duplicated(bar$rep),]
    rowindex <- (i-1)*2
    summary[j+rowindex,3] <- nrow(bar1)
    
    #false positives
    qux <- baz[baz$K == K[j],]
    
    temp <- vector(mode="integer", length=10)
    
    for (k in 1:length(repl)) {
      rux <- qux[qux$rep == repl[k],]
      rux1 <- rux[!duplicated(rux$locus),] # removes any loci detected on >1 axis in same replicate
      temp[k] <- nrow(rux1)
    }
    
    summary[j+rowindex,4] <- sum(temp)
    summary[j+rowindex,5] <- sd(temp)
  }
}

### CHANGE OUTPUT FILE NAME: 
fname <- paste("/path/to/lfmm/output/", sim[s], "_50percent_500ind_Summary.csv", sep="")
write.csv(summary, file=fname, row.names = F)
