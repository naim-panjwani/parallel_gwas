args=(commandArgs(TRUE))
infile <- as.character(args[1])
outfile <- as.character(args[2])

library(data.table)
library(dplyr)
#library(RNOmni)
#library(quantreg)
#library(QRank)
#source("code/gJLS_NP.R")
#library(statmod)
#library(car)
#source("code/rankscore_test.R")

dosage_filename <- infile
# dosage_filename <- "tempfiles2012243/imputed_data_topmed_chr21_split_00_withHeader"
fam_filename <- "data/intermediate_files/02-king_unrelated.fam"
covar_filename <- "data/intermediate_files/covar_file_noAge561.txt"
outfilename <- outfile


# newphenos <- fread("data/raw_data/21_06_18_sz_freq_Delnaz.csv")[,1:4] # has the myoclonus or absence sz. frequency
# myo_abs_freq <- as.integer(apply(newphenos, 1, 
#                                  function(x) ifelse(!is.na(x['myo_freq_cut']), max(c(x['myo_freq_cut'], x['abs_freq_cut_inc_none_ever']),na.rm=T), NA)))
# newphenos <- cbind(newphenos, myo_abs_freq)
# nrow(newphenos)
# head(newphenos)

# length(which(!is.na(newphenos$myo_abs_freq)))

fam <- fread(fam_filename, header=F, stringsAsFactors = F)

colslist <- c("FID","IID"
              ,"PC1","PC2","PC3", "PC4"
              #,"genocohort2","genocohort3","genocohort4"
              #,"intervention_drug", "intervention_seizurefree"
              #,"age_consent", "eur"
             )
covar <- fread(covar_filename, header=T, stringsAsFactors = F, sep="\t", quote="")[, ..colslist]
# the other covariates will be added later; running simpler regression first

print(paste0("Loading data for chr ", dosage_filename))
dosage <- fread(dosage_filename, header=T, stringsAsFactors = F, sep="\t", na.strings = c(".","NA"), skip="#CHROM")
isX <- gsub("chr","",tolower(as.character(dosage[1,"#CHROM"]))) %in% c("X","23")

print('Removing multi-allelic variants')
if(length(grep(",", dosage$ALT))>0) dosage <- dosage[-grep(",", dosage$ALT),] # remove multi-allelic variants

dt <- dosage[,10:ncol(dosage),with=F]

#print("Extracting dosage values")
#dosagemat <- t(apply(dt,1,function(y) sapply(y, function(x) as.numeric(sapply(tstrsplit(x, ":"),"[[",1)[2]))))
dosagemat <- dt

# Fixing ID's in imputation header file
x <- data.frame(FID=tstrsplit(colnames(dosagemat),"_",keep=1)[[1]], IID=tstrsplit(colnames(dosagemat),"_",keep=2)[[1]], stringsAsFactors = F)
x$IID <- gsub("-","_",x$IID)
x <- cbind(x, index=1:nrow(x))
x[x$FID=='125-Round2','FID'] <- '235-Round1' # this FID issue was later corrected in the genotype array data BEAGLE imputation, but not UK10K and TOPMED imputations

# Merge
x2 <- merge(x, covar, by=c("FID","IID"), sort=F) # the index position in dosagemat of samples to keep

x <- x2

x3 <- merge(x2, fam[,c(2,5,6)], by.x="IID",by.y="V2", sort=F)
setnames(x3,"V5","Sex")
setnames(x3,"V6","Photosens")
x3$Sex <- x3$Sex-1
x3$Photosens <- x3$Photosens-1
setcolorder(x3,"FID")

x <- x3

#dosagemat2 <- dosagemat[, x2$index]
dosagemat2 <- dosagemat[, x$index, with=F]
dosage2 <- cbind(dosage[,1:9], dosagemat2)

# add realmaf to info column
realmaf <- apply(dosagemat2, 1, function(x) sum(x)/(2*ncol(dosagemat2)))
realmaf <- ifelse(realmaf>0.5,1-realmaf,realmaf)
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k)) # for formatting number
dsginfo <- paste0(dosage2$INFO,";REALMAF=",specify_decimal(realmaf,5))
dosage2$INFO <- dsginfo

dosage <- dosage2
dosagemat <- dosagemat2

# MAF filter of 1%
dosage3 <- dosage2[realmaf>0.01]
dosagemat3 <- dosagemat2[realmaf>0.01]

dosage <- dosage3
dosagemat <- dosagemat3


x4 <- na.omit(x3[,c("FID","IID","Photosens","Sex",paste0("PC",1:4))])
x <- x4

model0 <- formula(paste("Photosens ~ ", "Sex + ", paste0("PC",1:4,collapse="+")))

logistic_assoc <- function(dosage_i, x, isX) {
    # x contains the phenotype (Y),
    # and index of sample VCF columns to include
    # _i in dosage_i refers to row i in dosage 
    # (basically a vcf file processed to only have dosage)
    
    tryCatch({
        results <- NULL
        g <- dosage_i[1,10:ncol(dosage_i)]
        snp_i <- dosage_i[,c(1:5,8)] # realmaf is accurate as it was subsetted earlier
        
        # merge g with x
        #print(paste0("dim(g): ",dim(g)))
        x5 <- data.frame(FID=tstrsplit(names(g),"_",keep=1)[[1]], IID=tstrsplit(names(g),"_",keep=2)[[1]], g=as.numeric(g), stringsAsFactors = F)
        x5$IID <- gsub("-","_",x5$IID)
        x5[x5$FID=='125-Round2','FID'] <- '235-Round1'
        #print(paste0("dim(x5)",dim(x5)))
        x6 <- merge(x5, x, on=c("FID","IID"), sort=F)
        x6 <- na.omit(x6)
        #print(paste0("dim(x6)",dim(x6)))
        
        n_i <- nrow(x6)
        #print(paste0("n_i: ",n_i))
        
        # GWAS
        model <- update(model0, ~. + g)
        fit <- glm(model, family="binomial", data=x6)
        
        coefs <- coef(summary(fit))
        result_i <- data.table(
          TEST="G"
          ,OBS_CT=n_i
          ,BETA=coefs['g','Estimate']
          ,SE=coefs['g','Std. Error']
          ,T_STAT=coefs['g', 'z value']
          ,P=coefs['g','Pr(>|z|)']
        )
        results <- rbind(results, cbind(snp_i, result_i))
        covar_results_mat <- coefs[2:(nrow(coefs)-1),]
        covar_results_i <- data.table(
          TEST=rownames(covar_results_mat)
          ,OBS_CT=n_i
          ,BETA=covar_results_mat[,'Estimate']
          ,SE=covar_results_mat[,'Std. Error']
          ,T_STAT=covar_results_mat[,'z value']
          ,P=covar_results_mat[,'Pr(>|z|)']
        )
        results <- rbind(results, cbind(snp_i, covar_results_i))
        out <- rbind(cbind(snp_i, result_i)
                    ,cbind(snp_i, covar_results_i))
        if(file.exists(outfilename)){
            fwrite(out
                   , outfilename
                   , append=T
                   , quote=F
                   , row.names=F
                   , col.names=F, sep="\t", na="NA")
        } else {
            fwrite(out
                   , outfilename
                   , append=T
                   , quote=F
                   , row.names=F
                   , col.names=T, sep="\t", na="NA")
        }
    })
    return(results)
}

print("Running association analysis")
print(paste("Writing output to", outfilename))
#print(paste0("nrow(x): ",nrow(x)))
assoc <- bind_rows(lapply(1:nrow(dosage), function(i) logistic_assoc(dosage[i,], x, isX)))

#print("Saving results")
#fwrite(assoc, outfilename, quote=F, row.names=F, col.names=T, sep="\t", na="NA")
