
### Functions for 'analysis.Rmd' ###



### Data handling functions ----

getMSigDB <- function(database, select = NULL){
  
  # Get gene sets from MSigDB
  
  MSigDB <- getGmt(database)
  genesets <- geneIds(MSigDB)
  
  if (!is.null(select)){
    genesets <- genesets[select]
  }
  
  return(genesets)
}



pvalStr <- function(p, digits = 4, min_p = 0.0001){
  
  # Pval used in plots
  
  if (p < min_p){
    str <- paste0("pval < ", min_p)
  } else {
    str <- paste0("pval = ", format.pval(p, digits))
  }
  return(str)
}



cbindIDs <- function(data1, data2){
  ids <- intersect(rownames(data1), rownames(data2))
  cdata <- data.frame(row.names = ids, data1[ids,], data2[ids,])
  return(cdata)
}



naf <- function(v){
  v[is.na(v)] <- FALSE
  return(v)
}



pss2date <- function(x) as.Date(x/86400, origin = "1582-10-14") # Convert PSS to data



getGeneSets <- function(){
  
  # Gene sets of interest
  
  genesets <- c("GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
                "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
                "GO_REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION",
                "GO_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY",
                "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN",
                "GO_RESPONSE_TO_TYPE_I_INTERFERON",
                "GO_ANTIGEN_PROCESSING_AND_PRESENTATION",
                "GO_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
                "GO_NEGATIVE_REGULATION_OF_T_CELL_DIFFERENTIATION",
                "GO_REGULATION_OF_ALPHA_BETA_T_CELL_PROLIFERATION",
                "GO_POSITIVE_REGULATION_OF_INTERLEUKIN_1_SECRETION",
                "GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                "GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
                "GO_RESPONSE_TO_GAMMA_RADIATION",
                "GO_REGULATION_OF_MACROPHAGE_ACTIVATION",
                "GO_REGULATION_OF_INTERLEUKIN_8_SECRETION",
                "GO_REGULATION_OF_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                "GO_CELL_KILLING",
                "GO_REGULATION_OF_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
                "GO_PEPTIDE_ANTIGEN_BINDING",
                "GO_REGULATION_OF_RESPONSE_TO_INTERFERON_GAMMA",
                "GO_NEGATIVE_REGULATION_OF_LYMPHOCYTE_DIFFERENTIATION",
                "GO_NEGATIVE_REGULATION_OF_INNATE_IMMUNE_RESPONSE",
                "GO_RESPONSE_TO_INTERFERON_GAMMA",
                "GO_T_CELL_RECEPTOR_COMPLEX",
                "GO_REGULATION_OF_LYMPHOCYTE_APOPTOTIC_PROCESS",
                "GO_REGULATION_OF_T_CELL_PROLIFERATION",
                "GO_IMMUNOLOGICAL_SYNAPSE",
                "GO_POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION",
                "GO_POSITIVE_REGULATION_OF_T_CELL_PROLIFERATION",
                "GO_REGULATION_OF_INTERLEUKIN_10_PRODUCTION",
                "GO_PURINERGIC_NUCLEOTIDE_RECEPTOR_SIGNALING_PATHWAY",
                "GO_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS",
                "GO_POSITIVE_REGULATION_OF_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                "GO_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY",
                "GO_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY",
                "GO_POSITIVE_REGULATION_OF_LYMPHOCYTE_APOPTOTIC_PROCESS",
                "GO_POSITIVE_REGULATION_OF_CELL_KILLING",
                "GO_POSITIVE_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE",
                "GO_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE",
                "GO_REGULATION_OF_INTERLEUKIN_1_SECRETION",
                "GO_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY",
                "GO_REGULATION_OF_INTERLEUKIN_2_BIOSYNTHETIC_PROCESS",
                "GO_POSITIVE_REGULATION_OF_INTERLEUKIN_12_PRODUCTION",
                "GO_POSITIVE_REGULATION_OF_INTERLEUKIN_10_PRODUCTION",
                "GO_NEGATIVE_REGULATION_OF_INTERLEUKIN_10_PRODUCTION",
                "GO_TOLL_LIKE_RECEPTOR_4_SIGNALING_PATHWAY",
                "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                "HALLMARK_INFLAMMATORY_RESPONSE",
                "HALLMARK_IL2_STAT5_SIGNALING",
                "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                "HALLMARK_TGF_BETA_SIGNALING",
                "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  
  return(genesets)
}



getCoxNames1 <- function(){
  
  # More readable variable names for report
  
  df <- data.frame("Name" = c("Age",
                              "ATG",
                              "Chemotherapy: YES vs NO",
                              "ER-Status: Positive vs Negative",
                              "HER2-status measured by SNP6: GAIN vs NEUTRAL",
                              "HER2-status measured by SNP6: LOSS vs NEUTRAL",
                              "HER2-Status: Positive vs Negative",
                              "Histologic Grade: II vs I",
                              "Histologic Grade: III vs I",
                              "Histologic Subtype: Lobular vs Ductal/NST",
                              "Histologic Subtype: Medullary vs Ductal/NST",
                              "Histologic Subtype: Mixed vs Ductal/NST",
                              "Histologic Subtype: Mucinous vs Ductal/NST",
                              "Histologic Subtype: Other vs Ductal/NST",
                              "Histologic Subtype: Tubular/cribriform vs Ductal/NST",
                              "IFNA",                                                
                              "IFNG",
                              "Inferred Menopausal State: Post vs Pre",
                              "Oncotree Code: IDC vs BREAST",
                              "Oncotree Code: ILC vs BREAST",
                              "Oncotree Code: IMMC vs BREAST",
                              "Oncotree Code: MDLC vs BREAST",
                              "PR-Status: Positive vs Negative",
                              "Stage",
                              "Tumor Size",
                              "Stage (I-IV): II vs I",
                              "Stage (I-IV): III vs I",
                              "Stage (I-IV): IV vs I"),
                   "Var" = c("Age",
                             "ATG",
                             "ChemotherapyYES",
                             "ER.StatusPositive",
                             "HER2.status.measured.by.SNP6GAIN",
                             "HER2.status.measured.by.SNP6LOSS",
                             "HER2.StatusPositive",
                             "Histologic.Grade2",
                             "Histologic.Grade3",
                             "Histologic.SubtypeLobular",
                             "Histologic.SubtypeMedullary",
                             "Histologic.SubtypeMixed",
                             "Histologic.SubtypeMucinous",
                             "Histologic.SubtypeOther",
                             "Histologic.SubtypeTubular/ cribriform",
                             "IFNA",
                             "IFNG",
                             "Inferred.Menopausal.StatePost",
                             "Oncotree.CodeIDC",
                             "Oncotree.CodeILC",
                             "Oncotree.CodeIMMC",
                             "Oncotree.CodeMDLC",
                             "PR.StatusPositive",
                             "Stage",
                             "Tumor.Size",
                             "Tumor.Stage2",
                             "Tumor.Stage3",
                             "Tumor.Stage4"),
                   stringsAsFactors = FALSE)
  
  rownames(df) <- df$Var
  return(df)
}



getCoxNames2 <- function(oldnames){
  
  # More readable variable names for report
  
  df <- data.frame("Name" = c("Age",
                              "Chemotherapy",
                              "ER-Status",
                              "HER2-Status",
                              "HER2-status measured by SNP6",
                              "Histologic Grade",
                              "Histologic Subtype",
                              "Inferred Menopausal State",
                              "PR-Status",
                              "Stage",
                              "Tumor Size",
                              "Stage (I-IV)"),
                   "Var" = c("Age",
                             "Chemotherapy",
                             "ER.Status",
                             "HER2.Status",
                             "HER2.status.measured.by.SNP6",
                             "Histologic.Grade",
                             "Histologic.Subtype",
                             "Inferred.Menopausal.State",
                             "PR.Status",
                             "Stage",
                             "Tumor.Size",
                             "Tumor.Stage"),
                   stringsAsFactors = FALSE)
  
    rownames(df) <- df$Var
    return(df[oldnames,1])
}



getReport <- function(data){
  
  report <- data.frame(data[,c("Variable", "H", "pval", "CI95_lower", "CI95_upper")])
  colnames(report) <- c("Variable", "HR", "pval", "CI95 lower", "CI95 upper")
  
  report$Variable[report$Variable == "prim.mtDNAratio"] <- "mtDNA/gDNA primary tumor"
  report$Variable[report$Variable == "Histologic.Grade"] <- "Histologic Grade: III vs I+II"
  report$Variable[report$Variable == "Histologic.Subtype"] <- "Histologic Subtype: Lobular vs Ductal/NST"
  report <- report[rowSums(is.infinite(data.matrix(report[,-1])) != 0) == 0,]
  
  return(report)
}




### Differential expression analysis functions ----


runLimma <- function(data, HiLo, alpha=0.05) {
  
  reflev <- unique(HiLo)
  reflev <- as.character(reflev[grep("Lo", reflev)])
  group <- relevel(HiLo,ref=reflev)
  design <- model.matrix(~group)
  
  fit <- lmFit(data, design)
  fit <- eBayes(fit)
  DEGs <- topTable(fit, coef=colnames(design)[2], adjust="BH", number=nrow(data))
  
  DEGs$DE <- NA
  DEGs$DE[intersect(which(DEGs$logFC>0), which(DEGs$adj.P.Val<alpha))] <- "up"
  DEGs$DE[intersect(which(DEGs$logFC<0), which(DEGs$adj.P.Val<alpha))] <- "down"
  
  return(DEGs)
  
}




### Score calculation functions ----


calcGSVAscores <- function(data, database, maxcores = 20){

  # Calculate IFNA/IFNG scores
  
  # Load gene sets
  allgs <- read.gmt(database)
  gs <- allgs[grep("HALLMARK_INTERFERON", names(allgs))] # select gene sets
  names(gs)[names(gs) == "HALLMARK_INTERFERON_ALPHA_RESPONSE"] <- "IFNalpha"
  names(gs)[names(gs) == "HALLMARK_INTERFERON_GAMMA_RESPONSE"] <- "IFNgamma"
  
  # Run GSVA
  gsva.results <- gsva(data.matrix(data), gs, method="gsva", kcdf="Gaussian", parallel.sz = min(maxcores, detectCores()))
  return(t(gsva.results))
}



getATGscore <- function(data, genes){
  
  # Calculate ATG scores
  
  AutophagyScore <- apply(data[genes,],2,function(tmp){
      exp(median(tmp))
    })
  
  return(AutophagyScore)
}



getMitoScore <- function(data, file, mitogenes = "MitoOnly"){
  
  # Calculate Mito scores
  
  mito <- read.csv(file, sep="\t")
  
  mito$mitoonly <- mito$Subcellular.location=="Mitochondria"
  mito.genes <- strsplit(as.character(mito$Gene.synonym), ", ")
  names(mito.genes) <- mito$Gene
  
  if (mitogenes == "MitoAll") {
    genes<-as.character(mito$Gene)
    
  } else if (mitogenes == "MitoOnly") {
    genes <- as.character(mito$Gene[mito$mitoonly])
    
  }
  
  for (j in 1:length(genes)) {
    
    cgene<-genes[j]
    if (! cgene %in% rownames(data)) {
      
      syn <- as.character(mito.genes[[cgene]])
      ind <- which(syn %in% rownames(data))
      if (length(ind)>0) genes[j] <- syn[ind[1]]
    }
  }
  
  genes <- intersect(genes, rownames(data))
  
  MitoScore <- apply(data[genes,],2,function(tmp){exp(mean(tmp))})
  
  return(MitoScore)
  
}




### Survival analysis functions ----


getHiLo <- function(x, str = NULL) {
  
  # Stratify into Hi/Lo groups
  
  classes <- rep("Lo", length(x))
  classes[which(x>median(x))] <- "Hi"
  
  if (!is.null(str)) classes <- paste(str, classes, sep="_")
  
  names(classes) <- names(x)
  classes <- as.factor(classes)
  lowclass <- levels(classes)[grep("Lo", levels(classes))]
  classes <- relevel(classes, ref=lowclass)
  
  return(classes)
  
}



getStrata <- function(scores, n = 2){
  
  # Stratify into any number of groups
  
  if (n < 2){
    return(NULL)
  }
  
  qs <- quantile(scores, seq_along(1:(n-1))/n, na.rm = TRUE)
  
  # class names
  if (n == 2){
    classes.names <- c("Lo", "Hi")
  } else if (n == 3) {
    classes.names <- c("Lo", "Mid", "Hi")
  } else {
    names(qs) <- paste0( as.character(round( 100*seq_along(1:(n-1))/n, 1)), "%") 
    tmpnames <- c()
    if (n > 2){
      for (i in 1:(length(qs)-1)){
        tmpnames[i] <- paste0(names(qs)[i], "-", names(qs)[i+1])
      }
    }
    classes.names <- c(paste0("<", names(qs)[1]), tmpnames, paste0(">", names(qs)[length(qs)]))
  }
  
  classes <- rep(classes.names[1], length(scores))
  for (i in 1:length(qs)){
    classes[ scores >= qs[i] ] <- classes.names[i+1]
  }
  
  classes <- factor(classes, levels = classes.names, ordered = TRUE)
  names(classes) <- names(scores)
  classes[is.na(scores)] <- NA
  
  return(classes)
}



getDiseaseSpecificSurvival <- function(clinicalData){
  
  status <- seq_along(clinicalData$Patient.s.Vital.Status)*NA
  names(status) <- clinicalData$Patient.ID
  
  status[ clinicalData$Patient.s.Vital.Status == "Living" ] <- 0
  status[ clinicalData$Patient.s.Vital.Status == "Died of Other Causes" ] <- NA
  status[ clinicalData$Patient.s.Vital.Status == "Died of Disease" ] <- 1
  
  status <- status[!is.na(status)]
  
  survData <- data.frame(status = status, time = clinicalData$Overall.Survival..Months.[ clinicalData$Patient.ID %in% names(status)])
  return(survData)
}



getCovData <- function(clinicalData, minObs = 5, addStage = TRUE){
  
  # Get covariates from clinical data for survival analysis
  
  covNames <- c("Age.at.Diagnosis",
                "Tumor.Stage",
                "Chemotherapy",
                "Hormone.Therapy",
                "Radio.Therapy",
                "Neoplasm.Histologic.Grade",
                "Tumor.Size",
                "ER.Status",
                "HER2.Status",
                "PR.Status",
                "HER2.status.measured.by.SNP6",
                "Oncotree.Code",
                "Tumor.Other.Histologic.Subtype",
                "Inferred.Menopausal.State")
  
  covdata <- clinicalData
  rownames(covdata) <- clinicalData$Patient.ID
  
  covNames <- covNames[covNames %in% colnames(covdata)]
  
  str <- paste(covNames[!covNames %in% colnames(covdata)], collapse = ", ")
  if (nchar(str)>0){
    message(paste0("Removed ", str))
  }
  
  covdata <- covdata[, covNames]
  
  
  # Stage
  if (addStage == TRUE){
    tmp <- covdata[,"Tumor.Stage"]
    covdata$Stage <- as.numeric(as.character(tmp))
    message("Stage added.")
  }
  
  # Chemotherapy
  if ("Chemotherapy" %in% covNames){
    tmp <- covdata[,"Chemotherapy"]
    tmp <- as.factor(tmp)
    tmp <- relevel(tmp, ref = "NO")
    covdata[,"Chemotherapy"] <- tmp
  }
  
  # ER.Status
  if ("ER.Status" %in% covNames){
    tmp <- covdata[,"ER.Status"]
    tmp <- as.factor(tmp)
    tmp <- relevel(tmp, ref = "Negative")
    covdata[,"ER.Status"] <- tmp
  }
  
  # HER2.Status
  if ("HER2.Status" %in% covNames){
    tmp <- covdata[,"HER2.Status"]
    tmp <- as.factor(tmp)
    tmp <- relevel(tmp, ref = "Negative")
    covdata[,"HER2.Status"] <- tmp
  }
  
  # PR.Status
  if ("PR.Status" %in% covNames){
    tmp <- covdata[,"PR.Status"]
    tmp <- as.factor(tmp)
    tmp <- relevel(tmp, ref = "Negative")
    covdata[,"PR.Status"] <- tmp
  }
  
  
  # Tumor.Stage
  if ("Tumor.Stage" %in% covNames){
    tmp <- covdata[,"Tumor.Stage"]
    tmp[ as.numeric( as.character(tmp)) == 0] <- NA
    tmp <- as.factor(tmp)
    tmp <- relevel(tmp, ref = "1")
    covdata[,"Tumor.Stage"] <- tmp
  }
  
  # Tumor.Other.Histologic.Subtype
  if ("Tumor.Other.Histologic.Subtype" %in% covNames){
    tmp <- covdata[,"Tumor.Other.Histologic.Subtype"]
    tmp <- as.character(tmp)
    tmp[ !tmp %in% names(table(tmp))[table(tmp) > minObs]] <- NA
    tmp <- relevel(factor(tmp), ref = "Ductal/NST")
    covdata[,"Tumor.Other.Histologic.Subtype"] <- tmp
  }
  
  
  # HER2.status.measured.by.SNP6
  if ("HER2.status.measured.by.SNP6" %in% covNames){
    tmp <- covdata[,"HER2.status.measured.by.SNP6"]
    tmp <- as.character(tmp)
    tmp[tmp == "UNDEF"] <- NA
    tmp <- relevel(factor(tmp), ref = "NEUTRAL")
    covdata[,"HER2.status.measured.by.SNP6"] <- tmp
  }
  
  
  # Oncotree.Code
  if ("Oncotree.Code" %in% covNames){
    tmp <- covdata[,"Oncotree.Code"]
    tmp <- as.character(tmp)
    tmp[tmp == "MBC"] <- NA
    tmp <- relevel(factor(tmp), ref = "BREAST")
    covdata[,"Oncotree.Code"] <- tmp
  }
  
  
  # Neoplasm.Histologic.Grade
  if ("Neoplasm.Histologic.Grade" %in% covNames){
    tmp <- covdata[,"Neoplasm.Histologic.Grade"]
    tmp <- as.factor(tmp)
    covdata[,"Neoplasm.Histologic.Grade"] <- tmp
  }
  
  
  # Inferred.Menopausal.State
  if ("Inferred.Menopausal.State" %in% covNames){
    tmp <- covdata[,"Inferred.Menopausal.State"]
    tmp <- tmp <- relevel(factor(tmp), ref = "Pre")
    covdata[,"Inferred.Menopausal.State"] <- tmp
  }
  
  # Simplify names
  colnames(covdata)[colnames(covdata) == "Age.at.Diagnosis"] <- "Age"
  colnames(covdata)[colnames(covdata) == "Radio.Therapy"] <- "Radiotherapy"
  colnames(covdata)[colnames(covdata) == "Neoplasm.Histologic.Grade"] <- "Histologic.Grade"
  colnames(covdata)[colnames(covdata) == "Tumor.Other.Histologic.Subtype"] <- "Histologic.Subtype"
  
  return(covdata)
}



orderLevels <- function(data){
  data$HER2.status.measured.by.SNP6 <- ordered(factor(data$HER2.status.measured.by.SNP6, levels = c("LOSS", "NEUTRAL", "GAIN")))
  data$Histologic.Subtype <- ordered(factor(data$Histologic.Subtype, levels = c("Other", "Tubular/ cribriform", "Mucinous",  "Medullary", "Lobular", "Ductal/NST", "Mixed")))
  return(data)
}



KaplanMeier <- function(data, classes, returnData = FALSE){
  
  ids <- intersect(names(classes), rownames(data))
  survdata <- data.frame(classes = classes[ids], status = data[ids,"status"], time = data[ids,"time"])
  
  surv <- Surv(time = survdata$time, event = survdata$status)
  KM <- survfit( surv ~ classes, data = survdata)
  
  if (returnData == FALSE){
    return(KM)
  } else {
    res.survdiff <- survdiff(surv ~ classes, data = survdata)
    pval <- 1-pchisq(res.survdiff$chisq, length(res.survdiff$n)-1)
    return(list(KM = KM, surv = surv, survdata = data.frame(surv = surv, survdata), pval = pval))
  }
}



UniCoxRegression <- function(survdata, data, vars = NULL, ...){
  
  # Prepare data
  ids <- intersect(rownames(data), rownames(survdata))
  coxdata <- data.frame(data[ids,])
  surv <- Surv(time = survdata[ids,]$time, event = survdata[ids,]$status)
  
  # Check levels
  coxdata <- coxdata[,minLevElem(coxdata) > 0 & minLev(coxdata)]
  
  if (is.null(vars)){
    coxvars <- colnames(coxdata)
  } else {
    coxvars <- vars[vars %in% colnames(coxdata)]
    if (any(vars[!vars %in% colnames(coxdata)])){
      message( paste("Removed variable(s) ", vars[!vars %in% colnames(coxdata)], collapse = "") )
    }
  }
  
  # COX regression
  COX.list <- lapply(coxvars, function(var){ coxph(as.formula(paste('surv ~', var)), data = coxdata, ...) })
  names(COX.list) <- coxvars
  
  return(COX.list)
  
}



MultiCoxRegression <- function(survdata, data, vars = NULL, interactions = NULL, returnData = FALSE, ...){
  
  # Prepare data
  ids <- intersect(rownames(data), rownames(survdata))
  coxdata <- data.frame(data[ids,])
  surv <- Surv(time = survdata[ids,]$time, event = survdata[ids,]$status)
  
  # Check levels
  coxdata <- coxdata[,minLevElem(coxdata) > 0 & minLev(coxdata)]
  
  # COX regression
  if (is.null(vars)){
    print( paste('surv ~', paste( colnames(coxdata), collapse = " + ")) )
    COX <- coxph( surv ~ ., data = coxdata, ...)
  } else {
    
    coxvars <- vars[vars %in% colnames(coxdata)]
    if (length((vars[!vars %in% colnames(coxdata)])) > 0){
      message( paste0("Removed variable(s) ", paste(vars[!vars %in% colnames(coxdata)], collapse = ", ")) )
    }
    print( paste('surv ~', paste( coxvars, collapse = " + ")) )
    COX <- coxph(as.formula( paste('surv ~', paste( coxvars, collapse = " + ") )), data = coxdata, ...)
    
  }
  
  if (returnData == TRUE){
    return(list(Cox = COX, Coxdata = coxdata))
  } else {
    return(COX)
  }
  
}



getCox <- function(r){
  tmp <- summary(r)
  tmp.H <- tmp$coefficients[,"exp(coef)"]
  tmp.pval <- tmp$coefficients[,"Pr(>|z|)"]
  tmp.CI_95_lower <- tmp$conf.int[,"lower .95"]
  tmp.CI_95_upper <- tmp$conf.int[,"upper .95"]
  n <- tmp$n
  LRT <- tmp$logtest["pvalue"]
  W <- tmp$waldtest["pvalue"]
  Sctest <- tmp$sctest["pvalue"]
  
  res.df <- data.frame(H = tmp.H, pval = tmp.pval, CI95_lower = tmp.CI_95_lower, CI95_upper = tmp.CI_95_upper, n = n, LRT.model = LRT, Wald.model = W, Sctest.model = Sctest)
  rownames(res.df) <- rownames(tmp$coefficients)
  return(res.df)
}



getCoxres <- function(res){
  
  for (i in 1:length(res)){
    tmp.res <- getCox(res[[i]])
    tmp.res$Variable <- names(res)[i]
    if (i == 1){
      results <- tmp.res
    } else {
      results <- rbind(results, tmp.res)
    }
  }
  return(results)
}



minLevElem <- function(data){
  apply(data, 2, function(tmp){
    min(as.data.frame(table(tmp))$Freq)
  })
}



minLev <- function(data){
  res <- vector(length = ncol(data))
  for (i in 1:ncol(data)){
    if (is.numeric(data[,i])){
      res[i] <- TRUE
    } else {
      res[i] <- nlevels(data[,i]) > 1
    }
  }
  return(res)
}



CovCorr <- function(testvars, survData.csos, coxData.csos){
  
  # Get ATG-covariate correlations and their impacts on the prognostic value of ATG
  
  res <- data.frame(matrix(nrow = length(testvars), ncol = 8))
  coxData.csos.ex <- orderLevels(coxData.csos)
  
  for (i in 1:length(testvars)){
    message(testvars[-i])
    cox.multi.ex <- MultiCoxRegression(survData.csos, coxData.csos.ex, vars = testvars[-i])
    res.cox.multi.ex <- getCox(cox.multi.ex)
    res[i,] <- res.cox.multi.ex["ATG",]
  }
  colnames(res) <- gsub("_", " ", colnames(res.cox.multi.ex), fixed = TRUE)
  rownames(res) <- testvars
  res <- res[rownames(res) != "ATG",]
  
  res$Correlation <- sapply(rownames(res), function(var){
    tmpdat <- coxData.csos.ex[,var]
    if (is.numeric(tmpdat)){
      cor.test(coxData.csos.ex$ATG, tmpdat, method = "spearman")$estimate
    } else {
      hetcor(coxData.csos.ex$ATG, tmpdat, ML = TRUE)$correlations[2,1]}
  })
  
  res$"Corr p value" <- sapply(rownames(res), function(var){
    tmpdat <- coxData.csos.ex[,var]
    if (is.numeric(tmpdat)){
      cor.test(coxData.csos.ex$ATG, tmpdat, method = "spearman")$p.value
    } else {
      pr <- polyserial(coxData.csos.ex$ATG, tmpdat, ML = TRUE, std.err = TRUE)
      2*pnorm(-abs(pr$rho / sqrt(pr$var[1,1])))
    }
  })
  
  res <- res[, !colnames(res) %in% c("n", "LRT.model", "Wald.model", "Sctest.model")]
  colnames(res)[1:2] <- paste0("ATG ", colnames(res)[1:2])
  res <- res[order(rownames(res)),]
  report.exclude <- cbind("Excluded variable" = rownames(res), res)
  report.exclude$`Excluded variable` <- getCoxNames2(rownames(report.exclude))
  
  return(report.exclude)
}




### GSEA input ----


prepareGSEAinput <- function(data, groups, basename, outputdir = "./"){
  
  GCTfilename <- paste0(outputdir, basename, "_expression.gct")
  CLSfilename <- paste0(outputdir, basename, "_phenotypes.cls")
  
  ### 1. Prepare expression data
  use.data <- data[,names(groups)]
  export.data <- cbind(rownames(use.data), "na", use.data)
  colnames(export.data)[1:2] <- c("NAME", "Description")
  l1 <- "#1.2"
  l2 <- nrow(export.data)
  l3 <- rbind(colnames(export.data), as.matrix(export.data))
  gct.list <- list(l1, data.frame(l2, length(groups)), l3)
  
  ### 2. Prepare phenotype data
  l1 <- paste(length(groups), 2, 1, sep = " ")
  l2 <- paste("#", paste(unique(groups), collapse = " "), sep = " ")
  l3 <- paste(groups, collapse = " ")
  cls.list <- list(l1, l2, l3)
  
  ### 3. Write output files
  if(file.exists(GCTfilename)){file.remove(GCTfilename)}
  if(file.exists(CLSfilename)){file.remove(CLSfilename)}
  
  tmp <- lapply(gct.list, function(x) write.table( x, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE, sep = "\t", file = GCTfilename))
  tmp <- lapply(cls.list, function(x) write.table( x, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE, sep = "\t", file = CLSfilename))

  return(list(GCT = GCTfilename, CLS = CLSfilename))
}



prepareGeneSets <- function(GeneSets, basename, outputdir = "./"){
  # Get gene sets of interests from gene set database and write to file
  GeneSetsCol <- GeneSetCollection(lapply(names(GeneSets), function(gs){ GeneSet(GeneSets[[gs]], setName = gs) }))
  GMTfilename <- paste0(outputdir, basename, "_genesets.gmt")
  toGmt(GeneSetsCol, GMTfilename) # write to file
  return(GMTfilename)
}



getGSEAoutput <- function(resultsdir){
  
  # get GSEA results files from this dir
  filenames <- list.files(path = resultsdir, pattern = ".xls")
  files <- filenames[ grep("gsea_report_for_", filenames)]
  groupnames <- sapply(  strsplit(files, split = "_"), function(tmp) tmp[4])
  
  res1 <- read.delim(paste0(resultsdir, files[1])) # up in group 1
  res2 <- read.delim(paste0(resultsdir, files[2])) # up in group 2
  res <- rbind(res1, res2)
  results <- list(res)
  names(results) <- paste(groupnames, collapse = "_vs_")
  return(results)
}



getNES <- function(allres){
  NES <- sapply(allres, function(tmp){
    tmp[,"NES"]
  })
  rownames(NES) <- allres[[1]]$NAME
  return(NES)
}



getQval <- function(allres){
  qval <- sapply(allres, function(tmp){
    tmp[,"FDR.q.val"]
  })
  rownames(qval) <- allres[[1]]$NAME
  return(qval)
}


getGSEAreport <- function(resultsdir){
  gsea.all <- getGSEAoutput(resultsdir)
  gsea.results <- cbindIDs(getNES(gsea.all), getQval(gsea.all))
  colnames(gsea.results) <- c("NES", "pval")
  gsea.results <- gsea.results[order(gsea.results$pval),]
  return(gsea.results)
}



### Plotting functions ----


HeatMap <- function(heatdata, titlestr, annDF, annColors, file, dendrogram = NULL){
  
  mycols <- colorRamp2(breaks = c(-2, 0, 2), colors=c("#313695", "white", "#d73027"))
  
  ha <- HeatmapAnnotation(df=annDF, col=annColors)
  
  if (!is.null(dendrogram)){
    ccol <- dendrogram
  } else {
    ccol <- TRUE
  }
  
  hm <- Heatmap(data.matrix(heatdata), 
                name="", 
                col=mycols,
                show_row_names=FALSE,
                show_column_names=FALSE,
                column_title=titlestr, 
                row_title="Top DE genes",
                top_annotation=ha,
                cluster_rows = FALSE,
                cluster_columns = ccol)
  
  
  pdf(file, onefile=TRUE)
    draw(hm)
  dev.off()
  
  
}



plotKaplanMeier <- function(KMobj, xlabel = "time", ylabel = "survival", colors = NULL, CI = FALSE, axislabsize = 15, ticklabsize = 13){

  if (is.null(colors)){
    colors <- c("orangered", "dodgerblue", "forestgreen", "darkmagenta", topo.colors(20))[1:nlevels(KM$survdata$classes)]
  } else {
    colors <- colors[levels(KM$survdata$classes)]
    names(colors) <- NULL
  }

  gp <- ggsurvplot( KMobj$KM,
                    data = KMobj$survdata,
                    risk.table = FALSE,
                    pval = pvalStr(KMobj$pval),
                    palette = colors,
                    ggtheme=theme_classic(),
                    risk.table.y.text.col = TRUE,
                    risk.table.y.text = FALSE,
                    xlab = xlabel,
                    ylab = ylabel,
                    font.x = axislabsize,
                    font.y = axislabsize,
                    font.tickslab = ticklabsize,
                    conf.int = ifelse(CI, TRUE, FALSE),
                    legend.title = "",
                    legend.labs =  levels(KM$survdata$classes))

  return(gp)
}



ViolinPlot <- function(data, vars, cols, ptext){
  
  groups <- data[,vars[1]]
  
  ggplot(data, aes(x=data[,vars[1]], y=data[,vars[2]], fill=groups)) + 
    geom_violin(trim=FALSE, show.legend=FALSE)+
    labs(title="", x = vars[1], y = vars[2], caption = ptext)+
    geom_boxplot(width=0.1, show.legend=FALSE,  outlier.size = 0.8)+
    theme_classic() + theme(text = element_text(size=20, colour=rgb(0,0,0)), axis.text = element_text(size=20, colour=rgb(0,0,0))) + scale_fill_manual(values=cols)
  
}



plotViolin <- function(scoreData, filename){
  
  minp <- 10^-10
  
  pdf(filename, onefile=TRUE, width = 8, height = 8)
  
  pvals <- matrix(nrow = 3, ncol = ncol(scoreData)-1)
  rownames(pvals) <- paste0("strat", 2:4)
  colnames(pvals) <- colnames(scoreData)[-1]
  for (i in 2:length(colnames(scoreData))){
    
    for (q in 2:4){
      
      if (q <= 2){
        res <- kruskal.test(scoreData[,i] ~ getStrata(scoreData$ATG, q))
        ptext <- paste0("pval=", format.pval(res$p.value, digits = 4, eps = minp))
        plotcols <- ATGcolors
      } else {
        res <- kruskal.test(scoreData[,i] ~ getStrata(scoreData$ATG, q))
        ptext <- paste0("pval=", format.pval(res$p.value, digits = 4, eps = minp))
        plotcols <- rev(colorRampPalette(ATGcolors)(q))
        names(plotcols) <- levels(getStrata(scoreData$ATG, q))
      }
      
      plotData <- cbind(ATG = getStrata(scoreData$ATG, q), scoreData[,-1,drop=FALSE])
      print(ViolinPlot(plotData, c("ATG", colnames(scoreData)[i]), cols = plotcols, ptext))
      
      pvals[q-1,i-1] <- res$p.value
    }
    
  }
  
  dev.off()
  
  return(pvals)
}


