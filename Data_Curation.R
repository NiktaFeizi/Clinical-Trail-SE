
# Curating SE objects from clinical trial data

options(encoding = "UTF-8" , timeout = 1000)
work_dir = "~/Roche/Task4"
setwd(work_dir)
load("~/Data/Ensembl.v99.annotation.RData")

library(stringr)
library(GEOquery)
library(SummarizedExperiment)

#### Curation Functions ####
geo_rma<- function (raw_tar, goe_path , brain_path){
  
  # Download CEL files
  cel.file <- paste(work_dir, raw_tar, sep="/")
  download.file(goe_path, cel.file)
  
  # Unpack the CEL files
  exit_dir <- sub("\\..*", "", raw_tar)
  untar(raw_tar, exdir=exit_dir)
  cel.files<-list.files(paste(exit_dir, "/", sep=""), pattern = "gz")
  sapply(paste(exit_dir, cel.files, sep="/"), GEOquery::gunzip)
  
  # Install the cdf brainarray package
  pack_name <- tail(str_split(brain_path, "/")[[1]], n=1) 
  cfd.file <- paste(work_dir, pack_name, sep="/")
  download.file( brain_path, cfd.file)
  install.packages(pack_name, repos = NULL, type = "source")
  
  # RMA normalization
  cdf <- sub("_.*","",pack_name)
  cels <- affy::list.celfiles(paste(exit_dir, "/", sep=""), full.names = TRUE)
  rma.norm <- affy::justRMA(filenames = cels, verbose = TRUE, cdfname = cdf)
  
  return(rma.norm)
  
} 

# ExpressionSet 
ESet <-  function(rma.norm, GSE , delim, GPL, delim2=NA, cancer_type=NA, treatment = NA, response = NA, survival_type = NA,
                  survival_time = NA,  survival_time_unit=NA, event_occured = NA){
  
  # Assay data
  assay_data <- as.data.frame(exprs(rma.norm)) 
  rownames(assay_data) <-sub("_at","",rownames(assay_data)) # ask how to improve
  colnames(assay_data) <-sub(delim ,"",colnames(assay_data)) # ask how to improve
  if (!identical(c(grep("AFFX", rownames(assay_data))) , integer(0))) { 
    assay_data <- assay_data[-c(grep("AFFX", rownames(assay_data))),] } # Removing control genes (which start with "AFFX")
  
  # Pheno data
  gset <- GEOquery::getGEO(GSE, GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset) > 1) idx <- grep(GPL, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  
  if (all(colnames(assay_data) %in% rownames(gset@phenoData@data))){
    pheno_data <- gset@phenoData@data[colnames(assay_data),]
  } else {colnames(assay_data) <-sub(delim2 ,"",colnames(assay_data))
  pheno_data <- gset@phenoData@data[colnames(assay_data),]}
  
  # adding standard columns:
  pheno_data$unique_patient_ID <- rownames(pheno_data)
  prior <- c("cancer_type","treatment", "response", "survival_type", "survival_time", "survival_time_unit", "event_occured")
  
  for (c in prior){
    if (is.na(get(c))){
      pheno_data[[c]] <- NA
    } else {
      pheno_data[[c]] <- pheno_data [,  get(c)] 
    }
  }
  
  # Ordering columns based on priority
  pheno_data <- pheno_data[, c("unique_patient_ID", prior, colnames(pheno_data)[!colnames(pheno_data) %in% prior])]
  
  # Replacing "NA" with NA_character_ entries
  pheno_data[pheno_data=="NA"]=NA
  
  
  # Feature data
  feat_data <- merge(features_gene, data.frame(gene_id = rownames(assay_data)), all.y=T, by = "gene_id")
  rownames(feat_data) <- feat_data$gene_id
  
  eSet <- ExpressionSet(assayData = as.matrix(assay_data), 
                        phenoData = AnnotatedDataFrame(pheno_data), 
                        featureData = AnnotatedDataFrame(feat_data[rownames(assay_data),])) 
}

# ExpressionSet to SummarizedExperiment
eSetToSE <- function(eSet) {
  
  BiocGenerics::annotation(eSet) <- "rna"
  stopifnot(all(rownames(fData(eSet)) == rownames(exprs(eSet))))
  stopifnot(all(rownames(pData(eSet)) == colnames(exprs(eSet))))
  
  # Build summarized experiment from eSet
  SE <- SummarizedExperiment::SummarizedExperiment(
    assays=SimpleList(as.list(Biobase::assayData(eSet))),
    # Switch rearrange columns so that IDs are first, probes second
    rowData=S4Vectors::DataFrame(Biobase::fData(eSet)),
    colData=S4Vectors::DataFrame(Biobase::pData(eSet)),
    metadata=list("experimentData" = eSet@experimentData, 
                  "annotation" = Biobase::annotation(eSet), 
                  "protocolData" = Biobase::protocolData(eSet)))
  
  # Extract names from expression set                  
  SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)
  
  stopifnot(all(rownames(colData(SE)) == rownames(pData(eSet))))
  stopifnot(all(rownames(rowData(SE)) == rownames(fData(eSet))))
  stopifnot(all(rownames(colData(SE)) == colnames(assay(SE))))
  
  return(SE)
}

# Run all the functions
run_all <- function(raw_tar, goe_path, brain_path, delim, GPL, cancer_type ,treatment, response, survival_type,
                    survival_time, survival_time_unit, event_occured, delim2=delim2){
  
  rma <- geo_rma(raw_tar = raw_tar, 
                 goe_path = goe_path , 
                 brain_path= brain_path)
  
  GSE <- sub("_.*", "",raw_tar)
  eset <- ESet(rma.norm = rma, GSE = GSE, delim = delim, GPL = GPL, delim2=delim2, cancer_type=cancer_type ,
               treatment =treatment, response= response, survival_type= survival_type,
               survival_time= survival_time, survival_time_unit=survival_time_unit,  event_occured= event_occured)
  SE <- eSetToSE(eset)
  
  return(SE)
}


##### Curating summerized experiments #####
#==== 1 ====
GSE41998_SE <- run_all(raw_tar = "GSE41998_RAW.tar" , 
                       goe_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE41998&format=file" , 
                       brain_path= "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133a2hsensgcdf_24.0.0.tar.gz",
                       delim = "_.*" , GPL = "GPL571", survival_type = NA, survival_time = NA, survival_time_unit=NA, 
                       event_occured = NA, cancer_type = NA, treatment = "treatment arm:ch1", response = "ac response:ch1")

colData(GSE41998_SE) [["cancer_type"]] <- "Breast cancer"
colData(GSE41998_SE) $treatment <- colData(GSE41998_SE)$treatment.arm.ch1
colData(GSE41998_SE) $response <- colData(GSE41998_SE)$ac.response.ch1 

#==== 2 ====
GSE14671_SE <- run_all(raw_tar = "GSE14671_RAW.tar" , 
                       goe_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE14671&format=file" , 
                       brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133plus2hsensgcdf_24.0.0.tar.gz",
                       delim = ".CEL.*" , GPL = "GPL570",cancer_type=NA, treatment = NA, response = NA, survival_type = NA,
                       survival_time = NA, survival_time_unit=NA,  event_occured = NA)

colData(GSE14671_SE) [["cancer_type"]] <- "Chronic myelogenous leukemia (CML)"
colData(GSE14671_SE)[["treatment"]] <- "Imatinib"
colData(GSE14671_SE)[["response"]] <- ifelse(grepl("NonResponder", colData(GSE14671_SE)$title), "NonResponder", "Responder")

#==== 3 ====
GSE23554_SE <- run_all(raw_tar = "GSE23554_RAW.tar" , 
                       goe_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE23554&format=file" , 
                       brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133ahsensgcdf_24.0.0.tar.gz",
                       delim = ".CEL.*" , GPL = "GPL96", cancer_type=NA,treatment = NA, survival_type = NA,
                       event_occured = NA, response = "cisplatin response (complete response or incomplete response):ch1",
                       survival_time_unit= NA, survival_time = "overall survival in days:ch1")

colData(GSE23554_SE)[["cancer_type"]] <- "Ovarian Cancer"
colData(GSE23554_SE)[["treatment"]] <- "Cisplatin"
colData(GSE23554_SE)[["survival_type"]] <- "Overall_Survival"
colData(GSE23554_SE)[["survival_time_unit"]] <- "day"

#==== 4 ====
# R and NR based on the overal survival time.
GSE14764_SE <- run_all(raw_tar = "GSE14764_RAW.tar" , 
                       goe_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE14764&format=file" , 
                       brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133ahsensgcdf_24.0.0.tar.gz",
                       delim = ".cel.*" , GPL = "GPL96",cancer_type=NA, treatment = NA, response = NA, survival_type = NA,
                       survival_time = "overall survival time:ch1", survival_time_unit= NA, event_occured = "overall survival event:ch1")

colData(GSE14764_SE) [["cancer_type"]] <- "Ovarian Cancer"
colData(GSE14764_SE) [["treatment"]] <- "Carboplatin+Paclitaxel (adjuvant)"
colData(GSE14764_SE) [["survival_type"]] <- "Overall_Survival"
colData(GSE14764_SE) [["survival_time_unit"]]  <- "month"

#==== 5 ====
GSE20194_SE <- run_all(raw_tar = "GSE20194_RAW.tar" , 
                       goe_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE20194&format=file" , 
                       brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133ahsensgcdf_24.0.0.tar.gz",
                       delim = "_.*" , GPL = "GPL96",cancer_type=NA, treatment = "treatment code:ch1" , response = "pcr_vs_rd:ch1",
                       survival_type = NA, survival_time = NA, survival_time_unit= NA,  event_occured = NA)

colData(GSE20194_SE)[["cancer_type"]] <- "Breast cancer"
colData(GSE20194_SE)$treatment[is.na(colData(GSE20194_SE)$treatment) & colData(GSE20194_SE)$treatments.comments.ch1 == "Taxol x 12 FAC x 4"] <- "TFAC"
colData(GSE20194_SE)$treatment[is.na(colData(GSE20194_SE)$treatment) & colData(GSE20194_SE)$treatments.comments.ch1 == "Taxol x 12 FEC x 4"] <- "TFEC"
colData(GSE20194_SE)$treatment[is.na(colData(GSE20194_SE)$treatment) & colData(GSE20194_SE)$treatments.comments.ch1 == "Taxol x 12  FEC x 4"] <- "TFEC"
colData(GSE20194_SE)$treatment[is.na(colData(GSE20194_SE)$treatment) & colData(GSE20194_SE)$treatments.comments.ch1 == "FEC x 3 ( d/t N.C)Taxol x 12"] <- "FECT" # Double check
colData(GSE20194_SE)$treatment[is.na(colData(GSE20194_SE)$treatment)] <- colData(GSE20194_SE)$treatments.comments.ch1[is.na(colData(GSE20194_SE)$treatment)] 


#==== 6 ====
GSE50948_SE <- run_all(raw_tar = "GSE50948_RAW.tar" , 
                       goe_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE50948&format=file" , 
                       brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133plus2hsensgcdf_24.0.0.tar.gz",
                       delim = "_.*" , GPL = "GPL570",  cancer_type=NA, treatment = "treatment:ch1", response = "pcr:ch1", 
                       survival_type = NA, survival_time = NA, survival_time_unit= NA, event_occured = NA)

colData(GSE50948_SE)[["cancer_type"]] <- "Breast cancer"

#==== 7 ====
GSE33072_SE <- run_all(raw_tar = "GSE33072_RAW.tar" , 
                       goe_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE33072&format=file" , 
                       brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hugene10sthsensgcdf_24.0.0.tar.gz",
                       delim = "_.*" , GPL = "GPL6244", delim2=".CEL.*",cancer_type=NA, treatment = "treatment:ch1", response = NA, 
                       survival_type = NA, survival_time = NA, survival_time_unit= NA, event_occured = NA)

colData(GSE33072_SE)[["cancer_type"]] <- "Lung cancer"
colData(GSE33072_SE) [["survival_type"]] <- "progression_free_survival"
colData(GSE33072_SE)$survival_time <- ifelse(is.na(colData(GSE33072_SE)$progression.free.survival.time..months..ch1), colData(GSE33072_SE)$pfsm..month..ch1, colData(GSE33072_SE)$progression.free.survival.time..months..ch1)
colData(GSE33072_SE)[["survival_time_unit"]] <- "month"

#==== 8 ====
GSE25066_SE <- run_all(raw_tar = "GSE25066_RAW.tar" , 
                       goe_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE25066&format=file" , 
                       brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133ahsensgcdf_24.0.0.tar.gz",
                       delim = "_.*" , GPL = "GPL96", cancer_type=NA,treatment = NA, response = "pathologic_response_pcr_rd:ch1" , 
                       survival_type = NA, survival_time = "drfs_even_time_years:ch1",  
                       survival_time_unit= NA, event_occured = "drfs_1_event_0_censored:ch1")

colData(GSE25066_SE)[["cancer_type"]] <- "Breast cancer"
colData(GSE25066_SE)[["treatment"]] <- "taxane-anthracycline"
colData(GSE25066_SE)[["survival_type"]] <- "Distant recurrence-free survival (DRFS)"
colData(GSE25066_SE)[["survival_time_unit"]] <- "year"

#=====Save the 6 objects =======
qsave(GSE41998_SE, "GSE41998_SE.qs")
qsave(GSE14671_SE, "GSE14671_SE.qs")
qsave(GSE23554_SE, "GSE23554_SE.qs")
qsave(GSE14764_SE, "GSE14764_SE.qs")
qsave(GSE20194_SE, "GSE20194_SE.qs")
qsave(GSE50948_SE, "GSE50948_SE.qs")
qsave(GSE33072_SE, "GSE33072_SE.qs")
qsave(GSE25066_SE, "GSE25066_SE.qs")





