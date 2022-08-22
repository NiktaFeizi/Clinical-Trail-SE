## ----setup, include=FALSE-----------------------------------------------------
nthread <- 2
knitr::opts_chunk$set(echo = TRUE)


## -----------------------------------------------------------------------------
#' configure download options for this session
options(encoding = "UTF-8", timeout = 1000)
work_dir = setwd(paste0(getwd()))


## -----------------------------------------------------------------------------
#' fetch Gencode v33 annotations from BHKLAB-Pachyderm/Annotations
gencode_url <- "https://github.com/BHKLAB-Pachyderm/Annotations/raw/master/Gencode.v33.annotation.RData"
gencode_file <- file.path(work_dir, "Ensembl.v99.annotation.RData")
download.file(gencode_url, destfile=gencode_file)
load("Ensembl.v99.annotation.RData", verbose = TRUE)


## -----------------------------------------------------------------------------
#' required libraries
library(stringr)
library(GEOquery)
library(SummarizedExperiment)
library(qs)
library(affy)


## -----------------------------------------------------------------------------
#'
#'
#' @description
#' Warning: this function installs packages and thus modifies the global
#' state of your R environmnt.
#'
#' @param raw_tar Name of the tar file that contains raw expression values
#' @param geo_path Path to the expression tar file
#' @param brain_path Path to the microarray source CDF package which gets
#'   installed when this functon runs.
#'
#' @return RMA normalized expression data
#' @md
#' @export
geo_rma_norm <- function (raw_tar, geo_path , brain_path){
  
  # capture user options
  opts <- options()
  on.exit(options(opts))
  # set curl flags for donwload.file
  options(download.file.method="curl", download.file.extra="-k -L")
  
  #' Download CEL files
  cel.file <- file.path(work_dir, raw_tar)
  download.file(geo_path, cel.file)
  
  #' Unpack the CEL files
  exit_dir <- sub("\\..*", "", raw_tar)
  untar(raw_tar, exdir=exit_dir)
  cel.files<-list.files(file.path(exit_dir), pattern = "gz")
  sapply(file.path(exit_dir, cel.files), GEOquery::gunzip)
  
  #' Install the cdf brainarray package
  pack_name <- tail(stringr::str_split(brain_path, "/")[[1]], n=1)
  cdf.file <- file.path(work_dir, pack_name)
  download.file( brain_path, cdf.file)
  install.packages(pack_name, repos = NULL, type = "source")
  
  #' RMA normalization
  cdf <- sub("_.*", "", pack_name)
  cels <- affy::list.celfiles(file.path(exit_dir), full.names = TRUE)
  rma.norm <- affy::justRMA(filenames = cels, verbose = TRUE, cdfname = cdf)
  
  
  return(rma.norm)
  
}


## -----------------------------------------------------------------------------
#' ExpressionSet
modify_norm_eset <-  function(rma.norm, GSE , delim, GPL, delim2=NA, tissueid=NA, treatmentid = NA, response = NA, survival_type = NA,
                              survival_time = NA,  survival_time_unit=NA, event_occured = NA){
  
  #' Assay data
  assay_data <- as.data.frame(exprs(rma.norm))
  rownames(assay_data) <- sub("_at", "", rownames(assay_data))
  colnames(assay_data) <- sub(delim, "", colnames(assay_data))
  if (!identical(c(grep("AFFX", rownames(assay_data))), integer(0))) {
    assay_data <- assay_data[-c(grep("AFFX", rownames(assay_data))), ]
  } # Removing control genes (which start with "AFFX")
  
  #' Pheno data
  gset <- GEOquery::getGEO(GSE, GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset) > 1) idx <- grep(GPL, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  
  if (all(colnames(assay_data) %in% rownames(gset@phenoData@data))){
    pheno_data <- gset@phenoData@data[colnames(assay_data),]
  } else {colnames(assay_data) <-sub(delim2 ,"",colnames(assay_data))
  pheno_data <- gset@phenoData@data[colnames(assay_data),]}
  
  #' adding standard columns:
  pheno_data$patientid <- rownames(pheno_data)
  prior <- c("tissueid","treatmentid", "response", "survival_type", "survival_time", "survival_time_unit", "event_occured")
  
  for (c in prior){
    if (is.na(get(c))){
      pheno_data[[c]] <- NA
    } else {
      pheno_data[[c]] <- pheno_data [,  get(c)]
    }
  }
  
  #' Ordering columns based on priority
  pheno_data <- pheno_data[, c("patientid", prior, colnames(pheno_data)[!colnames(pheno_data) %in% prior])]
  
  #' Replacing "NA" with NA_character_ entries
  pheno_data[pheno_data=="NA"]=NA
  
  
  #' Feature data
  feat_data <- merge(features_gene, data.frame(gene_id = rownames(assay_data)), all.y=TRUE, by = "gene_id")
  rownames(feat_data) <- feat_data$gene_id
  
  eSet <- ExpressionSet(assayData = as.matrix(assay_data),
                        phenoData = AnnotatedDataFrame(pheno_data),
                        featureData = AnnotatedDataFrame(feat_data[rownames(assay_data),]))
}


## -----------------------------------------------------------------------------
#' ExpressionSet to SummarizedExperiment
eSetToSE <- function(eSet) {
  
  BiocGenerics::annotation(eSet) <- "rna"
  stopifnot(all(rownames(fData(eSet)) == rownames(exprs(eSet))))
  stopifnot(all(rownames(pData(eSet)) == colnames(exprs(eSet))))
  
  #' Build summarized experiment from eSet
  SE <- SummarizedExperiment::SummarizedExperiment(
    assays=SimpleList(as.list(Biobase::assayData(eSet))),
    #' Switch rearrange columns so that IDs are first, probes second
    rowData=S4Vectors::DataFrame(Biobase::fData(eSet)),
    colData=S4Vectors::DataFrame(Biobase::pData(eSet)),
    metadata=list("experimentData" = eSet@experimentData,
                  "annotation" = Biobase::annotation(eSet),
                  "protocolData" = Biobase::protocolData(eSet)))
  
  #' Extract names from expression set
  SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)
  
  stopifnot(all(rownames(colData(SE)) == rownames(pData(eSet))))
  stopifnot(all(rownames(rowData(SE)) == rownames(fData(eSet))))
  stopifnot(all(rownames(colData(SE)) == colnames(assay(SE))))
  
  return(SE)
}


## -----------------------------------------------------------------------------
#' Run all the functions
createSE <- function(raw_tar, geo_path, brain_path, delim, GPL, tissueid ,treatmentid, response, survival_type,
                     survival_time, survival_time_unit, event_occured, delim2=delim2){
  
  rma <- geo_rma_norm(raw_tar = raw_tar,
                      geo_path = geo_path,
                      brain_path = brain_path)
  
  GSE <- sub("_.*", "",raw_tar)
  eset <- modify_norm_eset(rma.norm = rma, GSE = GSE, delim = delim, GPL = GPL, delim2=delim2, tissueid=tissueid ,
                           treatmentid =treatmentid, response= response, survival_type= survival_type,
                           survival_time= survival_time, survival_time_unit=survival_time_unit,  event_occured= event_occured)
  SE <- eSetToSE(eset)
  
  return(SE)
}



## -----------------------------------------------------------------------------
#' ==== 1 ====
GSE41998_SE <- createSE(raw_tar = "GSE41998_RAW.tar",
                        geo_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE41998&format=file" ,
                        brain_path= "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133a2hsensgcdf_24.0.0.tar.gz",
                        delim = "_.*" , GPL = "GPL571", survival_type = NA, survival_time = NA, survival_time_unit=NA,
                        event_occured = NA, tissueid = NA, treatmentid = "treatment arm:ch1", response = "ac response:ch1")

colData(GSE41998_SE) [["tissueid"]] <- "Breast"


## -----------------------------------------------------------------------------
#' ==== 2 ====
GSE14671_SE <- createSE(raw_tar = "GSE14671_RAW.tar" ,
                        geo_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE14671&format=file" ,
                        brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133plus2hsensgcdf_24.0.0.tar.gz",
                        delim = ".CEL.*" , GPL = "GPL570",tissueid=NA, treatmentid = NA, response = NA, survival_type = NA,
                        survival_time = NA, survival_time_unit=NA,  event_occured = NA)

colData(GSE14671_SE) [["tissueid"]] <- "Myeloid"
colData(GSE14671_SE)[["treatmentid"]] <- "Imatinib"
colData(GSE14671_SE)[["response"]] <- ifelse(grepl("NonResponder", colData(GSE14671_SE)$title), "NonResponder", "Responder")



## -----------------------------------------------------------------------------
#' ==== 3 ====
GSE23554_SE <- createSE(raw_tar = "GSE23554_RAW.tar" ,
                        geo_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE23554&format=file" ,
                        brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133ahsensgcdf_24.0.0.tar.gz",
                        delim = ".CEL.*" , GPL = "GPL96", tissueid=NA,treatmentid = NA, survival_type = NA,
                        event_occured = NA, response = "cisplatin response (complete response or incomplete response):ch1",
                        survival_time_unit= NA, survival_time = "overall survival in days:ch1")

colData(GSE23554_SE)[["tissueid"]] <- "Ovary/Fallopian Tube"
colData(GSE23554_SE)[["treatmentid"]] <- "Cisplatin"
colData(GSE23554_SE)[["survival_type"]] <- "Overall_Survival"
colData(GSE23554_SE)[["survival_time_unit"]] <- "day"



## -----------------------------------------------------------------------------
#' ==== 4 ====
#' R and NR based on the overal survival time.
GSE14764_SE <- createSE(raw_tar = "GSE14764_RAW.tar" ,
                        geo_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE14764&format=file" ,
                        brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133ahsensgcdf_24.0.0.tar.gz",
                        delim = ".cel.*" , GPL = "GPL96",tissueid=NA, treatmentid = NA, response = "residual tumor:ch1", survival_type = NA,
                        survival_time = "overall survival time:ch1", survival_time_unit= NA, event_occured = "overall survival event:ch1")

colData(GSE14764_SE) [["tissueid"]] <- "Ovary/Fallopian Tube"
colData(GSE14764_SE) [["treatmentid"]] <- "Carboplatin:Paclitaxel"
colData(GSE14764_SE) [["survival_type"]] <- "Overall_Survival"
colData(GSE14764_SE) [["survival_time_unit"]]  <- "month"



## -----------------------------------------------------------------------------
#' ==== 5 ====
GSE20194_SE <- createSE(raw_tar = "GSE20194_RAW.tar" ,
                        geo_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE20194&format=file" ,
                        brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133ahsensgcdf_24.0.0.tar.gz",
                        delim = "_.*" , GPL = "GPL96",tissueid=NA, treatmentid = "treatment code:ch1" , response = "pcr_vs_rd:ch1",
                        survival_type = NA, survival_time = NA, survival_time_unit= NA,  event_occured = NA)

colData(GSE20194_SE)[["tissueid"]] <- "Breast"
colData(GSE20194_SE)$treatmentid[is.na(colData(GSE20194_SE)$treatmentid) & colData(GSE20194_SE)$treatments.comments.ch1 == "Taxol q.3 week x 4 FAC x 4"] <- "TFAC" 
colData(GSE20194_SE)$treatmentid[is.na(colData(GSE20194_SE)$treatmentid) & colData(GSE20194_SE)$treatments.comments.ch1 == "Taxol x 12 FAC x 4"] <- "TFAC"
colData(GSE20194_SE)$treatmentid[is.na(colData(GSE20194_SE)$treatmentid) & colData(GSE20194_SE)$treatments.comments.ch1 == "Taxol x 12 FEC x 4"] <- "TFEC"
colData(GSE20194_SE)$treatmentid[is.na(colData(GSE20194_SE)$treatmentid) & colData(GSE20194_SE)$treatments.comments.ch1 == "Taxol x 12  FEC x 4"] <- "TFEC"
colData(GSE20194_SE)$treatmentid[is.na(colData(GSE20194_SE)$treatmentid) & colData(GSE20194_SE)$treatments.comments.ch1 == "FEC x 3 ( d/t N.C)Taxol x 12"] <- "FECT" 
colData(GSE20194_SE)$treatmentid[is.na(colData(GSE20194_SE)$treatmentid)] <- colData(GSE20194_SE)$treatments.comments.ch1[is.na(colData(GSE20194_SE)$treatmentid)]

#' Reivsing drug names

t_code <- c("TFAC", "Tonly",
            "TH/FAC" , "TH/FEC", 
            "FAC" , "FEC",
            "TXFAC" , NA, 
            "TFEC", 
            "FECT" , "TFAC/HT", 
            "FACT", "FACT+XRT/X") # XRT is loco-regional radiation

t_name <- c( "Paclitaxel:5-Fluorouracil:Doxorubicin:Cyclophosphamide" , "Paclitaxel", 
             "Paclitaxel:Herceptin:5-Fluorouracil:Doxorubicin:Cyclophosphamide" , "Paclitaxel:Herceptin:5-Fluorouracil:Epirubicin:Cyclophosphamide",
             "5-Fluorouracil:Doxorubicin:Cyclophosphamide" , "5-Fluorouracil:Epirubicin:Cyclophosphamide",
             "Paclitaxel:Xeloda:5-Fluorouracil:Doxorubicin:Cyclophosphamide" , NA,
             "Paclitaxel:5-Fluorouracil:Epirubicin:Cyclophosphamide",
             "5-Fluorouracil:Epirubicin:Cyclophosphamide:Paclitaxel" , "Paclitaxel:5-Fluorouracil:Doxorubicin:Cyclophosphamide:Herceptin",
             "5-Fluorouracil:Doxorubicin:Cyclophosphamide:Paclitaxel" , "5-Fluorouracil:Doxorubicin:Cyclophosphamide:Paclitaxel:Xeloda:XRT")


for (i in seq(length(t_code))){
  
  colData(GSE20194_SE)$treatmentid[colData(GSE20194_SE)$treatmentid == t_code[i]] <- t_name[i] 
  
}


## -----------------------------------------------------------------------------
#' ==== 6 ====
GSE50948_SE <- createSE(raw_tar = "GSE50948_RAW.tar" ,
                        geo_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE50948&format=file" ,
                        brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133plus2hsensgcdf_24.0.0.tar.gz",
                        delim = "_.*" , GPL = "GPL570",  tissueid=NA, treatmentid = "treatment:ch1", response = "pcr:ch1",
                        survival_type = NA, survival_time = NA, survival_time_unit= NA, event_occured = NA)


colData(GSE50948_SE)[["tissueid"]] <- "Breast"

colData(GSE50948_SE)$treatmentid <- ifelse(colData(GSE50948_SE)$treatmentid == "neoadjuvant doxorubicin/paclitaxel (AT) followed by cyclophosphamide/methotrexate/fluorouracil (CMF)",
                                           "Doxorubicin:Paclitaxel:Cyclophosphamide:Methotrexate:Fluorouracil",
                                           "Doxorubicin:Paclitaxel:Cyclophosphamide:Methotrexate:Fluorouracil:Trastuzumab")

## -----------------------------------------------------------------------------
#' ==== 7 ====
GSE33072_SE <- createSE(raw_tar = "GSE33072_RAW.tar" ,
                        geo_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE33072&format=file" ,
                        brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hugene10sthsensgcdf_24.0.0.tar.gz",
                        delim = "_.*" , GPL = "GPL6244", delim2=".CEL.*",tissueid=NA, treatmentid = "treatment:ch1", response = NA,
                        survival_type = NA, survival_time = NA, survival_time_unit= NA, event_occured = NA)

colData(GSE33072_SE)[["tissueid"]] <- "Lung"
colData(GSE33072_SE) [["survival_type"]] <- "progression_free_survival"
colData(GSE33072_SE)$survival_time <- ifelse(is.na(colData(GSE33072_SE)$progression.free.survival.time..months..ch1), colData(GSE33072_SE)$pfsm..month..ch1, colData(GSE33072_SE)$progression.free.survival.time..months..ch1)
colData(GSE33072_SE)$event_occured <- ifelse(is.na(colData(GSE33072_SE)$progression.free.survival.status.ch1), colData(GSE33072_SE)$pfsc..1.progressed..0.not.progressed..ch1, colData(GSE33072_SE)$progression.free.survival.status.ch1)
colData(GSE33072_SE)[["survival_time_unit"]] <- "month"

# Drug name revision:
colData(GSE33072_SE)$treatmentid <- str_to_title(colData(GSE33072_SE)$treatmentid)
colData(GSE33072_SE)$treatmentid[colData(GSE33072_SE)$treatmentid == "Erlotinib+Bexarotene"] <- "Erlotinib:Bexarotene"


## -----------------------------------------------------------------------------
#' ==== 8 ====
GSE25066_SE <- createSE(raw_tar = "GSE25066_RAW.tar" ,
                        geo_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE25066&format=file" ,
                        brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133ahsensgcdf_24.0.0.tar.gz",
                        delim = "_.*" , GPL = "GPL96", tissueid=NA,treatmentid = NA, response = "pathologic_response_pcr_rd:ch1" ,
                        survival_type = NA, survival_time = "drfs_even_time_years:ch1",
                        survival_time_unit= NA, event_occured = "drfs_1_event_0_censored:ch1")

colData(GSE25066_SE)[["tissueid"]] <- "Breast"
colData(GSE25066_SE)[["treatmentid"]] <- "Taxane:Anthracycline"
colData(GSE25066_SE)[["survival_type"]] <- "Distant recurrence-free survival (DRFS)"
colData(GSE25066_SE)[["survival_time_unit"]] <- "year"



## -----------------------------------------------------------------------------
#' ==== 9 ====
GSE15622_SE <- createSE(raw_tar = "GSE15622_RAW.tar" ,
                        geo_path = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE15622&format=file" ,
                        brain_path = "http://mbni.org/customcdf/24.0.0/ensg.download/hgu133a2hsensgcdf_24.0.0.tar.gz",
                        delim = ".CEL" , GPL = "GPL8414", tissueid=NA,treatmentid = "treatment:ch1", response = "response:ch1" ,
                        survival_type = NA, survival_time = NA,
                        survival_time_unit= NA, event_occured = NA)

colData(GSE15622_SE)[["tissueid"]] <- "Ovary/Fallopian Tube"
colData(GSE15622_SE)$treatmentid[colData(GSE15622_SE)$treatmentid == "Both"] <- "Paclitaxel:Carboplatin"


## -----------------------------------------------------------------------------
#' ===== Save objects =======
SE_objects <- c("GSE41998_SE", "GSE14671_SE", "GSE23554_SE", "GSE14764_SE", "GSE20194_SE", "GSE50948_SE", "GSE33072_SE", "GSE25066_SE", "GSE15622_SE")
for (dataset in SE_objects) {qsave(get(dataset), file=paste0(dataset, ".qs"), nthread=nthread)}



## -----------------------------------------------------------------------------
#' ===== Automatically upload SE objects to Zenodo =======

#' To prevent reinstalling the package every time this script is run
if (!require(AnnotationGx)) remotes::install_github("bhklab/AnnotationGx@development")
library(AnnotationGx)

#' input the GSE id or name of the dataset
gse_id <- "GSE50948"
fill_meta <- zenodoMetadata(title = paste("Clinical Dataset SE objects - ", gse_id, sep = ""),
                            description = paste("Curated clinical datasets for Roche collaboration, Task 4a. Data source - ", gse_id, sep = ""))

#' Please make sure to create and store Zenodo access token in the environment. See ?depositZenodo
depo_zen <- depositZenodo(file_path = paste(gse_id,"_SE.qs",sep = ""), metadata = fill_meta, publish = TRUE)

