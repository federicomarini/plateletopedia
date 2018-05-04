




#' Fetch SRA RunInfo
#' 
#' Fetch the information on a RunInfo object from SRA
#' 
#' This function initializes the samplesinfo-like list as a standard component to
#' store all the necessary information.
#' 
#' The fundamental elements are a runinfo data.frame, storing info like a SraRunInfo 
#' file (allowing also unpublished/in-house data to be included); a datasetID as 
#' unique identifier ("firstauthor_last-paper_topic-year"); and a data_dir, which
#' defaults to "_publicdata", as the overarching folder for storing all datasets
#'
#' @param SRAid The SRA identifier, normally the SRP id for the relevant project
#' @param datasetID  A string, ideally informative about the project/dataset it refers to.
#' An example could be "firstauthor_last-paper_topic-year". Defaults to the SRAid value, 
#' if not provided
#' @param data_dir A string, with the name of the main folder where the subfolder called 
#' \code{datasetID} will be created
#' @param outname A string, with the filename where the SraRunInfo file will be stored
#' @param quietDownload Logical, will define the behavior during the download via 
#' \code{download.file()}
#' @param force Logical. If the file to be created is already available, it can be 
#' overwritten by setting to TRUE. Defaults to FALSE
#'
#' @return A samplesinfo-like list to store all the required info and steps. This object
#' can be passed to the subsequent steps and further updated, thus containing all the 
#' information on the processing steps
#' 
#' @export
#'
#' @examples
#' samplesinfo <- fetch_sraruninfo("SRP098699",datasetID = "__TEST__mills_ingolia-pelo_decay-2017")
fetch_sraruninfo <- function(SRAid,
                             datasetID = SRAid,
                             data_dir = "_publicdata",
                             outname = paste0(SRAid,"_SraRunInfo.csv"),
                             quietDownload = FALSE,
                             force = FALSE) {
  
  if(is.null(datasetID))
    stop("You need to provide an ID for the dataset")
  
  if(file.exists(file.path(data_dir,datasetID,outname)) & !force) {
    message("The file for this SRAid and/or datasetID is already available. Reading in...")
    sri_file <- list.files(file.path(data_dir,datasetID),
                           pattern = "SraRunInfo.csv$", full.names = TRUE)
    sri <- read.delim(sri_file, header = TRUE, sep =",", as.is = TRUE)
    message(sri_file)
    
  } else {
    if(!dir.exists(file.path(data_dir,datasetID))) {
      dir.create(file.path(data_dir,datasetID),recursive = TRUE,showWarnings = FALSE)
      message("Directory created at ",file.path(data_dir,datasetID) )
    } else {
      message("Using directory ",file.path(data_dir,datasetID) )
    }
    # https://www.ncbi.nlm.nih.gov/books/NBK242621/
    # wget -O <file_name.csv> 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=<query>'
    download.file(url = paste0("http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=",
                               SRAid),
                  destfile = file.path(data_dir,datasetID,outname),quiet = quietDownload)
    message("SraRunInfo for SRA id ", SRAid, " downloaded at ",file.path(data_dir,datasetID,outname))
    sri_file <- list.files(file.path(data_dir,datasetID),
                           pattern = "SraRunInfo.csv$", full.names = TRUE)
    sri <- read.delim(sri_file, header = TRUE, sep =",", as.is = TRUE)
    
  }
  
  samplesinfo <- list(runinfo = sri,
                      datasetID = datasetID,
                      data_dir = data_dir,         # these first three arguments are mandatory - allow private data!
                      sri_file = sri_file)
  return(samplesinfo)
}



#' Fetch GEO info
#' 
#' Fetch the information on a series matrix object from GEO
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param GEOid The GEO identifier, normally the GSE id for the relevant project
#' @param datasetID A string, ideally informative about the project/dataset it refers to.
#' An example could be "firstauthor_last-paper_topic-year". Defaults to the datasetID 
#' element of the provided samplesinfo object
#' @param data_dir A string, with the name of the main folder where the subfolder called 
#' \code{datasetID} was created. Defaults to the data_dir element of the provided 
#' samplesinfo object
#' @param outname A string, with the filename where the series_matrix.txt.gz file 
#' will be stored
#' @param force Logical. If the file to be created is already available, it can be 
#' overwritten by setting to TRUE. Defaults to FALSE
#'
#' @return A samplesinfo-like list to store all the required info and steps. This object
#' can be passed to the subsequent steps and further updated, thus containing all the 
#' information on the processing steps
#' 
#' @export
#'
#' @examples
#' samplesinfo <- fetch_sraruninfo("SRP098699",datasetID = "__TEST__mills_ingolia-pelo_decay-2017")
#' samplesinfo <- fetch_geoinfo("GSE94384",samplesinfo = samplesinfo)
fetch_geoinfo <- function(samplesinfo,
                          GEOid,
                          datasetID = samplesinfo$datasetID,
                          data_dir = samplesinfo$data_dir,
                          outname = paste0(GEOid,"_series_matrix.txt.gz"),
                          force = FALSE) {
  
  if(is.null(datasetID))
    stop("You need to provide an ID for the dataset")
  if(file.exists(file.path(data_dir,datasetID,outname)) & !force) {
    message("The file for this GEOid and/or datasetID is already available. Reading in...")
    mygq <- GEOquery::getGEO(filename = file.path(data_dir,datasetID,outname),
                             GSEMatrix = TRUE,getGPL = FALSE)
  } else {
    if(!dir.exists(file.path(data_dir,datasetID))) {
      dir.create(file.path(data_dir,datasetID),recursive = TRUE,showWarnings = FALSE)
      message("Directory created at ",file.path(data_dir,datasetID) )
    } else {
      message("Using directory ",file.path(data_dir,datasetID) )
    }
    mygq <- GEOquery::getGEO(GEOid,
                             destdir = file.path(data_dir,datasetID),
                             GSEMatrix = TRUE,getGPL = FALSE)[[1]]
    message("Series matrix for GEO id ", GEOid, " downloaded at ",file.path(data_dir,datasetID,outname))
  }
  
  # extend the existing samplesinfo object
  samplesinfo$geoinfo <- mygq
  samplesinfo$geo_file <- file.path(data_dir,datasetID,outname)
  return(samplesinfo)
}




#' Create analysis folders
#' 
#' Creates all the subfolders in the specific datasetID folder
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param datasetID A string, ideally informative about the project/dataset it refers to.
#' An example could be "firstauthor_last-paper_topic-year". Defaults to the datasetID 
#' element of the provided samplesinfo object
#' @param data_dir A string, with the name of the main folder where the subfolder called 
#' \code{datasetID} was created. Defaults to the data_dir element of the provided 
#' samplesinfo object
#' 
#' @return A samplesinfo-like list to store all the required info and steps. This object
#' can be passed to the subsequent steps and further updated, thus containing all the 
#' information on the processing steps
#' 
#' @export
#'
#' @examples
#' samplesinfo <- fetch_sraruninfo("SRP098699",datasetID = "__TEST__mills_ingolia-pelo_decay-2017")
#' samplesinfo <- create_analysisfolders(samplesinfo)
create_analysisfolders <- function(samplesinfo,
                                   datasetID = samplesinfo$datasetID,
                                   data_dir = samplesinfo$data_dir) {
  
  message("Creating/checking the subfolder structure for ",datasetID," in the ",data_dir," main folder...")
  if(!dir.exists(file.path(data_dir,datasetID,"_sra")))
    dir.create(file.path(data_dir,datasetID,"_sra"))        # for SRA raw data
  if(!dir.exists(file.path(data_dir,datasetID,"_fastq")))
    dir.create(file.path(data_dir,datasetID,"_fastq"))      # for dumped/downloaded files
  if(!dir.exists(file.path(data_dir,datasetID,"_qc")))
    dir.create(file.path(data_dir,datasetID,"_qc"))         # for quality controls
  if(!dir.exists(file.path(data_dir,datasetID,"_aligned")))
    dir.create(file.path(data_dir,datasetID,"_aligned"))    # for aligned sam/bam/bai files
  if(!dir.exists(file.path(data_dir,datasetID,"_quants")))
    dir.create(file.path(data_dir,datasetID,"_quants"))     # for quantifications of expression
  if(!dir.exists(file.path(data_dir,datasetID,"_report")))
    dir.create(file.path(data_dir,datasetID,"_report"))     # well, for writing out reports
  
  # returns the unmodified object, so that the functions are ideally pipeable?
  return(samplesinfo)
}




#' Get data from SRA
#' 
#' Create a script to retrieve all the data from SRA for a specific project
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param datasetID A string, ideally informative about the project/dataset it refers to.
#' An example could be "firstauthor_last-paper_topic-year". Defaults to the datasetID 
#' element of the provided samplesinfo object
#' @param data_dir A string, with the name of the main folder where the subfolder called 
#' \code{datasetID} was created. Defaults to the data_dir element of the provided 
#' samplesinfo object
#' @param sra_runinfo A data frame with the info from a SRARunInfo file, defaults to the
#' runinfo element of the samplesinfo object
#' @param create_script Logical, whether to create the script - which needs to be run
#' from the terminal at a later moment
#' @param force Logical. If the file to be created is already available, it can be 
#' overwritten by setting to TRUE. Defaults to FALSE
#'
#' @return A samplesinfo-like list to store all the required info and steps. This object
#' can be passed to the subsequent steps and further updated, thus containing all the 
#' information on the processing steps
#' 
#' @export
#'
#' @examples
#' samplesinfo <- fetch_sraruninfo("SRP098699",datasetID = "__TEST__mills_ingolia-pelo_decay-2017")
#' samplesinfo <- create_analysisfolders(samplesinfo)
#' samplesinfo <- get_sradata(samplesinfo)
get_sradata <- function(samplesinfo,
                        datasetID = samplesinfo$datasetID,
                        data_dir = samplesinfo$data_dir,
                        sra_runinfo = samplesinfo$runinfo, # to override?
                        # samplenames_col = "SampleName",
                        create_script = TRUE,
                        force = FALSE) {
  
  # sri_file <- list.files(file.path(data_dir,datasetID),
  # pattern = "SraRunInfo.csv$", full.names = TRUE)
  sri <- sra_runinfo
  
  # grab all the samples - can be that some are listed more than once and have different SRR
  # samples <- unique(sri[, samplenames_col])
  library_types <- sri$LibraryLayout
  runs <- sri$Run
  runs_dlpath <- sri$download_path
  runs_files <- file.path(data_dir,datasetID,"_sra",paste0(runs,".sra"))
  cmd <- paste0("wget -c -O ",runs_files," ",runs_dlpath)
  out_bashscript <- file.path(data_dir,datasetID,paste0("cmd_batchretrieve_",unique(sri$SRAStudy),".sh"))
  # maybe prepend with shebang?
  if(create_script) {
    if(file.exists(out_bashscript) & !force) {
      message("Bash script for data retrieval has already been created at ", out_bashscript)
    } else {
      writeLines(cmd, con = out_bashscript)
      message("Bash script for data retrieval generated at ", out_bashscript)
    }
  }
  
  samplesinfo$files_sra <- runs_files
  samplesinfo$cmd_sraretrieve <- cmd
  
  # sradata_out <- list(runs_files = runs_files, library_types = library_types, cmd = cmd, sri = sri,
  # datasetID = datasetID)        
  
  return(samplesinfo)
  # system(cmd)
  
}




