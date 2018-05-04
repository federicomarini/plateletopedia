




#' Fetch SRA RunInfo
#' 
#' Fetch the information on a RunInfo object from SRA
#' 
#' This function initializes the samplesinfo-like list as a standard component to
#' store all the necessary information.
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



