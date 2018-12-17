




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





#' Convert the sra files to fastq
#' 
#' Convert the sra archive files to compressed fastq files
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param datasetID A string, ideally informative about the project/dataset it refers to.
#' An example could be "firstauthor_last-paper_topic-year". Defaults to the datasetID 
#' element of the provided samplesinfo object
#' @param data_dir A string, with the name of the main folder where the subfolder called 
#' \code{datasetID} was created. Defaults to the data_dir element of the provided 
#' samplesinfo object
#' @param ncores A numeric value, with the number of cores to be used in the call to
#' fastq-dump, via the parallel-fastq-dump wrapper
#' @param create_script Logical, whether to create the script - which needs to be run
#' from the terminal at a later moment
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
#' samplesinfo <- sra_to_fastq(samplesinfo)
sra_to_fastq <- function(samplesinfo,
                         datasetID = samplesinfo$datasetID,
                         data_dir = samplesinfo$data_dir,
                         ncores = 12,
                         create_script = TRUE,
                         force = FALSE
) {
  sra_files <- samplesinfo$files_sra
  sra_libtypes <- samplesinfo$runinfo$LibraryLayout
  
  param_set <- "--gzip --skip-technical"

  outdir <- file.path(data_dir,datasetID,"_fastq")
  
  # fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR_ID
  # fastq-dump

  # parallel-fastq-dump, highly recommended
  # see https://github.com/rvalieris/parallel-fastq-dump
  cmd <- paste("parallel-fastq-dump --threads",ncores,
               param_set,
               ifelse(sra_libtypes=="SINGLE","","--split-files"), # make it depend on the single libraries
               "--outdir",outdir,"--sra-id",sra_files)
  
  out_bashscript <- file.path(data_dir,datasetID,paste0("cmd_batchfastqdump_",unique(samplesinfo$runinfo$SRAStudy),".sh"))
  # maybe prepend with shebang?
  if(create_script) {
    if(file.exists(out_bashscript) & !force) {
      message("Bash script for dumping sra to fastq.gz has already been created at ", out_bashscript)
    } else {
      writeLines(cmd, con = out_bashscript)
      message("Bash script for dumping sra to fastq.gz generated at ", out_bashscript)
    }
  }                
  # sratofastq_out <- list(sra_files = sra_files, library_types = sra_libtypes, cmd = cmd, sri = samplesinfo$sri,
  # datasetID = datasetID)                               
  samplesinfo$cmd_fastqdump <- cmd
  return(samplesinfo)
  
}





#' Match fastq to sra files
#' 
#' Match the fastq files to the original sra files
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
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
#' samplesinfo <- get_sradata(samplesinfo)
#' samplesinfo <- sra_to_fastq(samplesinfo)
#' samplesinfo <- match_fastq(samplesinfo)
match_fastq <- function(samplesinfo,
                        data_dir = samplesinfo$data_dir) {
  files_sra <- samplesinfo$files_sra
  sra_libtypes <- samplesinfo$runinfo$LibraryLayout
  
  allfastqs <- list.files(file.path(data_dir,samplesinfo$datasetID,"_fastq"),
                          pattern = ".fastq.gz$",
                          full.names = TRUE)
  
  if(length(allfastqs)==0) {
    message("No fastq.gz files were found in the corresponding folders. ",
            "Did you run the scripts to retrieve and convert the sra files?")
    return(samplesinfo)
  }
  
  # if libtype is single, one file is expected; otherwise it should be two
  list_matchedfastqs <- lapply(seq_len(length(files_sra)), function(arg){
    thissrafile <- basename(files_sra[arg])
    thislib <- sra_libtypes[arg]
    if(thislib=="SINGLE") {
      # should be only one
      single_fastq <- allfastqs[grep(gsub(".sra","",thissrafile),allfastqs)]
      thisfastqfiles <- single_fastq
    } else {
      # should be two
      paired_fastq <- allfastqs[grep(gsub(".sra","",thissrafile),allfastqs)]
      thisfastqfiles <- list(r1 = paired_fastq[1], r2 = paired_fastq[2])
    }
    return(thisfastqfiles)
    
  }) 
  names(list_matchedfastqs) <- gsub(".sra","",basename(files_sra))
  # return(list_matchedfastqs)
  samplesinfo[["files_fastq"]] <- list_matchedfastqs
  message("Performed assignment of fastq files to corresponding sra archives. ")
  return(samplesinfo)
}


#' Merge fastq files of the same sample
#' 
#' Merge fastq files of the same sample, concatenating the different fastq of each run
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param data_dir A string, with the name of the main folder where the subfolder called 
#' \code{datasetID} was created. Defaults to the data_dir element of the provided 
#' samplesinfo object
#' @param run_commands Logical, whether to run the commands to merge together fastq files 
#' belonging to the same sample name, but sequenced in different runs (SRR)
#'
#' @return A samplesinfo-like list to store all the required info and steps. This object
#' can be passed to the subsequent steps and further updated, thus containing all the 
#' information on the processing steps
#' @export
#'
#' @examples
mergetech_fastq <- function(samplesinfo,
                            data_dir = samplesinfo$data_dir,
                            run_commands = TRUE) {
  
  samplenames <- samplesinfo$runinfo$SampleName
  message("Found ", length(unique(samplenames)), " different samples in the ", nrow(samplesinfo$runinfo), " runs...")
  
  aftermerging <- vector("list", length(unique(samplenames)))
  aftermerging_names <- character(length(unique(samplenames)))
  for(i in 1:length(unique(samplenames))){
    corresp_fastq <- samplesinfo$files_fastq[which(samplenames == unique(samplenames)[i])]
    if(length(corresp_fastq) >= 2) {
      merged_file <- file.path(samplesinfo$data_dir,samplesinfo$datasetID,"_fastq",
                               paste0(paste0(names(corresp_fastq),collapse = "-"),".fastq.gz"))
      merge_cmd <- paste("cat", 
                         paste(corresp_fastq, collapse = " "),
                         ">", merged_file
      )
      message(merge_cmd)
      file.exists(merged_file)
      # if(run_commands) {
      # system(merge_cmd)
      # }
      # once done, remove originals?
      
      aftermerge_samplename <- file.path(samplesinfo$data_dir,samplesinfo$datasetID,"_fastq",
                                         paste0(paste0(names(corresp_fastq),collapse = "-"),".fastq.gz"))
      names(aftermerge_samplename) <- paste0(names(corresp_fastq),collapse = "-")
    } else {
      aftermerge_samplename <- corresp_fastq
      names(aftermerge_samplename) <- names(corresp_fastq)
    }
    aftermerging[i] <- aftermerge_samplename
    aftermerging_names[i] <- names(aftermerge_samplename)
  }
  names(aftermerging) <- aftermerging_names
  
  samplesinfo[["files_fastq_proc"]] <- aftermerging
  return(samplesinfo)
}



## TODO: will need a function for alternative entry point with own data
## ideally: fetch_sraruninfo - like
##          create_analysisfolders
##          match_fastq - like, to populate the samplesinfo$files_fastq slot



#' Create a samplesinfo object from local data
#' 
#' Initialize a samplesinfo object, from a non-public data repository.
#'
#' @param datasetID  A string, ideally informative about the project/dataset it refers to.
#' An example could be "firstauthor_last-paper_topic-year".
#' @param data_dir A string, with the name of the main folder where the subfolder called 
#' \code{datasetID} is supposed to be found, plus where the rest of the computations would
#' be performed
#' @param samples_metadata The file containing information in a format similar to a SRARunInfo 
#' file. This is expected to be in csv format
#' @param files_fastq A vector containing the full names of the fastq(.gz) files. Typically
#' the output of a call to \code{list.files(..., full.names = TRUE)} 
#' @param study_ID A string identifier, thought to mimicry the SRP id normally provided 
#' from SRA. Recommended to be short, e.g. 'F07', if related to Project folder F07
#'
#' @return
#' @export
#'
#' @examples
initialize_samplesinfo <- function(datasetID = NULL, 
                                   data_dir = NULL,
                                   samples_metadata = NULL,
                                   files_fastq = NULL,
                                   study_ID = NULL) {
  # checks:
  stopifnot(!is.null(datasetID))
  stopifnot(!is.null(data_dir))
  stopifnot(!is.null(samples_metadata))
  
  samplesinfo <- list()
  
  stopifnot(file.exists(samples_metadata))
  
  message("Reading in the sample metadata provided...")
  
  runinfo <- read.delim(samples_metadata, header = TRUE, sep =",", as.is = TRUE)
  # need to provide at least the LibraryLayout, as in the SRA run info files
  if(!"LibraryLayout" %in% names(runinfo))
    stop("No column with name 'LibraryLayout' provided! Subsequent functions might not work correctly, please provide this information!")
  # checking the content of the LibraryLayout: must be 'SINGLE' or 'PAIRED'
  
  if(!"Run" %in% names(runinfo))
    stop("No column with name 'Run' provided! Please provide such information!")
  
  if(!"ScientificName" %in% names(runinfo))
    stop("No column with name 'ScientificName' provided! Please provide the species information!")
  if(!"TaxID" %in% names(runinfo))
    stop("No column with name 'TaxID' provided! Please provide the species information!")
  
  if(length(runinfo$Run) != length(unique(runinfo$Run)))
    stop("No unique names provided in the 'Run' column")
  
  # adding also a fake "SRAStudy" column for compatibility with the other functions
  runinfo$Study <- study_ID  
  runinfo$SRAStudy <- study_ID
    
  samplesinfo[["runinfo"]] <- runinfo
  samplesinfo[["datasetID"]] <- datasetID
  samplesinfo[["data_dir"]] <- data_dir
  samplesinfo[["srifile"]] <- samples_metadata  
  
  # check that all data provided actually exist
  if(!all(file.exists(files_fastq)))
    stop("Error: not all files provided seem to exist. Possible reason: you did not specify the full path? Retry providing the output of 'list.files(..., full.names = TRUE)'")
  
  # if libtype is single, one file is expected; otherwise it should be two
  list_matchedfastqs <- lapply(seq_len(length(runinfo$Run)), function(arg){
    thislib <- runinfo$LibraryLayout[arg]
    thisrun <- runinfo$Run[arg]
    if(thislib=="SINGLE") {
      # should be only one
      single_fastq <- files_fastq[grep(thisrun,files_fastq)]
      thisfastqfiles <- single_fastq
    } else {
      # should be two
      paired_fastq <- files_fastq[grep(thisrun,files_fastq)]
      thisfastqfiles <- list(r1 = paired_fastq[1], r2 = paired_fastq[2])
    }
    return(thisfastqfiles)
  })
  
  names(list_matchedfastqs) <- runinfo$Run
  
  samplesinfo[["files_fastq"]] <- list_matchedfastqs
 
  return(samplesinfo)
}





#' Run FastQC for quality control
#' 
#' Creates a script for running FastQC on the fastq.gz files
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param fastqc_bin A string, with the location of the FastQC binary. This needs to be 
#' provided according to the local installation. Defaults to fastqc, as installed via conda
#' @param fastqc_ncores A numeric value, with the number of cores to be used in the call
#' to FastQC
#' @param fastqc_dir A string, where to store the output of the tool. Defaults to _qc,
#' created via create_analysisfolders
#' @param create_script Logical, whether to create the script - which needs to be run
#' from the terminal at a later moment
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
#' samplesinfo <- sra_to_fastq(samplesinfo)
#' samplesinfo <- match_fastq(samplesinfo)
#' samplesinfo <- run_fastqc(samplesinfo)
run_fastqc <- function(samplesinfo, # contains the locations of each file/file pair
                       fastqc_bin = "fastqc",
                       fastqc_ncores = 4,
                       fastqc_dir = "_qc",
                       create_script = TRUE,
                       force = FALSE
) {
  all_fastqsamples <- unlist(samplesinfo$files_fastq)
  fastqc_calls <- paste(fastqc_bin,
                        "--threads", fastqc_ncores,
                        "--outdir ",file.path(samplesinfo$data_dir,samplesinfo$datasetID,fastqc_dir),
                        all_fastqsamples
  )
  names(fastqc_calls) <- names(all_fastqsamples)
  
  out_bashscript <- file.path(samplesinfo$data_dir,
                              samplesinfo$datasetID,
                              paste0("cmd_batchFastQCrun_",unique(samplesinfo$runinfo$SRAStudy),".sh"))
  # maybe prepend with shebang?
  cmd <- (fastqc_calls)
  
  if(create_script) {
    if(file.exists(out_bashscript) & !force) {
      message("Bash script for running FastQC generated at ", out_bashscript)
    } else {
      writeLines(cmd, con = out_bashscript)
      message("Bash script for running FastQC generated at ", out_bashscript)
    }
  } 
  
  samplesinfo$cmd_fastqc <- cmd
  
  return(samplesinfo)
}


run_fastqc_multiT <- function(samplesinfo, # contains the locations of each file/file pair
                       fastqc_bin = "fastqc",
                       fastqc_ncores = 4,
                       fastqc_dir = "_qc",
                       create_script = TRUE,
                       force = FALSE
) {
  # check that files are there and do exist
  if(is.null(samplesinfo$files_fastq)) {
    warning("No fastq files provided in the samplesinfo object!")
    return(samplesinfo)
  }
  if(!all(file.exists(unlist(samplesinfo$files_fastq)))) {
    warning("Not all fastq files of the samplesinfo object are actually existing!")
    return(samplesinfo)
  }
  
  ## TODO: check that the output is not already there
  
  all_fastqsamples <- unlist(samplesinfo$files_fastq)
  # fastqc_calls <- paste(fastqc_bin,
  #                       "\\ \n  --threads", fastqc_ncores,
  #                       "\\ \n  --outdir",file.path(samplesinfo$data_dir,samplesinfo$datasetID,fastqc_dir),
  #                       "\\ \n ",paste(all_fastqsamples, collapse = " \\ \n  ")
  # )
  fastqc_calls <- paste(fastqc_bin,"--threads", fastqc_ncores,
                        "--outdir",file.path(samplesinfo$data_dir,samplesinfo$datasetID,fastqc_dir),
                        paste(all_fastqsamples, collapse = " ")
  )
  # names(fastqc_calls) <- names(all_fastqsamples)
  
  out_bashscript <- file.path(samplesinfo$data_dir,
                              samplesinfo$datasetID,
                              paste0("cmd_batchFastQCrun_M_",unique(samplesinfo$runinfo$SRAStudy),".sh"))
  # maybe prepend with shebang?
  cmd <- (fastqc_calls)
  
  if(create_script) {
    if(file.exists(out_bashscript) & !force) {
      message("Bash script for running FastQC already available at ", out_bashscript)
    } else {
      writeLines(cmd, con = out_bashscript)
      message("Bash script for running FastQC generated at ", out_bashscript)
    }
  } 
  
  samplesinfo$cmd_fastqc <- cmd
  
  return(samplesinfo)
}




#' Run salmon
#' 
#' Creates a script for running salmon and generating transcript-level quantifications
#' 
#' Providing two indices allows the handling of data sets where more than one species
#' is available, without the need to split the set into smaller subsets
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param salmon_bin A string, with the location of the salmon binary. This needs to be 
#' provided according to the local installation. Defaults to salmon, as installed via conda
#' @param salmon_index_human A string, the location of the salmon index for human samples
#' @param salmon_index_mouse A string, the location of the salmon index for mouse samples
#' @param salmon_dir A string, where to store the output of the tool. Defaults to _quants,
#' created via create_analysisfolders
#' @param salmon_ncores A numeric value, with the number of cores to be used in the call
#' to salmon
#' @param salmon_nbootstraps Numeric, the number of bootstrap samples to generate
#' @param salmon_gcbias Logical, perform sequence-specific bias correction
#' @param create_script Logical, whether to create the script - which needs to be run
#' from the terminal at a later moment
#' @param force 
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
#' samplesinfo <- sra_to_fastq(samplesinfo)
#' samplesinfo <- match_fastq(samplesinfo)
#' samplesinfo <- run_fastqc(samplesinfo)
#' samplesinfo <- run_salmon(samplesinfo)
run_salmon <- function(samplesinfo, # contains the locations of each file/file pair
                       salmon_bin = "salmon",
                       salmon_index_human,
                       salmon_index_mouse, # both, to correctly handle sets with >1 species
                       salmon_dir = "_quants",
                       salmon_ncores = 12,
                       salmon_nbootstraps = 100,
                       salmon_gcbias = TRUE,
                       create_script = TRUE,
                       force = FALSE
) {
  # check that files are there and do exist
  if(is.null(samplesinfo$files_fastq)) {
    warning("No fastq files provided in the samplesinfo object!")
    return(samplesinfo)
  }
  if(!all(file.exists(unlist(samplesinfo$files_fastq)))) {
    warning("Not all fastq files of the samplesinfo object are actually existing!")
    return(samplesinfo)
  }
  
  # check existence of indices
  if(!dir.exists(salmon_index_human))
    stop("salmon index for human not found!")
  if(!dir.exists(salmon_index_mouse))
    stop("salmon index for mouse not found!")
  
  
  N <- nrow(samplesinfo$runinfo)
  salmon_calls <- lapply(seq_len(N), function (i){
    this_libtype <- samplesinfo$runinfo$LibraryLayout[i]
    this_species <- samplesinfo$runinfo$ScientificName[i]
    this_index <- if(this_species == "Homo sapiens") {
      salmon_index_human
    } else if(this_species == "Mus musculus") {
      salmon_index_mouse
    } else {
      message("Found a dataset with a species currently not supported...")
      return(paste("#",samplesinfo$runinfo$Run[i],"Found a dataset with a species currently not supported..."))
    }
    # this_fastqset <- 
    if (this_libtype == "SINGLE") {
      this_fastqset <- samplesinfo$files_fastq[[i]]
      salmon_call <- paste(salmon_bin,"quant",
                           "--threads",salmon_ncores,
                           "--numBootstraps",salmon_nbootstraps,
                           ifelse(salmon_gcbias, "--seqBias", ""),
                           "--libType A",
                           "--index", this_index,
                           "-r",this_fastqset,
                           "-o",file.path(samplesinfo$data_dir,samplesinfo$datasetID,salmon_dir,
                                          paste0(names(samplesinfo$files_fastq)[i],"_salmon")))
    } else if (this_libtype == "PAIRED"){      
      this_fastqset <- samplesinfo$files_fastq[[i]]
      salmon_call <- paste(salmon_bin,"quant",
                           "--threads",salmon_ncores,
                           "--numBootstraps",salmon_nbootstraps,
                           ifelse(salmon_gcbias, "--seqBias", ""),
                           "--libType A",
                           "--index", this_index,
                           "-1",this_fastqset$r1,"-2",this_fastqset$r2,
                           "-o",file.path(samplesinfo$data_dir,samplesinfo$datasetID,salmon_dir,
                                          paste0(names(samplesinfo$files_fastq)[i],"_salmon")))
    }
    
    
    return(salmon_call)
  })
  names(salmon_calls) <- names(samplesinfo$files_fastq)
  
  out_bashscript <- file.path(samplesinfo$data_dir,
                              samplesinfo$datasetID,
                              paste0("cmd_batchsalmonrun_",unique(samplesinfo$runinfo$SRAStudy),".sh"))
  # maybe prepend with shebang?
  cmd <- unlist(salmon_calls)
  
  if(create_script) {
    if(file.exists(out_bashscript) & !force) {
      message("Bash script for running salmon already available at ", out_bashscript)
    } else {
      writeLines(cmd, con = out_bashscript)
      message("Bash script for running salmon generated at ", out_bashscript)
    }
  } 
  
  samplesinfo$cmd_salmon <- cmd
  
  
  return(samplesinfo)
}




#' Run kallisto
#' 
#' Creates a script for running kallisto and generating transcript-level quantifications
#' 
#' Providing two indices allows the handling of data sets where more than one species
#' is available, without the need to split the set into smaller subsets
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param kallisto_bin A string, with the location of the kallisto binary. This needs to be 
#' provided according to the local installation. Defaults to kallisto, as installed via conda
#' @param kallisto_index_human A string, the location of the kallisto index for human samples
#' @param kallisto_index_mouse A string, the location of the kallisto index for mouse samples
#' @param kallisto_dir A string, where to store the output of the tool. Defaults to _quants,
#' created via create_analysisfolders
#' @param kallisto_ncores A numeric value, with the number of cores to be used in the call
#' to kallisto
#' @param kallisto_nbootstraps Numeric, the number of bootstrap samples to generate
#' @param kallisto_gcbias Logical, perform sequence-specific bias correction
#' @param create_script Logical, whether to create the script - which needs to be run
#' from the terminal at a later moment
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
#' samplesinfo <- sra_to_fastq(samplesinfo)
#' samplesinfo <- match_fastq(samplesinfo)
#' samplesinfo <- run_fastqc(samplesinfo)
#' samplesinfo <- run_kallisto(samplesinfo)
run_kallisto <- function(samplesinfo, # contains the locations of each file/file pair
                         kallisto_bin = "kallisto",
                         kallisto_index_human,
                         kallisto_index_mouse,
                         kallisto_dir = "_quants",
                         kallisto_ncores = 12,
                         kallisto_nbootstraps = 100,
                         kallisto_gcbias = TRUE,
                         create_script = TRUE,
                         force = TRUE
) {
  # check that files are there and do exist
  if(is.null(samplesinfo$files_fastq)) {
    warning("No fastq files provided in the samplesinfo object!")
    return(samplesinfo)
  }
  if(!all(file.exists(unlist(samplesinfo$files_fastq)))) {
    warning("Not all fastq files of the samplesinfo object are actually existing!")
    return(samplesinfo)
  }
  
  # check existence of indices
  if(!file.exists(kallisto_index_human))
    stop("kallisto index for human not found!")
  if(!file.exists(kallisto_index_mouse))
    stop("kallisto index for mouse not found!")
  
  N <- length(samplesinfo$files_sra)
  kallisto_calls <- lapply(seq_len(N), function (i){
    this_libtype <- samplesinfo$runinfo$LibraryLayout[i]
    this_species <- samplesinfo$runinfo$ScientificName[i]
    this_index <- if(this_species == "Homo sapiens") {
      kallisto_index_human
    } else if(this_species == "Mus musculus") {
      kallisto_index_mouse
    } else {
      message("Found a dataset with a species currently not supported...")
      return(paste("#",samplesinfo$runinfo$Run[i],"Found a dataset with a species currently not supported..."))
    }
    if (this_libtype == "SINGLE") {
      this_fastqset <- samplesinfo$files_fastq[[i]]
      kallisto_call <- paste(kallisto_bin,"quant",
                             "--threads",kallisto_ncores,
                             "--bootstrap-samples",kallisto_nbootstraps,
                             ifelse(kallisto_gcbias, "--bias", ""),
                             "--index", this_index,
                             "--single -l 200 -s 30", # for single end this needs to be specified
                             # think of having evtl. pseudobam too?
                             "--output-dir",file.path(samplesinfo$data_dir,samplesinfo$datasetID,kallisto_dir,
                                                      paste0(names(samplesinfo$files_fastq)[i],"_kallisto")),
                             this_fastqset
      )
    } else if (this_libtype == "PAIRED"){
      this_fastqset <- samplesinfo$files_fastq[[i]]
      kallisto_call <- paste(kallisto_bin,"quant",
                             "--threads",kallisto_ncores,
                             "--bootstrap-samples",kallisto_nbootstraps,
                             ifelse(kallisto_gcbias, "--bias", ""),
                             "--index", this_index,
                             # think of having evtl. pseudobam too?
                             "--output-dir",file.path(samplesinfo$data_dir,samplesinfo$datasetID,kallisto_dir,
                                                      paste0(names(samplesinfo$files_fastq)[i],"_kallisto")),
                             this_fastqset$r1,this_fastqset$r2
      )
    }
    return(kallisto_call)
  })
  names(kallisto_calls) <- names(samplesinfo$files_fastq)
  
  out_bashscript <- file.path(samplesinfo$data_dir,
                              samplesinfo$datasetID,
                              paste0("cmd_batchkallistorun_",unique(samplesinfo$runinfo$SRAStudy),".sh"))
  # maybe prepend with shebang?
  cmd <- unlist(kallisto_calls)
  
  if(create_script) {
    if(file.exists(out_bashscript) & !force) {
      message("Bash script for running kallisto already available at ", out_bashscript)
    } else {
      writeLines(cmd, con = out_bashscript)
      message("Bash script for running kallisto generated at ", out_bashscript)
    }
  } 
  
  samplesinfo$cmd_kallisto <- cmd
  
  return(samplesinfo)
}



#' Run STAR
#' 
#' Creates a script for running STAR to align fastq raw files
#'
#' Providing two indices allows the handling of data sets where more than one species
#' is available, without the need to split the set into smaller subsets
#' 
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param star_bin A string, with the location of the STAR binary. This needs to be 
#' provided according to the local installation. Defaults to star, as installed via conda
#' @param star_index_human A string, the location of the STAR index for human samples
#' @param star_index_mouse A string, the location of the STAR index for mouse samples
#' @param star_gtffile_human A string, the location of the gtf annotation file to use in STAR
#' for human samples
#' @param star_gtffile_mouse A string, the location of the gtf annotation file to use in STAR
#' for mouse samples
#' @param star_dir A string, where to store the output of the tool. Defaults to _aligned,
#' created via create_analysisfolders
#' @param star_ncores A numeric value, with the number of cores to be used in the call
#' to STAR
#' @param create_script Logical, whether to create the script - which needs to be run
#' from the terminal at a later moment
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
#' samplesinfo <- sra_to_fastq(samplesinfo)
#' samplesinfo <- match_fastq(samplesinfo)
#' samplesinfo <- run_fastqc(samplesinfo)
#' samplesinfo <- run_STAR(samplesinfo)
run_STAR <- function(samplesinfo, # contains the locations of each file/file pair
                     star_bin = "STAR",
                     star_index_human,
                     star_index_mouse,
                     star_gtffile_human,
                     star_gtffile_mouse,
                     star_dir = "_aligned",
                     star_ncores = 12,
                     create_script = TRUE,
                     force = FALSE) {
  # check that files are there and do exist
  if(is.null(samplesinfo$files_fastq)) {
    warning("No fastq files provided in the samplesinfo object!")
    return(samplesinfo)
  }
  if(!all(file.exists(unlist(samplesinfo$files_fastq)))) {
    warning("Not all fastq files of the samplesinfo object are actually existing!")
    return(samplesinfo)
  }
  
  # check existence of indices
  if(!dir.exists(star_index_human))
    stop("STAR index for human not found!")
  if(!dir.exists(star_index_mouse))
    stop("STAR index for mouse not found!")
  
  
  N <- length(samplesinfo$files_sra)
  star_calls <- lapply(seq_len(N), function (i){
    this_libtype <- samplesinfo$runinfo$LibraryLayout[i]
    this_species <- samplesinfo$runinfo$ScientificName[i]
    this_index <- if(this_species == "Homo sapiens") {
      star_index_human
    } else if(this_species == "Mus musculus") {
      star_index_mouse
    } else {
      message("Found a dataset with a species currently not supported...")
      return(paste("#",samplesinfo$runinfo$Run[i],"Found a dataset with a species currently not supported..."))
    }
    
    this_gtffile <- if(this_species == "Homo sapiens") {
      star_gtffile_human
    } else if(this_species == "Mus musculus") {
      star_gtffile_mouse
    } else {
      # message("Found a dataset with a species currently not supported...")
      return(paste("#",samplesinfo$runinfo$Run[i],"Found a dataset with a species currently not supported..."))
    }
    # this_fastqset <- 
    if (this_libtype == "SINGLE") {
      this_fastqset <- samplesinfo$files_fastq[[i]]
      star_call <- paste0(star_bin,
                          " --runThreadN ",star_ncores,
                          " --genomeDir ",this_index,
                          " --sjdbGTFfile ",this_gtffile,
                          " --readFilesIn <(zcat ",this_fastqset,")",
                          " --outFileNamePrefix ",file.path(samplesinfo$data_dir,samplesinfo$datasetID,star_dir,
                                                            paste0(samplesinfo$runinfo$Run[i],"_")),
                          " --outSAMtype BAM SortedByCoordinate")
      
    } else if (this_libtype == "PAIRED"){
      this_fastqset <- samplesinfo$files_fastq[[i]]
      star_call <- paste0(star_bin,
                          " --runThreadN ",star_ncores,
                          " --genomeDir ",this_index,
                          " --sjdbGTFfile ",this_gtffile,
                          " --readFilesIn <(zcat ",this_fastqset$r1,")", " <(zcat ",this_fastqset$r2,")",
                          " --outFileNamePrefix ",file.path(samplesinfo$data_dir,samplesinfo$datasetID,star_dir,
                                                            paste0(samplesinfo$runinfo$Run[i],"_")),
                          " --outSAMtype BAM SortedByCoordinate")
    }
    return(star_call)
  })
  names(star_calls) <- names(samplesinfo$files_fastq)
  
  out_bashscript <- file.path(samplesinfo$data_dir,
                              samplesinfo$datasetID,
                              paste0("cmd_batchSTARrun_",unique(samplesinfo$runinfo$SRAStudy),".sh"))
  # maybe prepend with shebang?
  cmd <- unlist(star_calls)
  
  if(create_script) {
    if(file.exists(out_bashscript) & !force) {
      message("Bash script for running STAR already available at ", out_bashscript)
    } else {
      writeLines(cmd, con = out_bashscript)
      message("Bash script for running STAR generated at ", out_bashscript)
    }
  } 
  
  samplesinfo$cmd_STAR <- cmd
  
  return(samplesinfo)
}



# we can/could/should check the md5 hash of the downloaded sra files

# one way of doing this: using vdb-validate from sra-toolkit
# https://www.biostars.org/p/147148/

#' Validate the sra files
#' 
#' Validate sra files via vdb-validate from the sra toolkit
#' 
#' A working installation of the sra-toolkit has to be present.
#' 
#' Ideally, after the sra files are downloaded and checked, they could be deleted
#' after dumping the fastq files
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param nthreads The number of cores to use if parallel calls are to be used. The 
#' BiocParallel package is used to handle the multicore parameters
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
#' samplesinfo <- validate_sra(samplesinfo)
validate_sra <- function(samplesinfo,
                         nthreads = 4,
                         force = FALSE) {
  stopifnot(!is.null(samplesinfo$runinfo))
  stopifnot(!is.null(samplesinfo$files_sra))
  
  st <- Sys.time()
  timestamp()
  
  if(!is.null(samplesinfo$validated_sra) & !force){
    message("sra files have been already validated")
    if(all(samplesinfo$validated_sra == 0))
      message("All files validated, YAY!")
    else
      message("Following files were not ok: ",
              names(samplesinfo$validated_sra[samplesinfo$validated_sra != 0]))
    
    return(samplesinfo)
    
  } else {
    
    message("Checking that all sra files exist...")
    allthere <- all(file.exists(samplesinfo$files_sra))
    message("All there... ", allthere)
    if(!allthere)
      message("Missing datasets:", samplesinfo$files_sra[!file.exists(samplesinfo$files_sra)])
    
    N <- length(samplesinfo$files_sra)
    validate_calls <- lapply(seq_len(N), function (i) {
      mycall <- paste0("vdb-validate ",samplesinfo$files_sra[i])
      return(mycall)  
    })
    names(validate_calls) <- samplesinfo$runinfo$Run
    
    if(nthreads > 1){
      library(BiocParallel)
      register(MulticoreParam(workers = nthreads))
      myret <- unlist(bplapply(validate_calls,system))
    } else {
      myret <- unlist(lapply(validate_calls,system))
    }
    
    if(all(myret == 0))
      message("All files validated, YAY!")
    else
      message("Following files were not ok: ",
              names(myret[myret != 0]))
    
    # add some info to the samplesinfo that we did perform validation
    samplesinfo$validated_sra <- myret
    
    et <- Sys.time()
    elapsed_time <- et - st
    message("Elapsed time: ", round(elapsed_time,3), " ", attr(elapsed_time,"units"))
    
    return(samplesinfo)
  }
}




#' Check the fastq.gz files
#' 
#' Checks the integrity of the fastq.gz files, based on a system call to gunzip -t
#' 
#' Note: no parallelization is expected in this case, as it normally clutters the 
#' I/O rate of the storage system, thus slowing everything down.
#' 
#' This function might require some time to run, especially for sets with many
#' samples or big sized data. Just sit tight and grab a coffee (or two)
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps. 
#' @param force Logical. If the file to be created is already available, it can be 
#' overwritten by setting to TRUE. Defaults to FALSE
#'
#' @return
#' @export
#'
#' @examples
check_fastq <- function(samplesinfo,
                        force = FALSE) {
  stopifnot(!is.null(samplesinfo$runinfo))
  stopifnot(!is.null(samplesinfo$files_sra))
  
  # fastq files need to be there already
  stopifnot(!is.null(samplesinfo[["files_fastq"]]))
  
  nrsamples <- length(samplesinfo$files_fastq)
  nrsamples_runinfo <- nrow(samplesinfo$runinfo)
  # if there is a discrepancy, flag it?
  
  # somewhat inspired by the approach in https://www.biostars.org/p/147148/
  
  st <- Sys.time()
  timestamp()
  
  if(!is.null(samplesinfo$checked_fastq) & !force){
    message("fastq files have been already checked for integrity")
    if(all(samplesinfo$checked_fastq == 0))
      message("All files checked, YAY!")
    else
      message("Following files were not ok: ",
              names(samplesinfo$checked_fastq[samplesinfo$checked_fastq != 0]))
    return(samplesinfo)
  } else {
    # no need to replicate the single/paired structure of the data. yet...
    # samplesinfo$checked_fastq <- samplesinfo$files_fastq
    
    check_calls <- lapply(seq_len(nrsamples), function (i) {
      this_libtype <- samplesinfo$runinfo$LibraryLayout[i]
      if (this_libtype == "SINGLE") {
        this_fastqset <- samplesinfo$files_fastq[[i]]
        mycall <- paste0("gunzip -t ",this_fastqset)
        return(mycall)  
      } else if (this_libtype == "PAIRED"){
        this_fastqset <- samplesinfo$files_fastq[[i]]
        mycall_r1 <- paste0("gunzip -t ",this_fastqset$r1)
        mycall_r2 <- paste0("gunzip -t ",this_fastqset$r2)
        mycall <- c(mycall_r1,mycall_r2)
        return(mycall)  
      }
    })
    names(check_calls) <- samplesinfo$runinfo$Run
    
    # gracefully transform into vector and keep matched names
    unlisted_check_calls <- unlist(check_calls)
    
    allfastqs <- unlist(lapply(seq_len(nrsamples), function (i) {
      this_libtype <- samplesinfo$runinfo$LibraryLayout[i]
      if (this_libtype == "SINGLE")
        return(names(samplesinfo$files_fastq)[i])
      else if (this_libtype == "PAIRED")
        return(c(paste0(names(samplesinfo$files_fastq)[i],"_1"),
                 paste0(names(samplesinfo$files_fastq)[i],"_2"))
              )
    }))
    
    names(unlisted_check_calls) <- allfastqs
    
    # this one takes some time
    ## would be nice to have some kind of progress messages while running?
    ### something like R.utils::gunzip with remove=FALSE
    ### probably a for cycle with messages delivered is good enough
    
    myret <- rep(NA,length(unlisted_check_calls))
    names(myret) <- names(unlisted_check_calls)
    
    for(i in seq_len(length(unlisted_check_calls))) {
      message("Checking file ",i," of ",length(unlisted_check_calls),
              " - filename: ", names(unlisted_check_calls)[i],".fastq.gz...")
      myret[i] <- system(unlisted_check_calls[i])
    }
    # myret <- unlist(lapply(unlisted_check_calls,system))
    
    if(all(myret == 0))
      message("All files checked, YAY!")
    else
      message("Following files were not ok: ",
              names(myret[myret != 0]))
    
    # add some info to the samplesinfo that we did perform validation
    samplesinfo$checked_fastq <- myret
    
    et <- Sys.time()
    elapsed_time <- et - st
    message("Elapsed time: ", round(elapsed_time,3), " ", attr(elapsed_time,"units"))
    return(samplesinfo)
  }
}







