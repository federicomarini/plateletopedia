




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
  sra_libtypes <- samplesinfo$library_types
  
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
  
  return(samplesinfo)
}


## TODO: will need a function for alternative entry point with own data
## ideally: fetch_sraruninfo - like
##          create_analysisfolders
##          match_fastq - like, to populate the samplesinfo$files_fastq slot



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





#' Run salmon
#' 
#' Creates a script for running salmon and generating transcript-level quantifications
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param salmon_bin A string, with the location of the salmon binary. This needs to be 
#' provided according to the local installation. Defaults to salmon, as installed via conda
#' @param salmon_index A string, the location of the salmon index
#' @param salmon_dir A string, where to store the output of the tool. Defaults to _quants,
#' created via create_analysisfolders
#' @param salmon_ncores A numeric value, with the number of cores to be used in the call
#' to salmon
#' @param salmon_nbootstraps Numeric, the number of bootstrap samples to generate
#' @param salmon_gcbias Logical, perform sequence-specific bias correction
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
#' samplesinfo <- run_salmon(samplesinfo)
run_salmon <- function(samplesinfo, # contains the locations of each file/file pair
                       salmon_bin = "salmon",
                       salmon_index,
                       salmon_dir = "_quants",
                       salmon_ncores = 12,
                       salmon_nbootstraps = 100,
                       salmon_gcbias = TRUE,
                       create_script = TRUE,
                       force = FALSE
) {
  N <- length(samplesinfo$files_sra)
  salmon_calls <- lapply(seq_len(N), function (i){
    this_libtype <- samplesinfo$runinfo$LibraryLayout[i]
    # this_fastqset <- 
    if (this_libtype == "SINGLE") {
      this_fastqset <- samplesinfo$files_fastq[[i]]
      salmon_call <- paste(salmon_bin,"quant",
                           "--threads",salmon_ncores,
                           "--numBootstraps",salmon_nbootstraps,
                           "--seqBias",salmon_gcbias,
                           "--libType A",
                           "--index", salmon_index,
                           "-r",this_fastqset,
                           "-o",file.path(samplesinfo$data_dir,samplesinfo$datasetID,salmon_dir,
                                          paste0(names(samplesinfo$files_fastq)[i],"_salmon")))
    } else if (this_libtype == "PAIRED"){
      salmon_call <- paste(salmon_bin,"quant",
                           "--threads",salmon_ncores,
                           "--numBootstraps",salmon_nbootstraps,
                           "--seqBias",salmon_gcbias,
                           "--libType A",
                           "--index", salmon_index,
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
      message("Bash script for running salmon generated at ", out_bashscript)
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
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param kallisto_bin A string, with the location of the kallisto binary. This needs to be 
#' provided according to the local installation. Defaults to kallisto, as installed via conda
#' @param kallisto_index A string, the location of the kallisto index
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
                         kallisto_index,
                         kallisto_dir = "_quants",
                         kallisto_ncores = 12,
                         kallisto_nbootstraps = 100,
                         kallisto_gcbias = TRUE,
                         create_script = TRUE,
                         force = TRUE
) {
  N <- length(samplesinfo$files_sra)
  kallisto_calls <- lapply(seq_len(N), function (i){
    this_libtype <- samplesinfo$runinfo$LibraryLayout[i]
    if (this_libtype == "SINGLE") {
      this_fastqset <- samplesinfo$files_fastq[[i]]
      kallisto_call <- paste(kallisto_bin,"quant",
                             "--threads",kallisto_ncores,
                             "--bootstrap-samples",kallisto_nbootstraps,
                             "--bias",kallisto_gcbias,
                             "--index", kallisto_index,
                             "--single -l 200 -s 30", # for single end this needs to be specified
                             # think of having evtl. pseudobam too?
                             "--output-dir",file.path(samplesinfo$data_dir,samplesinfo$datasetID,kallisto_dir,
                                                      paste0(names(samplesinfo$files_fastq)[i],"_kallisto")),
                             this_fastqset
      )
    } else if (this_libtype == "PAIRED"){
      kallisto_call <- paste(kallisto_bin,"quant",
                             "--threads",kallisto_ncores,
                             "--bootstrap-samples",kallisto_nbootstraps,
                             "--bias",kallisto_gcbias,
                             "--index", kallisto_index,
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
      message("Bash script for running kallisto generated at ", out_bashscript)
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
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param star_bin A string, with the location of the STAR binary. This needs to be 
#' provided according to the local installation. Defaults to star, as installed via conda
#' @param star_index A string, the location of the STAR index
#' @param star_gtffile A string, the location of the gtf annotation file to use in STAR
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
                     star_bin = "star",
                     star_index,
                     star_gtffile,
                     star_dir = "_aligned",
                     star_ncores = 12,
                     create_script = TRUE,
                     force = FALSE) {
  N <- length(samplesinfo$files_sra)
  star_calls <- lapply(seq_len(N), function (i){
    this_libtype <- samplesinfo$runinfo$LibraryLayout[i]
    # this_fastqset <- 
    if (this_libtype == "SINGLE") {
      this_fastqset <- samplesinfo$files_fastq[[i]]
      star_call <- paste0(star_bin,
                          " --runThreadN ",star_ncores,
                          " --genomeDir ",star_index,
                          " --sjdbGTFfile ",star_gtffile,
                          " --readFilesIn <(zcat ",this_fastqset,")",
                          " --outFileNamePrefix ",file.path(samplesinfo$data_dir,samplesinfo$datasetID,star_dir,
                                                            paste0(samplesinfo$runinfo$Run[i],"_")),
                          " --outSAMtype BAM SortedByCoordinate")
      
    } else if (this_libtype == "PAIRED"){
      this_fastqset <- samplesinfo$files_fastq[[i]]
      star_call <- paste0(star_bin,
                          " --runThreadN ",star_ncores,
                          " --genomeDir ",star_index,
                          " --sjdbGTFfile ",star_gtffile,
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
      message("Bash script for running STAR generated at ", out_bashscript)
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

#' Title
#'
#' @param samplesinfo A samplesinfo-like list to store all the required info and steps.
#' @param nthreads 
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
    message(allthere)
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
    return(samplesinfo)
  }
}



