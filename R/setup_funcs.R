#' setup_folders
#'
#' Creating folder structure for the main reosurce
#'
#' @param project_name
#'
#' @return
#' @export
#'
#' @examples
setup_folders <- function(project_name = "plateletopedia") {
  dir.create(paste0(project_name, "_ref"))         # for reference genomes, annotations, and so on
  dir.create(paste0(project_name, "_db"))          # for sql(ite) dbs, or tabular dbs retrieved
  dir.create(paste0(project_name, "_samplesinfo"))
  dir.create(paste0(project_name, "_publicdata"))
}

# ## -------------------------------------------------------------------------- ##
# Retrieving via SRAdb the whole dump of the SRA repository
# ## -------------------------------------------------------------------------- ##
# library(SRAdb)
# sqlfile <- '_db/SRAmetadb.sqlite'
# if(!file.exists('_db/SRAmetadb.sqlite'))
#   sqlfile <<- getSRAdbFile()
#
# # Then, create a connection for later queries. The standard DBI functionality as implemented
# # in RSQLite function dbConnect makes the connection to the database. The dbDisconnect
# # function disconnects the connection.
# sra_con <- dbConnect(SQLite(),sqlfile)
#
# # some examples on how this could be used:
# rs <- getSRAinfo ( c("SRX000122"), sra_con, sraType = "sra" )
# rs[1:3,]
#
# getFASTQinfo( c("SRR000648","SRR000657"), sra_con, srcType = 'ftp' )
# getSRAfile( c("SRR000648","SRR000657"), sra_con, fileType = 'fastq' )
#
#
# # # Additionally, to obtain the SRAid correspondences as in https://www.biostars.org/p/244150/
# # accessions_SRA_file <- "ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab"
# # download.file(accessions_SRA_file,paste0("_db/SRA_Accessions_SNAPSHOT_",Sys.Date(),".tab"))
# #
# # # accessions_SRA <- readr::read_tsv(accessions_SRA_file)
# #
# # grep '^SRR' _db/SRA_Accessions_SNAPSHOT_2018-01-31.tab | grep 'GSM' > _db/SRA_Accessions_RUNSonly_SNAPSHOT_2018-01-31.tab
# #
# # # insert header as first line
# # # Accession       Submission      Status  Updated Published       Received        Type    Center  Visibility      Alias   Experiment      Sample  Study   Loaded  Spots   Bases   Md5sum  BioSample     BioProject      ReplacedBy
# # sed -i "1i\
# # $(head -n 1 _db/SRA_Accessions_SNAPSHOT_2018-01-31.tab)" _db/SRA_Accessions_RUNSonly_SNAPSHOT_2018-01-31.tab
# #
# # mydf <- data.table::fread("_db/SRA_Accessions_RUNSonly_SNAPSHOT_2018-01-31.tab")
# # # colnames(mydf) <- dbListFields(sraacc_con,name = "SRAaccs")
# #
# #
# # # to be done once
# # library(DBI)
# # mySRAaccessiondb <- dbConnect(RSQLite::SQLite(), paste0("_db/SRA_Accessions_SNAPSHOT_",Sys.Date(),".sqlite"))
# # dbWriteTable(mySRAaccessiondb, "SRAaccs", mydf)
# #
# # # to use it:
# # sraacc_con <- dbConnect(SQLite(),paste0("_db/SRA_Accessions_SNAPSHOT_",Sys.Date(),".sqlite"))
# # dbListTables(sraacc_con)
# #
# # dbListFields(sraacc_con,name = "SRAaccs")
# #
# # dbGetQuery(sraacc_con, 'SELECT * FROM SRAaccs WHERE "Alias" == :x',
# #            params = list(x = "GSM2474455"))
#
#
# rs <- getSRA( search_terms = "platelets",out_types = c('run','study'), sra_con )
# # retrieve info programmatically
# unique(rs$study)   ## as a starter



# Retrieving for human and mouse all the relevant references, annotation, and so on

#' Title
#'
#' @param version_number
#' @param species
#' @param unzip_files
#' @param ref_dir
#' @param quiet_dl
#' @param just_check
#'
#' @return
#' @export
#'
#' @examples
retrieve_ensembl_refs <- function(version_number = 102,
                                  species = c("human","mouse"),
                                  # maybe parametrize species, and use different folder structure?
                                  unzip_files = FALSE,
                                  ref_dir = "plateletopedia_ref",
                                  quiet_dl = FALSE,
                                  just_check = TRUE) {

  cur_species <- match.arg(species)

  all_organisms <- setNames(c("Homo_sapiens","Mus_musculus"),c("human","mouse"))
  all_builds <- setNames(c("GRCh38","GRCm38"),c("human","mouse"))
  # all_prefixes <- setNames(c("Homo_sapiens.GRCh38","Mus_musculus.GRCm38"),c("human","mouse"))
  cur_organism <- all_organisms[cur_species]
  cur_build <- all_builds[cur_species]

  if(just_check) {
    if(cur_species == "human")
      browseURL("https://www.ensembl.org/Homo_sapiens/Info/Annotation")
    if(cur_species == "mouse")
      browseURL("https://www.ensembl.org/Mus_musculus/Info/Annotation")
    return(invisible())
  }

  message("Retrieving all required files from ENSEMBL version ", version_number)
  this_dir <- file.path(ref_dir,paste0("Ensembl.",cur_build,".",version_number))

  if(!dir.exists(this_dir)) {
    dir.create(this_dir)
    message(Sys.time(), " --- ", this_dir, " folder created...")
  } else {
    message(Sys.time(), " --- ", this_dir, " folder is already existing")
  }

  message("\nRetrieving references for ", paste(cur_organism, cur_build))

  base_ftp <- paste0("ftp://ftp.ensembl.org/pub/release-",version_number)
  url_genomesequence <- paste0(base_ftp,"/fasta/",tolower(cur_organism),"/dna/",cur_organism,".",cur_build,".dna.primary_assembly.fa.gz")
  url_cdnasequence <-   paste0(base_ftp,"/fasta/",tolower(cur_organism),"/cdna/",cur_organism,".",cur_build,".cdna.all.fa.gz")
  url_codingseq <-      paste0(base_ftp,"/fasta/",tolower(cur_organism),"/cds/",cur_organism,".",cur_build,".cds.all.fa.gz")
  url_ncrna <-          paste0(base_ftp,"/fasta/",tolower(cur_organism),"/ncrna/",cur_organism,".",cur_build,".ncrna.fa.gz")
  url_annogtf <-        paste0(base_ftp,"/gtf/",tolower(cur_organism),"/",cur_organism,".",cur_build,".",version_number,".gtf.gz")
  url_annogff3 <-       paste0(base_ftp,"/gff3/",tolower(cur_organism),"/",cur_organism,".",cur_build,".",version_number,".gff3.gz")

  dest_genomesequence <- file.path(this_dir,paste0(cur_organism,".",cur_build,".dna.primary_assembly.fa.gz"))
  dest_cdnasequence <-   file.path(this_dir,paste0(cur_organism,".",cur_build,".cdna.all.fa.gz"))
  dest_codingseq <-      file.path(this_dir,paste0(cur_organism,".",cur_build,".cds.all.fa.gz"))
  dest_ncrna <-          file.path(this_dir,paste0(cur_organism,".",cur_build,".ncrna.fa.gz"))
  dest_annogtf <-        file.path(this_dir,paste0(cur_organism,".",cur_build,".",version_number,".gtf.gz"))
  dest_annogff3 <-       file.path(this_dir,paste0(cur_organism,".",cur_build,".",version_number,".gff3.gz"))

  if(!file.exists(dest_genomesequence)){
    message("Fetching genome sequence in fasta format")
    download.file(url = url_genomesequence,destfile = dest_genomesequence, quiet = quiet_dl)
  } else {
    message("Already found genome sequence in fasta format")
  }

  if(!file.exists(dest_cdnasequence)){
    message("Fetching cdna sequence in fasta format")
    download.file(url = url_cdnasequence,destfile = dest_cdnasequence, quiet = quiet_dl)
  } else {
    message("Already found cdna sequence in fasta format")
  }

  if(!file.exists(dest_codingseq)){
    message("Fetching coding sequences in fasta format")
    download.file(url = url_codingseq,destfile = dest_codingseq, quiet = quiet_dl)
  } else {
    message("Already found coding sequences in fasta format")
  }

  if(!file.exists(dest_ncrna)){
    message("Fetching noncoding rna sequences in fasta format")
    download.file(url = url_ncrna,destfile = dest_ncrna, quiet = quiet_dl)
  } else {
    message("Already found noncoding rna sequences in fasta format")
  }

  if(!file.exists(dest_annogtf)){
    message("Fetching gene annotation in gtf format")
    download.file(url = url_annogtf,destfile = dest_annogtf, quiet = quiet_dl)
  } else {
    message("Already found gene annotation in gtf format")
  }

  if(!file.exists(dest_annogff3)){
    message("Fetching gene annotation in gff3 format")
    download.file(url = url_annogff3,destfile = dest_annogff3, quiet = quiet_dl)
  } else {
    message("Already found gene annotation in gff3 format")
  }

  if(unzip_files) {
    # unzipping all files so that they are directly accessible for the program to generate indices and co.
    message("\nUnzipping the downloaded resources... (this might take a while)")
    message("Unzipping references in ", this_dir)
    system(paste0("gunzip ", file.path(this_dir,"*gz")))
    message("Done unzipping")
  }

  message("\n",Sys.time(), " --- Done retrieving all files for ENSEMBL version ", version_number, " for ", paste0(cur_organism, " - build ", cur_build))
  return(invisible())
}

retrieve_ensembl_refs(version_number = 102, species = "mouse", just_check = FALSE)
retrieve_ensembl_refs(version_number = 102, species = "human", just_check = FALSE)


#' Title
#'
#' @param version_number
#' @param unzip_files
#' @param ref_dir
#' @param quiet_dl
#' @param just_check
#'
#' @return
#' @export
#'
#' @examples
retrieve_gencode_refs <- function(version_number = 30, # or "M21"
                                  unzip_files = FALSE,
                                  ref_dir = "plateletopedia_ref",
                                  quiet_dl = FALSE,
                                  just_check = TRUE){
  ## downloading the essential from the GENCODE project as well (https://www.gencodegenes.org/releases/current.html)
  if(just_check) {
    browseURL("https://www.gencodegenes.org")
    return(invisible())
  }
  message("Retrieving all required files from GENCODE version ", version_number)
  this_dir <- file.path(ref_dir,paste0("GENCODE_",version_number))

  if(!dir.exists(file.path(ref_dir,paste0("GENCODE_",version_number)))) {
    dir.create(this_dir)
    message(Sys.time(), " --- ", this_dir, " folder created...")
  } else {
    message(Sys.time(), " --- ", this_dir, " folder is already existing")
  }

  gencode_organism <- ifelse(grepl("^M",version_number),"Gencode_mouse","Gencode_human")
  message("\nRetrieving references for GENCODE version ",version_number)
  message("Organism: ", gencode_organism)


  if(grepl("^M",version_number)) {
    # in mouse
    vers <- gsub("^M", "", version_number)
    if(vers <= 25) {
      gencode_fastaref <- "GRCm38.primary_assembly.genome.fa.gz"
    } else {
      gencode_fastaref <- "GRCm39.primary_assembly.genome.fa.gz"
    }
  } else {
    # in human
    if(version_number <= 38) {
      gencode_fastaref <- "GRCh38.primary_assembly.genome.fa.gz"
    } else {
      # in case it is changed with 39
      gencode_fastaref <- "GRCh38.primary_assembly.genome.fa.gz"
    }
  }

  base_ftp <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/"
  url_genomesequence <- paste0(base_ftp ,gencode_organism,"/release_",version_number,"/",gencode_fastaref)
  url_cdnasequence <-   paste0(base_ftp, gencode_organism,"/release_",version_number,"/gencode.v",version_number,".transcripts.fa.gz")
  url_annogtf <-        paste0(base_ftp, gencode_organism,"/release_",version_number,"/gencode.v",version_number,".annotation.gtf.gz")
  url_annogff3 <-       paste0(base_ftp,gencode_organism, "/release_",version_number,"/gencode.v",version_number,".annotation.gff3.gz")

  dest_genomesequence <- file.path(this_dir,gencode_fastaref)
  dest_cdnasequence <-   file.path(this_dir,paste0("gencode.v",version_number,".transcripts.fa.gz"))
  dest_annogtf <-        file.path(this_dir,paste0("gencode.v",version_number,".annotation.gtf.gz"))
  dest_annogff3 <-       file.path(this_dir,paste0("gencode.v",version_number,".annotation.gff3.gz"))

  if(!file.exists(dest_genomesequence)){
    message("Fetching genome sequence in fasta format")
    download.file(url = url_genomesequence,destfile = dest_genomesequence, quiet = quiet_dl)
  } else {
    message("Already found genome sequence in fasta format")
  }

  if(!file.exists(dest_cdnasequence)){
    message("Fetching cdna sequence in fasta format")
    download.file(url = url_cdnasequence,destfile = dest_cdnasequence, quiet = quiet_dl)
  } else {
    message("Already found cdna sequence in fasta format")
  }

  if(!file.exists(dest_annogtf)){
    message("Fetching gene annotation in gtf format")
    download.file(url = url_annogtf,destfile = dest_annogtf, quiet = quiet_dl)
  } else {
    message("Already found gene annotation in gtf format")
  }

  if(!file.exists(dest_annogff3)){
    message("Fetching gene annotation in gff3 format")
    download.file(url = url_annogff3,destfile = dest_annogff3, quiet = quiet_dl)
  } else {
    message("Already found gene annotation in gff3 format")
  }

  if(unzip_files) {
    # unzipping all files so that they are directly accessible for the program to generate indices and co.
    message("\nUnzipping the downloaded resources... (this might take a while)")
    message("Unzipping references for ",gencode_organism)
    system(paste0("gunzip ", file.path(this_dir,"*gz")))
    message("Done unzipping")
  }

  message("\n",Sys.time(), " --- Done retrieving all files for GENCODE version ", version_number)

  return(invisible())
}


