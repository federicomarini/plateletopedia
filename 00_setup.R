## -------------------------------------------------------------------------- ##
##  Creating folder structure for the main reosurce -------------------------
## -------------------------------------------------------------------------- ##

dir.create("_ref")        # for reference genomes, annotations, and so on
dir.create("_db")         # for sql(ite) dbs, or tabular dbs retrieved

## -------------------------------------------------------------------------- ##
##  Retrieving via SRAdb the whole dump of the SRA repository ---------------
## -------------------------------------------------------------------------- ##
library(SRAdb)
sqlfile <- '_db/SRAmetadb.sqlite'
if(!file.exists('_db/SRAmetadb.sqlite')) 
  sqlfile <<- getSRAdbFile()

# Then, create a connection for later queries. The standard DBI functionality as implemented
# in RSQLite function dbConnect makes the connection to the database. The dbDisconnect 
# function disconnects the connection.
sra_con <- dbConnect(SQLite(),sqlfile)

# some examples on how this could be used:
rs <- getSRAinfo ( c("SRX000122"), sra_con, sraType = "sra" )
rs[1:3,]

getFASTQinfo( c("SRR000648","SRR000657"), sra_con, srcType = 'ftp' ) 
getSRAfile( c("SRR000648","SRR000657"), sra_con, fileType = 'fastq' )


# # Additionally, to obtain the SRAid correspondences as in https://www.biostars.org/p/244150/
# accessions_SRA_file <- "ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab"
# download.file(accessions_SRA_file,paste0("_db/SRA_Accessions_SNAPSHOT_",Sys.Date(),".tab"))
# 
# # accessions_SRA <- readr::read_tsv(accessions_SRA_file)
# 
# grep '^SRR' _db/SRA_Accessions_SNAPSHOT_2018-01-31.tab | grep 'GSM' > _db/SRA_Accessions_RUNSonly_SNAPSHOT_2018-01-31.tab
# 
# # insert header as first line
# # Accession       Submission      Status  Updated Published       Received        Type    Center  Visibility      Alias   Experiment      Sample  Study   Loaded  Spots   Bases   Md5sum  BioSample     BioProject      ReplacedBy
# sed -i "1i\
# $(head -n 1 _db/SRA_Accessions_SNAPSHOT_2018-01-31.tab)" _db/SRA_Accessions_RUNSonly_SNAPSHOT_2018-01-31.tab
# 
# mydf <- data.table::fread("_db/SRA_Accessions_RUNSonly_SNAPSHOT_2018-01-31.tab")
# # colnames(mydf) <- dbListFields(sraacc_con,name = "SRAaccs")
# 
# 
# # to be done once
# library(DBI)
# mySRAaccessiondb <- dbConnect(RSQLite::SQLite(), paste0("_db/SRA_Accessions_SNAPSHOT_",Sys.Date(),".sqlite"))
# dbWriteTable(mySRAaccessiondb, "SRAaccs", mydf)
# 
# # to use it:
# sraacc_con <- dbConnect(SQLite(),paste0("_db/SRA_Accessions_SNAPSHOT_",Sys.Date(),".sqlite"))
# dbListTables(sraacc_con)
# 
# dbListFields(sraacc_con,name = "SRAaccs")
# 
# dbGetQuery(sraacc_con, 'SELECT * FROM SRAaccs WHERE "Alias" == :x', 
#            params = list(x = "GSM2474455"))


rs <- getSRA( search_terms = "platelets",out_types = c('run','study'), sra_con )
# retrieve info programmatically
unique(rs$study)   ## as a starter


## -------------------------------------------------------------------------- ##
##  Retrieve the whole refs and annotations ---------------------------------
## -------------------------------------------------------------------------- ##



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
retrieve_ensembl_refs <- function(version_number = 96, 
                                  species = c("human","mouse"),
                                  # maybe parametrize species, and use different folder structure?
                                  unzip_files = FALSE, 
                                  ref_dir = "_ref", 
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

retrieve_ensembl_refs(version_number = 96, species = "mouse", just_check = FALSE)
retrieve_ensembl_refs(version_number = 96, species = "human", just_check = FALSE)


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
                                  ref_dir = "_ref", 
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
  
  gencode_fastaref <- ifelse(
    grepl("^M",version_number),"GRCm38.p6.genome.fa.gz","GRCh38.p12.genome.fa.gz")
  
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

# retrieve_gencode_refs(version_number = "30",just_check = FALSE)
# retrieve_gencode_refs(version_number = "M21",just_check = FALSE)


## -------------------------------------------------------------------------- ##
##  Installing the binaries of the tools that are required ------------------
## -------------------------------------------------------------------------- ##

# nice solution: via bioconda!
conda install -c bioconda ucsc-gtftogenepred  ucsc-genepredtobed
conda install -c bioconda ucsc-bedgraphtobigwig ucsc-bedtobigbed ucsc-bedgraphtobigwig ucsc-bigwiginfo ucsc-bigwigtobedgraph ucsc-bedtobigbed ucsc-bigbedtobed ucsc-bigwigtowig ucsc-bigwigsummary
conda install sra-tools

conda install parallel-fastqdump

conda install suppa

# Do something like an environment.yaml file to be read and fed directly to conda?


## -------------------------------------------------------------------------- ##
##  Creating the relevant indices -------------------------------------------
## -------------------------------------------------------------------------- ##


# one fasta file for both coding and noncoding genes
cat _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.all.fa \
    _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.ncrna.fa > \
    _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.ncrna.fa
cat _ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.all.fa \
    _ref/ENSEMBL_release92/Mus_musculus.GRCm38.ncrna.fa > \
    _ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.ncrna.fa



## STAR indices

mkdir _ref/index_STAR_hs_GRCh38.92
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir _ref/index_STAR_hs_GRCh38.92 \
--genomeFastaFiles _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.gtf --sjdbOverhang 100 # as default value, works for range of cases

mkdir _ref/index_STAR_mm_GRCm38.92
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir _ref/index_STAR_mm_GRCm38.92 \
--genomeFastaFiles _ref/ENSEMBL_release92/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--sjdbGTFfile _ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.gtf --sjdbOverhang 100 # as default value, works for range of cases

# salmon indices - this was done with version 0.9.1

salmon index -t _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.ncrna.fa -i _ref/index_salmon_hs_GRCh38.92_cdna.ncrna.sidx --type quasi -k 31
salmon index -t _ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.ncrna.fa -i _ref/index_salmon_mm_GRCm38.92_cdna.ncrna.sidx --type quasi -k 31

# also with cdna only...
salmon index -t _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.all.fa -i _ref/index_salmon_hs_GRCh38.92_cdna.sidx --type quasi -k 31
salmon index -t _ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.all.fa -i _ref/index_salmon_mm_GRCm38.92_cdna.sidx --type quasi -k 31

# and why not, also with gencode transcripts
salmon index -t _ref/GENCODE_human_release28/gencode.v28.transcripts.fa -i _ref/index_salmon_hs_gencode_28.sidx --type quasi -k 31
salmon index -t _ref/GENCODE_mouse_release17/gencode.vM17.transcripts.fa -i _ref/index_salmon_mm_gencode_17.sidx --type quasi -k 31

## salmon indices with recent version - according to the algorithm change...

salmon index -t _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.ncrna.fa -i _ref/index_salmon_v0.11.2_hs_GRCh38.92_cdna.ncrna.sidx --type quasi -k 31
salmon index -t _ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.ncrna.fa -i _ref/index_salmon_v0.11.2_mm_GRCm38.92_cdna.ncrna.sidx --type quasi -k 31

# also with cdna only...
salmon index -t _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.all.fa -i _ref/index_salmon_v0.11.2_hs_GRCh38.92_cdna.sidx --type quasi -k 31
salmon index -t _ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.all.fa -i _ref/index_salmon_v0.11.2_mm_GRCm38.92_cdna.sidx --type quasi -k 31

# and why not, also with gencode transcripts
salmon index -t _ref/GENCODE_human_release28/gencode.v28.transcripts.fa -i _ref/index_salmon_v0.11.2_hs_gencode_28.sidx --type quasi -k 31
salmon index -t _ref/GENCODE_mouse_release18/gencode.vM18.transcripts.fa -i _ref/index_salmon_v0.11.2_mm_gencode_18.sidx --type quasi -k 31




# kallisto indices

kallisto index --index=_ref/index_kallisto_hs_GRCh38.92_cdna.ncrna.kidx -k 31 _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.ncrna.fa
kallisto index --index=_ref/index_kallisto_mm_GRCm38.92_cdna.ncrna.kidx -k 31 _ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.ncrna.fa


## -------------------------------------------------------------------------- ##
##  Converting the reference and annotation files to formats used by other tools ----------
## -------------------------------------------------------------------------- ##

# converting fasta to 2bit format
## require installation of kent tools faToTwoBit
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x faToTwoBit
./faToTwoBit _ref _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.dna.primary_assembly.fa _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.dna.primary_assembly.2bit
./faToTwoBit _ref _ref/ENSEMBL_release92/Mus_musculus.GRCm38.dna.primary_assembly.fa _ref/ENSEMBL_release92/Mus_musculus.GRCm38.dna.primary_assembly.2bit


# gtf to genePred formats, as well as to refFlat, and also to bed
gtfToGenePred -genePredExt -geneNameAsName2 _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.gtf refFlat.tmp.txt
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.refFlat.txt
rm refFlat.tmp.txt

gtfToGenePred -genePredExt -geneNameAsName2 _ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.gtf refFlat.tmp.txt
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > _ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.refFlat.txt
rm refFlat.tmp.txt


gtfToGenePred _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.gtf _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.genePred
gtfToGenePred _ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.gtf _ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.genePred

genePredToBed _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.genePred _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.bed
genePredToBed _ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.genePred _ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.bed


## for SUPPA - generate events files
# generate local AS events
mkdir _ref/index_suppa_hs_GRCh38.92
suppa.py generateEvents -i _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.gtf \
-o suppa_hs_GRCh38.92_events -f ioe -e SE SS MX RI FL
#Put all the ioe events in the same file:
awk '
FNR==1 && NR!=1 { while (/^<header>/) getline; }
1 {print}
' suppa_hs_GRCh38.92_events*.ioe > suppa_hs_GRCh38.92_events.ioe
# moving all the files to _ref location (cannot generate them directly, raises error?)
mv suppa_hs_GRCh38.92_events* _ref/index_suppa_hs_GRCh38.92/
  
mkdir _ref/index_suppa_mm_GRCm38.92
suppa.py generateEvents -i _ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.gtf \
-o suppa_mm_GRCm38.92_events -f ioe -e SE SS MX RI FL
awk '
FNR==1 && NR!=1 { while (/^<header>/) getline; }
1 {print}
' suppa_mm_GRCm38.92_events*.ioe > suppa_mm_GRCm38.92_events.ioe
mv suppa_mm_GRCm38.92_events* _ref/index_suppa_mm_GRCm38.92/

# generate the transcript "events"  
suppa.py generateEvents -i _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.gtf \
-o _ref/index_suppa_hs_GRCh38.92/suppa_hs_GRCh38.92_txevents -f ioi 
suppa.py generateEvents -i _ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.gtf \
-o _ref/index_suppa_mm_GRCm38.92/suppa_mm_GRCm38.92_txevents -f ioi 

## -------------------------------------------------------------------------- ##
##  Retrieving the annotations from AnnotationHub and making them pkgs ------
## -------------------------------------------------------------------------- ##


library(AnnotationHub)
## Load the annotation resource.
ah <- AnnotationHub()

## Query for all available EnsDb databases
query(ah, c("92", "EnsDb","Homo"))
query(ah, c("92", "EnsDb","Mus"))

edb_v92_human <- ah[["AH60977"]]
edb_v92_mouse <- ah[["AH60992"]]

makeEnsembldbPackage(ensdb = dbfile(dbconn(edb_v92_human)),
                     version = "0.0.0.9000", 
                     maintainer = "Federico Marini <marinif@uni-mainz.de>", 
                     author = "F Marini",
                     destDir = "_ref")

makeEnsembldbPackage(ensdb = dbfile(dbconn(edb_v92_mouse)),
                     version = "0.0.0.9000", 
                     maintainer = "Federico Marini <marinif@uni-mainz.de>",
                     author = "F Marini",
                     destDir = "_ref")

# these should be then built/installed as usual if one wants to avoid using the (still very handy) AnnotationHub






## -------------------------------------------------------------------------- ##
##  Generating tx2gene mappings and GRanges objects                    ------
## -------------------------------------------------------------------------- ##

# First, the tx2gene objects

# ercc <- readDNAStringSet(ercc_fa)
# ercc <- data.frame(tx = names(ercc), gene = names(ercc), stringsAsFactors = FALSE)

# using the excellent code from https://github.com/markrobinsonuzh/conquer/blob/master/01_build_index.R
library("Biostrings")

## Human
cdna <- readDNAStringSet("_ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.all.fa")
ncrna <- readDNAStringSet("_ref/ENSEMBL_release92/Homo_sapiens.GRCh38.ncrna.fa")
cdna_df <- data.frame(t(sapply(as.character(names(cdna)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  gene_symbol <- gsub("^gene_symbol:", "", a[grep("^gene_symbol:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), 
    gene = ifelse(length(gene) != 0, gene, NA), 
    gene_symbol = ifelse(length(gene_symbol) != 0, gene_symbol, NA))
})), stringsAsFactors = FALSE)
ncrna_df <- data.frame(t(sapply(as.character(names(ncrna)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  gene_symbol <- gsub("^gene_symbol:", "", a[grep("^gene_symbol:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), 
    gene = ifelse(length(gene) != 0, gene, NA), 
    gene_symbol = ifelse(length(gene_symbol) != 0, gene_symbol, NA))
})), stringsAsFactors = FALSE)
txgenemap <- rbind(cdna_df, ncrna_df) # , ercc)
colnames(txgenemap) <- c("TXNAME","GENEID","GENENAME")
saveRDS(txgenemap, 
        file = "_ref/tx2gene_Homo_sapiens.GRCh38.92.cdna.ncrna.txgenemap.rds")

## Mouse
cdna <- readDNAStringSet("_ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.all.fa")
ncrna <- readDNAStringSet("_ref/ENSEMBL_release92/Mus_musculus.GRCm38.ncrna.fa")
cdna_df <- data.frame(t(sapply(as.character(names(cdna)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  gene_symbol <- gsub("^gene_symbol:", "", a[grep("^gene_symbol:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), 
    gene = ifelse(length(gene) != 0, gene, NA), 
    gene_symbol = ifelse(length(gene_symbol) != 0, gene_symbol, NA))
})), stringsAsFactors = FALSE)
ncrna_df <- data.frame(t(sapply(as.character(names(ncrna)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  gene_symbol <- gsub("^gene_symbol:", "", a[grep("^gene_symbol:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), 
    gene = ifelse(length(gene) != 0, gene, NA), 
    gene_symbol = ifelse(length(gene_symbol) != 0, gene_symbol, NA))
})), stringsAsFactors = FALSE)
txgenemap <- rbind(cdna_df, ncrna_df) # , ercc)
colnames(txgenemap) <- c("TXNAME","GENEID","GENENAME")
saveRDS(txgenemap, 
        file = "_ref/tx2gene_Mus_musculus.GRCm38.92.cdna.ncrna.txgenemap.rds")


# Probably even more elegant, creating TxDb objects

library(GenomicFeatures)

gtf <- "_ref/GENCODE_human_release28/gencode.v28.annotation.gtf"
txdb.filename <- "_ref/TxDb_gencode_human_v28_annotation.sqlite"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)

gtf <- "_ref/GENCODE_mouse_release17/gencode.vM17.annotation.gtf"
txdb.filename <- "_ref/TxDb_gencode_mouse_v17_annotation.sqlite"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)

gtf <- "_ref/GENCODE_mouse_release18/gencode.vM18.annotation.gtf"
txdb.filename <- "_ref/TxDb_gencode_mouse_v18_annotation.sqlite"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)

gtf <- "_ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.gtf"
txdb.filename <- "_ref/TxDb_ensembl_human_v92_annotation.sqlite"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)

gtf <- "_ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.gtf"
txdb.filename <- "_ref/TxDb_ensembl_mouse_v92_annotation.sqlite"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)





# Now, the GRanges objects

library(GenomicRanges)
library(dplyr)

## Human
cdna <- readDNAStringSet("_ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.all.fa")
ncrna <- readDNAStringSet("_ref/ENSEMBL_release92/Homo_sapiens.GRCh38.ncrna.fa")

nm <- c(names(cdna), names(ncrna))
info <- data.frame(t(sapply(nm, function(w) {
  w <- strsplit(w, " ")[[1]]
  transcript <- w[1]
  gene <- gsub("gene:", "", w[grep("^gene:", w)])
  position <- gsub("chromosome:", "", w[grep("^chromosome:", w)])
  if (length(position) == 0) position <- gsub("scaffold:", "", w[grep("^scaffold:", w)])
  symbol <- gsub("gene_symbol:", "", w[grep("^gene_symbol", w)])
  c(transcript = transcript, 
    gene = gene,
    genome = strsplit(position, ":")[[1]][1], 
    chromosome = strsplit(position, ":")[[1]][2],
    start = strsplit(position, ":")[[1]][3],
    end = strsplit(position, ":")[[1]][4],
    strand = ifelse(strsplit(position, ":")[[1]][5] == "1", "+", "-"),
    symbol = symbol
  )
})), stringsAsFactors = FALSE)
rownames(info) <- NULL
info$start <- as.numeric(info$start)
info$end <- as.numeric(info$end)

txgr <- GRanges(seqnames = info$chromosome, 
                ranges = IRanges(start = info$start, end = info$end), 
                strand = info$strand)
mcols(txgr) <- info[, c("transcript", "gene", "genome", "symbol")]

# erccgtf <- import(ercc_gtf, format = "gtf")
# ercctx <- erccgtf
# mcols(ercctx) <- data.frame(transcript = ercctx$gene_id,
#                             gene = ercctx$gene_id,
#                             genome = ercctx$source,
#                             symbol = ercctx$gene_id)
# 
# txgr <- suppressWarnings(c(txgr, ercctx))
names(txgr) <- txgr$transcript

geneinfo <- info %>% group_by(gene) %>% 
  summarise(genome = unique(genome),
            chromosome = unique(chromosome),
            start = min(start),
            end = max(end),
            strand = unique(strand),
            symbol = unique(symbol)) %>% 
  ungroup()

ggr <- GRanges(seqnames = geneinfo$chromosome,
               ranges = IRanges(start = geneinfo$start, end = geneinfo$end),
               strand = geneinfo$strand)
mcols(ggr) <- geneinfo[, c("gene", "genome", "symbol")]

# erccgene <- erccgtf
# mcols(erccgene) <- data.frame(gene = erccgene$gene_id,
#                               genome = erccgene$source,
#                               symbol = erccgene$gene_id)
# ggr <- suppressWarnings(c(ggr, erccgene))
# names(ggr) <- ggr$gene

gene_granges <- ggr
tx_granges <- txgr
saveRDS(list(gene_granges = gene_granges, tx_granges = tx_granges),
        file = "_ref/granges_Homo_sapiens.GRCh38.92.cdna.ncrna.granges.gene_tx.rds")

## Mouse
cdna <- readDNAStringSet("_ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.all.fa")
ncrna <- readDNAStringSet("_ref/ENSEMBL_release92/Mus_musculus.GRCm38.ncrna.fa")

nm <- c(names(cdna), names(ncrna))
info <- data.frame(t(sapply(nm, function(w) {
  w <- strsplit(w, " ")[[1]]
  transcript <- w[1]
  gene <- gsub("gene:", "", w[grep("^gene:", w)])
  position <- gsub("chromosome:", "", w[grep("^chromosome:", w)])
  if (length(position) == 0) position <- gsub("scaffold:", "", w[grep("^scaffold:", w)])
  symbol <- gsub("gene_symbol:", "", w[grep("^gene_symbol", w)])
  c(transcript = transcript, 
    gene = gene,
    genome = strsplit(position, ":")[[1]][1], 
    chromosome = strsplit(position, ":")[[1]][2],
    start = strsplit(position, ":")[[1]][3],
    end = strsplit(position, ":")[[1]][4],
    strand = ifelse(strsplit(position, ":")[[1]][5] == "1", "+", "-"),
    symbol = symbol
  )
})), stringsAsFactors = FALSE)
rownames(info) <- NULL
info$start <- as.numeric(info$start)
info$end <- as.numeric(info$end)

txgr <- GRanges(seqnames = info$chromosome, 
                ranges = IRanges(start = info$start, end = info$end), 
                strand = info$strand)
mcols(txgr) <- info[, c("transcript", "gene", "genome", "symbol")]

# erccgtf <- import(ercc_gtf, format = "gtf")
# ercctx <- erccgtf
# mcols(ercctx) <- data.frame(transcript = ercctx$gene_id,
#                             gene = ercctx$gene_id,
#                             genome = ercctx$source,
#                             symbol = ercctx$gene_id)
# 
# txgr <- suppressWarnings(c(txgr, ercctx))
names(txgr) <- txgr$transcript

geneinfo <- info %>% group_by(gene) %>% 
  summarise(genome = unique(genome),
            chromosome = unique(chromosome),
            start = min(start),
            end = max(end),
            strand = unique(strand),
            symbol = unique(symbol)) %>% 
  ungroup()

ggr <- GRanges(seqnames = geneinfo$chromosome,
               ranges = IRanges(start = geneinfo$start, end = geneinfo$end),
               strand = geneinfo$strand)
mcols(ggr) <- geneinfo[, c("gene", "genome", "symbol")]

# erccgene <- erccgtf
# mcols(erccgene) <- data.frame(gene = erccgene$gene_id,
#                               genome = erccgene$source,
#                               symbol = erccgene$gene_id)
# ggr <- suppressWarnings(c(ggr, erccgene))
# names(ggr) <- ggr$gene

gene_granges <- ggr
tx_granges <- txgr
saveRDS(list(gene_granges = gene_granges, tx_granges = tx_granges),
        file = "_ref/granges_Mus_musculus.GRCm38.92.cdna.ncrna.granges.gene_tx.rds")





