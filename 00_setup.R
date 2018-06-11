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

dir.create("_ref/ENSEMBL_release91")
dir.create("_ref/ENSEMBL_release92")

# release 91 of ENSEMBL
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",destfile = "_ref/ENSEMBL_release91/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz",destfile = "_ref/ENSEMBL_release91/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",destfile = "_ref/ENSEMBL_release91/Homo_sapiens.GRCh38.cdna.all.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz",destfile = "_ref/ENSEMBL_release91/Mus_musculus.GRCm38.cdna.all.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz",destfile = "_ref/ENSEMBL_release91/Homo_sapiens.GRCh38.cds.all.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",destfile = "_ref/ENSEMBL_release91/Mus_musculus.GRCm38.cds.all.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz",destfile = "_ref/ENSEMBL_release91/Homo_sapiens.GRCh38.ncrna.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz",destfile = "_ref/ENSEMBL_release91/Mus_musculus.GRCm38.ncrna.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz",destfile = "_ref/ENSEMBL_release91/Homo_sapiens.GRCh38.91.gtf.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/gtf/mus_musculus/Mus_musculus.GRCm38.91.gtf.gz",destfile = "_ref/ENSEMBL_release91/Mus_musculus.GRCm38.91.gtf.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/gff3/homo_sapiens/Homo_sapiens.GRCh38.91.gff3.gz",destfile = "_ref/ENSEMBL_release91/Homo_sapiens.GRCh38.91.gff3.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-91/gff3/mus_musculus/Mus_musculus.GRCm38.91.gff3.gz",destfile = "_ref/ENSEMBL_release91/Mus_musculus.GRCm38.91.gff3.gz")

# release 92 of ENSEMBL - april 2018
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",destfile = "_ref/ENSEMBL_release92/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz",destfile = "_ref/ENSEMBL_release92/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",destfile = "_ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.all.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz",destfile = "_ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.all.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz",destfile = "_ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cds.all.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz",destfile = "_ref/ENSEMBL_release92/Mus_musculus.GRCm38.cds.all.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz",destfile = "_ref/ENSEMBL_release92/Homo_sapiens.GRCh38.ncrna.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz",destfile = "_ref/ENSEMBL_release92/Mus_musculus.GRCm38.ncrna.fa.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz",destfile = "_ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.gtf.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/Mus_musculus.GRCm38.92.gtf.gz",destfile = "_ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.gtf.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/gff3/homo_sapiens/Homo_sapiens.GRCh38.92.gff3.gz",destfile = "_ref/ENSEMBL_release92/Homo_sapiens.GRCh38.92.gff3.gz")
download.file(url = "ftp://ftp.ensembl.org/pub/release-92/gff3/mus_musculus/Mus_musculus.GRCm38.92.gff3.gz",destfile = "_ref/ENSEMBL_release92/Mus_musculus.GRCm38.92.gff3.gz")


# unzipping all files so that they are directly accessible for the program to generate indices and co.
system("gunzip _ref/ENSEMBL_release91/*")
system("gunzip _ref/ENSEMBL_release92/*")



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

# salmon indices

salmon index -t _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.ncrna.fa -i _ref/index_salmon_hs_GRCh38.92_cdna.ncrna.sidx --type quasi -k 31
salmon index -t _ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.ncrna.fa -i _ref/index_salmon_mm_GRCm38.92_cdna.ncrna.sidx --type quasi -k 31

# also with cdna only...
salmon index -t _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.all.fa -i _ref/index_salmon_hs_GRCh38.92_cdna.sidx --type quasi -k 31
salmon index -t _ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.all.fa -i _ref/index_salmon_mm_GRCm38.92_cdna.sidx --type quasi -k 31


# kallisto indices

kallisto index --index=_ref/index_kallisto_hs_GRCh38.92_cdna.ncrna.kidx -k 31 _ref/ENSEMBL_release92/Homo_sapiens.GRCh38.cdna.ncrna.fa
kallisto index --index=_ref/index_kallisto_mm_GRCm38.92_cdna.ncrna.kidx -k 31 _ref/ENSEMBL_release92/Mus_musculus.GRCm38.cdna.ncrna.fa


## -------------------------------------------------------------------------- ##
##  Converting the annotation files to formats used by other tools ----------
## -------------------------------------------------------------------------- ##



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











