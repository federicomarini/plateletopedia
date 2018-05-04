# Creating folder structure for the main reosurce

dir.create("_ref")        # for reference genomes, annotations, and so on
dir.create("_db")         # for sql(ite) dbs, or tabular dbs retrieved

# Retrieving via SRAdb the whole dump of the SRA repository
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



# retrieve the whole refs and annotations



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



## installing the binaries of the tools that are required

# nice solution: via bioconda!
conda install -c bioconda ucsc-gtftogenepred  ucsc-genepredtobed
conda install -c bioconda ucsc-bedgraphtobigwig ucsc-bedtobigbed ucsc-bedgraphtobigwig ucsc-bigwiginfo ucsc-bigwigtobedgraph ucsc-bedtobigbed ucsc-bigbedtobed ucsc-bigwigtowig ucsc-bigwigsummary
conda install sra-tools

conda install parallel-fastqdump

# Do something like an environment.yaml file to be read and fed directly to conda?

## STAR


mkdir _ref/index_STAR_hs_GRCh38.91
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir _ref/index_STAR_hs_GRCh38.91 \
--genomeFastaFiles _ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile _ref/Homo_sapiens.GRCh38.91.gtf --sjdbOverhang 100 # as default value, works for range of cases

mkdir _ref/index_STAR_mm_GRCm38.91
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir _ref/index_STAR_mm_GRCm38.91 \
--genomeFastaFiles _ref/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--sjdbGTFfile _ref/Mus_musculus.GRCm38.91.gtf --sjdbOverhang 100 # as default value, works for range of cases

# salmon indices
cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.ncrna.fa > Homo_sapiens.GRCh38.cdna.ncrna.fa
cat Mus_musculus.GRCm38.cdna.all.fa Mus_musculus.GRCm38.ncrna.fa > Mus_musculus.GRCm38.cdna.ncrna.fa

salmon index -t Homo_sapiens.GRCh38.cdna.ncrna.fa -i index_salmon_hs_GRCh38.91_cdna.ncrna.sidx --type quasi -k 31
salmon index -t Mus_musculus.GRCm38.cdna.ncrna.fa -i index_salmon_mm_GRCm38.91_cdna.ncrna.sidx --type quasi -k 31

# kallisto indices

kallisto index --index=index_kallisto_hs_GRCh38.91_cdna.ncrna.kidx -k 31 Homo_sapiens.GRCh38.cdna.ncrna.fa
kallisto index --index=index_kallisto_mm_GRCm38.91_cdna.ncrna.kidx -k 31 Mus_musculus.GRCm38.cdna.ncrna.fa



# gtf to genePred formats, as well as to refFlat, and also to bed
gtfToGenePred -genePredExt -geneNameAsName2 _ref/Homo_sapiens.GRCh38.91.gtf refFlat.tmp.txt
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > _ref/Homo_sapiens.GRCh38.91.refFlat.txt

gtfToGenePred -genePredExt -geneNameAsName2 _ref/Mus_musculus.GRCm38.91.gtf refFlat.tmp.txt
paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 refFlat.tmp.txt) > _ref/Mus_musculus.GRCm38.91.refFlat.txt

gtfToGenePred _ref/Homo_sapiens.GRCh38.91.gtf _ref/Homo_sapiens.GRCh38.91.genePred
gtfToGenePred _ref/Mus_musculus.GRCm38.91.gtf _ref/Mus_musculus.GRCm38.91.genePred

genePredToBed _ref/Homo_sapiens.GRCh38.91.genePred _ref/Homo_sapiens.GRCh38.91.bed
genePredToBed _ref/Mus_musculus.GRCm38.91.genePred _ref/Mus_musculus.GRCm38.91.bed



## Path to the kallisto binary
kallistobin <- "/usr/local/software/kallisto_linux-v0.44.0/kallisto"
kallistoversion <- "0.44"

## Human cDNA & ncRNA + ERCC
cmd <- sprintf("bash -c '%s index -i %s %s'",
               kallistobin,
               gsub("cdna.all.fa$", paste0(human_ensembl_version,
                                           ".cdna.ncrna.ercc92.", kallistoversion,
                                           ".kidx"), human_cdna_fa),
               paste0("<(cat ", human_cdna_fa, " ", human_ncrna_fa, " ", ercc_fa, ")"))
cmd
system(cmd)
