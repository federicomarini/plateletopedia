
# general setup and so
library(magrittr)
source("01_helper_functions.R")






## -------------------------------------------------------------------------- ##
##  alhasan_jackson-circ_degradation-2016  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_alhasan_jackson_2016 <- fetch_sraruninfo("SRP058654",datasetID = "alhasan_jackson-circ_degradation-2016") %>% 
  fetch_geoinfo("GSE69192") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_alhasan_jackson_2016 <- validate_sra(samplesinfo_alhasan_jackson_2016)
samplesinfo_alhasan_jackson_2016 <- check_fastq(samplesinfo_alhasan_jackson_2016)
save(samplesinfo_alhasan_jackson_2016, file = "_samplesinfo/samplesinfo_alhasan_jackson_2016.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_alhasan_jackson_2016$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_alhasan_jackson_2016$files_sra)
}






## -------------------------------------------------------------------------- ##
##  an_gallagher-erythroid_diff-2014  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_an_gallagher_2014 <- fetch_sraruninfo("SRP035312",datasetID = "an_gallagher-erythroid_diff-2014") %>% 
  fetch_geoinfo("GSE53983") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_an_gallagher_2014 <- validate_sra(samplesinfo_an_gallagher_2014)
samplesinfo_an_gallagher_2014 <- check_fastq(samplesinfo_an_gallagher_2014)
save(samplesinfo_an_gallagher_2014, file = "_samplesinfo/samplesinfo_an_gallagher_2014.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_an_gallagher_2014$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_an_gallagher_2014$files_sra)
}






## -------------------------------------------------------------------------- ##
##  beauchemin_moroy-megs_gfi1b-2017  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_beauchemin_moroy_2017 <- fetch_sraruninfo("SRP061548",datasetID = "beauchemin_moroy-megs_gfi1b-2017") %>% 
  fetch_geoinfo("GSE71310") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_beauchemin_moroy_2017 <- validate_sra(samplesinfo_beauchemin_moroy_2017)
samplesinfo_beauchemin_moroy_2017 <- check_fastq(samplesinfo_beauchemin_moroy_2017)
save(samplesinfo_beauchemin_moroy_2017, file = "_samplesinfo/samplesinfo_beauchemin_moroy_2017.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_beauchemin_moroy_2017$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_beauchemin_moroy_2017$files_sra)
}






## -------------------------------------------------------------------------- ##
##  best_wurdinger-TEPs-2015  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_best_wurdinger_2015 <- fetch_sraruninfo("SRP057500",datasetID = "best_wurdinger-TEPs-2015") %>% 
  fetch_geoinfo("GSE68086") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_best_wurdinger_2015 <- validate_sra(samplesinfo_best_wurdinger_2015)
samplesinfo_best_wurdinger_2015 <- check_fastq(samplesinfo_best_wurdinger_2015)
save(samplesinfo_best_wurdinger_2015, file = "_samplesinfo/samplesinfo_best_wurdinger_2015.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_best_wurdinger_2015$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_best_wurdinger_2015$files_sra)
}






## -------------------------------------------------------------------------- ##
##  best_wurdinger-TEPs_swarm-2017  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_best_wurdinger_2017 <- fetch_sraruninfo("SRP093349",datasetID = "best_wurdinger-TEPs_swarm-2017") %>% 
  fetch_geoinfo("GSE85864") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_best_wurdinger_2017 <- validate_sra(samplesinfo_best_wurdinger_2017)
samplesinfo_best_wurdinger_2017 <- check_fastq(samplesinfo_best_wurdinger_2017)
save(samplesinfo_best_wurdinger_2017, file = "_samplesinfo/samplesinfo_best_wurdinger_2017.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_best_wurdinger_2017$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_best_wurdinger_2017$files_sra)
}









samplesinfo_best_wurdinger_2017 <- fetch_sraruninfo("SRP093349",datasetID = "best_wurdinger-TEPs_swarm-2017") %>% 
  fetch_geoinfo("GSE85864") %>% 
  create_analysisfolders() %>% 
  get_sradata() %>% 
  sra_to_fastq() %>% 
  match_fastq() 






## -------------------------------------------------------------------------- ##
##  boisset_vanoudenaarden-network_bonemarrow-2016  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_boisset_vanoudenaarden_2016 <- fetch_sraruninfo("SRP092389",datasetID = "boisset_vanoudenaarden-network_bonemarrow-2016") %>% 
  fetch_geoinfo("GSE89378") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_boisset_vanoudenaarden_2016 <- validate_sra(samplesinfo_boisset_vanoudenaarden_2016)
samplesinfo_boisset_vanoudenaarden_2016 <- check_fastq(samplesinfo_boisset_vanoudenaarden_2016)
save(samplesinfo_boisset_vanoudenaarden_2016, file = "_samplesinfo/samplesinfo_boisset_vanoudenaarden_2016.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_boisset_vanoudenaarden_2016$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_boisset_vanoudenaarden_2016$files_sra)
}






## -------------------------------------------------------------------------- ##
##  bray_rigoutsos-complex_landscape-2013  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_bray_rigoutsos_2013 <- fetch_sraruninfo("SRA062032",datasetID = "bray_rigoutsos-complex_landscape-2013") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_bray_rigoutsos_2013 <- validate_sra(samplesinfo_bray_rigoutsos_2013)
samplesinfo_bray_rigoutsos_2013 <- check_fastq(samplesinfo_bray_rigoutsos_2013)
save(samplesinfo_bray_rigoutsos_2013, file = "_samplesinfo/samplesinfo_bray_rigoutsos_2013.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_bray_rigoutsos_2013$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_bray_rigoutsos_2013$files_sra)
}








## -------------------------------------------------------------------------- ##
##  campbell_rondina-granzymea_platelets-2017  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_campbell_rondina_2017 <- fetch_sraruninfo("SRP114983",datasetID = "campbell_rondina-granzymea_platelets-2017") %>%
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_campbell_rondina_2017 <- validate_sra(samplesinfo_campbell_rondina_2017)
samplesinfo_campbell_rondina_2017 <- check_fastq(samplesinfo_campbell_rondina_2017)
save(samplesinfo_campbell_rondina_2017, file = "_samplesinfo/samplesinfo_campbell_rondina_2017.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_campbell_rondina_2017$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_campbell_rondina_2017$files_sra)
}






## -------------------------------------------------------------------------- ##
##  cimmino_golino-mirna_modulation-2015  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_cimmino_golino_2015 <- fetch_sraruninfo("ERP004316",datasetID = "cimmino_golino-mirna_modulation-2015") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_cimmino_golino_2015 <- validate_sra(samplesinfo_cimmino_golino_2015)
samplesinfo_cimmino_golino_2015 <- check_fastq(samplesinfo_cimmino_golino_2015)
save(samplesinfo_cimmino_golino_2015, file = "_samplesinfo/samplesinfo_cimmino_golino_2015.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_cimmino_golino_2015$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_cimmino_golino_2015$files_sra)
}







## -------------------------------------------------------------------------- ##
##  corces_chang-lineagespecific_hematopoiesis-2016  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_corces_chang_2016 <- fetch_sraruninfo("SRP065216",datasetID = "corces_chang-lineagespecific_hematopoiesis-2016") %>% 
  fetch_geoinfo("GSE74246") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_corces_chang_2016 <- validate_sra(samplesinfo_corces_chang_2016)
samplesinfo_corces_chang_2016 <- check_fastq(samplesinfo_corces_chang_2016)
save(samplesinfo_corces_chang_2016, file = "_samplesinfo/samplesinfo_corces_chang_2016.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_corces_chang_2016$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_corces_chang_2016$files_sra)
}






## -------------------------------------------------------------------------- ##
##  delbridge_grabow-puma_ttp-2016  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_delbridge_grabow_2016 <- fetch_sraruninfo("SRP067232",datasetID = "delbridge_grabow-puma_ttp-2016") %>% 
  fetch_geoinfo("GSE75896") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_delbridge_grabow_2016 <- validate_sra(samplesinfo_delbridge_grabow_2016)
samplesinfo_delbridge_grabow_2016 <- check_fastq(samplesinfo_delbridge_grabow_2016)
save(samplesinfo_delbridge_grabow_2016, file = "_samplesinfo/samplesinfo_delbridge_grabow_2016.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_delbridge_grabow_2016$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_delbridge_grabow_2016$files_sra)
}






## -------------------------------------------------------------------------- ##
##  duff_graveley-rnaseq_20humantissues-2015  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_duff_graveley_2015 <- fetch_sraruninfo("SRP056969",datasetID = "duff_graveley-rnaseq_20humantissues-2015") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_duff_graveley_2015 <- validate_sra(samplesinfo_duff_graveley_2015)
samplesinfo_duff_graveley_2015 <- check_fastq(samplesinfo_duff_graveley_2015)
save(samplesinfo_duff_graveley_2015, file = "_samplesinfo/samplesinfo_duff_graveley_2015.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_duff_graveley_2015$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_duff_graveley_2015$files_sra)
}







## -------------------------------------------------------------------------- ##
##  eicher_johnson-acute_mi-2016  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_eicher_johnson_2016 <- fetch_sraruninfo("SRP053296",datasetID = "eicher_johnson-acute_mi-2016") %>% 
  fetch_geoinfo("GSE65705") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_eicher_johnson_2016 <- validate_sra(samplesinfo_eicher_johnson_2016)
samplesinfo_eicher_johnson_2016 <- check_fastq(samplesinfo_eicher_johnson_2016)
save(samplesinfo_eicher_johnson_2016, file = "_samplesinfo/samplesinfo_eicher_johnson_2016.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_eicher_johnson_2016$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_eicher_johnson_2016$files_sra)
}






## -------------------------------------------------------------------------- ##
##  fagerberg_uhlen-HPA_normaltissues-2016  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_fagerberg_uhlen_2016 <- fetch_sraruninfo("ERP003613",datasetID = "fagerberg_uhlen-HPA_normaltissues-2016") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_fagerberg_uhlen_2016 <- validate_sra(samplesinfo_fagerberg_uhlen_2016)
samplesinfo_fagerberg_uhlen_2016 <- check_fastq(samplesinfo_fagerberg_uhlen_2016)
save(samplesinfo_fagerberg_uhlen_2016, file = "_samplesinfo/samplesinfo_fagerberg_uhlen_2016.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_fagerberg_uhlen_2016$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_fagerberg_uhlen_2016$files_sra)
}







## -------------------------------------------------------------------------- ##
##  grover_nerlov-singlecell_hsc-2016  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_grover_nerlov_2016 <- fetch_sraruninfo("SRP060557",datasetID = "grover_nerlov-singlecell_hsc-2016") %>% 
  fetch_geoinfo("GSE70657") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_grover_nerlov_2016 <- validate_sra(samplesinfo_grover_nerlov_2016)
samplesinfo_grover_nerlov_2016 <- check_fastq(samplesinfo_grover_nerlov_2016)
save(samplesinfo_grover_nerlov_2016, file = "_samplesinfo/samplesinfo_grover_nerlov_2016.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_grover_nerlov_2016$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_grover_nerlov_2016$files_sra)
}






## -------------------------------------------------------------------------- ##
##  illumina_bodymap2-2013  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_illumina_bodymap2_2013 <- fetch_sraruninfo("ERP000546",datasetID = "illumina_bodymap2-2013") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_illumina_bodymap2_2013 <- validate_sra(samplesinfo_illumina_bodymap2_2013)
samplesinfo_illumina_bodymap2_2013 <- check_fastq(samplesinfo_illumina_bodymap2_2013)
save(samplesinfo_illumina_bodymap2_2013, file = "_samplesinfo/samplesinfo_illumina_bodymap2_2013.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_illumina_bodymap2_2013$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_illumina_bodymap2_2013$files_sra)
}







## -------------------------------------------------------------------------- ##
##  kissopoulou_osman-polyA-2013  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_kissopoulou_osman_2013 <- fetch_sraruninfo("ERP000803",datasetID = "kissopoulou_osman-polyA-2013") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_kissopoulou_osman_2013 <- validate_sra(samplesinfo_kissopoulou_osman_2013)
samplesinfo_kissopoulou_osman_2013 <- check_fastq(samplesinfo_kissopoulou_osman_2013)
save(samplesinfo_kissopoulou_osman_2013, file = "_samplesinfo/samplesinfo_kissopoulou_osman_2013.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_kissopoulou_osman_2013$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_kissopoulou_osman_2013$files_sra)
}







## -------------------------------------------------------------------------- ##
##  kissopoulou_osman-ribo_depl-2013  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_kissopoulou_osman2_2013 <- fetch_sraruninfo("ERP003815",datasetID = "kissopoulou_osman-ribo_depl-2013") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_kissopoulou_osman2_2013 <- validate_sra(samplesinfo_kissopoulou_osman2_2013)
samplesinfo_kissopoulou_osman2_2013 <- check_fastq(samplesinfo_kissopoulou_osman2_2013)
save(samplesinfo_kissopoulou_osman2_2013, file = "_samplesinfo/samplesinfo_kissopoulou_osman2_2013.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_kissopoulou_osman2_2013$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_kissopoulou_osman2_2013$files_sra)
}







## -------------------------------------------------------------------------- ##
##  lefrancais_looney-lungs_bm-2017  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_lefrancais_looney_2017 <- fetch_sraruninfo("SRP097794",datasetID = "lefrancais_looney-lungs_bm-2017") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_lefrancais_looney_2017 <- validate_sra(samplesinfo_lefrancais_looney_2017)
samplesinfo_lefrancais_looney_2017 <- check_fastq(samplesinfo_lefrancais_looney_2017)
save(samplesinfo_lefrancais_looney_2017, file = "_samplesinfo/samplesinfo_lefrancais_looney_2017.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_lefrancais_looney_2017$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_lefrancais_looney_2017$files_sra)
}







## -------------------------------------------------------------------------- ##
##  londin_rigoutsos-txome_proteome-2014  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_londin_rigoutsos_2014 <- fetch_sraruninfo("SRP028846",datasetID = "londin_rigoutsos-txome_proteome-2014") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_londin_rigoutsos_2014 <- validate_sra(samplesinfo_londin_rigoutsos_2014)
samplesinfo_londin_rigoutsos_2014 <- check_fastq(samplesinfo_londin_rigoutsos_2014)
save(samplesinfo_londin_rigoutsos_2014, file = "_samplesinfo/samplesinfo_londin_rigoutsos_2014.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_londin_rigoutsos_2014$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_londin_rigoutsos_2014$files_sra)
}







## -------------------------------------------------------------------------- ##
##  maass_rajewsky-human_circ-2017  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_maass_rajewsky_2017 <- fetch_sraruninfo("SRP109805",datasetID = "maass_rajewsky-human_circ-2017") %>% 
  fetch_geoinfo("GSE100242") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_maass_rajewsky_2017 <- validate_sra(samplesinfo_maass_rajewsky_2017)
samplesinfo_maass_rajewsky_2017 <- check_fastq(samplesinfo_maass_rajewsky_2017)
save(samplesinfo_maass_rajewsky_2017, file = "_samplesinfo/samplesinfo_maass_rajewsky_2017.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_maass_rajewsky_2017$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_maass_rajewsky_2017$files_sra)
}






## -------------------------------------------------------------------------- ##
##  meinders_philipsen-sp1sp3_hallmarks-2015  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_meinders_philipsen_2015 <- fetch_sraruninfo("SRP043469",datasetID = "meinders_philipsen-sp1sp3_hallmarks-2015") %>% 
  fetch_geoinfo("GSE58707") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_meinders_philipsen_2015 <- validate_sra(samplesinfo_meinders_philipsen_2015)
samplesinfo_meinders_philipsen_2015 <- check_fastq(samplesinfo_meinders_philipsen_2015)
save(samplesinfo_meinders_philipsen_2015, file = "_samplesinfo/samplesinfo_meinders_philipsen_2015.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_meinders_philipsen_2015$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_meinders_philipsen_2015$files_sra)
}






## -------------------------------------------------------------------------- ##
##  mills_ingolia-pelo_decay-2017  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_mills_ingolia_2017 <- fetch_sraruninfo("SRP098699",datasetID = "mills_ingolia-pelo_decay-2017") %>% 
  fetch_geoinfo("GSE94384") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_mills_ingolia_2017 <- validate_sra(samplesinfo_mills_ingolia_2017)
samplesinfo_mills_ingolia_2017 <- check_fastq(samplesinfo_mills_ingolia_2017)
save(samplesinfo_mills_ingolia_2017, file = "_samplesinfo/samplesinfo_mills_ingolia_2017.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_mills_ingolia_2017$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_mills_ingolia_2017$files_sra)
}






## -------------------------------------------------------------------------- ##
##  mills_ingolia-riboprofiling-2016  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_mills_ingolia_2016 <- fetch_sraruninfo("SRP082436",datasetID = "mills_ingolia-riboprofiling-2016") %>% 
  fetch_geoinfo("GSE85864") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_mills_ingolia_2016 <- validate_sra(samplesinfo_mills_ingolia_2016)
samplesinfo_mills_ingolia_2016 <- check_fastq(samplesinfo_mills_ingolia_2016)
save(samplesinfo_mills_ingolia_2016, file = "_samplesinfo/samplesinfo_mills_ingolia_2016.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_mills_ingolia_2016$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_mills_ingolia_2016$files_sra)
}






## -------------------------------------------------------------------------- ##
##  nassa_tarallo-splicing_proteome-2018  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_nassa_tarallo_2018 <- fetch_sraruninfo("ERP104860",datasetID = "nassa_tarallo-splicing_proteome-2018") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_nassa_tarallo_2018 <- validate_sra(samplesinfo_nassa_tarallo_2018)
samplesinfo_nassa_tarallo_2018 <- check_fastq(samplesinfo_nassa_tarallo_2018)
save(samplesinfo_nassa_tarallo_2018, file = "_samplesinfo/samplesinfo_nassa_tarallo_2018.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_nassa_tarallo_2018$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_nassa_tarallo_2018$files_sra)
}







## -------------------------------------------------------------------------- ##
##  nurnberg_ouwehand-invitro_megs-2012  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_nurnberg_ouwehand_2012 <- fetch_sraruninfo("ERP001115",datasetID = "nurnberg_ouwehand-invitro_megs-2012") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_nurnberg_ouwehand_2012 <- validate_sra(samplesinfo_nurnberg_ouwehand_2012)
samplesinfo_nurnberg_ouwehand_2012 <- check_fastq(samplesinfo_nurnberg_ouwehand_2012)
save(samplesinfo_nurnberg_ouwehand_2012, file = "_samplesinfo/samplesinfo_nurnberg_ouwehand_2012.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_nurnberg_ouwehand_2012$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_nurnberg_ouwehand_2012$files_sra)
}







## -------------------------------------------------------------------------- ##
##  osman_provost-pathogen_reduction-2015  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_osman_provost_2015 <- fetch_sraruninfo("ERP009260",datasetID = "osman_provost-pathogen_reduction-2015") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_osman_provost_2015 <- validate_sra(samplesinfo_osman_provost_2015)
samplesinfo_osman_provost_2015 <- check_fastq(samplesinfo_osman_provost_2015)
save(samplesinfo_osman_provost_2015, file = "_samplesinfo/samplesinfo_osman_provost_2015.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_osman_provost_2015$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_osman_provost_2015$files_sra)
}







## -------------------------------------------------------------------------- ##
##  pontes_burbano-mirna_celldamage-2015  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_pontes_burbano_2015 <- fetch_sraruninfo("SRP048290",datasetID = "pontes_burbano-mirna_celldamage-2015") %>% 
  fetch_geoinfo("GSE61856") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_pontes_burbano_2015 <- validate_sra(samplesinfo_pontes_burbano_2015)
samplesinfo_pontes_burbano_2015 <- check_fastq(samplesinfo_pontes_burbano_2015)
save(samplesinfo_pontes_burbano_2015, file = "_samplesinfo/samplesinfo_pontes_burbano_2015.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_pontes_burbano_2015$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_pontes_burbano_2015$files_sra)
}






## -------------------------------------------------------------------------- ##
##  preusser_bindereif-release_circ-2018  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_preusser_bindereif_2018 <- fetch_sraruninfo("SRP118609",datasetID = "preusser_bindereif-release_circ-2018") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_preusser_bindereif_2018 <- validate_sra(samplesinfo_preusser_bindereif_2018)
samplesinfo_preusser_bindereif_2018 <- check_fastq(samplesinfo_preusser_bindereif_2018)
save(samplesinfo_preusser_bindereif_2018, file = "_samplesinfo/samplesinfo_preusser_bindereif_2018.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_preusser_bindereif_2018$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_preusser_bindereif_2018$files_sra)
}







## -------------------------------------------------------------------------- ##
##  ramirez_mortazavi-dynamic_myeloid-2017  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_ramirez_mortazavi_2017 <- fetch_sraruninfo("SRP071547",datasetID = "ramirez_mortazavi-dynamic_myeloid-2017") %>% 
  fetch_geoinfo("GSE79044") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_ramirez_mortazavi_2017 <- validate_sra(samplesinfo_ramirez_mortazavi_2017)
samplesinfo_ramirez_mortazavi_2017 <- check_fastq(samplesinfo_ramirez_mortazavi_2017)
save(samplesinfo_ramirez_mortazavi_2017, file = "_samplesinfo/samplesinfo_ramirez_mortazavi_2017.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_ramirez_mortazavi_2017$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_ramirez_mortazavi_2017$files_sra)
}






## -------------------------------------------------------------------------- ##
##  rowley_weyrich-human_mouse-2011  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_rowley_weyrich_2011 <- fetch_sraruninfo("SRP119431",datasetID = "rowley_weyrich-human_mouse-2011") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_rowley_weyrich_2011 <- validate_sra(samplesinfo_rowley_weyrich_2011)
samplesinfo_rowley_weyrich_2011 <- check_fastq(samplesinfo_rowley_weyrich_2011)
save(samplesinfo_rowley_weyrich_2011, file = "_samplesinfo/samplesinfo_rowley_weyrich_2011.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_rowley_weyrich_2011$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_rowley_weyrich_2011$files_sra)
}







## -------------------------------------------------------------------------- ##
##  sawai_reizis-hsc_multilineage-2016  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_sawai_reizis_2016 <- fetch_sraruninfo("SRP071090",datasetID = "sawai_reizis-hsc_multilineage-2016") %>% 
  fetch_geoinfo("GSE78855") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_sawai_reizis_2016 <- validate_sra(samplesinfo_sawai_reizis_2016)
samplesinfo_sawai_reizis_2016 <- check_fastq(samplesinfo_sawai_reizis_2016)
save(samplesinfo_sawai_reizis_2016, file = "_samplesinfo/samplesinfo_sawai_reizis_2016.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_sawai_reizis_2016$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_sawai_reizis_2016$files_sra)
}






## -------------------------------------------------------------------------- ##
##  shi_weyrich-proteasome_platelets-2014  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_shi_weyrich_2014 <- fetch_sraruninfo("SRP062023",datasetID = "shi_weyrich-proteasome_platelets-2014") %>% 
  fetch_geoinfo("GSE58202") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_shi_weyrich_2014 <- validate_sra(samplesinfo_shi_weyrich_2014)
samplesinfo_shi_weyrich_2014 <- check_fastq(samplesinfo_shi_weyrich_2014)
save(samplesinfo_shi_weyrich_2014, file = "_samplesinfo/samplesinfo_shi_weyrich_2014.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_shi_weyrich_2014$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_shi_weyrich_2014$files_sra)
}






## -------------------------------------------------------------------------- ##
##  soellner_simon-mouserat_atlas-2017  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_soellner_simon_2017 <- fetch_sraruninfo("ERP104395",datasetID = "soellner_simon-mouserat_atlas-2017") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_soellner_simon_2017 <- validate_sra(samplesinfo_soellner_simon_2017)
samplesinfo_soellner_simon_2017 <- check_fastq(samplesinfo_soellner_simon_2017)
save(samplesinfo_soellner_simon_2017, file = "_samplesinfo/samplesinfo_soellner_simon_2017.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_soellner_simon_2017$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_soellner_simon_2017$files_sra)
}







## -------------------------------------------------------------------------- ##
##  szabo_salzman-human_fetaldevel-2015  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_szabo_salzman_2015 <- fetch_sraruninfo("SRP051249",datasetID = "szabo_salzman-human_fetaldevel-2015") %>% 
  fetch_geoinfo("GSE64283") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_szabo_salzman_2015 <- validate_sra(samplesinfo_szabo_salzman_2015)
samplesinfo_szabo_salzman_2015 <- check_fastq(samplesinfo_szabo_salzman_2015)
save(samplesinfo_szabo_salzman_2015, file = "_samplesinfo/samplesinfo_szabo_salzman_2015.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_szabo_salzman_2015$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_szabo_salzman_2015$files_sra)
}






## -------------------------------------------------------------------------- ##
##  unigiessen-plt_activation-unp  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_unigiessen_unp <- fetch_sraruninfo("SRP118609",datasetID = "unigiessen-plt_activation-unp") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_unigiessen_unp <- validate_sra(samplesinfo_unigiessen_unp)
samplesinfo_unigiessen_unp <- check_fastq(samplesinfo_unigiessen_unp)
save(samplesinfo_unigiessen_unp, file = "_samplesinfo/samplesinfo_unigiessen_unp.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_unigiessen_unp$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_unigiessen_unp$files_sra)
}







## -------------------------------------------------------------------------- ##
##  UNK_jefferson-human_plts-2013  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_UNK_jefferson_2013 <- fetch_sraruninfo("SRP034558",datasetID = "UNK_jefferson-human_plts-2013") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_UNK_jefferson_2013 <- validate_sra(samplesinfo_UNK_jefferson_2013)
samplesinfo_UNK_jefferson_2013 <- check_fastq(samplesinfo_UNK_jefferson_2013)
save(samplesinfo_UNK_jefferson_2013, file = "_samplesinfo/samplesinfo_UNK_jefferson_2013.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_UNK_jefferson_2013$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_UNK_jefferson_2013$files_sra)
}







## -------------------------------------------------------------------------- ##
##  yu_maheswaran-circulating_tumor-2013  ----------------------------------
## -------------------------------------------------------------------------- ##

samplesinfo_yu_maheswaran_2013 <- fetch_sraruninfo("SRP015945",datasetID = "yu_maheswaran-circulating_tumor-2013") %>% 
  fetch_geoinfo("GSE41245") %>% 
  create_analysisfolders() %>% 
  get_sradata(force = TRUE) %>% 
  sra_to_fastq(force = TRUE) %>% 
  match_fastq() 
samplesinfo_yu_maheswaran_2013 <- validate_sra(samplesinfo_yu_maheswaran_2013)
samplesinfo_yu_maheswaran_2013 <- check_fastq(samplesinfo_yu_maheswaran_2013)
save(samplesinfo_yu_maheswaran_2013, file = "_samplesinfo/samplesinfo_yu_maheswaran_2013.RData")

# if sra are all validated, they can be deleted
if(all(samplesinfo_yu_maheswaran_2013$validated_sra == 0)) {
  # delete from the R environment the sra files
  unlink(samplesinfo_yu_maheswaran_2013$files_sra)
}






























