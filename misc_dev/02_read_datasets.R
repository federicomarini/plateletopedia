# access the googledocs sheet
# install.packages("googlesheets")
library("googlesheets")
gs_auth(new_user = TRUE)

library(dplyr)
gs_ls() %>% View

plt_sets <- gs_read(gs_title("Platelets_datasets_collection"))
# open up the different Pubmed IDs, to retrieve the doi
## not required after the first time :)

# once the doi is there, use rcrossref to get the full textual citation
# devtools::install_github("ropensci/rcrossref")
library('rcrossref')
View(plt_sets)
## see here some documentation: https://github.com/ropensci/rcrossref
# cr_cn(dois = plt_sets$doi[31], format = "text", style = "apa")
cr_cn(dois = plt_sets$doi[36], format = "text", style = "apa")
View(get_styles()) # and there search for something with 'no-et-al' to list all authors
# cr_cn(dois = plt_sets$DOI[36], format = "text", style = "american-medical-association-no-et-al")
all_crossref_infos <- cr_cn(dois = plt_sets$doi, format = "text", style = "american-medical-association-no-et-al")

# store this info in a multi-line text document - even better if the exported file is in markdown, it is directly previewable
all_crossref_infos %>% unlist() %>% writeLines()

# writeLines(paste0("pubstatus_",Sys.Date(),"_exported.md"))

# this one gives some bibtex
cat(cr_cn(dois = plt_sets$doi[36], format = "bibtex"))

cr_abstract(plt_sets$doi[3])
# does not work...

