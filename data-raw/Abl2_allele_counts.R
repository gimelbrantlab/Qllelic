## code to prepare `Abl2_allele_counts` dataset goes here
Abl2_allele_counts = GetGatkPipelineTabs("/home/asya/Dropbox (Partners HealthCare)/data/Abl/Abl2_processed_gene.v3.1.txt", c(2), 1:19)

usethis::use_data(Abl2_allele_counts)
