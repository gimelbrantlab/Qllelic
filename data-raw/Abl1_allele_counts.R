## code to prepare `Abl1_allele_counts` dataset goes here
Abl1_allele_counts = GetGatkPipelineTabs("/home/asya/Dropbox (Partners HealthCare)/data/Abl/Abl1_processed_gene.v3.1.txt", c(2), 1:19)

usethis::use_data(Abl1_allele_counts)
