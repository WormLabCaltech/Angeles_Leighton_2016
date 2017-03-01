source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
biocLite("biomaRt")

setwd('~/WormFiles/agingRNAseq/rdocs')
#gene info for sleuth
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "celegans_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

#install.packages("devtools")

#say no to binary files
#devtools::install_github("pachterlab/sleuth")

library("sleuth")

#point to your directory+
base_dir <- "~/WormFiles/agingRNAseq/rdocs/sleuth"

#get ids
sample_id <- dir(file.path(base_dir, "results"))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id, "kallisto"))
kal_dirs


s2c <- read.table(file.path(base_dir, "aging_rnaseq_info.txt"), header = TRUE, stringsAsFactors= FALSE)
s2c <- dplyr::select(s2c, sample = experiment, genotype, age)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)


#prepend and make object, state maximum model here
so <- sleuth_prep(s2c, ~ genotype+age+genotype*age, target_mapping= t2g)

#fit the models
so <- sleuth_fit(so,~ genotype + age + genotype*age, fit_name = 'interaction')
so <- sleuth_fit(so,~ genotype + age, fit_name = 'full')

#Wald test implementations
so <- sleuth_wt(so, which_beta = 'ageold', which_model = 'full')
so <- sleuth_wt(so, which_beta = 'genotypemt', which_model = 'full')
#model with interaction
so <- sleuth_wt(so, which_beta = 'ageold', which_model = 'interaction')
so <- sleuth_wt(so, which_beta = 'genotypemt', which_model = 'interaction')
so <- sleuth_wt(so, which_beta = 'genotypemt:ageold', which_model = 'interaction')
#likelihood test
so <- sleuth_lrt(so, 'full', 'interaction')

#if you want to look at shiny
sleuth_live(so)

#write results to tables
results_table <- sleuth_results(so, 'ageold','interaction', test_type= 'wt')
write.csv(results_table, "~/WormFiles/agingRNAseq/input/si2_aging_analysis.csv")

results_table <- sleuth_results(so, 'genotypemt','interaction', test_type= 'wt')
write.csv(results_table, "~/WormFiles/agingRNAseq/input/si3_genotype_analysis.csv")

results_table <- sleuth_results(so, 'genotypemt:ageold','interaction', test_type= 'wt')
write.csv(results_table, "~/WormFiles/agingRNAseq/input/si4_interaction_analysis.csv")

results_table <- sleuth_results(so, 'full:interaction', test_type='lrt')
write.csv(results_table, "~/WormFiles/agingRNAseq/input/lrt.csv")


