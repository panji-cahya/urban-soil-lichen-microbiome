### R Analysis of Bacterial Fungal Interactions in Soil and Lichen Samples
#For the paper published in npj Biofilms and Microbiomes
#Unveiling the ecological processes driving soil and lichen microbiome assembly along an urbanization gradient
#Authored by: Panji Cahya Mawarda, Rens van der Kaaij, Francisco Dini-Andreote, Deniz Duijker, Michael Stech, Arjen Speksnijder

#Created by: Panji Cahya Mawarda 04-04-2024

#call necessary library packages

#--------- Input data for tempzone one ----

#In Excel, combine temperaturezone1.tsv from bacteria and fungal communities

temperaturezone1_combined=read.csv("temperaturezone1_combined.tsv", row.names=1, header = T, sep = '\t')

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table temperaturezone1_combined.tsv --correlation temperaturezone1_cor.tsv --covariance temperaturezone1_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table temperaturezone1_combined.tsv --number 1000 --prefix bootstrap_counts/temperaturezone1 -t 4
# mkdir bootstrap_correlation
# parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
# fastspar_pvalues --otu_table temperaturezone1.tsv --correlation temperaturezone1_cor.tsv --prefix bootstrap_correlation/cor_temperaturezone1_ --permutations 1000 --outfile temperaturezone1_pval.tsv
# fastspar_reduce -r temperaturezone1_cor.tsv -p temperaturezone1_pval.tsv -o temperaturezone1
# rm -r bootstrap_counts
# rm -r bootstrap_correlation

fs_Rtempone_cor <- read.table('temperaturezone1_filtered_correlation.tsv', 
                              header = T, sep = '\t')

fs_Rtempone_cor <- read.table('temperaturezone1_filtered_correlation.tsv', 
                              header = T, sep = '\t')
fs_Rtempone_p <- read.table('temperaturezone1_filtered_pvalue.tsv', 
                            header = T, sep = '\t')

fs_Rtempone_combine <- cbind.data.frame(fs_Rtempone_cor, fs_Rtempone_p[,3])
colnames(fs_Rtempone_combine)[3:4] = c('cor','pval')

fs_Rtempone_combine %>%
  write.table('./temperaturezone1_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

# Filter based on r and p values ----
fs_Rtempone_filt_nine <- fs_Rtempone_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
fs_Rtempone_filt_eight <- fs_Rtempone_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rtempone_filt_bf <- fs_Rtempone_combine %>% filter(abs(cor) > 0.7 & pval < 0.05)

fs_Rtempone_filt_bf %>%
  write.table('./fs_Rtempone_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Choose bacterial fungal having R2 more than 0.07
fs_Rtempone_filt_bf <- read.table('Temperature1_bacterial_fungal.csv', 
                                  header = T, sep = ",")

# Export to Cytoscape ----
fs_Rtempone_filt_nine %>%
  write.table('./temperaturezone1_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rtempone_filt_eight %>%
  write.table('./temperaturezone1_net2Cys_eight.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_tempone <- temperaturezone1_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_tempone[,-1])

otu_annot_df_tempone_nine <- otu_df_tempone %>%
  filter(OTU %in% unique(c(fs_Rtempone_filt_nine $otu_1, fs_Rtempone_filt_nine $otu2)))

otu_annot_df_tempone_eight <- otu_df_tempone %>%
  filter(OTU %in% unique(c(fs_Rtempone_filt_eight $otu_1, fs_Rtempone_filt_eight $otu2)))

otu_annot_df_tempone_bf <- otu_df_tempone %>%
  filter(OTU %in% unique(c(fs_Rtempone_filt_bf $otu_1, fs_Rtempone_filt_bf $otu2)))

otu_annot_df_tempone_nine$mean_abun <- rowMeans(otu_annot_df_tempone_nine[, -1])
otu_annot_df_mean_tempone_nine <- otu_annot_df_tempone_nine[,c(1,ncol(otu_annot_df_tempone_nine))]

otu_annot_df_tempone_eight$mean_abun <- rowMeans(otu_annot_df_tempone_eight[, -1])
otu_annot_df_mean_tempone_eight <- otu_annot_df_tempone_eight[,c(1,ncol(otu_annot_df_tempone_eight))]

otu_annot_df_tempone_bf$mean_abun <- rowMeans(otu_annot_df_tempone_bf[, -1])
otu_annot_df_mean_tempone_bf <- otu_annot_df_tempone_bf[,c(1,ncol(otu_annot_df_tempone_bf))]

tax_df=read.csv("taxonomy_combined.csv", header = T)

tax_annot_df_tempone_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtempone_filt_nine $otu_1, fs_Rtempone_filt_nine $otu2)))

tax_annot_df_tempone_eight <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtempone_filt_eight $otu_1, fs_Rtempone_filt_eight $otu2)))

tax_annot_df_tempone_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtempone_filt_bf $otu_1, fs_Rtempone_filt_bf $otu2)))

tax_annot_df_tempone_nine %>%
  left_join(otu_annot_df_mean_tempone_nine) %>%
  write.table('temperaturezone1_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_tempone_eight %>%
  left_join(otu_annot_df_mean_tempone_eight) %>%
  write.table('temperaturezone1_eight_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_tempone_bf %>%
  left_join(otu_annot_df_mean_tempone_bf) %>%
  write.table('temperaturezone1_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#--------- Input data for tempzone two ----

#In Excel, combine temperaturezone2.tsv from bacteria and fungal communities

temperaturezone2_combined=read.csv("temperaturezone2_combined.tsv", row.names=1, header = T, sep = '\t')

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table temperaturezone2_combined.tsv --correlation temperaturezone2_cor.tsv --covariance temperaturezone2_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table temperaturezone2_combined.tsv --number 1000 --prefix bootstrap_counts/temperaturezone2 -t 4
# mkdir bootstrap_correlation
# parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
# fastspar_pvalues --otu_table temperaturezone2.tsv --correlation temperaturezone2_cor.tsv --prefix bootstrap_correlation/cor_temperaturezone2_ --permutations 1000 --outfile temperaturezone2_pval.tsv
# fastspar_reduce -r temperaturezone2_cor.tsv -p temperaturezone2_pval.tsv -o temperaturezone2
# rm -r bootstrap_counts
# rm -r bootstrap_correlation

fs_Rtemptwo_cor <- read.table('temperaturezone2_filtered_correlation.tsv', 
                              header = T, sep = '\t')
fs_Rtemptwo_p <- read.table('temperaturezone2_filtered_pvalue.tsv', 
                            header = T, sep = '\t')

fs_Rtemptwo_combine <- cbind.data.frame(fs_Rtemptwo_cor, fs_Rtemptwo_p[,3])
colnames(fs_Rtemptwo_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rtemptwo_filt_nine <- fs_Rtemptwo_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
fs_Rtemptwo_filt_eight <- fs_Rtemptwo_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rtemptwo_filt_bf <- fs_Rtemptwo_combine %>% filter(abs(cor) > 0.7 & pval < 0.05)

fs_Rtemptwo_filt_bf %>%
  write.table('./fs_Rtemptwo_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Choose bacterial fungal having R2 more than 0.07
fs_Rtemptwo_filt_bf <- read.table('Temperaturezone2_bacterial_fungal.csv', 
                                  header = T, sep = ",")

# Export to Cytoscape ----
fs_Rtemptwo_filt_nine %>%
  write.table('./temperaturezone2_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rtemptwo_filt_eight %>%
  write.table('./temperaturezone2_net2Cys_eight.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_temptwo <- temperaturezone2_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_temptwo[,-1])

otu_annot_df_temptwo_nine <- otu_df_temptwo %>%
  filter(OTU %in% unique(c(fs_Rtemptwo_filt_nine $otu_1, fs_Rtemptwo_filt_nine $otu2)))

otu_annot_df_temptwo_eight <- otu_df_temptwo %>%
  filter(OTU %in% unique(c(fs_Rtemptwo_filt_eight $otu_1, fs_Rtemptwo_filt_eight $otu2)))

otu_annot_df_temptwo_bf <- otu_df_temptwo %>%
  filter(OTU %in% unique(c(fs_Rtemptwo_filt_bf $otu_1, fs_Rtemptwo_filt_bf $otu2)))

otu_annot_df_temptwo_nine$mean_abun <- rowMeans(otu_annot_df_temptwo_nine[, -1])
otu_annot_df_mean_temptwo_nine <- otu_annot_df_temptwo_nine[,c(1,ncol(otu_annot_df_temptwo_nine))]

otu_annot_df_temptwo_eight$mean_abun <- rowMeans(otu_annot_df_temptwo_eight[, -1])
otu_annot_df_mean_temptwo_eight <- otu_annot_df_temptwo_eight[,c(1,ncol(otu_annot_df_temptwo_eight))]

otu_annot_df_temptwo_bf$mean_abun <- rowMeans(otu_annot_df_temptwo_bf[, -1])
otu_annot_df_mean_temptwo_bf <- otu_annot_df_temptwo_bf[,c(1,ncol(otu_annot_df_temptwo_bf))]

tax_df=read.csv("taxonomy_combined.csv", header = T)

tax_annot_df_temptwo_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtemptwo_filt_nine $otu_1, fs_Rtemptwo_filt_nine $otu2)))

tax_annot_df_temptwo_eight <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtemptwo_filt_eight $otu_1, fs_Rtemptwo_filt_eight $otu2)))

tax_annot_df_temptwo_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtemptwo_filt_bf $otu_1, fs_Rtemptwo_filt_bf $otu2)))

tax_annot_df_temptwo_nine %>%
  left_join(otu_annot_df_mean_temptwo_nine) %>%
  write.table('temperaturezone2_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_temptwo_eight %>%
  left_join(otu_annot_df_mean_temptwo_eight) %>%
  write.table('temperaturezone2_eight_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_temptwo_bf %>%
  left_join(otu_annot_df_mean_temptwo_bf) %>%
  write.table('temperaturezone2_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#--------- Input data for tempzone three ----
#In Excel, combine temperaturezone3.tsv from bacteria and fungal communities

temperaturezone3_combined=read.csv("temperaturezone3_combined.tsv", row.names=1, header = T, sep = '\t')

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table temperaturezone3_combined.tsv --correlation temperaturezone3_cor.tsv --covariance temperaturezone3_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table temperaturezone3_combined.tsv --number 1000 --prefix bootstrap_counts/temperaturezone3 -t 4
# mkdir bootstrap_correlation
# parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
# fastspar_pvalues --otu_table temperaturezone3.tsv --correlation temperaturezone3_cor.tsv --prefix bootstrap_correlation/cor_temperaturezone3_ --permutations 1000 --outfile temperaturezone3_pval.tsv
# fastspar_reduce -r temperaturezone3_cor.tsv -p temperaturezone3_pval.tsv -o temperaturezone3
# rm -r bootstrap_counts
# rm -r bootstrap_correlation

fs_Rtempthree_cor <- read.table('temperaturezone3_filtered_correlation.tsv', 
                                header = T, sep = '\t')
fs_Rtempthree_p <- read.table('temperaturezone3_filtered_pvalue.tsv', 
                              header = T, sep = '\t')

fs_Rtempthree_combine <- cbind.data.frame(fs_Rtempthree_cor, fs_Rtempthree_p[,3])
colnames(fs_Rtempthree_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rtempthree_filt_nine <- fs_Rtempthree_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
fs_Rtempthree_filt_eight <- fs_Rtempthree_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rtempthree_filt_bf <- fs_Rtempthree_combine %>% filter(abs(cor) > 0.7 & pval < 0.05)

fs_Rtempthree_filt_bf %>%
  write.table('./fs_Rtempthree_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Choose bacterial fungal having R2 more than 0.07
fs_Rtempthree_filt_bf <- read.table('Temperature3_bacterial_fungal.csv', 
                                    header = T, sep = ",")

# Export to Cytoscape ----
fs_Rtempthree_filt_nine %>%
  write.table('./temperaturezone3_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rtempthree_filt_eight %>%
  write.table('./temperaturezone3_net2Cys_eight.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_tempthree <- temperaturezone3_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_tempthree[,-1])

otu_annot_df_tempthree_nine <- otu_df_tempthree %>%
  filter(OTU %in% unique(c(fs_Rtempthree_filt_nine $otu_1, fs_Rtempthree_filt_nine $otu2)))

otu_annot_df_tempthree_eight <- otu_df_tempthree %>%
  filter(OTU %in% unique(c(fs_Rtempthree_filt_eight $otu_1, fs_Rtempthree_filt_eight $otu2)))

otu_annot_df_tempthree_bf <- otu_df_tempthree %>%
  filter(OTU %in% unique(c(fs_Rtempthree_filt_bf $otu_1, fs_Rtempthree_filt_bf $otu2)))

otu_annot_df_tempthree_nine$mean_abun <- rowMeans(otu_annot_df_tempthree_nine[, -1])
otu_annot_df_mean_tempthree_nine <- otu_annot_df_tempthree_nine[,c(1,ncol(otu_annot_df_tempthree_nine))]

otu_annot_df_tempthree_eight$mean_abun <- rowMeans(otu_annot_df_tempthree_eight[, -1])
otu_annot_df_mean_tempthree_eight <- otu_annot_df_tempthree_eight[,c(1,ncol(otu_annot_df_tempthree_eight))]

otu_annot_df_tempthree_bf$mean_abun <- rowMeans(otu_annot_df_tempthree_bf[, -1])
otu_annot_df_mean_tempthree_bf <- otu_annot_df_tempthree_bf[,c(1,ncol(otu_annot_df_tempthree_bf))]

tax_df=read.csv("taxonomy_combined.csv", header = T)

tax_annot_df_tempthree_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtempthree_filt_nine $otu_1, fs_Rtempthree_filt_nine $otu2)))

tax_annot_df_tempthree_eight <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtempthree_filt_eight $otu_1, fs_Rtempthree_filt_eight $otu2)))

tax_annot_df_tempthree_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtempthree_filt_bf $otu_1, fs_Rtempthree_filt_bf $otu2)))

tax_annot_df_tempthree_nine %>%
  left_join(otu_annot_df_mean_tempthree_nine) %>%
  write.table('temperaturezone3_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_tempthree_eight %>%
  left_join(otu_annot_df_mean_tempthree_eight) %>%
  write.table('temperaturezone3_eight_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_tempthree_bf %>%
  left_join(otu_annot_df_mean_tempthree_bf) %>%
  write.table('temperaturezone3_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#--------- Input data for Candelaria concolor ----
#In Excel, combine candelariaconcolor.tsv from bacteria and fungal communities

candelariaconcolor_combined=read.csv("candelariaconcolor_combined.tsv", row.names=1, header = T, sep = '\t')
candelariaconcolor_tempone <- candelariaconcolor_combined[c("T0130_1103", "T0130_1104", "T0130_1105", "T0130_1117")]
candelariaconcolor_temptwo <- candelariaconcolor_combined[c("T0130_1106", "T0130_1114", "T0130_1115")]
candelariaconcolor_tempthree <- candelariaconcolor_combined[c("T0130_1108", "T0130_1109", "T0130_1111", "T0130_1112")]

candelariaconcolor_tempone %>%
  write.table('./candelaria_tempone.tsv', 
              col.names = T, row.names = T, quote = F, sep = '\t')

candelariaconcolor_temptwo %>%
  write.table('./candelaria_temptwo.tsv', 
              col.names = T, row.names = T, quote = F, sep = '\t')

candelariaconcolor_tempthree %>%
  write.table('./candelaria_tempthree.tsv', 
              col.names = T, row.names = T, quote = F, sep = '\t')

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table candelariaconcolor_combined.tsv --correlation candelariaconcolor_cor.tsv --covariance candelariaconcolor_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table candelariaconcolor_combined.tsv --number 1000 --prefix bootstrap_counts/candelariaconcolor -t 4
# mkdir bootstrap_correlation
# parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
# fastspar_pvalues --otu_table candelariaconcolor.tsv --correlation candelariaconcolor_cor.tsv --prefix bootstrap_correlation/cor_candelariaconcolor_ --permutations 1000 --outfile candelariaconcolor_pval.tsv
# fastspar_reduce -r candelariaconcolor_cor.tsv -p candelariaconcolor_pval.tsv -o candelariaconcolor
# rm -r bootstrap_counts
# rm -r bootstrap_correlation

fs_Rcandelaria_cor <- read.table('candelariaconcolor_filtered_correlation.tsv', 
                                 header = T, sep = '\t')
fs_Rcandelaria_p <- read.table('candelariaconcolor_filtered_pvalue.tsv', 
                               header = T, sep = '\t')

fs_Rcandelaria_combine <- cbind.data.frame(fs_Rcandelaria_cor, fs_Rcandelaria_p[,3])
colnames(fs_Rcandelaria_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rcandelaria_filt_nine <- fs_Rcandelaria_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
fs_Rcandelaria_filt_eight <- fs_Rcandelaria_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rcandelaria_filt_bf <- fs_Rcandelaria_combine %>% filter(abs(cor) > 0.8 & pval < 0.01)

fs_Rcandelaria_filt_bf %>%
  write.table('./fs_Rcandelaria_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Choose bacterial fungal having R2 more than 0.07
fs_Rcandelaria_filt_bf <- read.table('candelaria_bacterial_fungal.csv', 
                                     header = T, sep = ",")

# Export to Cytoscape ----
fs_Rcandelaria_filt_nine %>%
  write.table('./candelariaconcolor_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rcandelaria_filt_eight %>%
  write.table('./candelariaconcolor_net2Cys_eight.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_candelaria <- candelariaconcolor_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_candelaria[,-1])

otu_annot_df_candelaria_nine <- otu_df_candelaria %>%
  filter(OTU %in% unique(c(fs_Rcandelaria_filt_nine $otu_1, fs_Rcandelaria_filt_nine $otu2)))

otu_annot_df_candelaria_eight <- otu_df_candelaria %>%
  filter(OTU %in% unique(c(fs_Rcandelaria_filt_eight $otu_1, fs_Rcandelaria_filt_eight $otu2)))

otu_annot_df_candelaria_bf <- otu_df_candelaria %>%
  filter(OTU %in% unique(c(fs_Rcandelaria_filt_bf $otu_1, fs_Rcandelaria_filt_bf $otu2)))

otu_annot_df_candelaria_nine$mean_abun <- rowMeans(otu_annot_df_candelaria_nine[, -1])
otu_annot_df_mean_candelaria_nine <- otu_annot_df_candelaria_nine[,c(1,ncol(otu_annot_df_candelaria_nine))]

otu_annot_df_candelaria_eight$mean_abun <- rowMeans(otu_annot_df_candelaria_eight[, -1])
otu_annot_df_mean_candelaria_eight <- otu_annot_df_candelaria_eight[,c(1,ncol(otu_annot_df_candelaria_eight))]

otu_annot_df_candelaria_bf$mean_abun <- rowMeans(otu_annot_df_candelaria_bf[, -1])
otu_annot_df_mean_candelaria_bf <- otu_annot_df_candelaria_bf[,c(1,ncol(otu_annot_df_candelaria_bf))]

tax_df=read.csv("taxonomy_lichen_combined.csv", header = T)

tax_annot_df_candelaria_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rcandelaria_filt_nine $otu_1, fs_Rcandelaria_filt_nine $otu2)))

tax_annot_df_candelaria_eight <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rcandelaria_filt_eight $otu_1, fs_Rcandelaria_filt_eight $otu2)))

tax_annot_df_candelaria_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rcandelaria_filt_bf $otu_1, fs_Rcandelaria_filt_bf $otu2)))

tax_annot_df_candelaria_nine %>%
  left_join(otu_annot_df_mean_candelaria_nine) %>%
  write.table('candelariaconcolor_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_candelaria_eight %>%
  left_join(otu_annot_df_mean_candelaria_eight) %>%
  write.table('candelariaconcolor_eight_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_candelaria_bf %>%
  left_join(otu_annot_df_mean_candelaria_bf) %>%
  write.table('candelariaconcolor_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#--------- Input data for Physcia adscendens ----
#In Excel, combine physciaadscendens.tsv from bacteria and fungal communities

physciaadscendens_combined=read.csv("physciaadscendens_combined.tsv", row.names=1, header = T, sep = '\t')
physciaadscendens_tempone <- physciaadscendens_combined[c("T0130_1088", "T0130_1092", "T0130_1101", "T0130_1102")]
physciaadscendens_temptwo <- physciaadscendens_combined[c("T0130_1089", "T0130_1090", "T0130_1091", "T0130_1096", "T0130_1099")]
physciaadscendens_tempthree <- physciaadscendens_combined[c("T0130_1093", "T0130_1094", "T0130_1095", "T0130_1097", "T0130_1098", "T0130_1111", "T0130_1112")]

physciaadscendens_tempone %>%
  write.table('./physcia_tempone.tsv', 
              col.names = T, row.names = T, quote = F, sep = '\t')

physciaadscendens_temptwo %>%
  write.table('./physcia_temptwo.tsv', 
              col.names = T, row.names = T, quote = F, sep = '\t')

physciaadscendens_tempthree %>%
  write.table('./physcia_tempthree.tsv', 
              col.names = T, row.names = T, quote = F, sep = '\t')

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table physciaadscendens_combined.tsv --correlation physciaadscendens_cor.tsv --covariance physciaadscendens_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table physciaadscendens_combined.tsv --number 1000 --prefix bootstrap_counts/physciaadscendens -t 4
# mkdir bootstrap_correlation
# parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
# fastspar_pvalues --otu_table physciaadscendens.tsv --correlation physciaadscendens_cor.tsv --prefix bootstrap_correlation/cor_physciaadscendens_ --permutations 1000 --outfile physciaadscendens_pval.tsv
# fastspar_reduce -r physciaadscendens_cor.tsv -p physciaadscendens_pval.tsv -o physciaadscendens
# rm -r bootstrap_counts
# rm -r bootstrap_correlation

fs_Rphyscia_cor <- read.table('physciaadscendens_filtered_correlation.tsv', 
                              header = T, sep = '\t')
fs_Rphyscia_p <- read.table('physciaadscendens_filtered_pvalue.tsv', 
                            header = T, sep = '\t')

fs_Rphyscia_combine <- cbind.data.frame(fs_Rphyscia_cor, fs_Rphyscia_p[,3])
colnames(fs_Rphyscia_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rphyscia_filt_nine <- fs_Rphyscia_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
fs_Rphyscia_filt_eight <- fs_Rphyscia_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rphyscia_filt_bf <- fs_Rphyscia_combine %>% filter(abs(cor) > 0.8 & pval < 0.01)

fs_Rphyscia_filt_bf %>%
  write.table('./fs_Rphyscia_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Choose bacterial fungal having R2 more than 0.07
fs_Rphyscia_filt_bf <- read.table('physcia_bacterial_fungal.csv', 
                                  header = T, sep = ",")

# Export to Cytoscape ----
fs_Rphyscia_filt_nine %>%
  write.table('./physciaadscendens_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rphyscia_filt_eight %>%
  write.table('./physciaadscendens_net2Cys_eight.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_physcia <- physciaadscendens_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_physcia[,-1])

otu_annot_df_physcia_nine <- otu_df_physcia %>%
  filter(OTU %in% unique(c(fs_Rphyscia_filt_nine $otu_1, fs_Rphyscia_filt_nine $otu2)))

otu_annot_df_physcia_eight <- otu_df_physcia %>%
  filter(OTU %in% unique(c(fs_Rphyscia_filt_eight $otu_1, fs_Rphyscia_filt_eight $otu2)))

otu_annot_df_physcia_bf <- otu_df_physcia %>%
  filter(OTU %in% unique(c(fs_Rphyscia_filt_bf $otu_1, fs_Rphyscia_filt_bf $otu2)))

otu_annot_df_physcia_nine$mean_abun <- rowMeans(otu_annot_df_physcia_nine[, -1])
otu_annot_df_mean_physcia_nine <- otu_annot_df_physcia_nine[,c(1,ncol(otu_annot_df_physcia_nine))]

otu_annot_df_physcia_eight$mean_abun <- rowMeans(otu_annot_df_physcia_eight[, -1])
otu_annot_df_mean_physcia_eight <- otu_annot_df_physcia_eight[,c(1,ncol(otu_annot_df_physcia_eight))]

otu_annot_df_physcia_bf$mean_abun <- rowMeans(otu_annot_df_physcia_bf[, -1])
otu_annot_df_mean_physcia_bf <- otu_annot_df_physcia_bf[,c(1,ncol(otu_annot_df_physcia_bf))]

tax_df=read.csv("taxonomy_lichen_combined.csv", header = T)

tax_annot_df_physcia_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rphyscia_filt_nine $otu_1, fs_Rphyscia_filt_nine $otu2)))

tax_annot_df_physcia_eight <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rphyscia_filt_eight $otu_1, fs_Rphyscia_filt_eight $otu2)))

tax_annot_df_physcia_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rphyscia_filt_bf $otu_1, fs_Rphyscia_filt_bf $otu2)))

tax_annot_df_physcia_nine %>%
  left_join(otu_annot_df_mean_physcia_nine) %>%
  write.table('physciaadscendens_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_physcia_eight %>%
  left_join(otu_annot_df_mean_physcia_eight) %>%
  write.table('physciaadscendens_eight_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_physcia_bf %>%
  left_join(otu_annot_df_mean_physcia_bf) %>%
  write.table('physciaadscendens_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#--------- Input data for Xanthoria parietina ----
#In Excel, combine xanthoriaparietina.tsv from bacteria and fungal communities

xanthoriaparietina_combined=read.csv("xanthoriaparietina_combined.tsv", row.names=1, header = T, sep = '\t')
xanthoriaparietina_tempone <- xanthoriaparietina_combined[c("T0130_1237", "T0130_1238", "T0130_1239", "T0130_1241", "T0130_1242")]
xanthoriaparietina_temptwo <- xanthoriaparietina_combined[c("T0130_1240", "T0130_1243", "T0130_1244", "T0130_1245", "T0130_1246")]
xanthoriaparietina_tempthree <- xanthoriaparietina_combined[c("T0130_1247", "T0130_1248", "T0130_1249", "T0130_1250", "T0130_1251")]


xanthoriaparietina_tempone %>%
  write.table('./xanthoria_tempone.tsv', 
              col.names = T, row.names = T, quote = F, sep = '\t')

xanthoriaparietina_temptwo %>%
  write.table('./xanthoria_temptwo.tsv', 
              col.names = T, row.names = T, quote = F, sep = '\t')

xanthoriaparietina_tempthree %>%
  write.table('./xanthoria_tempthree.tsv', 
              col.names = T, row.names = T, quote = F, sep = '\t')

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table xanthoriaparietina_combined.tsv --correlation xanthoriaparietina_cor.tsv --covariance xanthoriaparietina_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table xanthoriaparietina_combined.tsv --number 1000 --prefix bootstrap_counts/xanthoriaparietina -t 4
# mkdir bootstrap_correlation
# parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
# fastspar_pvalues --otu_table xanthoriaparietina.tsv --correlation xanthoriaparietina_cor.tsv --prefix bootstrap_correlation/cor_xanthoriaparietina_ --permutations 1000 --outfile xanthoriaparietina_pval.tsv
# fastspar_reduce -r xanthoriaparietina_cor.tsv -p xanthoriaparietina_pval.tsv -o xanthoriaparietina
# rm -r bootstrap_counts
# rm -r bootstrap_correlation

fs_Rxanthoria_cor <- read.table('xanthoriaparietina_filtered_correlation.tsv', 
                                header = T, sep = '\t')
fs_Rxanthoria_p <- read.table('xanthoriaparietina_filtered_pvalue.tsv', 
                              header = T, sep = '\t')

fs_Rxanthoria_combine <- cbind.data.frame(fs_Rxanthoria_cor, fs_Rxanthoria_p[,3])
colnames(fs_Rxanthoria_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rxanthoria_filt_nine <- fs_Rxanthoria_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
fs_Rxanthoria_filt_eight <- fs_Rxanthoria_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rxanthoria_filt_bf <- fs_Rxanthoria_combine %>% filter(abs(cor) > 0.8 & pval < 0.01)

fs_Rxanthoria_filt_bf %>%
  write.table('./fs_Rxanthoria_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Choose bacterial fungal having R2 more than 0.07
fs_Rxanthoria_filt_bf <- read.table('xanthoria_bacterial_fungal.csv', 
                                    header = T, sep = ",")

# Export to Cytoscape ----
fs_Rxanthoria_filt_nine %>%
  write.table('./xanthoriaparietina_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rxanthoria_filt_eight %>%
  write.table('./xanthoriaparietina_net2Cys_eight.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_xanthoria <- xanthoriaparietina_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_xanthoria[,-1])

otu_annot_df_xanthoria_nine <- otu_df_xanthoria %>%
  filter(OTU %in% unique(c(fs_Rxanthoria_filt_nine $otu_1, fs_Rxanthoria_filt_nine $otu2)))

otu_annot_df_xanthoria_eight <- otu_df_xanthoria %>%
  filter(OTU %in% unique(c(fs_Rxanthoria_filt_eight $otu_1, fs_Rxanthoria_filt_eight $otu2)))

otu_annot_df_xanthoria_bf <- otu_df_xanthoria %>%
  filter(OTU %in% unique(c(fs_Rxanthoria_filt_bf $otu_1, fs_Rxanthoria_filt_bf $otu2)))

otu_annot_df_xanthoria_nine$mean_abun <- rowMeans(otu_annot_df_xanthoria_nine[, -1])
otu_annot_df_mean_xanthoria_nine <- otu_annot_df_xanthoria_nine[,c(1,ncol(otu_annot_df_xanthoria_nine))]

otu_annot_df_xanthoria_eight$mean_abun <- rowMeans(otu_annot_df_xanthoria_eight[, -1])
otu_annot_df_mean_xanthoria_eight <- otu_annot_df_xanthoria_eight[,c(1,ncol(otu_annot_df_xanthoria_eight))]

otu_annot_df_xanthoria_bf$mean_abun <- rowMeans(otu_annot_df_xanthoria_bf[, -1])
otu_annot_df_mean_xanthoria_bf <- otu_annot_df_xanthoria_bf[,c(1,ncol(otu_annot_df_xanthoria_bf))]

tax_df=read.csv("taxonomy_lichen_combined.csv", header = T)

tax_annot_df_xanthoria_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rxanthoria_filt_nine $otu_1, fs_Rxanthoria_filt_nine $otu2)))

tax_annot_df_xanthoria_eight <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rxanthoria_filt_eight $otu_1, fs_Rxanthoria_filt_eight $otu2)))

tax_annot_df_xanthoria_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rxanthoria_filt_bf $otu_1, fs_Rxanthoria_filt_bf $otu2)))

tax_annot_df_xanthoria_nine %>%
  left_join(otu_annot_df_mean_xanthoria_nine) %>%
  write.table('xanthoriaparietina_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_xanthoria_eight %>%
  left_join(otu_annot_df_mean_xanthoria_eight) %>%
  write.table('xanthoriaparietina_eight_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_xanthoria_bf %>%
  left_join(otu_annot_df_mean_xanthoria_bf) %>%
  write.table('xanthoriaparietina_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')


#Candelaria Temperature One --------------

fs_Rcandelariatempone_cor <- read.table('candelariatempone_filtered_correlation.tsv', 
                                        header = T, sep = '\t')
fs_Rcandelariatempone_p <- read.table('candelariatempone_filtered_pvalue.tsv', 
                                      header = T, sep = '\t')

fs_Rcandelariatempone_combine <- cbind.data.frame(fs_Rcandelariatempone_cor, fs_Rcandelariatempone_p[,3])
colnames(fs_Rcandelariatempone_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rcandelariatempone_filt_nine <- fs_Rcandelariatempone_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
#fs_Rcandelariatempone_filt_eight <- fs_Rcandelariatempone_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rcandelariatempone_filt_bf <- fs_Rcandelariatempone_combine %>% filter(abs(cor) > 0.8 & pval < 0.01)


# Export to Cytoscape ----
fs_Rcandelariatempone_filt_nine %>%
  write.table('./candelariaconcolortempone_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rcandelariatempone_filt_bf %>%
  write.table('./fs_Rcandelariatempone_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Create file in excel combining the first and second _nine and _bf
fs_Rcandelariatempone_filt_bf <- read.table('mycobiont_photobiont_candelaria_tempone.tsv', 
                                            header = T, sep = '\t')

#fs_Rcandelariatempone_filt_eight %>%
#  write.table('./candelariaconcolortempone_net2Cys_eight.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_candelaria <- candelariaconcolor_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_candelaria[,-1])

otu_annot_df_candelariatempone_nine <- otu_df_candelaria %>%
  filter(OTU %in% unique(c(fs_Rcandelariatempone_filt_nine $otu_1, fs_Rcandelariatempone_filt_nine $otu2)))

#otu_annot_df_candelaria_eight <- otu_df_candelaria %>%
#  filter(OTU %in% unique(c(fs_Rcandelariatempone_filt_eight $otu_1, fs_Rcandelariatempone_filt_eight $otu2)))

otu_annot_df_candelariatempone_bf <- otu_df_candelaria %>%
  filter(OTU %in% unique(c(fs_Rcandelariatempone_filt_bf $otu_1, fs_Rcandelariatempone_filt_bf $otu2)))

otu_annot_df_candelariatempone_nine$mean_abun <- rowMeans(otu_annot_df_candelariatempone_nine[, -1])
otu_annot_df_mean_candelariatempone_nine <- otu_annot_df_candelariatempone_nine[,c(1,ncol(otu_annot_df_candelariatempone_nine))]

#otu_annot_df_candelaria_eight$mean_abun <- rowMeans(otu_annot_df_candelaria_eight[, -1])
#otu_annot_df_mean_candelaria_eight <- otu_annot_df_candelaria_eight[,c(1,ncol(otu_annot_df_candelaria_eight))]

otu_annot_df_candelariatempone_bf$mean_abun <- rowMeans(otu_annot_df_candelariatempone_bf[, -1])
otu_annot_df_mean_candelariatempone_bf <- otu_annot_df_candelariatempone_bf[,c(1,ncol(otu_annot_df_candelariatempone_bf))]

tax_df=read.csv("taxonomy_lichen_combined.csv", header = T)

tax_annot_df_candelariatempone_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rcandelariatempone_filt_nine $otu_1, fs_Rcandelariatempone_filt_nine $otu2)))

#tax_annot_df_candelaria_eight <- tax_df %>%
#  filter(OTU %in% unique(c(fs_Rcandelaria_filt_eight $otu_1, fs_Rcandelaria_filt_eight $otu2)))

tax_annot_df_candelariatempone_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rcandelariatempone_filt_bf $otu_1, fs_Rcandelariatempone_filt_bf $otu2)))

tax_annot_df_candelariatempone_nine %>%
  left_join(otu_annot_df_mean_candelariatempone_nine) %>%
  write.table('candelariaconcolortempone_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#tax_annot_df_candelaria_eight %>%
#  left_join(otu_annot_df_mean_candelaria_eight) %>%
#  write.table('candelariaconcolor_eight_annot_abundance.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_candelariatempone_bf %>%
  left_join(otu_annot_df_mean_candelariatempone_bf) %>%
  write.table('candelariaconcolortempone_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Candelaria Temperature Two --------------

fs_Rcandelariatemptwo_cor <- read.table('candelariatemptwo_filtered_correlation.tsv', 
                                        header = T, sep = '\t')
fs_Rcandelariatemptwo_p <- read.table('candelariatemptwo_filtered_pvalue.tsv', 
                                      header = T, sep = '\t')

fs_Rcandelariatemptwo_combine <- cbind.data.frame(fs_Rcandelariatemptwo_cor, fs_Rcandelariatemptwo_p[,3])
colnames(fs_Rcandelariatemptwo_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rcandelariatemptwo_filt_nine <- fs_Rcandelariatemptwo_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
#fs_Rcandelariatemptwo_filt_eight <- fs_Rcandelariatemptwo_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rcandelariatemptwo_filt_bf <- fs_Rcandelariatemptwo_combine %>% filter(abs(cor) > 0.8 & pval < 0.01)


# Export to Cytoscape ----
fs_Rcandelariatemptwo_filt_nine %>%
  write.table('./candelariaconcolortemptwo_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rcandelariatemptwo_filt_bf %>%
  write.table('./fs_Rcandelariatemptwo_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Create file in excel combining the first and second _nine and _bf
fs_Rcandelariatemptwo_filt_bf <- read.table('mycobiont_photobiont_candelaria_temptwo.tsv', 
                                            header = T, sep = '\t')

#fs_Rcandelariatemptwo_filt_eight %>%
#  write.table('./candelariaconcolortemptwo_net2Cys_eight.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_candelaria <- candelariaconcolor_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_candelaria[,-1])

otu_annot_df_candelariatemptwo_nine <- otu_df_candelaria %>%
  filter(OTU %in% unique(c(fs_Rcandelariatemptwo_filt_nine $otu_1, fs_Rcandelariatemptwo_filt_nine $otu2)))

#otu_annot_df_candelaria_eight <- otu_df_candelaria %>%
#  filter(OTU %in% unique(c(fs_Rcandelariatemptwo_filt_eight $otu_1, fs_Rcandelariatemptwo_filt_eight $otu2)))

otu_annot_df_candelariatemptwo_bf <- otu_df_candelaria %>%
  filter(OTU %in% unique(c(fs_Rcandelariatemptwo_filt_bf $otu_1, fs_Rcandelariatemptwo_filt_bf $otu2)))

otu_annot_df_candelariatemptwo_nine$mean_abun <- rowMeans(otu_annot_df_candelariatemptwo_nine[, -1])
otu_annot_df_mean_candelariatemptwo_nine <- otu_annot_df_candelariatemptwo_nine[,c(1,ncol(otu_annot_df_candelariatemptwo_nine))]

#otu_annot_df_candelaria_eight$mean_abun <- rowMeans(otu_annot_df_candelaria_eight[, -1])
#otu_annot_df_mean_candelaria_eight <- otu_annot_df_candelaria_eight[,c(1,ncol(otu_annot_df_candelaria_eight))]

otu_annot_df_candelariatemptwo_bf$mean_abun <- rowMeans(otu_annot_df_candelariatemptwo_bf[, -1])
otu_annot_df_mean_candelariatemptwo_bf <- otu_annot_df_candelariatemptwo_bf[,c(1,ncol(otu_annot_df_candelariatemptwo_bf))]

tax_df=read.csv("taxonomy_lichen_combined.csv", header = T)

tax_annot_df_candelariatemptwo_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rcandelariatemptwo_filt_nine $otu_1, fs_Rcandelariatemptwo_filt_nine $otu2)))

#tax_annot_df_candelaria_eight <- tax_df %>%
#  filter(OTU %in% unique(c(fs_Rcandelaria_filt_eight $otu_1, fs_Rcandelaria_filt_eight $otu2)))

tax_annot_df_candelariatemptwo_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rcandelariatemptwo_filt_bf $otu_1, fs_Rcandelariatemptwo_filt_bf $otu2)))

tax_annot_df_candelariatemptwo_nine %>%
  left_join(otu_annot_df_mean_candelariatemptwo_nine) %>%
  write.table('candelariaconcolortemptwo_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#tax_annot_df_candelaria_eight %>%
#  left_join(otu_annot_df_mean_candelaria_eight) %>%
#  write.table('candelariaconcolor_eight_annot_abundance.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_candelariatemptwo_bf %>%
  left_join(otu_annot_df_mean_candelariatemptwo_bf) %>%
  write.table('mycobiont_photobiont_candelariaconcolortemptwo_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Candelaria Temperature Three --------------

fs_Rcandelariatempthree_cor <- read.table('candelariatempthree_filtered_correlation.tsv', 
                                          header = T, sep = '\t')
fs_Rcandelariatempthree_p <- read.table('candelariatempthree_filtered_pvalue.tsv', 
                                        header = T, sep = '\t')

fs_Rcandelariatempthree_combine <- cbind.data.frame(fs_Rcandelariatempthree_cor, fs_Rcandelariatempthree_p[,3])
colnames(fs_Rcandelariatempthree_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rcandelariatempthree_filt_nine <- fs_Rcandelariatempthree_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
#fs_Rcandelariatempthree_filt_eight <- fs_Rcandelariatempthree_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rcandelariatempthree_filt_bf <- fs_Rcandelariatempthree_combine %>% filter(abs(cor) > 0.8 & pval < 0.01)


# Export to Cytoscape ----
fs_Rcandelariatempthree_filt_nine %>%
  write.table('./candelariaconcolortempthree_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rcandelariatempthree_filt_bf %>%
  write.table('./fs_Rcandelariatempthree_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Create file in excel combining the first and second _nine and _bf
fs_Rcandelariatempthree_filt_bf <- read.table('mycobiont_photobiont_candelaria_tempthree.tsv', 
                                              header = T, sep = '\t')

#fs_Rcandelariatempthree_filt_eight %>%
#  write.table('./candelariaconcolortempthree_net2Cys_eight.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_candelaria <- candelariaconcolor_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_candelaria[,-1])

otu_annot_df_candelariatempthree_nine <- otu_df_candelaria %>%
  filter(OTU %in% unique(c(fs_Rcandelariatempthree_filt_nine $otu_1, fs_Rcandelariatempthree_filt_nine $otu2)))

#otu_annot_df_candelaria_eight <- otu_df_candelaria %>%
#  filter(OTU %in% unique(c(fs_Rcandelariatempthree_filt_eight $otu_1, fs_Rcandelariatempthree_filt_eight $otu2)))

otu_annot_df_candelariatempthree_bf <- otu_df_candelaria %>%
  filter(OTU %in% unique(c(fs_Rcandelariatempthree_filt_bf $otu_1, fs_Rcandelariatempthree_filt_bf $otu2)))

otu_annot_df_candelariatempthree_nine$mean_abun <- rowMeans(otu_annot_df_candelariatempthree_nine[, -1])
otu_annot_df_mean_candelariatempthree_nine <- otu_annot_df_candelariatempthree_nine[,c(1,ncol(otu_annot_df_candelariatempthree_nine))]

#otu_annot_df_candelaria_eight$mean_abun <- rowMeans(otu_annot_df_candelaria_eight[, -1])
#otu_annot_df_mean_candelaria_eight <- otu_annot_df_candelaria_eight[,c(1,ncol(otu_annot_df_candelaria_eight))]

otu_annot_df_candelariatempthree_bf$mean_abun <- rowMeans(otu_annot_df_candelariatempthree_bf[, -1])
otu_annot_df_mean_candelariatempthree_bf <- otu_annot_df_candelariatempthree_bf[,c(1,ncol(otu_annot_df_candelariatempthree_bf))]

tax_df=read.csv("taxonomy_lichen_combined.csv", header = T)

tax_annot_df_candelariatempthree_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rcandelariatempthree_filt_nine $otu_1, fs_Rcandelariatempthree_filt_nine $otu2)))

#tax_annot_df_candelaria_eight <- tax_df %>%
#  filter(OTU %in% unique(c(fs_Rcandelaria_filt_eight $otu_1, fs_Rcandelaria_filt_eight $otu2)))

tax_annot_df_candelariatempthree_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rcandelariatempthree_filt_bf $otu_1, fs_Rcandelariatempthree_filt_bf $otu2)))

tax_annot_df_candelariatempthree_nine %>%
  left_join(otu_annot_df_mean_candelariatempthree_nine) %>%
  write.table('candelariaconcolortempthree_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#tax_annot_df_candelaria_eight %>%
#  left_join(otu_annot_df_mean_candelaria_eight) %>%
#  write.table('candelariaconcolor_eight_annot_abundance.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_candelariatempthree_bf %>%
  left_join(otu_annot_df_mean_candelariatempthree_bf) %>%
  write.table('mycobiont_photobiont_candelariaconcolortempthree_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#xanthoria Temperature One --------------

fs_Rxanthoriatempone_cor <- read.table('xanthoriatempone_filtered_correlation.tsv', 
                                       header = T, sep = '\t')
fs_Rxanthoriatempone_p <- read.table('xanthoriatempone_filtered_pvalue.tsv', 
                                     header = T, sep = '\t')

fs_Rxanthoriatempone_combine <- cbind.data.frame(fs_Rxanthoriatempone_cor, fs_Rxanthoriatempone_p[,3])
colnames(fs_Rxanthoriatempone_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rxanthoriatempone_filt_nine <- fs_Rxanthoriatempone_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
#fs_Rxanthoriatempone_filt_eight <- fs_Rxanthoriatempone_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rxanthoriatempone_filt_bf <- fs_Rxanthoriatempone_combine %>% filter(abs(cor) > 0.8 & pval < 0.01)


# Export to Cytoscape ----
fs_Rxanthoriatempone_filt_nine %>%
  write.table('./xanthoriaparietinatempone_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rxanthoriatempone_filt_bf %>%
  write.table('./fs_Rxanthoriatempone_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Create file in excel combining the first and second _nine and _bf
fs_Rxanthoriatempone_filt_bf <- read.table('mycobiont_photobiont_xanthoria_tempone.tsv', 
                                           header = T, sep = '\t')

#fs_Rxanthoriatempone_filt_eight %>%
#  write.table('./xanthoriaparietinatempone_net2Cys_eight.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_xanthoria <- xanthoriaparietina_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_xanthoria[,-1])

otu_annot_df_xanthoriatempone_nine <- otu_df_xanthoria %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatempone_filt_nine $otu_1, fs_Rxanthoriatempone_filt_nine $otu2)))

#otu_annot_df_xanthoria_eight <- otu_df_xanthoria %>%
#  filter(OTU %in% unique(c(fs_Rxanthoriatempone_filt_eight $otu_1, fs_Rxanthoriatempone_filt_eight $otu2)))

otu_annot_df_xanthoriatempone_bf <- otu_df_xanthoria %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatempone_filt_bf $otu_1, fs_Rxanthoriatempone_filt_bf $otu2)))

otu_annot_df_xanthoriatempone_nine$mean_abun <- rowMeans(otu_annot_df_xanthoriatempone_nine[, -1])
otu_annot_df_mean_xanthoriatempone_nine <- otu_annot_df_xanthoriatempone_nine[,c(1,ncol(otu_annot_df_xanthoriatempone_nine))]

#otu_annot_df_xanthoria_eight$mean_abun <- rowMeans(otu_annot_df_xanthoria_eight[, -1])
#otu_annot_df_mean_xanthoria_eight <- otu_annot_df_xanthoria_eight[,c(1,ncol(otu_annot_df_xanthoria_eight))]

otu_annot_df_xanthoriatempone_bf$mean_abun <- rowMeans(otu_annot_df_xanthoriatempone_bf[, -1])
otu_annot_df_mean_xanthoriatempone_bf <- otu_annot_df_xanthoriatempone_bf[,c(1,ncol(otu_annot_df_xanthoriatempone_bf))]

tax_df=read.csv("taxonomy_lichen_combined.csv", header = T)

tax_annot_df_xanthoriatempone_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatempone_filt_nine $otu_1, fs_Rxanthoriatempone_filt_nine $otu2)))

#tax_annot_df_xanthoria_eight <- tax_df %>%
#  filter(OTU %in% unique(c(fs_Rxanthoria_filt_eight $otu_1, fs_Rxanthoria_filt_eight $otu2)))

tax_annot_df_xanthoriatempone_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatempone_filt_bf $otu_1, fs_Rxanthoriatempone_filt_bf $otu2)))

tax_annot_df_xanthoriatempone_nine %>%
  left_join(otu_annot_df_mean_xanthoriatempone_nine) %>%
  write.table('xanthoriaparietinatempone_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#tax_annot_df_xanthoria_eight %>%
#  left_join(otu_annot_df_mean_xanthoria_eight) %>%
#  write.table('xanthoriaparietina_eight_annot_abundance.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_xanthoriatempone_bf %>%
  left_join(otu_annot_df_mean_xanthoriatempone_bf) %>%
  write.table('mycobiont_photobiont_xanthoriaparietinatempone_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#xanthoria Temperature Two --------------

fs_Rxanthoriatemptwo_cor <- read.table('xanthoriatemptwo_filtered_correlation.tsv', 
                                       header = T, sep = '\t')
fs_Rxanthoriatemptwo_p <- read.table('xanthoriatemptwo_filtered_pvalue.tsv', 
                                     header = T, sep = '\t')

fs_Rxanthoriatemptwo_combine <- cbind.data.frame(fs_Rxanthoriatemptwo_cor, fs_Rxanthoriatemptwo_p[,3])
colnames(fs_Rxanthoriatemptwo_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rxanthoriatemptwo_filt_nine <- fs_Rxanthoriatemptwo_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
#fs_Rxanthoriatemptwo_filt_eight <- fs_Rxanthoriatemptwo_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rxanthoriatemptwo_filt_bf <- fs_Rxanthoriatemptwo_combine %>% filter(abs(cor) > 0.8 & pval < 0.01)


# Export to Cytoscape ----
fs_Rxanthoriatemptwo_filt_nine %>%
  write.table('./xanthoriaparietinatemptwo_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rxanthoriatemptwo_filt_bf %>%
  write.table('./fs_Rxanthoriatemptwo_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Create file in excel combining the first and second _nine and _bf
fs_Rxanthoriatemptwo_filt_bf <- read.table('mycobiont_photobiont_xanthoria_temptwo.tsv', 
                                           header = T, sep = '\t')

#fs_Rxanthoriatemptwo_filt_eight %>%
#  write.table('./xanthoriaparietinatemptwo_net2Cys_eight.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_xanthoria <- xanthoriaparietina_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_xanthoria[,-1])

otu_annot_df_xanthoriatemptwo_nine <- otu_df_xanthoria %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatemptwo_filt_nine $otu_1, fs_Rxanthoriatemptwo_filt_nine $otu2)))

#otu_annot_df_xanthoria_eight <- otu_df_xanthoria %>%
#  filter(OTU %in% unique(c(fs_Rxanthoriatemptwo_filt_eight $otu_1, fs_Rxanthoriatemptwo_filt_eight $otu2)))

otu_annot_df_xanthoriatemptwo_bf <- otu_df_xanthoria %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatemptwo_filt_bf $otu_1, fs_Rxanthoriatemptwo_filt_bf $otu2)))

otu_annot_df_xanthoriatemptwo_nine$mean_abun <- rowMeans(otu_annot_df_xanthoriatemptwo_nine[, -1])
otu_annot_df_mean_xanthoriatemptwo_nine <- otu_annot_df_xanthoriatemptwo_nine[,c(1,ncol(otu_annot_df_xanthoriatemptwo_nine))]

#otu_annot_df_xanthoria_eight$mean_abun <- rowMeans(otu_annot_df_xanthoria_eight[, -1])
#otu_annot_df_mean_xanthoria_eight <- otu_annot_df_xanthoria_eight[,c(1,ncol(otu_annot_df_xanthoria_eight))]

otu_annot_df_xanthoriatemptwo_bf$mean_abun <- rowMeans(otu_annot_df_xanthoriatemptwo_bf[, -1])
otu_annot_df_mean_xanthoriatemptwo_bf <- otu_annot_df_xanthoriatemptwo_bf[,c(1,ncol(otu_annot_df_xanthoriatemptwo_bf))]

tax_df=read.csv("taxonomy_lichen_combined.csv", header = T)

tax_annot_df_xanthoriatemptwo_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatemptwo_filt_nine $otu_1, fs_Rxanthoriatemptwo_filt_nine $otu2)))

#tax_annot_df_xanthoria_eight <- tax_df %>%
#  filter(OTU %in% unique(c(fs_Rxanthoria_filt_eight $otu_1, fs_Rxanthoria_filt_eight $otu2)))

tax_annot_df_xanthoriatemptwo_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatemptwo_filt_bf $otu_1, fs_Rxanthoriatemptwo_filt_bf $otu2)))

tax_annot_df_xanthoriatemptwo_nine %>%
  left_join(otu_annot_df_mean_xanthoriatemptwo_nine) %>%
  write.table('xanthoriaparietinatemptwo_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#tax_annot_df_xanthoria_eight %>%
#  left_join(otu_annot_df_mean_xanthoria_eight) %>%
#  write.table('xanthoriaparietina_eight_annot_abundance.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_xanthoriatemptwo_bf %>%
  left_join(otu_annot_df_mean_xanthoriatemptwo_bf) %>%
  write.table('mycobiont_photobiont_xanthoriaparietinatemptwo_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#xanthoria Temperature Three --------------

fs_Rxanthoriatempthree_cor <- read.table('xanthoriatempthree_filtered_correlation.tsv', 
                                         header = T, sep = '\t')
fs_Rxanthoriatempthree_p <- read.table('xanthoriatempthree_filtered_pvalue.tsv', 
                                       header = T, sep = '\t')

fs_Rxanthoriatempthree_combine <- cbind.data.frame(fs_Rxanthoriatempthree_cor, fs_Rxanthoriatempthree_p[,3])
colnames(fs_Rxanthoriatempthree_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rxanthoriatempthree_filt_nine <- fs_Rxanthoriatempthree_combine %>% filter(abs(cor) > 0.9 & pval < 0.01)
#fs_Rxanthoriatempthree_filt_eight <- fs_Rxanthoriatempthree_combine %>% filter(abs(cor) > 0.85 & pval < 0.01)

fs_Rxanthoriatempthree_filt_bf <- fs_Rxanthoriatempthree_combine %>% filter(abs(cor) > 0.8 & pval < 0.01)


# Export to Cytoscape ----
fs_Rxanthoriatempthree_filt_nine %>%
  write.table('./xanthoriaparietinatempthree_net2Cys_nine.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

fs_Rxanthoriatempthree_filt_bf %>%
  write.table('./fs_Rxanthoriatempthree_filt_bf.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Create file in excel combining the first and second _nine and _bf
fs_Rxanthoriatempthree_filt_bf <- read.table('mycobiont_photobiont_xanthoria_tempthree.tsv', 
                                             header = T, sep = '\t')

#fs_Rxanthoriatempthree_filt_eight %>%
#  write.table('./xanthoriaparietinatempthree_net2Cys_eight.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_xanthoria <- xanthoriaparietina_combined %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_xanthoria[,-1])

otu_annot_df_xanthoriatempthree_nine <- otu_df_xanthoria %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatempthree_filt_nine $otu_1, fs_Rxanthoriatempthree_filt_nine $otu2)))

#otu_annot_df_xanthoria_eight <- otu_df_xanthoria %>%
#  filter(OTU %in% unique(c(fs_Rxanthoriatempthree_filt_eight $otu_1, fs_Rxanthoriatempthree_filt_eight $otu2)))

otu_annot_df_xanthoriatempthree_bf <- otu_df_xanthoria %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatempthree_filt_bf $otu_1, fs_Rxanthoriatempthree_filt_bf $otu2)))

otu_annot_df_xanthoriatempthree_nine$mean_abun <- rowMeans(otu_annot_df_xanthoriatempthree_nine[, -1])
otu_annot_df_mean_xanthoriatempthree_nine <- otu_annot_df_xanthoriatempthree_nine[,c(1,ncol(otu_annot_df_xanthoriatempthree_nine))]

#otu_annot_df_xanthoria_eight$mean_abun <- rowMeans(otu_annot_df_xanthoria_eight[, -1])
#otu_annot_df_mean_xanthoria_eight <- otu_annot_df_xanthoria_eight[,c(1,ncol(otu_annot_df_xanthoria_eight))]

otu_annot_df_xanthoriatempthree_bf$mean_abun <- rowMeans(otu_annot_df_xanthoriatempthree_bf[, -1])
otu_annot_df_mean_xanthoriatempthree_bf <- otu_annot_df_xanthoriatempthree_bf[,c(1,ncol(otu_annot_df_xanthoriatempthree_bf))]

tax_df=read.csv("taxonomy_lichen_combined.csv", header = T)

tax_annot_df_xanthoriatempthree_nine <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatempthree_filt_nine $otu_1, fs_Rxanthoriatempthree_filt_nine $otu2)))

#tax_annot_df_xanthoria_eight <- tax_df %>%
#  filter(OTU %in% unique(c(fs_Rxanthoria_filt_eight $otu_1, fs_Rxanthoria_filt_eight $otu2)))

tax_annot_df_xanthoriatempthree_bf <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rxanthoriatempthree_filt_bf $otu_1, fs_Rxanthoriatempthree_filt_bf $otu2)))

tax_annot_df_xanthoriatempthree_nine %>%
  left_join(otu_annot_df_mean_xanthoriatempthree_nine) %>%
  write.table('xanthoriaparietinatempthree_nine_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#tax_annot_df_xanthoria_eight %>%
#  left_join(otu_annot_df_mean_xanthoria_eight) %>%
#  write.table('xanthoriaparietina_eight_annot_abundance.tsv', 
#              col.names = T, row.names = F, quote = F, sep = '\t')

tax_annot_df_xanthoriatempthree_bf %>%
  left_join(otu_annot_df_mean_xanthoriatempthree_bf) %>%
  write.table('mycobiont_photobiont_xanthoriaparietinatempthree_bf_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

