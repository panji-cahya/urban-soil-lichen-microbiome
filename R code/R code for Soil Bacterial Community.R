### R Analysis of Soil Bacterial Community 
#For the paper published in npj Biofilms and Microbiomes
#Unveiling the ecological processes driving soil and lichen microbiome assembly along an urbanization gradient
#Authored by: Panji Cahya Mawarda, Rens van der Kaaij, Francisco Dini-Andreote, Deniz Duijker, Michael Stech, Arjen Speksnijder

#Created by: Panji Cahya Mawarda 04-04-2024

#call necessary library packages

#####################Start here##########################
biomfile=biomformat::read_biom("kraken 16S soil.biom2")

taxonomy_biom=observation_metadata(biomfile)
taxonomy_biom

write.table(taxonomy_biom, 'taxonomy_16S_Soil_New.csv')

#Taxonomy table for the tree----------
taxonomy_tree = as.data.frame(taxonomy_biom) %>%
  mutate_all(list(~str_replace_all(., regex("\\W+"), "_")))

taxonomy_tree["taxid"] = rownames(taxonomy_tree)

names (taxonomy_tree) <- NULL

write.table(taxonomy_tree, sep=",", quote=FALSE, row.names= FALSE, fileEncoding="UTF-8", 'taxonomy_16S_soil_tree.csv')

#perform the function below in powershell windows, inside the directory where you have both taxonomy_to_tree.pl and taxonomy_ITS_lichen_tree.csv
#perl .\taxonomy_to_tree.pl -b 60 -s "," -f .\taxonomy_ITS_Lichen_tree.csv > tree_ITS_lichen
#afterwards you have to open your file in notepad, and save as UTF-08

otu_biom=import_biom("kraken 16S soil.biom2")
write.csv(x = otu_table(otu_biom), file = "otu_table_16S_Soil_new.csv")

#reading files, you need four inputs: 
#otu table, taxonomy table, metadata

otu_table=read.csv("otu_table_16S_Soil_new.csv", row.names=1, header = T)
head(otu_table)
otu_table=as.matrix(otu_table)

taxonomy=read.csv("taxonomy_16S_Soil_New.csv", row.names=1, header = T)
replace(taxonomy, taxonomy == " ", NA) -> tax1
taxonomy1=as.matrix(taxonomy)

metadata=read.csv("metadata_16S_soil_new.csv", row.names=1, header=T)
metadata

OTU=otu_table(otu_table, taxa_are_rows=TRUE)

TAX=tax_table(taxonomy1)
head(TAX)

META=sample_data(metadata)
metadata

tre_fpath = "tree_soil.nwk"
tree <- ape::read.tree(tre_fpath)

#make sure files have the same sample names
sample_names(OTU)
sample_names(META)

#make the phyloseq object for downstream analysis
physeq=phyloseq(OTU, TAX, META)
physeq

sample_data(physeq)
sample_sums(physeq)
colnames(tax_table(physeq))
dim(tax_table(physeq))
tax_table(physeq)[1:5, 1:7]
sample_variables(physeq)
otu_table(physeq)[1:5, 1:54]
phy_tree(physeq)

############ Filter OTU table ############

#Remove Archea Mitochondria, Chloroplast, and unnecessary taxons 
physeqbac <- physeq %>% subset_taxa(Kingdom == "Bacteria" &
                                      Family != "Mitochondria" &
                                      Class != "Chloroplast" &
                                      Kingdom != "Archaea" &
                                      Kingdom != "Unassigned")

get_taxa_unique(physeqbac,"Kingdom")
get_taxa_unique(physeqbac,"Family")
get_taxa_unique(physeqbac,"Class")

colnames(tax_table(physeq))
otu_table(physeqbac)[1:5, 1:54]
tax_table(physeqbac)[1:5, 1:7]
physeqbac

tax_table(physeqbac)[1:5, 1:7]
any(taxa_sums(physeqbac) < 1) #false
any(is.na(otu_table(physeqbac))) #false

###discard taxons with phylum Unassigned (optional)

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(physeqbac))
sample_sum_df


# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  ylab("Number of samples") +
  theme(axis.title.y = element_blank())


##################################      Rarefaction         ######################################

# rarefy your dataset to lowest sequencing length, 711 should be kept constant:
phy_rare <- rarefy_even_depth(physeqbac, sample.size = 27194, rngseed = 711, replace = F, trimOTUs = TRUE, 	verbose = F)
phy_rare

ggrare(phy_rare, step = 100, color = "uhi_range", se = F)
sample_names(phy_rare)
sample_sums(phy_rare)

################################### Alpha diversity ############################################

estimate_richness(phy_rare, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
pAlpha = plot_richness(phy_rare,
                       shape = "uhi_range",
                       color = "temperature_zone",
                       measures = c("Observed", "Simpson", "Shannon"),
                       title = "Alpha Diversity")
pAlpha + geom_point(size = 3)

dalltreatments<-pAlpha$data

dalltreatments


################ Phyrare sub sample Based on Temperature/UHI Zone #######################

#temperaturezone1= low urbanized
#temperaturezone2= medium urbanized
#temperaturezone3= high urbanized

phyraresubtempone <- subset_samples(phy_rare,!temperature_zone %in% c("temperaturezone2", "temperaturezone3"))

phyraresubtempone

phyraresubtemptwo <- subset_samples(phy_rare,!temperature_zone %in% c("temperaturezone1", "temperaturezone3"))

phyraresubtemptwo

phyraresubtempthree <- subset_samples(phy_rare,!temperature_zone %in% c("temperaturezone1", "temperaturezone2"))

phyraresubtempthree


##function to calculate mean and standard error
####Determine the summarySE function###########################################################
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N = length2(xx[[col]], na.rm=na.rm),
                     mean = mean (xx[[col]], na.rm=na.rm),
                     sd = sd(xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N) # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#################### ALPHA DIVERSITY ##################
########### Species Richness #####################

alphadte = data.table(dalltreatments)
alphadte

# Subset to just the richness index
alphadte.richness <- alphadte[(variable == "Observed")]
alphadte.richness
str(alphadte.richness)

#Richness all--------------------------------------------------------------

#normal distribution test
shapiro.test(alphadte.richness$value)

#homogeneity of variance test
bartlett.test(value ~ temperature_zone, alphadte.richness)

#ANOVA Richness Based on Temperature Zone------------------
res.aov.richness.tempzone <- aov(value ~ uhi_range, data = alphadte.richness)
summary(res.aov.richness.tempzone)

TukeyHSD(res.aov.richness.tempzone)

#Kruskal Richness All Based on Temperature Zone -------------------
# IF the data is not normally distributed and the variance is not homogen, do Kruskal Wallis test
kruskal.test(value ~ uhi_range, data = alphadte.richness)

#Posthoc analysis with Wilcox test
pairwise.wilcox.test(alphadte.richness$value, alphadte.richness$temperature_zone,
                     p.adjust.method = "BH")


library(FSA)
dunnTest(value ~ temperature_zone,
         data=alphadte.richness,
         method="bh")

richness_all_tempzone <- ggplot(data = alphadte.richness,
                                aes(x=temperature_zone, y = value, color=uhi_range))+
  geom_boxplot(lwd=0.1)+
  geom_jitter(position=position_jitter(0.1), cex=3, alpha=0.7) +
  expand_limits(y=1600:4000) +
  scale_y_continuous(breaks=seq(1600,4000, by=400)) +
  labs(x="UHI range", y="Observed ASVs", title="ASV Richness")+
  theme_bw()+
  scale_color_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

richness_all_tempzone

richness_all <- ggarrange(richness_all_tempzone, richness_all_tempzone, richness_all_tempzone, richness_all_tempzone,
                          align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                          common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

richness_all

########### Shannon Diversity Index #####################

alphadte = data.table(dalltreatments)
alphadte

# Subset to just the richness index
alphadte.shannon <- alphadte[(variable == "Shannon")]
alphadte.shannon
str(alphadte.shannon)

shannon_all_tempzone <- ggplot(data = alphadte.shannon,
                               aes(x=temperature_zone, y = value, color=uhi_range))+
  geom_boxplot(lwd=0.1)+
  geom_jitter(position=position_jitter(0.1), cex=3, alpha=0.7) +
  expand_limits(y=2:8) +
  scale_y_continuous(breaks=seq(2,8, by=1)) +
  labs(x="UHI range", y="Shannon Diversity Index", title="Shannon")+
  theme_bw()+
  scale_color_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

shannon_all_tempzone

shannon_all <- ggarrange(shannon_all_tempzone, shannon_all_tempzone, shannon_all_tempzone, shannon_all_tempzone,
                         align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                         common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

shannon_all

#Shannon all--------------------------------------------------------------
#normal distribution test
shapiro.test(alphadte.shannon$value)

#homogeneity of variance test
bartlett.test(value ~ temperature_zone, alphadte.shannon)

#ANOVA shannon Based on Temperature Zone------------------
res.aov.shannon.tempzone <- aov(value ~ uhi_range, data = alphadte.shannon)
summary(res.aov.shannon.tempzone)

TukeyHSD(res.aov.shannon.tempzone)

#-----------------------BETA DIVERSITY --------------------------------

#Exporting OTU Table --------------------------------------------------
write.csv(x = otu_table(phy_rare), file = "phyraresuball_asv_16S_soil.csv")

# Based temperature zone----------------------------------------------------------------------------
#all
comall_tempzone <- read.csv("phyraresuball_asv_16S_soil.csv", header = TRUE, row.names = 1)
comall_tempzone <- t(comall_tempzone)
comall_tempzone [1:5, 1:3]

#all
group_info_all <- data.frame(row.names=rownames(comall_tempzone), t(as.data.frame(strsplit(rownames(comall_tempzone),"_"))))
head(group_info_all)

#Bray Curtis  Distance all-------------
dist_comall <- vegdist(comall_tempzone, method="bray", binary=FALSE,diag=1)
re_comall <- pcoa(dist_comall, correction="none",rn=NULL)
str(re_comall)

df_bray_comall <- data.frame(x=re_comall$vectors[,1], y=re_comall$vectors[,2],
                             tempzone=as.factor(group_info_all[,1]),
                             treespecies=as.factor(group_info_all[,2]),
                             replicates=as.factor(group_info_all[,3]))

str(df_bray_comall)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

pcoa_bray_all <- ggplot(df_bray_comall, aes(x,y, shape=tempzone, fill=tempzone)) +
  geom_point(size=5, alpha=0.7)+
  labs(x=paste("PCoA1: ", round(re_comall$values$Relative_eig[1] * 100, 2), "%", sep=""),
       y=paste("PCoA2: ", round(re_comall$values$Relative_eig[2] * 100, 2), "%", sep=""),
       title="Bray Curtis") +
  scale_fill_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme

pcoa_bray_all

#Average PCoA all------------

dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data1_all <- summaryBy(x~tempzone, data=df_bray_comall, FUN=dstats)
data2_all <- summaryBy(y~tempzone, data=df_bray_comall, FUN=dstats)

data_all <- cbind(data1_all[,3:5], data2_all)
str(data_all)
str(df_bray_comall)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

av_pcoa_bray_all_tempzone <- ggplot(data_all, aes(x=x.mean, y=y.mean, shape=tempzone, fill=tempzone)) +
  geom_point(size=5,alpha=0.7) +
  geom_errorbar(aes(ymin=y.mean-y.se, ymax=y.mean+y.se)) +
  geom_errorbar(aes(xmin=x.mean-x.se, xmax=x.mean+x.se)) +
  labs(x=paste("PCoA: ", round(re_comall$values$Relative_eig[1]*100, 2), "%", sep=""),
       y=paste("PCoA: ", round(re_comall$values$Relative_eig[2]*100, 2), "%", sep=""),
       title= "Bray Curtis") +
  scale_fill_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme
av_pcoa_bray_all_tempzone

bray_all <- ggarrange(av_pcoa_bray_all_tempzone, av_pcoa_bray_all_tempzone, av_pcoa_bray_all_tempzone, av_pcoa_bray_all_tempzone,
                      align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                      common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

bray_all

#Permanova Bray all------------------------
str(comall_tempzone)
comall_tempzone[1:6, 1:2]

df_all <- data.frame(row.names=rownames(comall_tempzone), t(as.data.frame(strsplit(rownames(comall_tempzone),"_"))))

df_all <- plyr::rename(df_all, replace = c("X1"="tempzone", "X2"="treespecies", "X3"="replicates"))
head(df_all)

set.seed(123)
Permanova_bray_all <- adonis2(comall_tempzone ~ tempzone, data=df_all, method="bray", permutation=9999)
Permanova_bray_all

Pairwise_Permanova_all <- pairwise.adonis(comall_tempzone, df_all$tempzone, sim.function = "vegdist", sim.method = "bray", perm = 999)
Pairwise_Permanova_all

############################## CO-OCCURENCE NETWORK #####################################

#--------- Input data for tempzone one ----
temperaturezone_1 <- subset_samples(phy_rare, temperature_zone == 'temperaturezone1')
temperaturezone_1

temperaturezone_1_df <- as.data.frame(otu_table(temperaturezone_1))
temperaturezone_1_df_filt <- temperaturezone_1_df[rowSums(temperaturezone_1_df) > 20 & rowSums(temperaturezone_1_df > 0) >= 3, ] #minimal abundance 20 ditemukan di 3 sampel, rule of thumb nya ambil ASV yang banyak ditemukan di banyak sampel untuk mempermudah perhitungan dan visualisasi
temperaturezone_1_df_filt <- temperaturezone_1_df_filt %>% rownames_to_column('#OTU ID')
write.table(temperaturezone_1_df_filt, 'temperaturezone1.tsv', quote = F, sep = '\t', col.names = T, row.names = F) 

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table temperaturezone1.tsv --correlation temperaturezone1_cor.tsv --covariance temperaturezone1_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table temperaturezone1.tsv --number 1000 --prefix bootstrap_counts/temperaturezone1 -t 4
# mkdir bootstrap_correlation
# parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
# fastspar_pvalues --otu_table temperaturezone1.tsv --correlation temperaturezone1_cor.tsv --prefix bootstrap_correlation/cor_temperaturezone1_ --permutations 1000 --outfile temperaturezone1_pval.tsv
# fastspar_reduce -r temperaturezone1_cor.tsv -p temperaturezone1_pval.tsv -o temperaturezone1
# rm -r bootstrap_counts
# rm -r bootstrap_correlation

fs_Rtempone_cor <- read.table('tempone_filtered_correlation.tsv', 
                              header = T, sep = '\t')
fs_Rtempone_p <- read.table('tempone_filtered_pvalue.tsv', 
                            header = T, sep = '\t')

fs_Rtempone_combine <- cbind.data.frame(fs_Rtempone_cor, fs_Rtempone_p[,3])
colnames(fs_Rtempone_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rtempone_filt <- fs_Rtempone_combine  %>% filter(abs(cor) > 0.8 & pval < 0.01)

# Export to Cytoscape ----
fs_Rtempone_filt %>%
  write.table('./temperaturezone1_net2Cys.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Open cytoscape

## Create annotation file ----
otu_df_tempone <- otu_table(temperaturezone_1 ) %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_tempone[,-1])

otu_annot_df_tempone <- otu_df_tempone %>%
  filter(OTU %in% unique(c(fs_Rtempone_filt $otu_1, fs_Rtempone_filt $otu2)))

otu_annot_df_tempone$mean_abun <- rowMeans(otu_annot_df_tempone[, -1])
otu_annot_df_mean_tempone <- otu_annot_df_tempone[,c(1,ncol(otu_annot_df_tempone))]

tax_df <- tax_table(phy_rare) %>% 
  data.frame() %>%
  rownames_to_column('OTU')

tax_annot_df_tempone <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtempone_filt $otu_1, fs_Rtempone_filt $otu2)))

tax_annot_df_tempone %>%
  left_join(otu_annot_df_mean_tempone) %>%
  write.table('temperaturezone1_net2Cys_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#--------- Input data for tempzone two ----
temperaturezone_2 <- subset_samples(phy_rare, temperature_zone == 'temperaturezone2')
temperaturezone_2

temperaturezone_2_df <- as.data.frame(otu_table(temperaturezone_2))
temperaturezone_2_df_filt <- temperaturezone_2_df[rowSums(temperaturezone_2_df) > 20 & rowSums(temperaturezone_2_df > 0) >= 3, ] #minimal abundance 20 ditemukan di 3 sampel, rule of thumb nya ambil ASV yang banyak ditemukan di banyak sampel untuk mempermudah perhitungan dan visualisasi
temperaturezone_2_df_filt <- temperaturezone_2_df_filt %>% rownames_to_column('#OTU ID')
write.table(temperaturezone_2_df_filt, 'temperaturezone2.tsv', quote = F, sep = '\t', col.names = T, row.names = F) 

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table temperaturezone2.tsv --correlation temperaturezone2_cor.tsv --covariance temperaturezone2_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table temperaturezone2.tsv --number 1000 --prefix bootstrap_counts/temperaturezone2 -t 4
# mkdir bootstrap_correlation
# parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
# fastspar_pvalues --otu_table temperaturezone2.tsv --correlation temperaturezone2_cor.tsv --prefix bootstrap_correlation/cor_temperaturezone2_ --permutations 1000 --outfile temperaturezone2_pval.tsv
# fastspar_reduce -r temperaturezone2_cor.tsv -p temperaturezone2_pval.tsv -o temperaturezone2
# rm -r bootstrap_counts
# rm -r bootstrap_correlation

fs_Rtemptwo_cor <- read.table('temptwo_filtered_correlation.tsv', 
                              header = T, sep = '\t')
fs_Rtemptwo_p <- read.table('temptwo_filtered_pvalue.tsv', 
                            header = T, sep = '\t')

fs_Rtemptwo_combine <- cbind.data.frame(fs_Rtemptwo_cor, fs_Rtemptwo_p[,3])
colnames(fs_Rtemptwo_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rtemptwo_filt <- fs_Rtemptwo_combine  %>% filter(abs(cor) > 0.8 & pval < 0.01)

# Export to Cytoscape ----
fs_Rtemptwo_filt %>%
  write.table('./temperaturezone2_net2Cys.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Open cytoscape

## Create annotation file ----
otu_df_temptwo <- otu_table(temperaturezone_2) %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_temptwo[,-1])

otu_annot_df_temptwo <- otu_df_temptwo %>%
  filter(OTU %in% unique(c(fs_Rtemptwo_filt $otu_1, fs_Rtemptwo_filt $otu2)))

otu_annot_df_temptwo$mean_abun <- rowMeans(otu_annot_df_temptwo[, -1])
otu_annot_df_mean_temptwo <- otu_annot_df_temptwo[,c(1,ncol(otu_annot_df_temptwo))]

tax_df <- tax_table(phy_rare) %>% 
  data.frame() %>%
  rownames_to_column('OTU')

tax_annot_df_temptwo <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtemptwo_filt $otu_1, fs_Rtemptwo_filt $otu2)))

tax_annot_df_temptwo %>%
  left_join(otu_annot_df_mean_temptwo) %>%
  write.table('temperaturezone2_net2Cys_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#--------- Input data for tempzone three ----
temperaturezone_3 <- subset_samples(phy_rare, temperature_zone == 'temperaturezone3')
temperaturezone_3

temperaturezone_3_df <- as.data.frame(otu_table(temperaturezone_3))
temperaturezone_3_df_filt <- temperaturezone_3_df[rowSums(temperaturezone_3_df) > 20 & rowSums(temperaturezone_3_df > 0) >= 3, ] #minimal abundance 20 ditemukan di 3 sampel, rule of thumb nya ambil ASV yang banyak ditemukan di banyak sampel untuk mempermudah perhitungan dan visualisasi
temperaturezone_3_df_filt <- temperaturezone_3_df_filt %>% rownames_to_column('#OTU ID')
write.table(temperaturezone_3_df_filt, 'temperaturezone3.tsv', quote = F, sep = '\t', col.names = T, row.names = F) 

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table temperaturezone3.tsv --correlation temperaturezone3_cor.tsv --covariance temperaturezone3_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table temperaturezone3.tsv --number 1000 --prefix bootstrap_counts/temperaturezone3 -t 4
# mkdir bootstrap_correlation
# parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
# fastspar_pvalues --otu_table temperaturezone3.tsv --correlation temperaturezone3_cor.tsv --prefix bootstrap_correlation/cor_temperaturezone3_ --permutations 1000 --outfile temperaturezone3_pval.tsv
# fastspar_reduce -r temperaturezone3_cor.tsv -p temperaturezone3_pval.tsv -o temperaturezone3
# rm -r bootstrap_counts
# rm -r bootstrap_correlation

fs_Rtempthree_cor <- read.table('tempthree_filtered_correlation.tsv', 
                                header = T, sep = '\t')
fs_Rtempthree_p <- read.table('tempthree_filtered_pvalue.tsv', 
                              header = T, sep = '\t')

fs_Rtempthree_combine <- cbind.data.frame(fs_Rtempthree_cor, fs_Rtempthree_p[,3])
colnames(fs_Rtempthree_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rtempthree_filt <- fs_Rtempthree_combine  %>% filter(abs(cor) > 0.8 & pval < 0.01)

# Export to Cytoscape ----
fs_Rtempthree_filt %>%
  write.table('./temperaturezone3_net2Cys.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#SOpen cytoscape

## Create annotation file ----
otu_df_tempthree <- otu_table(temperaturezone_3) %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_tempthree[,-1])

otu_annot_df_tempthree <- otu_df_tempthree %>%
  filter(OTU %in% unique(c(fs_Rtempthree_filt $otu_1, fs_Rtempthree_filt $otu2)))

otu_annot_df_tempthree$mean_abun <- rowMeans(otu_annot_df_tempthree[, -1])
otu_annot_df_mean_tempthree <- otu_annot_df_tempthree[,c(1,ncol(otu_annot_df_tempthree))]

tax_df <- tax_table(phy_rare) %>% 
  data.frame() %>%
  rownames_to_column('OTU')

tax_annot_df_tempthree <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rtempthree_filt $otu_1, fs_Rtempthree_filt $otu2)))

tax_annot_df_tempthree %>%
  left_join(otu_annot_df_mean_tempthree) %>%
  write.table('temperaturezone3_net2Cys_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Beta NTI------------------------
# Load packages ----
library(NST)
library(iCAMP)
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(phyloseq)
library(tibble)
library(vegan)

# 1. File and folder path and parameters ----

wd="."
save.wd="output/"
nworker=8 
memory_use=13
rand.time=9999 
prefix="Test"

# 2. Load data and match IDs ----

physeq_bnti=phyloseq(OTU, TAX, META, tree)
physeq_bnti

physeqbac_bnti <- physeq_bnti %>% subset_taxa(Kingdom == "Bacteria" &
                                                Family != "Mitochondria" &
                                                Class != "Chloroplast" &
                                                Kingdom != "Archaea" &
                                                Kingdom != "Unassigned")

get_taxa_unique(physeqbac_bnti,"Kingdom")
get_taxa_unique(physeqbac_bnti,"Family")
get_taxa_unique(physeqbac_bnti,"Class")

colnames(tax_table(physeqbac_bnti))
otu_table(physeqbac_bnti)[1:5, 1:54]
tax_table(physeqbac_bnti)[1:5, 1:7]
physeqbac_bnti

tax_table(physeqbac_bnti)[1:5, 1:7]
any(taxa_sums(physeqbac_bnti) < 1) #false
any(is.na(otu_table(physeqbac_bnti))) #false

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(physeqbac_bnti))
sample_sum_df

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  ylab("Number of samples") +
  theme(axis.title.y = element_blank())

# rarefy your dataset to lowest sequencing length, 711 should be kept constant:
phy_rare_bnti <- rarefy_even_depth(physeqbac_bnti, sample.size = 19999, rngseed = 711, replace = F, trimOTUs = TRUE, 	verbose = F)
phy_rare_bnti

ggrare(phy_rare_bnti, step = 100, color = "uhi_range", se = F)
sample_names(phy_rare_bnti)
sample_sums(phy_rare_bnti)

saveRDS(phy_rare_bnti, "phy_rare_bnti.rds")

otu_table_soil_filt <- otu_table(phy_rare_bnti)[rowSums (otu_table(phy_rare_bnti)) > 20 & rowSums(otu_table(phy_rare_bnti) > 0) >= 3, ]

comm=t(data.frame(otu_table_soil_filt))
group=data.frame(sample_data(phy_rare_bnti))
tree_bnti=phy_tree(phy_rare_bnti)

samp.ck=NST::match.name(rn.list=list(comm=comm,group=group))
comm=samp.ck$comm
comm=comm[,colSums(comm)>0,drop=FALSE]
group=samp.ck$group 

tax.ck=NST::match.name(cn.list = list(comm=comm), tree.list = list(tree=tree_bnti))
comm=tax.ck$comm
tre=tax.ck$tree_bnti

# 3. Calculate bNTI and RCI ----
## 3.1. Using iCAMP - Method 1 ----

set.seed(1234)
qpen <- qpen(comm = comm,
             pd = "./output/pdbig/path.rda", 
             pd.big.wd = NULL,
             pd.big.spname = NULL, 
             tree = tree_bnti,
             bNTI = NULL,
             RC = NULL,
             ab.weight = TRUE,
             meta.ab = NULL,
             exclude.conspecifics = FALSE,
             rand.time = rand.time, 
             sig.bNTI = 2, 
             sig.rc = 0.95,
             nworker = nworker, 
             memory.G = memory_use, 
             project = NA, 
             wd = save.wd,
             output.detail = FALSE, 
             save.bNTIRC = FALSE,
             taxo.metric = "bray",
             transform.method = NULL,
             logbase = 2,
             dirichlet = FALSE)

qpen <- readRDS('qpen.rds')
qpen

### Combine the results with group information ----
rc_data <- qpen$result %>%
  left_join(group %>% rownames_to_column('sample1'), 'sample1') %>% 
  left_join(group %>% rownames_to_column('sample2'), 'sample2') %>%
  filter(temperature_zone.x == temperature_zone.y & tree_species.x == tree_species.y)

plot_betanti_lichen <- rc_data %>%
  ggplot(aes(x = temperature_zone.x, y = bNTI, color = temperature_zone.x)) +
  geom_boxplot() +
  #facet_grid(~lichen_species.x, scales = 'free', space = 'free') +
  xlab('Temperature Zone') + ylab("betanti")+
  scale_color_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  scale_y_continuous(breaks=-10:10) +
  geom_hline(linewidth=1.5, yintercept=-2, linetype="dashed", color = "black")+
  geom_hline(linewidth=1.5, yintercept=2, linetype="dashed", color = "black")+
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_betanti_lichen

plot_rc_lichen <- rc_data %>%
  ggplot(aes(x = temperature_zone.x, y = RC, color = temperature_zone.x)) +
  geom_boxplot() +
  #geom_bar(aes(x = temperature_zone.x, y = RC, color = temperature_zone.x)) +
  #facet_grid(~lichen_species.x, scales = 'free', space = 'free') +
  xlab('Temperature Zone') + ylab("RC Bray")+
  scale_color_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  #scale_y_continuous(breaks=-0.95:0.95) +
  scale_y_continuous(breaks=seq(-1,1.35, by=0.35)) +
  geom_hline(linewidth=1.5, yintercept=-0.95, linetype="dashed", color = "black")+
  geom_hline(linewidth=1.5, yintercept=0.95, linetype="dashed", color = "black")+
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_rc_lichen

plot_betanti_tempzone <- rc_data %>%
  ggplot(aes(x = lichen_species.x, y = bNTI, color = lichen_species.x)) +
  geom_boxplot() +
  facet_grid(~temperature_zone.x, scales = 'free', space = 'free') +
  xlab('Day') + ylab("betanti")+
  scale_color_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  scale_y_continuous(breaks=-10:6) +
  geom_hline(linewidth=1.5, yintercept=-2, linetype="dashed", color = "black")+
  geom_hline(linewidth=1.5, yintercept=2, linetype="dashed", color = "black")+
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_betanti_tempzone

plot_rc_tempzone <- rc_data %>%
  ggplot(aes(x = lichen_species.x, y = RC, color = lichen_species.x)) +
  geom_boxplot() +
  facet_grid(~temperature_zone.x, scales = 'free', space = 'free') +
  xlab('Temperature Zone') + ylab("RC Bray")+
  scale_color_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  #scale_y_continuous(breaks=-0.95:0.95) +
  scale_y_continuous(breaks=seq(-1,1.35, by=0.35)) +
  geom_hline(linewidth=1.5, yintercept=-0.95, linetype="dashed", color = "black")+
  geom_hline(linewidth=1.5, yintercept=0.95, linetype="dashed", color = "black")+
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_rc_tempzone

plot_betanti_facetgrid <- ggarrange(plot_betanti_lichen, plot_rc_lichen, plot_betanti_lichen, plot_rc_lichen,
                                    align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                                    common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

plot_betanti_facetgrid

### Calculate relative importance (%) ----
rc_data_summary <- rc_data %>%
  group_by(process, temperature_zone.x) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(temperature_zone.x) %>%
  mutate(pct = count / sum(count) * 100) %>%
  ungroup()

### Plot the output
#rc_data_summary$Day.x <- factor(rc_data_summary$Day.x, 
#                                levels = c("Day 0","Day 35","Day 93","Day 103","Day 159")) #dibikin faktor biar berurutan

process_colors <- c("Homogeneous.Selection"="#DBFFD6",
                    "Heterogeneous.Selection"="#ECD4FF",
                    "Dispersal.Limitation"="orange",
                    "Undominated"="#ACE7FF",
                    "Homogenizing.Dispersal"="#FFABAB")

plot_process_tempzone <- rc_data_summary %>%
  ggplot(aes(x = temperature_zone.x, y = pct, fill = process)) +
  geom_bar(stat = 'identity') +
  #facet_grid(~lichen_species.x, scales = 'free', space = 'free') +
  xlab('Day') + ylab("Relative importance (%)") +
  scale_fill_manual(values = process_colors) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_process_tempzone

plot_process_methodone <- ggarrange(plot_process_tempzone, plot_process_tempzone, plot_process_tempzone, plot_process_tempzone,
                                    align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                                    common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

plot_process_methodone

### Combine the results with group information ----

rc_data_2 <- qpen_cm$result %>%
  left_join(group %>% rownames_to_column('sample1'), 'sample1') %>% 
  left_join(group %>% rownames_to_column('sample2'), 'sample2') %>%
  filter(temperature_zone.x == temperature_zone.y & lichen_species.x == lichen_species.y)

plot_betanti_lichen_2 <- rc_data_2 %>%
  ggplot(aes(x = temperature_zone.x, y = bNTI, color = temperature_zone.x)) +
  geom_boxplot() +
  facet_grid(~lichen_species.x, scales = 'free', space = 'free') +
  xlab('Day') + ylab("betanti")+
  scale_color_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  scale_y_continuous(breaks=-10:6) +
  geom_hline(linewidth=1.5, yintercept=-2, linetype="dashed", color = "black")+
  geom_hline(linewidth=1.5, yintercept=2, linetype="dashed", color = "black")+
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_betanti_lichen_2

plot_rc_lichen_2 <- rc_data_2 %>%
  ggplot(aes(x = temperature_zone.x, y = RC, color = temperature_zone.x)) +
  geom_boxplot() +
  facet_grid(~lichen_species.x, scales = 'free', space = 'free') +
  xlab('Temperature Zone') + ylab("RC Bray")+
  scale_color_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  #scale_y_continuous(breaks=-0.95:0.95) +
  scale_y_continuous(breaks=seq(-1,1.35, by=0.35)) +
  geom_hline(linewidth=1.5, yintercept=-0.95, linetype="dashed", color = "black")+
  geom_hline(linewidth=1.5, yintercept=0.95, linetype="dashed", color = "black")+
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_rc_lichen

plot_betanti_tempzone_2 <- rc_data_2 %>%
  ggplot(aes(x = lichen_species.x, y = bNTI, color = lichen_species.x)) +
  geom_boxplot() +
  facet_grid(~temperature_zone.x, scales = 'free', space = 'free') +
  xlab('Day') + ylab("betanti")+
  scale_color_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  scale_y_continuous(breaks=-10:6) +
  geom_hline(linewidth=1.5, yintercept=-2, linetype="dashed", color = "black")+
  geom_hline(linewidth=1.5, yintercept=2, linetype="dashed", color = "black")+
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_betanti_tempzone

plot_rc_tempzone_2 <- rc_data_2 %>%
  ggplot(aes(x = lichen_species.x, y = RC, color = lichen_species.x)) +
  geom_boxplot() +
  facet_grid(~temperature_zone.x, scales = 'free', space = 'free') +
  xlab('Temperature Zone') + ylab("RC Bray")+
  scale_color_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  #scale_y_continuous(breaks=-0.95:0.95) +
  scale_y_continuous(breaks=seq(-1,1.35, by=0.35)) +
  geom_hline(linewidth=1.5, yintercept=-0.95, linetype="dashed", color = "black")+
  geom_hline(linewidth=1.5, yintercept=0.95, linetype="dashed", color = "black")+
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_rc_tempzone_2

plot_betanti_facetgrid_2 <- ggarrange(plot_betanti_lichen_2, plot_betanti_tempzone_2, plot_rc_lichen_2, plot_rc_tempzone_2,
                                      align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                                      common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

plot_betanti_facetgrid_2

plot_betanti_facetgrid_2 <- ggarrange(plot_rc_lichen_2, plot_rc_tempzone, plot_rc_lichen, plot_rc_tempzone,
                                      align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                                      common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

plot_betanti_facetgrid_2

### Calculate relative importance (%) ----
rc_data_summary_2 <- rc_data_2 %>%
  group_by(process, temperature_zone.x, lichen_species.x) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(temperature_zone.x, lichen_species.x) %>%
  mutate(pct = count / sum(count) * 100) %>%
  ungroup()

### Plot the output ----
rc_data_summary_2$Day.x <- factor(rc_data_summary_2$Day.x, 
                                  levels = c("Day 0","Day 35","Day 93","Day 103","Day 159"))

plot_process_2_tempzone <- rc_data_summary_2 %>%
  ggplot(aes(x = temperature_zone.x, y = pct, fill = process)) +
  geom_bar(stat = 'identity') +
  facet_grid(~lichen_species.x, scales = 'free', space = 'free') +
  xlab('Temperature Zone') + ylab("Relative importance (%)") +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


plot_process_2_tempzone

plot_process_2_lichen <- rc_data_summary_2 %>%
  ggplot(aes(x = lichen_species.x, y = pct, fill = process)) +
  geom_bar(stat = 'identity') +
  facet_grid(~temperature_zone.x, scales = 'free', space = 'free') +
  xlab('Lichen Species') + ylab("Relative importance (%)") +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


plot_process_2_lichen

plot_process_facetgrid_2 <- ggarrange(plot_process_2_tempzone, plot_process_2_lichen, plot_process_2_tempzone, plot_process_2_lichen,
                                      align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                                      common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

plot_process_facetgrid_2

# Note ----
# bNTI < -2 = Homogeneous selection
# bNTI > +2 = Variable selection
# abs(bNTI) < 2 & RC < -0.95 = Homogenizing dispersal
# abs(bNTI) < 2 & RC < +0.95 = Undominated
# abs(bNTI) < 2 & RC > +0.95 = Dispersal limitation 

#---------------------------------Differential Abundance Across Urbanization Zone--------------------------------

#Converting Absolute Abundance to Relative Abundance

check_matrix = function(given_matrix, message = NULL) {
  
  if (!is.matrix(given_matrix) & !is(given_matrix, "sparseMatrix")) {
    stop("Provided ", message, " matrix is not a matrix")
  }
}

make_relative = function(abund_matrix) {
  
  # Check input type using internal function
  check_matrix(abund_matrix, "abundance")
  
  # Compute relative abundances matrix
  if (requireNamespace("Matrix", quietly = TRUE) &
      is(abund_matrix, "sparseMatrix")) {
    
    sites_abund = Matrix::rowSums(abund_matrix, na.rm = TRUE)
    
    rel_abund_matrix = abund_matrix / sites_abund
  } else {
    # Compute total site abundances
    sites_abund = rowSums(abund_matrix, na.rm = TRUE)
    
    # Divide each individual abundace by total site abundance
    rel_abund_matrix = sweep(abund_matrix, 1, sites_abund, "/")
  }
  
  return(rel_abund_matrix)
}

#------------------------------all-----------------------------------------

write.csv(x = otu_table(phy_rare), file = "phyrare_heatmap.csv")

#after you download, change the name of the column from code to temperature zone_tree_replicate

otu_table_all=read.csv("phyrare_heatmap.csv", row.names=1, header = 1, check.names = F)
head(otu_table_all)

otu_table_all <- otu_table_all[rowSums(otu_table_all) > 20 & rowSums(otu_table_all > 0) >= 10, ] #minimal abundance 20 ditemukan di 3 sampel, rule of thumb nya ambil ASV yang banyak ditemukan di banyak sampel untuk mempermudah perhitungan dan visualisasi


otu_table_all=as.matrix(otu_table_all)

otu_table_all <- make_relative(otu_table_all) 

otu_table_all <- na.omit(otu_table_all)

otu_table_all <- t(otu_table_all)

write.csv(x = (otu_table_all), file = "subset_significant_all.csv")

#Change OTUID in Excel Retain the name and match it later
#Make a new csv file in excel "subset_significant_all.csv" by transposing otu_table_all.csv
#Change Lalamusa_all_A to just Lalamusa
#Put "group" as a coloumn name of Day

#Subset Significant OTU Over Time
df_successionall_subset <- read.csv("subset_significant_all.csv", header=1, check.names = F)

for (i in 2:ncol(df_successionall_subset)) {
  res.aov.succession.all_subset <- aov(df_successionall_subset[,i] ~ df_successionall_subset$Group, data = df_successionall_subset)
  bb <- data.frame(summary(res.aov.succession.all_subset)[[1]])
  if (bb$Pr..F.[1] < 0.01) {
    print(colnames(df_successionall_subset)[i])
  }
}

df_successionall_significant <- df_successionall_subset %>% select(Group, OTU1,OTU2,OTU3,OTU4,OTU5,OTU6,OTU10,OTU11,OTU25,OTU39,OTU40,OTU41,OTU42,OTU45,OTU49,OTU50,OTU52,OTU59,OTU60,OTU63,OTU77,OTU106,OTU107,OTU109,OTU113,OTU116,OTU118,OTU121,OTU123,OTU127,OTU134,OTU148,OTU150,OTU162,OTU163,OTU168,OTU169,OTU172,OTU180,OTU216,OTU218,OTU227,OTU228,OTU230,OTU233,OTU243,OTU257,OTU258,OTU269,OTU270,OTU277,OTU279,OTU280,OTU306,OTU344,OTU351,OTU364,OTU365,OTU367,OTU368,OTU369,OTU370,OTU372,OTU373,OTU377,OTU378,OTU379,OTU391,OTU394,OTU405,OTU408,OTU409,OTU411,OTU420,OTU437,OTU442,OTU482,OTU555,OTU574,OTU587,OTU606,OTU610,OTU625,OTU662,OTU677,OTU716,OTU729,OTU737,OTU738,OTU740,OTU741,OTU746,OTU748,OTU752,OTU757,OTU759,OTU762,OTU771,OTU777,OTU818,OTU841,OTU851,OTU870,OTU881,OTU907,OTU908,OTU912,OTU914,OTU918,OTU919,OTU922,OTU924,OTU926,OTU927,OTU931,OTU932,OTU938,OTU939,OTU982,OTU985,OTU991,OTU994,OTU995,OTU996,OTU1002,OTU1006,OTU1009,OTU1013,OTU1019,OTU1020,OTU1021,OTU1035,OTU1049,OTU1052,OTU1056,OTU1064,OTU1081,OTU1106,OTU1111,OTU1119,OTU1136,OTU1147,OTU1148,OTU1161,OTU1177,OTU1184,OTU1190,OTU1196,OTU1198,OTU1206,OTU1244,OTU1248,OTU1250,OTU1251,OTU1254,OTU1256,OTU1280,OTU1285,OTU1302,OTU1310,OTU1312,OTU1315,OTU1316,OTU1320,OTU1321,OTU1322,OTU1323,OTU1329,OTU1331,OTU1338,OTU1354,OTU1356,OTU1369,OTU1376,OTU1382,OTU1384,OTU1385,OTU1386,OTU1388,OTU1389,OTU1391,OTU1396,OTU1399,OTU1400,OTU1401,OTU1402,OTU1404,OTU1410,OTU1426,OTU1430,OTU1435,OTU1436,OTU1453,OTU1455,OTU1461,OTU1469,OTU1476,OTU1481,OTU1482,OTU1484,OTU1493,OTU1518,OTU1535,OTU1550,OTU1562,OTU1563,OTU1565,OTU1566,OTU1613,OTU1614,OTU1621,OTU1640,OTU1641,OTU1642,OTU1643,OTU1658,OTU1679,OTU1752,OTU1755,OTU1760,OTU1762,OTU1782,OTU1790,OTU1791,OTU1794,OTU1829,OTU1839,OTU1842,OTU1851,OTU1861,OTU1869,OTU1885,OTU1886,OTU1898,OTU1907,OTU1909,OTU1920,OTU1922,OTU1934,OTU1947,OTU1961,OTU1963,OTU1967,OTU1974,OTU1984,OTU2007,OTU2023,OTU2024,OTU2029,OTU2043,OTU2051,OTU2056,OTU2107,OTU2133,OTU2223,OTU2313,OTU2314,OTU2321,OTU2323,OTU2327,OTU2343,OTU2362,OTU2377,OTU2515,OTU2521,OTU2522,OTU2547,OTU2585,OTU2653,OTU2663,OTU2762,OTU2798,OTU2802,OTU2804,OTU2805,OTU2815,OTU2838,OTU2839,OTU2878,OTU2891,OTU2892,OTU2903,OTU2914,OTU2915,OTU2921,OTU2966,OTU2968,OTU2989,OTU2994,OTU3048,OTU3097,OTU3139,OTU3141,OTU3145,OTU3148,OTU3153,OTU3164,OTU3166,OTU3168,OTU3171,OTU3174,OTU3176,OTU3178,OTU3186,OTU3190,OTU3196,OTU3197,OTU3226,OTU3241,OTU3268,OTU3350,OTU3351,OTU3357,OTU3359,OTU3376,OTU3403,OTU3593,OTU3596,OTU3608,OTU3611,OTU3660,OTU3692,OTU3697,OTU3762,OTU3770,OTU3799,OTU3828,OTU3856,OTU3872,OTU3873,OTU3891,OTU3902,OTU3933,OTU3973,OTU3998,OTU3999,OTU4025,OTU4094,OTU4098,OTU4128,OTU4130,OTU4137,OTU4138,OTU4139,OTU4185,OTU4224,OTU4237,OTU4241,OTU4274,OTU4293)
write.csv(x = (df_successionall_significant), file = "df_successionall_significant.csv")

#Make "df_successionall_significant.csv" by transposing 

library(OTUtable)

df_successionall_significant <- read.csv("df_successionall_significant.csv", row.names=1, header=1, check.names = F)
df_successionall_filter <- filter_taxa(df_successionall_significant, abundance = 1, persistence=12.34)
df_successionall_filter=as.matrix(df_successionall_filter)
df_successionall_filter=t(df_successionall_filter)

df_successionall_original <- read.csv("phyrare_heatmap.csv", row.names=1, header=1, check.names = F)
df_successionall_original_f <- filter_taxa(df_successionall_original, abundance = 0.00000000000000001, persistence=0)
df_successionall_original = as.matrix(df_successionall_original_f)
df_successionall_original = t(df_successionall_original)

#Procrustes Analysis
dist_relativeabundance_all <- vegdist(df_successionall_original, method="bray", binary=FALSE, diag=1)
library(ape)
re_relativeabundance_all <- pcoa(dist_relativeabundance_all, correction="none", rn=NULL)

#df_successionall_filter = t(df_successionall_filter)
dist_filter_all <- vegdist(df_successionall_filter, method="bray", binary=FALSE, diag=1)
re_dist_filter_all <- pcoa(dist_filter_all, correction="none", rn=NULL)

vare.proc.all <- procrustes(re_relativeabundance_all$vectors, re_dist_filter_all$vectors, symmetric = TRUE) 
protest(X = re_relativeabundance_all$vectors, Y = re_dist_filter_all$vectors, permutations = 9999)

plot(vare.proc.all)
plot(vare.proc.all, kind=2)
residuals(vare.proc.all)


#Put the taxonomy Phylum on excel file
#Create another excel file from filter table and replace the OTUCode to OTUID

df_successionall_filter=t(df_successionall_filter)
write.csv(x = (df_successionall_filter), file = "df_successionall_filter.csv")

# Heatmap all----

heatmap_PM1_data <- read.csv("df_heatmap_reduced_next.csv",
                             header=T, row.names=1)
Heatmap_PM1 <- heatmap_PM1_data[,1:28] # data
anno_one <- heatmap_PM1_data[,29] # annotation
anno_two <- heatmap_PM1_data[,30] # annotation
anno_one <- factor(anno_one, levels = c("temperaturezone1","temperaturezone2","temperaturezone3"))
annotation_col = data.frame(
  temperature_zone = factor(anno_one, levels = c("temperaturezone1","temperaturezone2","temperaturezone3")),
  uhi_level = anno_two)
heatmap_all <- pheatmap::pheatmap(t(Heatmap_PM1), 
                                  cluster_row=F,     
                                  cluster_col=F,         
                                  fontsize_col=8,        
                                  fontsize_row=8,         
                                  show_rownames = T,      
                                  show_colnames = F,     
                                  cellwidth=8,
                                  cellheight=8,
                                  color = colorRampPalette(rev(c("#4169E1","#00BFFF","#F8F8FF")))(100),
                                  border_color='grey60',
                                  scale = 'none',
                                  main = '', 
                                  fontsize = 8, 
                                  annotation_col = annotation_col,
                                  gaps_row = c(3,6,9),
                                  treeheight_row = 100
)

heatmap_all

pdf('Heatmap_all.pdf', width = 35, height = 150)
grid::grid.newpage() 
grid::grid.draw(heatmap_deseq$gtable) 
dev.off() 

otu_table_heatmap=read.csv("df_successionall_filter.csv", row.names=1, header = T)
head(otu_table_heatmap)
otu_table_heatmap=as.matrix(otu_table_heatmap)

OTU_heatmap=otu_table(otu_table_heatmap, taxa_are_rows=TRUE)

TAX=tax_table(taxonomy1)
head(TAX)

META=sample_data(metadata)
metadata

#make sure files have the same sample names
sample_names(OTU_heatmap)
sample_names(META)

#make the phyloseq object for downstream analysis
physeq_heatmap=phyloseq(OTU_heatmap, TAX, META)
physeq_heatmap

data.phyla.all <- phy_rare %>%
  tax_glom(taxrank = "Species") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>% 
  filter(Abundance > 0.01) %>%
  arrange(Phylum) 
data.phyla.all

#data.phyla.all$Treatment1 <- factor(data.phyla.all$Treatment1, levels=c("Soil_Rice_24D", "Soil_Rice", "Soil_24D", "Soil"))

#data.phyla.all[data.phyla.all == "Unidentified"] <- "Zunidentified"

phylum_colors <- c("Acidobacteria" = "pink", "Actinobacteria" = "blue", "Armatimonadetes" = "beige",
                   "Bacteroidetes" = "purple", "Chlamydiae" = "grey", "Chloroflexi" = "red",
                   "Cyanobacteria" = "yellow", "Firmicutes" = "cyan", "Gemmatimonadetes" = "lavender",
                   "Nitrospirae" = "gold", "OD1" = "salmon", "Planctomycetes" = "brown", 
                   "Proteobacteria" = "green", "TM7"= "skyblue", "Verrucomicrobia" = "violet",
                   "Deinococcus-Thermus" = "black")

class_colors <- c("Acidobacteriia" = "pink", "Acidimicrobiia" = "blue", "Alphaproteobacteria" = "lavender",
                  "Betaproteobacteria" = "purple", "Blastocatellia" = "grey", "Chitinophagia" = "red",
                  "Cytophagia" = "yellow", "Deinococci" = "beige", "Unassigned" = "cyan",
                  "Flavobacteriia" = "skyblue", "Gammaproteobacteria" = "salmon", "Sphingobacteriia" = "brown", 
                  "Thermoleophilia" = "green", "TM7"= "skyblue", "Verrucomicrobia" = "violet",
                  "Deltaproteobacteria" = "black", "Actinomycetia" = "orange")

class_colors <- c("Acidobacteriia" = "pink", "Actinomycetia" = "blue", "Alphaproteobacteria" = "lavender",
                  "Betaproteobacteria" = "purple", "Blastocatellia" = "cyan", "Bacilli" = "red",
                  "Holophagae" = "yellow", "Nitrospira" = "beige", "Oligoflexia" = "grey",
                  "Rubrobacteria" = "skyblue", "Gammaproteobacteria" = "salmon", "Spartobacteria" = "brown", 
                  "Thermoleophilia" = "green", "Terrimicrobia"= "black", "Vicinamibacteria" = "violet",
                  "Deltaproteobacteria" = "black", "Acidimicrobiia" = "orange")

write.csv(x = (data.phyla.all), file = "data.phyla.all.csv")
data.phyla.all=read.csv("data.phyla.all.csv", row.names=1, header = T)

m=ggplot(data.phyla.all, aes(x = temperature_zone, y = Abundance, fill = Class)) +
  facet_grid(.~temperature_zone, scales="free_x") +
  geom_bar(stat = "identity",position="fill") + # to equal 1
  scale_fill_manual(values = class_colors) +
  theme(text = element_text(size=15),axis.text.y = element_text(hjust=1)) +
  theme(strip.text = element_text(face="bold", size=18)) +
  scale_x_discrete(drop=TRUE) +
  theme(axis.title.x = element_text(size=12)) + # Remove x axis title
  ylab("Relative Abundance") +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m  

plot_taxa <- ggarrange(m, m, m, m,align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                       common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

plot_taxa

m_treatment=ggplot(data.phyla.all, aes(x = temperature_zone, y = Abundance, fill = Phylum)) +
  facet_grid(.~lichen_species, scales="free_x") +
  geom_bar(stat = "identity",position="fill") + # to equal 1
  scale_fill_manual(values = phylum_colors) +
  theme(text = element_text(size=15),axis.text.y = element_text(hjust=1)) +
  theme(strip.text = element_text(face="bold", size=18)) +
  scale_x_discrete(drop=TRUE) +
  theme(axis.title.x = element_text(size=12)) + # Remove x axis title
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.5)) + 
  theme(legend.position = "right")

m_treatment

#TURNOVER BASED ON BRAY CURTIS ------

#TURNOVER BASED ON BRAY CURTIS ------

bray_comall_turnover <- as.matrix(dist_comall)
bray_comall_turnover <- data.frame(as.table(bray_comall_turnover))[lower.tri(bray_comall_turnover, diag = FALSE), ]
cat("should got:", (33*33-33)/2, "pair-wise distance\n"); 
cat("Actually we got:", length(bray_comall_turnover$Freq), "pair-wise distance\n")


# extract samples from the same group
group1_comall<-data.frame(row.names=rownames(bray_comall_turnover), t(as.data.frame(strsplit(as.character(bray_comall_turnover$Var1), "_"))))
head(group1_comall)

group2_comall<-data.frame(row.names=rownames(bray_comall_turnover),t(as.data.frame(strsplit(as.character(bray_comall_turnover$Var2), "_"))))
head(group2_comall)

df_turnover_comall <- data.frame(urban1=group1_comall[,1],
                                 Treatment1=group1_comall[,2],
                                 Replicate1=group1_comall[,3],
                                 urban2=group2_comall[,1],
                                 Treatment2=group2_comall[,2],
                                 Replicate2=group2_comall[,3],
                                 bray_comall_turnover)

row.names(df_turnover_comall) <- paste(df_turnover_comall$Var1, df_turnover_comall$Var2, sep="/")

head(df_turnover_comall)
str(df_turnover_comall)

df_turnover_comall$urban_urban <- paste(df_turnover_comall$urban1, df_turnover_comall$urban2, sep="_")
head(df_turnover_comall)
str(df_turnover_comall)

# subset urban-urban distance
dflowmed_comall <- df_turnover_comall[df_turnover_comall$urban_urban=='temperaturezone1_temperaturezone2' | df_turnover_comall$urban_urban=='temperaturezone1_temperaturezone2',]; dflowmed_comall$urbangroup <- 'low_med'
dflowhigh_comall <- df_turnover_comall[df_turnover_comall$urban_urban=='temperaturezone1_temperaturezone3' | df_turnover_comall$urban_urban=='temperaturezone1_temperaturezone3',]; dflowhigh_comall$urbangroup <- 'low_high'
dfmedhigh_comall <- df_turnover_comall[df_turnover_comall$urban_urban=='temperaturezone2_temperaturezone3' | df_turnover_comall$urban_urban=='temperaturezone2_temperaturezone3',]; dfmedhigh_comall$urbangroup <- 'med_high'

dflowlow_comall <- df_turnover_comall[df_turnover_comall$urban_urban=='temperaturezone1_temperaturezone1' | df_turnover_comall$urban_urban=='temperaturezone1_temperaturezone1',]; dflowlow_comall$urbangroup <- 'low_low'
dfmedmed_comall <- df_turnover_comall[df_turnover_comall$urban_urban=='temperaturezone2_temperaturezone2' | df_turnover_comall$urban_urban=='temperaturezone2_temperaturezone2',]; dfmedmed_comall$urbangroup <- 'med_med'
dfhighhigh_comall <- df_turnover_comall[df_turnover_comall$urban_urban=='temperaturezone3_temperaturezone3' | df_turnover_comall$urban_urban=='temperaturezone3_temperaturezone3',]; dfhighhigh_comall$urbangroup <- 'high_high'

# combine subsets
new_bray.urban_comall_as <- rbind(dflowlow_comall, dfmedmed_comall, dfhighhigh_comall) 
new_bray.urban_comall_as

# reorder
new_bray.urban_comall_as$urbangroup <- factor(new_bray.urban_comall_as$urbangroup, levels=c("low_low", "med_med", "high_high"))

head(new_bray.urban_comall_as)
str(new_bray.urban_comall_as)
levels(new_bray.urban_comall_as$urbangroup)

p_comallnew_brayturnover_as <-ggplot(new_bray.urban_comall_as, aes(x=urbangroup, y=Freq))+ 
  geom_boxplot()+
  #scale_fill_brewer(palette = "Paired")+
  labs(x="Urbanization Level", y="Bray-Curtis Distance", title="Turnover Bacterial Community")+
  theme_bw()+ 
  theme(text = element_text(size=15),legend.position = "none",
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank())
p_comallnew_brayturnover_as

new_bray.urban_comall_se <- summarySE(new_bray.urban_comall, measurevar = "Freq", groupvars = c("urbangroup2"))
new_bray.urban_comall_se

new_bray.urban_comall_asse <- summarySE(new_bray.urban_comall_as, measurevar = "Freq", groupvars = c("urbangroup"))
new_bray.urban_comall_asse

turnoverplot_as <- ggplot(new_bray.urban_comall_asse, aes(x=urbangroup, y=Freq, fill=urbangroup)) +
  geom_col(colour = "black", width = 0.5, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin=Freq-se, ymax=Freq+se), width=.1, position = position_dodge(0.7)) +
  geom_text(
    aes(label= round(Freq, digits = 2)),
    colour = "black", size = 4,
    vjust = -2, position = position_dodge(0.7)
  )+
  scale_fill_manual(values = c('#D2691E', '#9370DB', '#DC143C')) +
  ggtitle("") +
  xlab("Urbanization Level") +
  ylab("Bray Curtis Distance") +
  #expand_limits(y=0:0.2) +                        
  #scale_y_continuous(breaks=0:0.2*0.02) +        
  theme_bw() +
  theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=12)) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  theme(legend.text = element_text(colour="black", size=12)) +
  theme(plot.title = element_text(face = "bold", size=14)) +
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=2)) +
  theme(legend.position="right", legend.box = "horizontal") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

turnoverplot_as

plot_turnover_facetgrid <- ggarrange(turnoverplot_as, turnoverplot_as, turnoverplot_as, turnoverplot_as,
                                     align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                                     common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

plot_turnover_facetgrid

#NormalityTest comall
shapiro.test(new_bray.urban_comall_as$Freq)

#HomogeneityofVariance turnover comall
bartlett.test(Freq ~ urbangroup, new_bray.urban_comall_as)

#Oneway Way Anova turnover comall
res.aov.turnover.Mycoides16scomall_new <- aov(Freq ~ urbangroup, data = new_bray.urban_comall_as)
summary(res.aov.turnover.Mycoides16scomall_new)
TukeyHSD(res.aov.turnover.Mycoides16scomall_new)

#Kruskal
kruskal.test(Freq ~ urbangroup2, data = new_bray.urban_comall_asse)

#Posthoc after Kruskal
pairwise.wilcox.test(new_bray.urban_comall_asse$Freq, new_bray.urban_comall_asse$urbangroup2,
                     p.adjust.method = "BH")

####### Pie chart soil bacteria ##########

function_bacteria=read.csv("function_bacteria_heatmap.csv", row.names=1, header = T, sep=";")

sorted_unique_names <- function_bacteria %>% distinct(Guild) %>% arrange(Guild)
print(sorted_unique_names)

function_bacteria_group1 <- subset(function_bacteria, Treatment == 'Group 1')
function_bacteria_group2 <- subset(function_bacteria, Treatment == 'Group 2')
function_bacteria_group3 <- subset(function_bacteria, Treatment == 'Group 3')

preference_counts_group1 <- table(function_bacteria_group1$Guild)

sorted_preferences_group1 <- sort(preference_counts_group1, decreasing = TRUE)
print(sorted_preferences_group1)

preference_counts_group2 <- table(function_bacteria_group2$Guild)

sorted_preferences_group2 <- sort(preference_counts_group2, decreasing = TRUE)
print(sorted_preferences_group2)

preference_counts_group3 <- table(function_bacteria_group3$Guild)

sorted_preferences_group3 <- sort(preference_counts_group3, decreasing = TRUE)
print(sorted_preferences_group3)

bacteriafunction_colors <- c("aerobic chemoheterotrophy_adipate degradation"="#DBFFD6",
                             "aerobic chemoheterotrophy_cellulolysis"="#ECD4FF",
                             "aerobic chemoheterotrophy_denitrification"="orange",
                             "aerobic chemoheterotrophy_fermentation"="#ACE7FF",
                             "aerobic chemoheterotrophy_hydrocarbon degradation"="#FFABAB",
                             "aerobic chemoheterotrophy_manganese oxidation" = "cyan",
                             "aerobic chemoheterotrophy_metalloid resistant" = "#ffa2c5",
                             "aerobic chemoheterotrophy_nitrification"= "#3498db",
                             "aerobic chemoheterotrophy_nitrogen fixation" = "#33ffbe",
                             "aerobic chemoheterotrophy_sulfur oxidation" = "#ff3333",
                             "aerobic chemoheterotrophy_unknown" = "#842704",
                             "aerobic chemoheterotrophy_ureolysis" = "#f4fb00",
                             "aerobic chemoheterotrophy_xenobiotic degradation" = "#a700fb",
                             "animal parasites or symbionts_human pathogens" = "#E7E8E8",
                             "animal parasites or symbionts_unknown" = "#dab821",
                             "chemoheterotrophy_adipate degradation" = "#FF4500",
                             "chemoheterotrophy_chitinolysis"= "#DA70D6",
                             "chemoheterotrophy_fermentation"= "#4682B4",
                             "chemoheterotrophy_hydrocarbon degradation"= "#006400",
                             "chemoheterotrophy_nitrogen fixation"= "#B8860B",
                             "chemoheterotrophy_sulfur oxidation"= "#FF1493",
                             "chemoheterotrophy_unknown"= "#7FFF00",
                             "chemoheterotrophy_ureolysis"= "#DC143C",
                             "methylotrophy_hydrocarbon degradation"= "#D2691E",
                             "methylotrophy_methanol oxidation"= "#A52A2A",
                             "methylotrophy_unknown"= "gold",
                             "oxygenic photoautotrophy_nitrogen fixation"= "navy",
                             "unknown_unknown" = "beige")

bacteriafunction_colors_part2 <- c("aerobic chemoheterotrophy_adipate degradation"="#DBFFD6",
                                   "aerobic chemoheterotrophy_cellulolysis"="#DBFFD6",
                                   "aerobic chemoheterotrophy_denitrification"="#DBFFD6",
                                   "aerobic chemoheterotrophy_fermentation"="#DBFFD6",
                                   "aerobic chemoheterotrophy_hydrocarbon degradation"="#DBFFD6",
                                   "aerobic chemoheterotrophy_manganese oxidation" = "#DBFFD6",
                                   "aerobic chemoheterotrophy_metalloid resistant" = "#DBFFD6",
                                   "aerobic chemoheterotrophy_nitrification"= "#DBFFD6",
                                   "aerobic chemoheterotrophy_nitrogen fixation" = "#DBFFD6",
                                   "aerobic chemoheterotrophy_sulfur oxidation" = "#DBFFD6",
                                   "aerobic chemoheterotrophy_unknown" = "#DBFFD6",
                                   "aerobic chemoheterotrophy_ureolysis" = "#DBFFD6",
                                   "aerobic chemoheterotrophy_xenobiotic degradation" = "#DBFFD6",
                                   "animal parasites or symbionts_human pathogens" = "#E7E8E8",
                                   "animal parasites or symbionts_unknown" = "#E7E8E8",
                                   "chemoheterotrophy_adipate degradation" = "#ECD4FF",
                                   "chemoheterotrophy_chitinolysis"= "#ECD4FF",
                                   "chemoheterotrophy_fermentation"= "#ECD4FF",
                                   "chemoheterotrophy_hydrocarbon degradation"= "#ECD4FF",
                                   "chemoheterotrophy_nitrogen fixation"= "#ECD4FF",
                                   "chemoheterotrophy_sulfur oxidation"= "#ECD4FF",
                                   "chemoheterotrophy_unknown"= "#ECD4FF",
                                   "chemoheterotrophy_ureolysis"= "#ECD4FF",
                                   "methylotrophy_hydrocarbon degradation"= "#ffa2c5",
                                   "methylotrophy_methanol oxidation"= "#ffa2c5",
                                   "methylotrophy_unknown"= "#ffa2c5",
                                   "oxygenic photoautotrophy_nitrogen fixation"= "gold",
                                   "unknown_unknown" = "beige")

piechart_bacteria=read.csv("piechart_heatmap_bacteria.csv", row.names=1, header = T, sep=";")

piechart_group1_bac <- subset(piechart_bacteria, treatment == 'Group1')
piechart_group2_bac <- subset(piechart_bacteria, treatment == 'Group2')
piechart_group3_bac <- subset(piechart_bacteria, treatment == 'Group3')

piechart_all <- ggplot(piechart_bacteria, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(label=abundance), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  facet_grid(.~treatment, scales="free_x") +
  scale_fill_manual(values = bacteriafunction_colors_part2) +
  guides(fill = guide_legend(title = "Guild")) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

piechart_all

plot_piechart_all <- ggarrange(piechart_all, piechart_all,
                               align = "v", font.label = list(size = 12, face = "bold", family="serif"), 
                               common.legend = TRUE, legend = "right", ncol = 1, nrow = 2)

plot_piechart_all

library(dplyr)
piechart_group1_try <- piechart_group1_bac %>% 
  mutate(end = 2 * pi * cumsum(abundance)/sum(abundance),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))


piechart_group2_try <- piechart_group2_bac %>% 
  mutate(end = 2 * pi * cumsum(abundance)/sum(abundance),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

piechart_group3_try <- piechart_group3_bac %>% 
  mutate(end = 2 * pi * cumsum(abundance)/sum(abundance),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

piechart_group1_try <- piechart_group1_try[order(piechart_group1_try$guild),]
piechart_group1_try

piechart_group2_try <- piechart_group2_try[order(piechart_group2_try$guild),]
piechart_group2_try

piechart_group3_try <- piechart_group3_try[order(piechart_group3_try$guild),]
piechart_group3_try

piechart_group1_try$guild <- factor(piechart_group1_try$guild, levels = sort(unique(piechart_group1_try$guild)))
piechart_group2_try$guild <- factor(piechart_group2_try$guild, levels = sort(unique(piechart_group2_try$guild)))
piechart_group3_try$guild <- factor(piechart_group3_try$guild, levels = sort(unique(piechart_group3_try$guild)))

library(ggforce) # for 'geom_arc_bar'
ggplot(piechart_group1_try) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, fill = guild), colour="white", linewidth=0.7) +
  geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = combine,
                hjust = hjust, vjust = vjust)) +
  scale_fill_manual(values = bacteriafunction_colors_part2) +
  coord_fixed () +
  scale_x_continuous(limits = c(-2, 4.5),  # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-3, 1.1),    # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        #axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        #axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Barplot
piechart_one_bac <- ggplot(piechart_group1_try, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(x = 1.5 * sin(middle), y = 1.5 * cos(middle), label = abundance,
                hjust = hjust, vjust = vjust), size = 5) +
  coord_polar(theta = "y")+
  scale_fill_manual(values = bacteriafunction_colors_part2) +
  guides(fill = guide_legend(title = "Guild")) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

piechart_one_bac

piechart_two_bac <- ggplot(piechart_group2_try, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(x = 1.5 * sin(middle), y = 1.5 * cos(middle), label = abundance,
                hjust = hjust, vjust = vjust), size = 6) +
  coord_polar(theta = "y")+
  scale_fill_manual(values = bacteriafunction_colors_part2) +
  guides(fill = guide_legend(title = "Guild")) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

piechart_two_bac

piechart_three_bac <- ggplot(piechart_group3_try, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(x = 1.5 * sin(middle), y = 1.5 * cos(middle), label = abundance,
                hjust = hjust, vjust = vjust), size = 6) +
  coord_polar(theta = "y")+
  scale_fill_manual(values = bacteriafunction_colors_part2) +
  guides(fill = guide_legend(title = "Guild")) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

piechart_three_bac

piechart_two_bac <- ggplot(piechart_group2_bac, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(label=abundance), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  scale_fill_manual(values = bacteriafunction_colors_part2) +
  guides(fill = guide_legend(title = "Guild")) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

piechart_two_bac

piechart_three_bac <- ggplot(piechart_group3_bac, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(label=abundance), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  scale_fill_manual(values = bacteriafunction_colors_part2) +
  guides(fill = guide_legend(title = "Guild")) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

piechart_three_bac

plot_piechart_bac <- ggarrange(piechart_one_bac, piechart_two_bac, piechart_three_bac, piechart_one_bac, piechart_two_bac, piechart_three_bac,
                               align = "v", font.label = list(size = 12, face = "bold", family="serif"), ncol = 3, nrow = 1, 
                               common.legend = TRUE, legend = "bottom" )

plot_piechart_bac
