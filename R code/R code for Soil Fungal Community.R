### R Analysis of Soil Bacterial Community 
#For the paper published in npj Biofilms and Microbiomes
#Unveiling the ecological processes driving soil and lichen microbiome assembly along an urbanization gradient
#Authored by: Panji Cahya Mawarda, Rens van der Kaaij, Francisco Dini-Andreote, Deniz Duijker, Michael Stech, Arjen Speksnijder

#Created by: Panji Cahya Mawarda 04-04-2024

#call necessary library packages

#####################Start here##########################

biomfile=biomformat::read_biom("ITS Soil.biom2")

taxonomy_biom=observation_metadata(biomfile)
taxonomy_biom

write.table(taxonomy_biom, 'taxonomy_ITS_soil_new.csv')

otu_biom=import_biom("ITS Soil.biom2")
write.csv(x = otu_table(otu_biom), file = "otu_table_ITS_soil_new.csv")

#reading files, you need four inputs: 
#otu table, taxonomy table, metadata

otu_table=read.csv("otu_table_ITS_soil_new.csv", row.names=1, header = T)
head(otu_table)
otu_table=as.matrix(otu_table)

taxonomy=read.csv("taxonomy_ITS_soil_new.csv", row.names=1, header = T)
replace(taxonomy, taxonomy == " ", NA) -> tax1
taxonomy1=as.matrix(taxonomy)

#Taxonomy table for the tree----------
taxonomy_tree = as.data.frame(taxonomy_biom)
taxonomy_tree["taxid"] = rownames(taxonomy_tree)
names (taxonomy_tree) <- NULL

write.table(taxonomy_tree, sep=",", quote=FALSE, row.names= FALSE, 'taxonomy_ITS_soil_tree.csv')

#perform the function below in powershell windows, inside the directory where you have both taxonomy_to_tree.pl and taxonomy_ITS_lichen_tree.csv
#perl .\taxonomy_to_tree.pl -b 60 -s "," -f .\taxonomy_ITS_Lichen_tree.csv > tree_ITS_lichen
#afterwards you have to open your file in notepad, and save as UTF-08

metadata=read.csv("metadata_ITS_soil_new.csv", row.names=1, header=T)
metadata

OTU=otu_table(otu_table, taxa_are_rows=TRUE)

TAX=tax_table(taxonomy1)
head(TAX)

META=sample_data(metadata)
metadata

tre_fpath = "tree_ITS_soil_01.nwk"
TRE <- ape::read.tree(tre_fpath)

#make sure files have the same sample names
sample_names(OTU)
sample_names(META)

#make the phyloseq object for downstream analysis
physeq=phyloseq(OTU, TAX, META, TRE)
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

#Remove Archaea, Mitochondria, Chloroplast, and unnecessary taxons 

physeqbac <- physeq %>% subset_taxa(Kingdom == "Fungi" &
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

ggrare(physeqbac, step = 100, color = "uhi_range", se = F)

# rarefy your dataset to lowest sequencing length, 711 should be kept constant:
phy_rare <- rarefy_even_depth(physeqbac, sample.size = 62836, rngseed = 711, replace = F, trimOTUs = TRUE, 	verbose = F)
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


################ Phyrare sub sample Based on Temperature Zone/UHI ZOne #######################

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
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
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

set.seed(123)
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
  expand_limits(y=0:2000) +
  scale_y_continuous(breaks=seq(0,2000, by=200)) +
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
  expand_limits(y=0:7) +
  scale_y_continuous(breaks=seq(0,7, by=1)) +
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

#Kruskal shannon All Based on Temperature Zone -------------------
# IF the data is not normally distributed and the variance is not homogen, do Kruskal Wallis test
kruskal.test(value ~ uhi_range, data = alphadte.shannon)

#Posthoc analysis with Wilcox test
pairwise.wilcox.test(alphadte.shannon$value, alphadte.shannon$temperature_zone,
                     p.adjust.method = "BH")

#-----------------------BETA DIVERSITY --------------------------------

#Exporting OTU Table --------------------------------------------------
write.csv(x = otu_table(phy_rare), file = "phyraresuball_asv_ITS_soil.csv")

# Based temperature zone----------------------------------------------------------------------------
#all
comall_tempzone <- read.csv("phyraresuball_asv_ITS_soil.csv", header = TRUE, row.names = 1)
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
rand.time=100
prefix="Test"

# 2. Load data and match IDs ----

comm=t(data.frame(otu_table(phy_rare)))
group=data.frame(sample_data(phy_rare))
tree=phy_tree(phy_rare)

samp.ck=NST::match.name(rn.list=list(comm=comm,group=group))
comm=samp.ck$comm
comm=comm[,colSums(comm)>0,drop=FALSE]
group=samp.ck$group 

tax.ck=NST::match.name(cn.list = list(comm=comm), tree.list = list(tree=tree))
comm=tax.ck$comm
tree=tax.ck$tree


# 3. Calculate bNTI and RCI ----
## 3.1. Using iCAMP - Method 1 ----

set.seed(1234)
qpen <- qpen(comm = comm,
             pd = NULL, 
             pd.big.wd = NULL,
             pd.big.spname = NULL, 
             tree = tree,
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

plot_betanti_facetgrid <- ggarrange(plot_betanti_tempzone, plot_rc_tempzone, plot_betanti_tempzone, plot_rc_tempzone,
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

# Note ----
# bNTI < -2 = Homogeneous selection
# bNTI > +2 = Variable selection
# abs(bNTI) < 2 & RC < -0.95 = Homogenizing dispersal
# abs(bNTI) < 2 & RC < +0.95 = Undominated
# abs(bNTI) < 2 & RC > +0.95 = Dispersal limitation 

#---------------------------------Succesional Pattern--------------------------------

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

write.csv(x = otu_table(phy_rare), file = "phyrare_heatmap_ok.csv")

#after you download, change the name of the column from code to temperature zone_tree_replicate

otu_table_all=read.csv("phyrare_heatmap.csv", row.names=1, header = 1, check.names = F)
head(otu_table_all)

otu_table_all <- otu_table_all[rowSums(otu_table_all) > 20 & rowSums(otu_table_all > 0) >= 10, ] #minimal abundance 20 ditemukan di 10 sampel, rule of thumb nya ambil ASV yang banyak ditemukan di banyak sampel untuk mempermudah perhitungan dan visualisasi


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

df_successionall_significant <- df_successionall_subset %>% select(Group, OTU2,OTU7,OTU8,OTU32,OTU33,OTU37,OTU41,OTU42,OTU45,OTU80,OTU82,OTU89,OTU119,OTU121,OTU122,OTU123,OTU124,OTU127,OTU130,OTU139,OTU140,OTU146,OTU160,OTU162,OTU230,OTU234,OTU261,OTU286,OTU316,OTU359,OTU372,OTU390,OTU399,OTU410,OTU423,OTU424,OTU425,OTU429,OTU441,OTU464,OTU465,OTU467,OTU472,OTU485,OTU507,OTU515,OTU622,OTU673,OTU675,OTU684,OTU700,OTU703,OTU704,OTU705,OTU712,OTU716,OTU740,OTU743,OTU771,OTU773,OTU776,OTU790,OTU791,OTU795,OTU801,OTU820,OTU830,OTU875,OTU876,OTU902,OTU903,OTU904,OTU912,OTU931,OTU934,OTU938,OTU947,OTU972,OTU990,OTU1007,OTU1016,OTU1047,OTU1088,OTU1094,OTU1099,OTU1122,OTU1137,OTU1138,OTU1151,OTU1237,OTU1244,OTU1276,OTU1281,OTU1286,OTU1287,OTU1295,OTU1312,OTU1317,OTU1326,OTU1332,OTU1334,OTU1347,OTU1354,OTU1360,OTU1378,OTU1403,OTU1453,OTU1458,OTU1471,OTU1479,OTU1480,OTU1515,OTU1535,OTU1545,OTU1553,OTU1568,OTU1578,OTU1609,OTU1613,OTU1676,OTU1692,OTU1694,OTU1714,OTU1718,OTU1851,OTU1853,OTU1856,OTU1882,OTU1908,OTU1916,OTU1923,OTU1968,OTU1976,OTU2014,OTU2061,OTU2063,OTU2067,OTU2158,OTU2162,OTU2182,OTU2189,OTU2198,OTU2207,OTU2208,OTU2212,OTU2213,OTU2218,OTU2219,OTU2220,OTU2226,OTU2233,OTU2240)
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

#Call:
#Procrustes Sum of Squares (m12 squared):        0.3048 
#Correlation in a symmetric Procrustes rotation: 0.8338 
#Significance:  1e-04 

#Permutation: free
#Number of permutations: 9999

plot(vare.proc.all)
plot(vare.proc.all, kind=2)
residuals(vare.proc.all)


#Put the taxonomy Phylum on excel file
#Create another excel file from filter table and replace the OTUCode to OTUID

df_successionall_filter=t(df_successionall_filter)
write.csv(x = (df_successionall_filter), file = "df_successionall_filter.csv")

# Heatmap all----

heatmap_PM1_data <- read.csv("df_heatmap_reduced.csv",
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

#Taxonomy Bar Plot ----------------------------------------------
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

data.phyla.all <- physeq_heatmap %>%
  tax_glom(taxrank = "Genus") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>% 
  filter(Abundance > 0.01) %>%
  arrange(Phylum) 
data.phyla.all

data.phyla.all[data.phyla.all == "Unassigned"] <- "Zunassigned"

phylum_colors <- c("Ascomycota" = "pink", "Basidiobolomycota" = "blue", "Armatimonadetes" = "beige",
                   "Basidiomycota" = "gold", "Chlamydiae" = "grey", "Glomeromycota" = "red",
                   "Mortierellomycota" = "purple", "Mucoromycota" = "cyan", "Gemmatimonadetes" = "purple",
                   "Nitrospirae" = "lavender", "OD1" = "salmon", "Planctomycetes" = "brown", 
                   "Zoopagomycota" = "darkgreen", "TM7"= "skyblue", "Verrucomicrobia" = "violet",
                   "Olpidiomycota" = "black")

class_colors <- c("Lecanoromycetes" = "pink", "Exobasidiomycetes" = "blue", "Eurotiomycetes" = "beige",
                  "Agaricomycetes" = "gold", "Mucoromycetes" = "grey", "Mortierellomycetes" = "red",
                  "Sordariomycetes" = "lavender", "Cystobasidiomycetes" = "cyan", "Dothideomycetes" = "purple",
                  "Archaeosporomycetes" = "maroon", "Leotiomycetes" = "salmon", "Tremellomycetes" = "brown", 
                  "Microbotryomycetes" = "darkgreen", "Olpidiomycota_cls_Incertae_sedis"= "skyblue", "Basidiomycota_cls_Incertae_sedis" = "violet",
                  "Zunassigned" = "black", "Zoopagomycetes" = "#808000", "Saccharomycetes" = "orange", 
                  "Pucciniomycetes" =  "#FF4500", "Pezizomycetes" = "#000080", "GS27" = "#ADFF2F" )

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

fs_Rtempone_cor <- read.table('temperaturezone1_filtered_correlation.tsv', 
                              header = T, sep = '\t')

fs_Rtempone_p <- read.table('temperaturezone1_filtered_pvalue.tsv', 
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

fs_Rtemptwo_cor <- read.table('temperaturezone2_filtered_correlation.tsv', 
                              header = T, sep = '\t')
fs_Rtemptwo_p <- read.table('temperaturezone2_filtered_pvalue.tsv', 
                            header = T, sep = '\t')

fs_Rtemptwo_combine <- cbind.data.frame(fs_Rtemptwo_cor, fs_Rtemptwo_p[,3])
colnames(fs_Rtemptwo_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rtemptwo_filt <- fs_Rtemptwo_combine  %>% filter(abs(cor) > 0.8 & pval < 0.01)

# Export to Cytoscape ----
fs_Rtemptwo_filt %>%
  write.table('./temperaturezone2_net2Cys.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

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

fs_Rtempthree_cor <- read.table('temperaturezone3_filtered_correlation.tsv', 
                                header = T, sep = '\t')
fs_Rtempthree_p <- read.table('temperaturezone3_filtered_pvalue.tsv', 
                              header = T, sep = '\t')

fs_Rtempthree_combine <- cbind.data.frame(fs_Rtempthree_cor, fs_Rtempthree_p[,3])
colnames(fs_Rtempthree_combine)[3:4] = c('cor','pval')

# Filter based on r and p values ----
fs_Rtempthree_filt <- fs_Rtempthree_combine  %>% filter(abs(cor) > 0.8 & pval < 0.01)

# Export to Cytoscape ----
fs_Rtempthree_filt %>%
  write.table('./temperaturezone3_net2Cys.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Open cytoscape

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

#Pie Chart Predicted Guild ---------------------

####### Pie chart soil fungi##########

library(ggplot2)
fungifunction_colors <- c("animal parasite"="#DBFFD6",
                          "litter saprotroph"="#ECD4FF",
                          "mycoparasite"="orange",
                          "soil saprotroph"="#ACE7FF",
                          "plant pathogen"="#FFABAB",
                          "unknown" = "beige",
                          "wood saprotroph" = "#ffa2c5",
                          "dung saprotroph"= "#3498db",
                          "ectomycorrhizal" = "#33ffbe",
                          "epiphyte" = "#ff3333",
                          "foliar_endophyte" = "#842704",
                          "lichenized" = "#f4fb00",
                          "sooty mold" = "#a700fb",
                          "nectar saprotroph" = "#E7E8E8",
                          "unspecified saprotroph" = "#dab821")

piechart_fungi=read.csv("piechart_fungi_function.csv", row.names=1, header = T, sep=";")

piechart_group1 <- subset(piechart_fungi, treatment == 'Group1')
piechart_group2 <- subset(piechart_fungi, treatment == 'Group2')
piechart_group3 <- subset(piechart_fungi, treatment == 'Group3')

piechart_all <- ggplot(piechart_fungi, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(label=abundance), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  facet_grid(.~treatment, scales="free_x") +
  scale_fill_manual(values = fungifunction_colors) +
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

# Barplot
piechart_one <- ggplot(piechart_group1, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(label=abundance), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  scale_fill_manual(values = fungifunction_colors) +
  guides(fill = guide_legend(title = "Guild")) +
  #scale_y_continuous(breaks = piechart_group1$abundance, labels = piechart_group1$guild) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

piechart_one

piechart_two <- ggplot(piechart_group2, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(label=abundance), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  scale_fill_manual(values = fungifunction_colors) +
  guides(fill = guide_legend(title = "Guild")) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

piechart_two

piechart_three <- ggplot(piechart_group3, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(label=abundance), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  scale_fill_manual(values = fungifunction_colors) +
  guides(fill = guide_legend(title = "Guild")) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

piechart_three

plot_piechart <- ggarrange(piechart_one, piechart_two, piechart_three, piechart_one, piechart_two, piechart_three,
                           align = "v", font.label = list(size = 12, face = "bold", family="serif"), 
                           common.legend = TRUE, legend = "right", ncol = 3, nrow = 2)

plot_piechart
