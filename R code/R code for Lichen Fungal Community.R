### R Analysis of Lichen Fungal Community 
#For the paper published in npj Biofilms and Microbiomes
#Unveiling the ecological processes driving soil and lichen microbiome assembly along an urbanization gradient
#Authored by: Panji Cahya Mawarda, Rens van der Kaaij, Francisco Dini-Andreote, Deniz Duijker, Michael Stech, Arjen Speksnijder

#Created by: Panji Cahya Mawarda 04-04-2024

#call necessary library packages

#####################Start here##########################

#reading files, you need three inputs: 
#otu table, taxonomy table, metadata
biomfile=biomformat::read_biom("ITS_Lichen.biom2")

taxonomy_biom=biomformat::observation_metadata(biomfile)
taxonomy_biom

write.table(taxonomy_tree, sep=",", 'taxonomy_ITS_lichen_tree.csv')

#Taxonomy table for the tree----------

taxonomy_tree = as.data.frame(taxonomy_biom)
taxonomy_tree["taxid"] = rownames(taxonomy_tree)
names (taxonomy_tree) <- NULL

write.table(taxonomy_tree, sep=",", quote=FALSE, row.names= FALSE, 'taxonomy_ITS_lichen_tree.csv')

#perform the function below in powershell windows, inside the directory where you have both taxonomy_to_tree.pl and taxonomy_ITS_lichen_tree.csv
#perl .\taxonomy_to_tree.pl -b 60 -s "," -f .\taxonomy_ITS_Lichen_tree.csv > tree_ITS_lichen
#afterwards you have to open your file in notepad, and save as UTF-08

otu_biom=phyloseq::import_biom("ITS_Lichen.biom2")
write.csv(x = otu_table(otu_biom), file = "otu_table_ITS_lichen_new.csv")

#reading files, you need four inputs: 
#otu table, taxonomy table, metadata

otu_table=read.csv("otu_table_ITS_lichen_new.csv", row.names=1, header = T)
head(otu_table)
otu_table=as.matrix(otu_table)

taxonomy=read.csv("taxonomy_ITS_lichen_new.csv", row.names=1, header = T)
replace(taxonomy, taxonomy == " ", NA) -> tax1
taxonomy1=as.matrix(taxonomy)

metadata=read.csv("metadata_ITS_Lichen_new.csv", row.names=1, header=T)
metadata

OTU=otu_table(otu_table, taxa_are_rows=TRUE)

TAX=tax_table(taxonomy1)
head(TAX)

META=sample_data(metadata)
metadata

tre_fpath = "tree_ITS_lichen_03.nwk"
TRE <- ape::read.tree(tre_fpath)

#make sure files have the same sample names
sample_names(OTU)
sample_names(META)

set.seed(123)
#make the phyloseq object for downstream analysis
physeq=phyloseq(OTU, TAX, META, TRE)
physeq

sample_data(physeq)
sample_sums(physeq)
colnames(tax_table(physeq))
dim(tax_table(physeq))
tax_table(physeq)[1:5, 1:7]
sample_variables(physeq)
otu_table(physeq)[1:5, 1:45]
phy_tree(physeq)

############ Filter OTU table ############

#Remove Mitochondria, Chloroplast, and unnecessary taxons 

physeqbac <- physeq %>% subset_taxa(Kingdom == "Fungi" &
                                      Family != "Mitochondria" &
                                      Class != "Chloroplast" &
                                      Kingdom != "Archaea" &
                                      Kingdom != "Unassigned")

get_taxa_unique(physeqbac,"Kingdom")
get_taxa_unique(physeqbac,"Family")
get_taxa_unique(physeqbac,"Class")

colnames(tax_table(physeq))
otu_table(physeqbac)[1:5, 1:45]
tax_table(physeqbac)[1:5, 1:7]
physeqbac

tax_table(physeqbac)[1:5, 1:7]
any(taxa_sums(physeqbac) < 1) #false
any(is.na(otu_table(physeqbac))) #false

###discard taxons with phylum Unidentified (optional)

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
phy_rare <- rarefy_even_depth(physeqbac, sample.size = 111227, rngseed = 711, replace = F, trimOTUs = TRUE, 	verbose = F)
phy_rare

library(phyloseq.extended)
ggrare(phy_rare, step = 100, color = "uhi_range", se = F)
sample_names(phy_rare)
sample_sums(phy_rare)

################################### Alpha diversity ############################################

## quantify species richness, shannon diversity index, Simpson, and Chao1,
##However, we usually just take species richness (observed), and shannon diversity index

#species richness refers to the number of different bacterial species or taxa present in a given microbiome sample
#shannon diversity index: It takes into account both the number of different species (species richness) and their 
#relative abundance. The index is calculated based on the natural logarithm of the number of individuals (or sequencing reads) belonging to each taxonomic group in the sample, 
#multiplied by the proportional abundance of that taxonomic group in the sample. 

estimate_richness(phy_rare, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
pAlpha = plot_richness(phy_rare,
                       shape = "uhi_range",
                       color = "lichen_species",
                       measures = c("Observed", "Simpson", "Shannon"),
                       title = "Alpha Diversity")
pAlpha + geom_point(size = 3)

dalltreatments<-pAlpha$data

dalltreatments


################ Phyrare sub sample Based on Temperature #######################

phyraresubtempone <- subset_samples(phy_rare,!temperature_zone %in% c("temperaturezone2", "temperaturezone3"))

phyraresubtempone

phyraresubtemptwo <- subset_samples(phy_rare,!temperature_zone %in% c("temperaturezone1", "temperaturezone3"))

phyraresubtemptwo

phyraresubtempthree <- subset_samples(phy_rare,!temperature_zone %in% c("temperaturezone1", "temperaturezone2"))

phyraresubtempthree

################ Phyrare sub sample Based on lichen species #######################

phyraresub_xanthoriaparietina <- subset_samples(phy_rare,lichen_species %in% c("Xanthoria parietina"))

phyraresub_xanthoriaparietina

phyraresub_physciaadscendens <- subset_samples(phy_rare,lichen_species %in% c("Physcia adscendens"))

phyraresub_physciaadscendens

phyraresub_candelariaconcolor <- subset_samples(phy_rare,lichen_species %in% c("Candelaria concolor"))

phyraresub_candelariaconcolor

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

richness_all_tempzone_facetgrid <- ggplot(data = alphadte.richness,
                                          aes(x=uhi_range, y = value, color=uhi_range))+
  facet_grid(.~lichen_species, scales="free_x") +
  geom_boxplot(lwd=0.1)+
  #geom_jitter(position=position_jitter(0.1), cex=5, alpha=0.7) +
  expand_limits(y=200:1000) +
  scale_y_continuous(breaks=seq(200,1000, by=100)) +
  labs(x="UHI range", y="ASV richness", title="ASV richness")+
  theme_bw()+
  scale_color_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

richness_all_tempzone_facetgrid

richness_all_lichen_facetgrid <- ggplot(data = alphadte.richness,
                                        aes(x=lichen_species, y = value, color=lichen_species))+
  facet_grid(.~uhi_range, scales="free_x") +
  geom_boxplot(lwd=0.1)+
  #geom_jitter(position=position_jitter(0.1), cex=5, alpha=0.7) +
  expand_limits(y=200:1000) +
  scale_y_continuous(breaks=seq(200,1000, by=100)) +
  labs(x="Lichen species", y="ASV richness", title="ASV richness")+
  theme_bw()+
  scale_color_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

richness_all_lichen_facetgrid

richness_all_tempzone <- ggplot(data = alphadte.richness,
                                aes(x=uhi_range, y = value, color=uhi_range))+
  geom_boxplot(lwd=0.1)+
  geom_jitter(position=position_jitter(0.1), cex=5, alpha=0.7) +
  expand_limits(y=200:1000) +
  scale_y_continuous(breaks=seq(200,1000, by=100)) +
  labs(x="UHI range", y="ASV richness", title="ASV richness")+
  theme_bw()+
  scale_color_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

richness_all_tempzone

richness_all_lichen <- ggplot(data = alphadte.richness,
                              aes(x=lichen_species, y = value, color=lichen_species))+
  geom_boxplot(lwd=0.1)+
  geom_jitter(position=position_jitter(0.1), cex=5, alpha=0.7) +
  expand_limits(y=200:1000) +
  scale_y_continuous(breaks=seq(200,1000, by=100)) +
  labs(x="Lichen species", y="ASV richness", title="ASV richness")+
  theme_bw()+
  scale_color_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

richness_all_lichen

richness_all_facetgrid <- ggarrange(richness_all_tempzone_facetgrid, richness_all_tempzone, richness_all_lichen_facetgrid, richness_all_lichen,
                                    align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                                    common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

richness_all_facetgrid

#Richness all--------------------------------------------------------------

alphadte.richness_tempzone1 <- alphadte.richness %>% filter(temperature_zone == "temperaturezone1")
alphadte.richness_tempzone2 <- alphadte.richness %>% filter(temperature_zone == "temperaturezone2")
alphadte.richness_tempzone3 <- alphadte.richness %>% filter(temperature_zone == "temperaturezone3")

alphadte.richness_candelaria <- alphadte.richness %>% filter(lichen_species == "Candelaria concolor")
alphadte.richness_physcia <- alphadte.richness %>% filter(lichen_species == "Physcia adscendens")
alphadte.richness_xanthoria <- alphadte.richness %>% filter(lichen_species == "Xanthoria parietina")


#normal distribution test
shapiro.test(alphadte.richness$value)
shapiro.test(alphadte.richness_tempzone1$value)
shapiro.test(alphadte.richness_tempzone2$value)
shapiro.test(alphadte.richness_tempzone3$value)

shapiro.test(alphadte.richness_candelaria$value)
shapiro.test(alphadte.richness_physcia$value)
shapiro.test(alphadte.richness_xanthoria$value)

#homogeneity of variance test
bartlett.test(value ~ temperature_zone, alphadte.richness)
bartlett.test(value ~ lichen_species, alphadte.richness_tempzone1)
bartlett.test(value ~ lichen_species, alphadte.richness_tempzone2)
bartlett.test(value ~ lichen_species, alphadte.richness_tempzone3)

bartlett.test(value ~ temperature_zone, alphadte.richness_candelaria)
bartlett.test(value ~ temperature_zone, alphadte.richness_physcia)
bartlett.test(value ~ temperature_zone, alphadte.richness_xanthoria)

#Two Way ANOVA Richness Based on Temperature Zone------------------
set.seed(123)
res.aov.richness.tempzone <- aov(value ~ uhi_range * lichen_species, data = alphadte.richness)
summary(res.aov.richness.tempzone)

#ANOVA Richness Based on Temperature Zone------------------
res.aov.richness.tempzone <- aov(value ~ uhi_range, data = alphadte.richness)
summary(res.aov.richness.tempzone)

res.aov.richness.tempzone_temp1 <- aov(value ~ lichen_species, data = alphadte.richness_tempzone1)
summary(res.aov.richness.tempzone_temp1)

TukeyHSD(res.aov.richness.tempzone_temp1)

set.seed(123)
res.aov.richness.tempzone_temp2 <- aov(value ~ lichen_species, data = alphadte.richness_tempzone2)
summary(res.aov.richness.tempzone_temp2)

TukeyHSD(res.aov.richness.tempzone_temp2)

res.aov.richness.tempzone_temp3 <- aov(value ~ lichen_species, data = alphadte.richness_tempzone3)
summary(res.aov.richness.tempzone_temp3)

TukeyHSD(res.aov.richness.tempzone_temp3)

#ANOVA Richness Based on Lichen Species------------------
res.aov.richness.lichen <- aov(value ~ lichen_species, data = alphadte.richness)
summary(res.aov.richness.lichen)

TukeyHSD(res.aov.richness.lichen)

res.aov.richness.lichen <- aov(value ~ uhi_range, data = alphadte.richness_candelaria)
summary(res.aov.richness.lichen)

res.aov.richness.lichen <- aov(value ~ uhi_range, data = alphadte.richness_physcia)
summary(res.aov.richness.lichen)

res.aov.richness.lichen <- aov(value ~ uhi_range, data = alphadte.richness_xanthoria)
summary(res.aov.richness.lichen)

#Kruskal Richness All Based on Temperature Zone -------------------
# IF the data is not normally distributed and the variance is not homogen, do Kruskal Wallis test
kruskal.test(value ~ uhi_range, data = alphadte.richness)

#Afterwards, you need to perform posthoc analysis, if you perform Kruskal Wallis, then do Conover test or Wilcox

#Posthoc analysis with Wilcox test
pairwise.wilcox.test(alphadte.richness$value, alphadte.richness$temperature_zone,
                     p.adjust.method = "BH")

library(FSA)
dunnTest(value ~ temperature_zone,
         data=alphadte.richness,
         method="bh")

#NA

library(DescTools)

NemenyiTest(x = alphadte.richness$value,
            g = alphadte.richness$temperature_zone,
            dist="tukey")

#Kruskal Richness All Based on lichen_species -------------------
# IF the data is not normally distributed and the variance is not homogen, do Kruskal Wallis test
kruskal.test(value ~ lichen_species, data = alphadte.richness)

#Afterwards, you need to perform posthoc analysis, if you perform Kruskal Wallis, then do Conover test or Wilcox

#Posthoc analysis with Wilcox test
pairwise.wilcox.test(alphadte.richness$value, alphadte.richness$lichen_species,
                     p.adjust.method = "BH")

#Two Way Kruskal---------------
library(rcompanion)

scheirerRayHare(value ~ temperature_zone + lichen_species, data = alphadte.richness)

########### Shannon Diversity Index #####################

alphadte = data.table(dalltreatments)
alphadte

# Subset to just the richness index
alphadte.shannon <- alphadte[(variable == "Shannon")]
alphadte.shannon
str(alphadte.shannon)

shannon_all_tempzone_facetgrid <- ggplot(data = alphadte.shannon,
                                         aes(x=temperature_zone, y = value, color=temperature_zone))+
  facet_grid(.~lichen_species, scales="free_x") +
  geom_boxplot(lwd=0.1)+
  #geom_jitter(position=position_jitter(0.1), cex=5, alpha=0.7) +
  expand_limits(y=1:4) +
  scale_y_continuous(breaks=seq(1,4, by=1)) +
  labs(x="Temperature Zone", y="Shannon Diversity Index", title="Shannon")+
  theme_bw()+
  scale_color_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

shannon_all_tempzone_facetgrid

shannon_all_lichen_facetgrid <- ggplot(data = alphadte.shannon,
                                       aes(x=lichen_species, y = value, color=lichen_species))+
  facet_grid(.~temperature_zone, scales="free_x") +
  geom_boxplot(lwd=0.1)+
  #geom_jitter(position=position_jitter(0.1), cex=5, alpha=0.7) +
  expand_limits(y=1:4) +
  scale_y_continuous(breaks=seq(1,4, by=1)) +
  labs(x="Lichen species", y="Shannon Diversity Index", title="Shannon")+
  theme_bw()+
  scale_color_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

shannon_all_lichen_facetgrid

shannon_all_tempzone <- ggplot(data = alphadte.shannon,
                               aes(x=temperature_zone, y = value, color=temperature_zone))+
  geom_boxplot(lwd=0.1)+
  geom_jitter(position=position_jitter(0.1), cex=5, alpha=0.7) +
  expand_limits(y=1:4) +
  scale_y_continuous(breaks=seq(1,4, by=1)) +
  labs(x="Temperature Zone", y="Shannon Diversity Index", title="Shannon")+
  theme_bw()+
  scale_color_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

shannon_all_tempzone

shannon_all_lichen <- ggplot(data = alphadte.shannon,
                             aes(x=lichen_species, y = value, color=lichen_species))+
  geom_boxplot(lwd=0.1)+
  geom_jitter(position=position_jitter(0.1), cex=5, alpha=0.7) +
  expand_limits(y=1:4) +
  scale_y_continuous(breaks=seq(1,4, by=1)) +
  labs(x="Lichen species", y="Shannon Diversity Index", title="Shannon")+
  theme_bw()+
  scale_color_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

shannon_all_lichen

shannon_all_facetgrid <- ggarrange(shannon_all_tempzone_facetgrid, shannon_all_tempzone, shannon_all_lichen_facetgrid, shannon_all_lichen,
                                   align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                                   common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

shannon_all_facetgrid

#Shannon all--------------------------------------------------------------
alphadte.shannon_tempzone1 <- alphadte.shannon %>% filter(temperature_zone == "temperaturezone1")
alphadte.shannon_tempzone2 <- alphadte.shannon %>% filter(temperature_zone == "temperaturezone2")
alphadte.shannon_tempzone3 <- alphadte.shannon %>% filter(temperature_zone == "temperaturezone3")

alphadte.shannon_candelaria <- alphadte.shannon %>% filter(lichen_species == "Candelaria concolor")
alphadte.shannon_physcia <- alphadte.shannon %>% filter(lichen_species == "Physcia adscendens")
alphadte.shannon_xanthoria <- alphadte.shannon %>% filter(lichen_species == "Xanthoria parietina")


#normal distribution test
shapiro.test(alphadte.shannon$value)
shapiro.test(alphadte.shannon_tempzone1$value)
shapiro.test(alphadte.shannon_tempzone2$value)
shapiro.test(alphadte.shannon_tempzone3$value)

shapiro.test(alphadte.shannon_candelaria$value)
shapiro.test(alphadte.shannon_physcia$value)
shapiro.test(alphadte.shannon_xanthoria$value)

#homogeneity of variance test
bartlett.test(value ~ temperature_zone, alphadte.shannon)
bartlett.test(value ~ lichen_species, alphadte.shannon_tempzone1)
bartlett.test(value ~ lichen_species, alphadte.shannon_tempzone2)
bartlett.test(value ~ lichen_species, alphadte.shannon_tempzone3)

bartlett.test(value ~ temperature_zone, alphadte.shannon_candelaria)
bartlett.test(value ~ temperature_zone, alphadte.shannon_physcia)
bartlett.test(value ~ temperature_zone, alphadte.shannon_xanthoria)

#Two Way ANOVA shannon Based on Temperature Zone------------------
res.aov.shannon.tempzone <- aov(value ~ uhi_range * lichen_species, data = alphadte.shannon)
summary(res.aov.shannon.tempzone)

#ANOVA shannon Based on Temperature Zone------------------
res.aov.shannon.tempzone <- aov(value ~ uhi_range, data = alphadte.shannon)
summary(res.aov.shannon.tempzone)

res.aov.shannon.tempzone_temp1 <- aov(value ~ lichen_species, data = alphadte.shannon_tempzone1)
summary(res.aov.shannon.tempzone_temp1)

TukeyHSD(res.aov.shannon.tempzone_temp1)

res.aov.shannon.tempzone_temp2 <- aov(value ~ lichen_species, data = alphadte.shannon_tempzone2)
summary(res.aov.shannon.tempzone_temp2)

TukeyHSD(res.aov.shannon.tempzone_temp2)

res.aov.shannon.tempzone_temp3 <- aov(value ~ lichen_species, data = alphadte.shannon_tempzone3)
summary(res.aov.shannon.tempzone_temp3)

TukeyHSD(res.aov.shannon.tempzone_temp3)

#ANOVA shannon Based on Lichen Species------------------
res.aov.shannon.lichen <- aov(value ~ lichen_species, data = alphadte.shannon)
summary(res.aov.shannon.lichen)

TukeyHSD(res.aov.shannon.lichen)

res.aov.shannon.lichen <- aov(value ~ uhi_range, data = alphadte.shannon_candelaria)
summary(res.aov.shannon.lichen)

res.aov.shannon.lichen <- aov(value ~ uhi_range, data = alphadte.shannon_physcia)
summary(res.aov.shannon.lichen)

res.aov.shannon.lichen <- aov(value ~ uhi_range, data = alphadte.shannon_xanthoria)
summary(res.aov.shannon.lichen)

#Kruskal shannon All Based on temperature zone -------------------
# IF the data is not normally distributed and the variance is not homogen, do Kruskal Wallis test
kruskal.test(value ~ temperature_zone, data = alphadte.shannon)

#Afterwards, you need to perform posthoc analysis, if you perform Kruskal Wallis, then do Conover test or Wilcox

#Posthoc analysis with Wilcox test
pairwise.wilcox.test(alphadte.shannon$value, alphadte.shannon$temperature_zone,
                     p.adjust.method = "BH")

#Kruskal shannon All Based on lichen_species -------------------
# IF the data is not normally distributed and the variance is not homogen, do Kruskal Wallis test
kruskal.test(value ~ lichen_species, data = alphadte.shannon)

#Afterwards, you need to perform posthoc analysis, if you perform Kruskal Wallis, then do Conover test or Wilcox

#Posthoc analysis with Wilcox test
pairwise.wilcox.test(alphadte.shannon$value, alphadte.shannon$lichen_species,
                     p.adjust.method = "BH")

#Two Way Kruskal-------------------------
library(rcompanion)

scheirerRayHare(value ~ temperature_zone + lichen_species, data = alphadte.shannon)


#-----------------------BETA DIVERSITY --------------------------------

#Exporting OTU Table --------------------------------------------------
write.csv(x = otu_table(phy_rare), file = "phyraresuball_asv_kraken.csv")
write.csv(x = otu_table(phyraresub_xanthoriaparietina), file = "phyraresub_xanthoriaparietina_kraken.csv")
write.csv(x = otu_table(phyraresub_physciaadscendens), file = "phyraresub_physciaadscendens_kraken.csv")
write.csv(x = otu_table(phyraresub_candelariaconcolor), file = "phyraresub_candelariaconcolor_kraken.csv")
write.csv(x = otu_table(phyraresubtempone), file = "phyraresubtempone_kraken.csv")
write.csv(x = otu_table(phyraresubtemptwo), file = "phyraresubtemptwo_kraken.csv")
write.csv(x = otu_table(phyraresubtempthree), file = "phyraresubtempthree_kraken.csv")

# Based temperature zone----------------------------------------------------------------------------
#all
comall_tempzone <- read.csv("phyraresuball_asv_kraken.csv", header = TRUE, row.names = 1)
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
                             lichenspecies=as.factor(group_info_all[,2]),
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

pcoa_bray_all <- ggplot(df_bray_comall, aes(x,y, shape=tempzone, fill=lichenspecies)) +
  geom_point(size=5, alpha=0.7)+
  labs(x=paste("PCoA1: ", round(re_comall$values$Relative_eig[1] * 100, 2), "%", sep=""),
       y=paste("PCoA2: ", round(re_comall$values$Relative_eig[2] * 100, 2), "%", sep=""),
       title="Bray Curtis") +
  scale_fill_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme

pcoa_bray_all

pcoa_bray_all_temp <- ggplot(df_bray_comall, aes(x,y, shape=lichenspecies, fill=tempzone)) +
  geom_point(size=5, alpha=0.7)+
  labs(x=paste("PCoA1: ", round(re_comall$values$Relative_eig[1] * 100, 2), "%", sep=""),
       y=paste("PCoA2: ", round(re_comall$values$Relative_eig[2] * 100, 2), "%", sep=""),
       title="Bray Curtis") +
  scale_fill_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme

pcoa_bray_all_temp

#Average PCoA all------------

dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data1_all <- summaryBy(x~tempzone, data=df_bray_comall, FUN=dstats)
data2_all <- summaryBy(y~tempzone, data=df_bray_comall, FUN=dstats)

data_all <- cbind(data1_all[,2:5], data2_all)
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

dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data1_alllichen <- summaryBy(x~lichenspecies, data=df_bray_comall, FUN=dstats)
data2_alllichen <- summaryBy(y~lichenspecies, data=df_bray_comall, FUN=dstats)

data_alllichen <- cbind(data1_alllichen[,3:5], data2_alllichen)
str(data_alllichen)
str(df_bray_comall)

av_pcoa_bray_all_lichen <- ggplot(data_alllichen, aes(x=x.mean, y=y.mean, fill=lichenspecies, shape=lichenspecies)) +
  geom_point(size=5,alpha=0.7) +
  geom_errorbar(aes(ymin=y.mean-y.se, ymax=y.mean+y.se)) +
  geom_errorbar(aes(xmin=x.mean-x.se, xmax=x.mean+x.se)) +
  labs(x=paste("PCoA: ", round(re_comall$values$Relative_eig[1]*100, 2), "%", sep=""),
       y=paste("PCoA: ", round(re_comall$values$Relative_eig[2]*100, 2), "%", sep=""),
       title= "Bray Curtis") +
  scale_fill_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  scale_shape_manual(values=c(21,22,23,21,21,21,21,21,21,21)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme
av_pcoa_bray_all_lichen

bray_temperature_lichen <- ggarrange(av_pcoa_bray_all_tempzone, av_pcoa_bray_all_lichen, av_pcoa_bray_all_tempzone, av_pcoa_bray_all_lichen,
                                     align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                                     common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

bray_temperature_lichen

#Permanova Bray all------------------------
str(comall_tempzone)
comall_tempzone[1:6, 1:2]

df_all <- data.frame(row.names=rownames(comall_tempzone), t(as.data.frame(strsplit(rownames(comall_tempzone),"_"))))

df_all <- plyr::rename(df_all, replace = c("X1"="tempzone", "X2"="lichenspecies", "X3"="replicates"))
head(df_all)

set.seed(123)
Permanova_bray_all <- adonis2(comall_tempzone ~ tempzone, data=df_all, method="bray", permutation=9999)
Permanova_bray_all

Pairwise_Permanova_all <- pairwise.adonis(comall_tempzone, df_all$tempzone, sim.function = "vegdist", sim.method = "bray", perm = 999)
Pairwise_Permanova_all

set.seed(123)
Permanova_bray_all <- adonis2(comall_tempzone ~ lichenspecies, data=df_all, method="bray", permutation=9999)
Permanova_bray_all

Pairwise_Permanova_all <- pairwise.adonis(comall_tempzone, df_all$lichenspecies, sim.function = "vegdist", sim.method = "bray", perm = 999)
Pairwise_Permanova_all

# Xanthoria parietina----------------------------------------------------------------------------

comxanthoriaparietina_xanthoriaparietina <- read.csv("phyraresub_xanthoriaparietina_kraken.csv", header = TRUE, row.names = 1)
comxanthoriaparietina_xanthoriaparietina <- t(comxanthoriaparietina_xanthoriaparietina)
comxanthoriaparietina_xanthoriaparietina [1:5, 1:3]


group_info_xanthoriaparietina <- data.frame(row.names=rownames(comxanthoriaparietina_xanthoriaparietina), t(as.data.frame(strsplit(rownames(comxanthoriaparietina_xanthoriaparietina),"_"))))
head(group_info_xanthoriaparietina)

#Bray Curtis  Distance xanthoriaparietina-------------
dist_comxanthoriaparietina <- vegdist(comxanthoriaparietina_xanthoriaparietina, method="bray", binary=FALSE,diag=1)
re_comxanthoriaparietina <- pcoa(dist_comxanthoriaparietina, correction="none",rn=NULL)
str(re_comxanthoriaparietina)

df_bray_comxanthoriaparietina <- data.frame(x=re_comxanthoriaparietina$vectors[,1], y=re_comxanthoriaparietina$vectors[,2],
                                            tempzone=as.factor(group_info_xanthoriaparietina[,1]),
                                            lichenspecies=as.factor(group_info_xanthoriaparietina[,2]),
                                            replicates=as.factor(group_info_xanthoriaparietina[,3]))

str(df_bray_comxanthoriaparietina)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

pcoa_bray_xanthoriaparietina <- ggplot(df_bray_comxanthoriaparietina, aes(x,y, shape=tempzone, fill=tempzone)) +
  geom_point(size=5, alpha=0.7)+
  labs(x=paste("PCoA1: ", round(re_comxanthoriaparietina$values$Relative_eig[1] * 100, 2), "%", sep=""),
       y=paste("PCoA2: ", round(re_comxanthoriaparietina$values$Relative_eig[2] * 100, 2), "%", sep=""),
       title="Bray Curtis Xanthoria parietina") +
  scale_fill_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme

pcoa_bray_xanthoriaparietina

#Average PCoA xanthoriaparietina------------

dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data1_xanthoriaparietina <- summaryBy(x~tempzone, data=df_bray_comxanthoriaparietina, FUN=dstats)
data2_xanthoriaparietina <- summaryBy(y~tempzone, data=df_bray_comxanthoriaparietina, FUN=dstats)

data_xanthoriaparietina <- cbind(data1_xanthoriaparietina[,2:5], data2_xanthoriaparietina)
str(data_xanthoriaparietina)
str(df_bray_comxanthoriaparietina)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

av_pcoa_bray_xanthoriaparietina <- ggplot(data_xanthoriaparietina, aes(x=x.mean, y=y.mean, shape=tempzone, fill=tempzone)) +
  geom_point(size=5,alpha=0.7) +
  geom_errorbar(aes(ymin=y.mean-y.se, ymax=y.mean+y.se)) +
  geom_errorbar(aes(xmin=x.mean-x.se, xmax=x.mean+x.se)) +
  labs(x=paste("PCoA: ", round(re_comxanthoriaparietina$values$Relative_eig[1]*100, 2), "%", sep=""),
       y=paste("PCoA: ", round(re_comxanthoriaparietina$values$Relative_eig[2]*100, 2), "%", sep=""),
       title= "Bray Curtis Xanthoria parietina") +
  scale_fill_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme
av_pcoa_bray_xanthoriaparietina

#Permanova Bray xanthoriaparietina------------------------
str(comxanthoriaparietina_xanthoriaparietina)
comxanthoriaparietina_xanthoriaparietina[1:6, 1:2]

df_xanthoriaparietina <- data.frame(row.names=rownames(comxanthoriaparietina_xanthoriaparietina), t(as.data.frame(strsplit(rownames(comxanthoriaparietina_xanthoriaparietina),"_"))))

df_xanthoriaparietina <- plyr::rename(df_xanthoriaparietina, replace = c("X1"="tempzone", "X2"="lichenspecies", "X3"="replicates"))
head(df_xanthoriaparietina)

set.seed(123)
Permanova_bray_xanthoriaparietina <- adonis2(comxanthoriaparietina_xanthoriaparietina ~ tempzone, data=df_xanthoriaparietina, method="bray", permutation=9999)
Permanova_bray_xanthoriaparietina

Pairwise_Permanova_xanthoriaparietina <- pairwise.adonis(comxanthoriaparietina_xanthoriaparietina, df_xanthoriaparietina$tempzone, sim.function = "vegdist", sim.method = "bray", perm = 999)
Pairwise_Permanova_xanthoriaparietina

# Physcia adscendens----------------------------------------------------------------------------

comphysciaadscendens_physciaadscendens <- read.csv("phyraresub_physciaadscendens_kraken.csv", header = TRUE, row.names = 1)
comphysciaadscendens_physciaadscendens <- t(comphysciaadscendens_physciaadscendens)
comphysciaadscendens_physciaadscendens [1:5, 1:3]


group_info_physciaadscendens <- data.frame(row.names=rownames(comphysciaadscendens_physciaadscendens), t(as.data.frame(strsplit(rownames(comphysciaadscendens_physciaadscendens),"_"))))
head(group_info_physciaadscendens)

#Bray Curtis  Distance physciaadscendens-------------
dist_comphysciaadscendens <- vegdist(comphysciaadscendens_physciaadscendens, method="bray", binary=FALSE,diag=1)
re_comphysciaadscendens <- pcoa(dist_comphysciaadscendens, correction="none",rn=NULL)
str(re_comphysciaadscendens)

df_bray_comphysciaadscendens <- data.frame(x=re_comphysciaadscendens$vectors[,1], y=re_comphysciaadscendens$vectors[,2],
                                           tempzone=as.factor(group_info_physciaadscendens[,1]),
                                           lichenspecies=as.factor(group_info_physciaadscendens[,2]),
                                           replicates=as.factor(group_info_physciaadscendens[,3]))

str(df_bray_comphysciaadscendens)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

pcoa_bray_physciaadscendens <- ggplot(df_bray_comphysciaadscendens, aes(x,y, shape=tempzone, fill=tempzone)) +
  geom_point(size=5, alpha=0.7)+
  labs(x=paste("PCoA1: ", round(re_comphysciaadscendens$values$Relative_eig[1] * 100, 2), "%", sep=""),
       y=paste("PCoA2: ", round(re_comphysciaadscendens$values$Relative_eig[2] * 100, 2), "%", sep=""),
       title="Bray Curtis Physcia adscendens") +
  scale_fill_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme

pcoa_bray_physciaadscendens

#Average PCoA physciaadscendens------------

dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data1_physciaadscendens <- summaryBy(x~tempzone, data=df_bray_comphysciaadscendens, FUN=dstats)
data2_physciaadscendens <- summaryBy(y~tempzone, data=df_bray_comphysciaadscendens, FUN=dstats)

data_physciaadscendens <- cbind(data1_physciaadscendens[,2:5], data2_physciaadscendens)
str(data_physciaadscendens)
str(df_bray_comphysciaadscendens)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

av_pcoa_bray_physciaadscendens <- ggplot(data_physciaadscendens, aes(x=x.mean, y=y.mean, shape=tempzone, fill=tempzone)) +
  geom_point(size=5,alpha=0.7) +
  geom_errorbar(aes(ymin=y.mean-y.se, ymax=y.mean+y.se)) +
  geom_errorbar(aes(xmin=x.mean-x.se, xmax=x.mean+x.se)) +
  labs(x=paste("PCoA: ", round(re_comphysciaadscendens$values$Relative_eig[1]*100, 2), "%", sep=""),
       y=paste("PCoA: ", round(re_comphysciaadscendens$values$Relative_eig[2]*100, 2), "%", sep=""),
       title= "Bray Curtis Physcia adscendens") +
  scale_fill_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme
av_pcoa_bray_physciaadscendens

#Permanova Bray physciaadscendens------------------------
str(comphysciaadscendens_physciaadscendens)
comphysciaadscendens_physciaadscendens[1:6, 1:2]

df_physciaadscendens <- data.frame(row.names=rownames(comphysciaadscendens_physciaadscendens), t(as.data.frame(strsplit(rownames(comphysciaadscendens_physciaadscendens),"_"))))

df_physciaadscendens <- plyr::rename(df_physciaadscendens, replace = c("X1"="tempzone", "X2"="lichenspecies", "X3"="replicates"))
head(df_physciaadscendens)

set.seed(123)
Permanova_bray_physciaadscendens <- adonis2(comphysciaadscendens_physciaadscendens ~ tempzone, data=df_physciaadscendens, method="bray", permutation=9999)
Permanova_bray_physciaadscendens

Pairwise_Permanova_physciaadscendens <- pairwise.adonis(comphysciaadscendens_tempzone, df_physciaadscendens$tempzone, sim.function = "vegdist", sim.method = "bray", perm = 999)
Pairwise_Permanova_physciaadscendens

# Candelaria concolor----------------------------------------------------------------------------

comcandelariaconcolor_candelariaconcolor <- read.csv("phyraresub_candelariaconcolor_kraken.csv", header = TRUE, row.names = 1)
comcandelariaconcolor_candelariaconcolor <- t(comcandelariaconcolor_candelariaconcolor)
comcandelariaconcolor_candelariaconcolor [1:5, 1:3]


group_info_candelariaconcolor <- data.frame(row.names=rownames(comcandelariaconcolor_candelariaconcolor), t(as.data.frame(strsplit(rownames(comcandelariaconcolor_candelariaconcolor),"_"))))
head(group_info_candelariaconcolor)

#Bray Curtis  Distance candelariaconcolor-------------
dist_comcandelariaconcolor <- vegdist(comcandelariaconcolor_candelariaconcolor, method="bray", binary=FALSE,diag=1)
re_comcandelariaconcolor <- pcoa(dist_comcandelariaconcolor, correction="none",rn=NULL)
str(re_comcandelariaconcolor)

df_bray_comcandelariaconcolor <- data.frame(x=re_comcandelariaconcolor$vectors[,1], y=re_comcandelariaconcolor$vectors[,2],
                                            tempzone=as.factor(group_info_candelariaconcolor[,1]),
                                            lichenspecies=as.factor(group_info_candelariaconcolor[,2]),
                                            replicates=as.factor(group_info_candelariaconcolor[,3]))

str(df_bray_comcandelariaconcolor)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

pcoa_bray_candelariaconcolor <- ggplot(df_bray_comcandelariaconcolor, aes(x,y, shape=tempzone, fill=tempzone)) +
  geom_point(size=5, alpha=0.7)+
  labs(x=paste("PCoA1: ", round(re_comcandelariaconcolor$values$Relative_eig[1] * 100, 2), "%", sep=""),
       y=paste("PCoA2: ", round(re_comcandelariaconcolor$values$Relative_eig[2] * 100, 2), "%", sep=""),
       title="Bray Curtis Candelaria concolor") +
  scale_fill_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme

pcoa_bray_candelariaconcolor

#Average PCoA candelariaconcolor------------

dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data1_candelariaconcolor <- summaryBy(x~tempzone, data=df_bray_comcandelariaconcolor, FUN=dstats)
data2_candelariaconcolor <- summaryBy(y~tempzone, data=df_bray_comcandelariaconcolor, FUN=dstats)

data_candelariaconcolor <- cbind(data1_candelariaconcolor[,2:5], data2_candelariaconcolor)
str(data_candelariaconcolor)
str(df_bray_comcandelariaconcolor)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

av_pcoa_bray_candelariaconcolor <- ggplot(data_candelariaconcolor, aes(x=x.mean, y=y.mean, shape=tempzone, fill=tempzone)) +
  geom_point(size=5,alpha=0.7) +
  geom_errorbar(aes(ymin=y.mean-y.se, ymax=y.mean+y.se)) +
  geom_errorbar(aes(xmin=x.mean-x.se, xmax=x.mean+x.se)) +
  labs(x=paste("PCoA: ", round(re_comcandelariaconcolor$values$Relative_eig[1]*100, 2), "%", sep=""),
       y=paste("PCoA: ", round(re_comcandelariaconcolor$values$Relative_eig[2]*100, 2), "%", sep=""),
       title= "Bray Curtis Candelaria concolor") +
  scale_fill_manual(values=c('#D2691E', '#9370DB', '#DC143C')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme
av_pcoa_bray_candelariaconcolor

#Permanova Bray candelariaconcolor------------------------
str(comcandelariaconcolor_candelariaconcolor)
comcandelariaconcolor_candelariaconcolor[1:6, 1:2]

df_candelariaconcolor <- data.frame(row.names=rownames(comcandelariaconcolor_candelariaconcolor), t(as.data.frame(strsplit(rownames(comcandelariaconcolor_candelariaconcolor),"_"))))

df_candelariaconcolor <- plyr::rename(df_candelariaconcolor, replace = c("X1"="tempzone", "X2"="lichenspecies", "X3"="replicates"))
head(df_candelariaconcolor)

set.seed(123)
Permanova_bray_candelariaconcolor <- adonis2(comcandelariaconcolor_candelariaconcolor ~ tempzone, data=df_candelariaconcolor, method="bray", permutation=9999)
Permanova_bray_candelariaconcolor

Pairwise_Permanova_candelariaconcolor <- pairwise.adonis(comcandelariaconcolor_tempzone, df_candelariaconcolor$tempzone, sim.function = "vegdist", sim.method = "bray", perm = 999)
Pairwise_Permanova_candelariaconcolor

bray_lichen <- ggarrange(av_pcoa_bray_xanthoriaparietina, av_pcoa_bray_physciaadscendens, av_pcoa_bray_candelariaconcolor, av_pcoa_bray_physciaadscendens, av_pcoa_bray_candelariaconcolor,
                         align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                         common.legend = TRUE, legend = "right", ncol = 3, nrow = 2)

bray_lichen

# Temperature zone 1----------------------------------------------------------------------------

comtempone_tempone <- read.csv("phyraresubtempone_kraken.csv", header = TRUE, row.names = 1)
comtempone_tempone <- t(comtempone_tempone)
comtempone_tempone [1:5, 1:3]


group_info_tempone <- data.frame(row.names=rownames(comtempone_tempone), t(as.data.frame(strsplit(rownames(comtempone_tempone),"_"))))
head(group_info_tempone)

#Bray Curtis  Distance tempone-------------
dist_comtempone <- vegdist(comtempone_tempone, method="bray", binary=FALSE,diag=1)
re_comtempone <- pcoa(dist_comtempone, correction="none",rn=NULL)
str(re_comtempone)

df_bray_comtempone <- data.frame(x=re_comtempone$vectors[,1], y=re_comtempone$vectors[,2],
                                 tempzone=as.factor(group_info_tempone[,1]),
                                 lichenspecies=as.factor(group_info_tempone[,2]),
                                 replicates=as.factor(group_info_tempone[,3]))

str(df_bray_comtempone)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

pcoa_bray_tempone <- ggplot(df_bray_comtempone, aes(x,y, shape=lichenspecies, fill=lichenspecies)) +
  geom_point(size=5, alpha=0.7)+
  labs(x=paste("PCoA1: ", round(re_comtempone$values$Relative_eig[1] * 100, 2), "%", sep=""),
       y=paste("PCoA2: ", round(re_comtempone$values$Relative_eig[2] * 100, 2), "%", sep=""),
       title="Bray Curtis temperature zone 1") +
  scale_fill_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme

pcoa_bray_tempone

#Average PCoA tempone------------

dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data1_tempone <- summaryBy(x~lichenspecies, data=df_bray_comtempone, FUN=dstats)
data2_tempone <- summaryBy(y~lichenspecies, data=df_bray_comtempone, FUN=dstats)

data_tempone <- cbind(data1_tempone[,2:5], data2_tempone)
str(data_tempone)
str(df_bray_comtempone)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

av_pcoa_bray_tempone <- ggplot(data_tempone, aes(x=x.mean, y=y.mean, shape=lichenspecies, fill=lichenspecies)) +
  geom_point(size=5,alpha=0.7) +
  geom_errorbar(aes(ymin=y.mean-y.se, ymax=y.mean+y.se)) +
  geom_errorbar(aes(xmin=x.mean-x.se, xmax=x.mean+x.se)) +
  labs(x=paste("PCoA: ", round(re_comtempone$values$Relative_eig[1]*100, 2), "%", sep=""),
       y=paste("PCoA: ", round(re_comtempone$values$Relative_eig[2]*100, 2), "%", sep=""),
       title= "Bray Curtis temperature zone 1") +
  scale_fill_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme
av_pcoa_bray_tempone

#Permanova Bray tempone------------------------
str(comtempone_tempone)
comtempone_tempone[1:9, 1:2]

df_tempone <- data.frame(row.names=rownames(comtempone_tempone), t(as.data.frame(strsplit(rownames(comtempone_tempone),"_"))))

df_tempone <- plyr::rename(df_tempone, replace = c("X1"="tempzone", "X2"="lichenspecies", "X3"="replicates"))
head(df_tempone)

set.seed(123)
Permanova_bray_tempone <- adonis2(comtempone_tempone ~ lichenspecies, data=df_tempone, method="bray", permutation=9999)
Permanova_bray_tempone

Pairwise_Permanova_tempone <- pairwise.adonis(comtempone_tempone, df_tempone$lichenspecies, sim.function = "vegdist", sim.method = "bray", perm = 999)
Pairwise_Permanova_tempone

# Temperature zone 2----------------------------------------------------------------------------

comtemptwo_temptwo <- read.csv("phyraresubtemptwo_kraken.csv", header = TRUE, row.names = 1)
comtemptwo_temptwo <- t(comtemptwo_temptwo)
comtemptwo_temptwo [1:9, 1:3]


group_info_temptwo <- data.frame(row.names=rownames(comtemptwo_temptwo), t(as.data.frame(strsplit(rownames(comtemptwo_temptwo),"_"))))
head(group_info_temptwo)

#Bray Curtis  Distance temptwo-------------
dist_comtemptwo <- vegdist(comtemptwo_temptwo, method="bray", binary=FALSE,diag=1)
re_comtemptwo <- pcoa(dist_comtemptwo, correction="none",rn=NULL)
str(re_comtemptwo)

df_bray_comtemptwo <- data.frame(x=re_comtemptwo$vectors[,1], y=re_comtemptwo$vectors[,2],
                                 tempzone=as.factor(group_info_temptwo[,1]),
                                 lichenspecies=as.factor(group_info_temptwo[,2]),
                                 replicates=as.factor(group_info_temptwo[,3]))

str(df_bray_comtemptwo)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

pcoa_bray_temptwo <- ggplot(df_bray_comtemptwo, aes(x,y, shape=lichenspecies, fill=lichenspecies)) +
  geom_point(size=5, alpha=0.7)+
  labs(x=paste("PCoA1: ", round(re_comtemptwo$values$Relative_eig[1] * 100, 2), "%", sep=""),
       y=paste("PCoA2: ", round(re_comtemptwo$values$Relative_eig[2] * 100, 2), "%", sep=""),
       title="Bray Curtis temperature zone 2") +
  scale_fill_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme

pcoa_bray_temptwo

#Average PCoA temptwo------------

dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data1_temptwo <- summaryBy(x~lichenspecies, data=df_bray_comtemptwo, FUN=dstats)
data2_temptwo <- summaryBy(y~lichenspecies, data=df_bray_comtemptwo, FUN=dstats)

data_temptwo <- cbind(data1_temptwo[,2:5], data2_temptwo)
str(data_temptwo)
str(df_bray_comtemptwo)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

av_pcoa_bray_temptwo <- ggplot(data_temptwo, aes(x=x.mean, y=y.mean, shape=lichenspecies, fill=lichenspecies)) +
  geom_point(size=5,alpha=0.7) +
  geom_errorbar(aes(ymin=y.mean-y.se, ymax=y.mean+y.se)) +
  geom_errorbar(aes(xmin=x.mean-x.se, xmax=x.mean+x.se)) +
  labs(x=paste("PCoA: ", round(re_comtemptwo$values$Relative_eig[1]*100, 2), "%", sep=""),
       y=paste("PCoA: ", round(re_comtemptwo$values$Relative_eig[2]*100, 2), "%", sep=""),
       title= "Bray Curtis temperature zone 2") +
  scale_fill_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme
av_pcoa_bray_temptwo

#Permanova Bray temptwo------------------------
str(comtemptwo_temptwo)
comtemptwo_temptwo[1:9, 1:2]

df_temptwo <- data.frame(row.names=rownames(comtemptwo_temptwo), t(as.data.frame(strsplit(rownames(comtemptwo_temptwo),"_"))))

df_temptwo <- plyr::rename(df_temptwo, replace = c("X1"="tempzone", "X2"="lichenspecies", "X3"="replicates"))
head(df_temptwo)

set.seed(123)
Permanova_bray_temptwo <- adonis2(comtemptwo_temptwo ~ lichenspecies, data=df_temptwo, method="bray", permutation=9999)
Permanova_bray_temptwo

Pairwise_Permanova_temptwo <- pairwise.adonis(comtemptwo_temptwo, df_temptwo$lichenspecies, sim.function = "vegdist", sim.method = "bray", perm = 999)
Pairwise_Permanova_temptwo

# Temperature zone 3----------------------------------------------------------------------------

comtempthree_tempthree <- read.csv("phyraresubtempthree_kraken.csv", header = TRUE, row.names = 1)
comtempthree_tempthree <- t(comtempthree_tempthree)
comtempthree_tempthree [1:9, 1:3]


group_info_tempthree <- data.frame(row.names=rownames(comtempthree_tempthree), t(as.data.frame(strsplit(rownames(comtempthree_tempthree),"_"))))
head(group_info_tempthree)

#Bray Curtis  Distance tempthree-------------
dist_comtempthree <- vegdist(comtempthree_tempthree, method="bray", binary=FALSE,diag=1)
re_comtempthree <- pcoa(dist_comtempthree, correction="none",rn=NULL)
str(re_comtempthree)

df_bray_comtempthree <- data.frame(x=re_comtempthree$vectors[,1], y=re_comtempthree$vectors[,2],
                                   tempzone=as.factor(group_info_tempthree[,1]),
                                   lichenspecies=as.factor(group_info_tempthree[,2]),
                                   replicates=as.factor(group_info_tempthree[,3]))

str(df_bray_comtempthree)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

pcoa_bray_tempthree <- ggplot(df_bray_comtempthree, aes(x,y, shape=lichenspecies, fill=lichenspecies)) +
  geom_point(size=5, alpha=0.7)+
  labs(x=paste("PCoA1: ", round(re_comtempthree$values$Relative_eig[1] * 100, 2), "%", sep=""),
       y=paste("PCoA2: ", round(re_comtempthree$values$Relative_eig[2] * 100, 2), "%", sep=""),
       title="Bray Curtis temperature zone 3") +
  scale_fill_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme

pcoa_bray_tempthree

#Average PCoA tempthree------------

dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data1_tempthree <- summaryBy(x~lichenspecies, data=df_bray_comtempthree, FUN=dstats)
data2_tempthree <- summaryBy(y~lichenspecies, data=df_bray_comtempthree, FUN=dstats)

data_tempthree <- cbind(data1_tempthree[,2:5], data2_tempthree)
str(data_tempthree)
str(df_bray_comtempthree)

mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

av_pcoa_bray_tempthree <- ggplot(data_tempthree, aes(x=x.mean, y=y.mean, shape=lichenspecies, fill=lichenspecies)) +
  geom_point(size=5,alpha=0.7) +
  geom_errorbar(aes(ymin=y.mean-y.se, ymax=y.mean+y.se)) +
  geom_errorbar(aes(xmin=x.mean-x.se, xmax=x.mean+x.se)) +
  labs(x=paste("PCoA: ", round(re_comtempthree$values$Relative_eig[1]*100, 2), "%", sep=""),
       y=paste("PCoA: ", round(re_comtempthree$values$Relative_eig[2]*100, 2), "%", sep=""),
       title= "Bray Curtis temperature zone 3") +
  scale_fill_manual(values=c('#009E73', '#56B4E9', '#E69F00')) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  mytheme
av_pcoa_bray_tempthree

#Permanova Bray tempthree------------------------
str(comtempthree_tempthree)
comtempthree_tempthree[1:9, 1:2]

df_tempthree <- data.frame(row.names=rownames(comtempthree_tempthree), t(as.data.frame(strsplit(rownames(comtempthree_tempthree),"_"))))

df_tempthree <- plyr::rename(df_tempthree, replace = c("X1"="tempzone", "X2"="lichenspecies", "X3"="replicates"))
head(df_tempthree)

set.seed(123)
Permanova_bray_tempthree <- adonis2(comtempthree_tempthree ~ lichenspecies, data=df_tempthree, method="bray", permutation=9999)
Permanova_bray_tempthree

library(pairwiseAdonis)
Pairwise_Permanova_tempthree <- pairwise.adonis(comtempthree_tempthree, df_tempthree$lichenspecies, sim.function = "vegdist", sim.method = "bray", perm = 999)
Pairwise_Permanova_tempthree

bray_temperature <- ggarrange(av_pcoa_bray_tempone, av_pcoa_bray_temptwo, av_pcoa_bray_tempthree, av_pcoa_bray_tempone, av_pcoa_bray_temptwo,
                              align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                              common.legend = TRUE, legend = "right", ncol = 3, nrow = 2)

bray_temperature

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

### Combine the results with group information ---- #group nya ini metadata
rc_data <- qpen$result %>%
  left_join(group %>% rownames_to_column('sample1'), 'sample1') %>% 
  left_join(group %>% rownames_to_column('sample2'), 'sample2') %>%
  filter(temperature_zone.x == temperature_zone.y & lichen_species.x == lichen_species.y)

rc_data_trial <- qpen$result %>%
  left_join(group %>% rownames_to_column('sample1'), 'sample1') %>% 
  left_join(group %>% rownames_to_column('sample2'), 'sample2') %>%
  filter(temperature_zone.x == temperature_zone.y)


plot_betanti_lichen <- rc_data %>%
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

plot_betanti_lichen

plot_rc_lichen <- rc_data %>%
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

plot_betanti_facetgrid <- ggarrange(plot_betanti_lichen, plot_betanti_tempzone, plot_rc_lichen, plot_rc_tempzone,
                                    align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                                    common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

plot_betanti_facetgrid

### Calculate relative importance (%) ----
rc_data_summary <- rc_data %>%
  group_by(process, temperature_zone.x, lichen_species.x) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(temperature_zone.x, lichen_species.x) %>%
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

# bNTI < -2 = Homogeneous selection
# bNTI > +2 = Variable selection
# abs(bNTI) < 2 & RC < -0.95 = Homogenizing dispersal
# abs(bNTI) < 2 & RC < +0.95 = Undominated
# abs(bNTI) < 2 & RC > +0.95 = Dispersal limitation

plot_process_lichen <- rc_data_summary %>%
  ggplot(aes(x = temperature_zone.x, y = pct, fill = process)) +
  geom_bar(stat = 'identity') +
  facet_grid(~lichen_species.x, scales = 'free', space = 'free') +
  xlab('Temperature Zone') + ylab("Relative importance (%)") +
  scale_fill_manual(values = process_colors) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_process_lichen

plot_process_tempzone <- rc_data_summary %>%
  ggplot(aes(x = lichen_species.x, y = pct, fill = process)) +
  geom_bar(stat = 'identity') +
  facet_grid(~temperature_zone.x, scales = 'free', space = 'free') +
  xlab('Temperature Zone') + ylab("Relative importance (%)") +
  scale_fill_manual(values = process_colors) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_process_tempzone

plot_process_methodone <- ggarrange(plot_process_lichen, plot_process_tempzone, plot_process_lichen, plot_process_tempzone,
                                    align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                                    common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

plot_process_methodone

### DESEQ ANALYSIS -------------------------------

#### candelaria vs physcia --------------------------------------

physeqbac_CP <- subset_samples(phy_rare, !is.na(lichen_species) & !lichen_species %in% c("Xanthoria parietina"))
physeqbac_CP

diagdds_CP = phyloseq_to_deseq2(physeqbac_CP, ~ lichen_species)

diagdds_CP = estimateSizeFactors(diagdds_CP, type = 'poscounts')

diagdds_CP = DESeq(diagdds_CP, test="Wald", fitType="parametric")

res_CP = results(diagdds_CP, cooksCutoff = FALSE)
head(res_CP)

alpha = 0.01
sigtab_CP = res_CP[which(res_CP$pvalue < alpha), ]
sigtab_CP = cbind(as(sigtab_CP, "data.frame"), as(tax_table(physeqbac_CP)[rownames(sigtab_CP), ], "matrix"))
head(sigtab_CP)
dim(sigtab_CP)

write.csv(x = (sigtab_CP), file = "Deseq2_CP.csv")

#### Candelaria vs Xanthoria --------------------------------------

physeqbac_CX <- subset_samples(phy_rare, !is.na(lichen_species) & !lichen_species %in% c("Physcia adscendens"))
physeqbac_CX

diagdds_CX = phyloseq_to_deseq2(physeqbac_CX, ~ lichen_species)

diagdds_CX = estimateSizeFactors(diagdds_CX, type = 'poscounts')

diagdds_CX = DESeq(diagdds_CX, test="Wald", fitType="parametric")

res_CX = results(diagdds_CX, cooksCutoff = FALSE)
head(res_CX)

alpha = 0.01
sigtab_CX = res_CX[which(res_CX$pvalue < alpha), ]
sigtab_CX = cbind(as(sigtab_CX, "data.frame"), as(tax_table(physeqbac_CX)[rownames(sigtab_CX), ], "matrix"))
head(sigtab_CX)
dim(sigtab_CX)

write.csv(x = (sigtab_CX), file = "Deseq2_CX.csv")

#### Physcia vs Xanthoria --------------------------------------

physeqbac_PX <- subset_samples(phy_rare, !is.na(lichen_species) & !lichen_species %in% c("Candelaria concolor"))
physeqbac_PX

diagdds_PX = phyloseq_to_deseq2(physeqbac_PX, ~ lichen_species)

diagdds_PX = estimateSizeFactors(diagdds_PX, type = 'poscounts')

diagdds_PX = DESeq(diagdds_PX, test="Wald", fitType="parametric")

res_PX = results(diagdds_PX, cooksCutoff = FALSE)
head(res_PX)

alpha = 0.01
sigtab_PX = res_PX[which(res_PX$pvalue < alpha), ]
sigtab_PX = cbind(as(sigtab_PX, "data.frame"), as(tax_table(physeqbac_PX)[rownames(sigtab_PX), ], "matrix"))
head(sigtab_PX)
dim(sigtab_PX)

write.csv(x = (sigtab_PX), file = "Deseq2_PX.csv")

# Circular heatmap ----

Deseq_CP <- read.csv("Deseq2_CP.csv", header = TRUE)
Deseq_CX <- read.csv("Deseq2_CX.csv", header = TRUE)
Deseq_PX <- read.csv("Deseq2_PX.csv", header = TRUE)

new_asv_code <- data.frame(Code=unique(c(Deseq_CP$Code,Deseq_CX$Code, Deseq_PX$Code)))

#Code is the name of ASV code coloumn
#there's 221 unique code

new_asv_code <- data.frame(Code=unique(c(Deseq_CP$Code,Deseq_CX$Code, Deseq_PX$Code)), 
                           ASV=paste0(paste0('B',seq(from = 1, to = 221)))) 

Deseq_CP <- Deseq_CP %>% left_join(new_asv_code, by = 'Code')
Deseq_CX <- Deseq_CX %>% left_join(new_asv_code, by = 'Code')
Deseq_PX <- Deseq_PX %>% left_join(new_asv_code, by = 'Code')

Deseq_CP <- Deseq_CP %>% 
  mutate(OTU=paste0(ASV,'-',Phylum,'-',Genus, '-',Species)) %>% 
  select(OTU, log2FoldChange) %>% 
  mutate(Sample='CP')
Deseq_CX <- Deseq_CX %>% 
  mutate(OTU=paste0(ASV,'-',Phylum,'-',Genus, '-',Species)) %>% 
  select(OTU, log2FoldChange) %>% 
  mutate(Sample='CX')
Deseq_PX <- Deseq_PX %>% 
  mutate(OTU=paste0(ASV,'-',Phylum,'-',Genus, '-',Species)) %>% 
  select(OTU, log2FoldChange) %>% 
  mutate(Sample='PX')

colnames(Deseq_CX)[2] <- "log2FoldChange"
colnames(Deseq_PX)[2] <- "log2FoldChange"

logFC_df <- Deseq_CP %>% rbind(Deseq_CX) %>% rbind(Deseq_PX) %>% 
  spread(key = 'Sample', value = 'log2FoldChange', fill = 0) %>%
  column_to_rownames('OTU')

head(logFC_df)
dim(logFC_df)

## Filter logFC <= 10

logFC_df_filt <- logFC_df[rowSums(abs(logFC_df)) > 5, ]
dim(logFC_df_filt)

write.csv(x = (logFC_df_filt), file = "logFC_df_filt.csv")
logFC_df_filt=read.csv("logFC_df_filt.csv", row.names=1, header = 1, check.names = F)

logFC_df_transpose <- t(logFC_df_filt)
mat_list <- logFC_df_transpose

dend_list <- as.dendrogram(hclust(dist(logFC_df_filt)))
dend_list <- raise.dendrogram(dend_list, 5)
plot(raise.dendrogram(dend_list, 5), horiz = F)

library(circlize)
pdf(file = 'Dendro.pdf', width = 30, height = 30)
col_fun = colorRamp2(breaks= c(-7, 0, 20), colors = c("blue","white", "red")) 

circos.clear()
circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0), gap.degree = 15) 
circos.initialize("a", xlim =c(0, nrow(logFC_df_filt)))

#adding new track for column label
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05,
             panel.fun = function(x, y) {
               for(i in seq_len(ncol(logFC_df_transpose))) {
                 circos.text(i-0.5, 0, 
                             colnames(logFC_df_transpose)[order.dendrogram(dend_list)][i], 
                             adj = c(0, 0.5), 
                             facing = "clockwise", niceFacing = TRUE,
                             cex = 1)                
               }
             })

circos.track(ylim = c(0, 3), bg.border = NA, panel.fun = function(x, y) {
  m = mat_list 
  dend = dend_list
  m2 = m[, order.dendrogram(dend)]
  col_mat = col_fun(m2)
  nr = nrow(m2)
  nc = ncol(m2)
  for(i in 1:nr) {
    circos.rect(1:nc - 1, rep(nr - i, nc), 
                1:nc, rep(nr - i + 1, nc), 
                border = col_mat[i, ], col = col_mat[i, ])
  }
  
  #adding row label
  circos.text(rep(1, 3), 1:3,
              rev(rownames(logFC_df_transpose)),
              facing = "downward", adj = c(1.1, 2), 
              cex = 1.6, font = 2) 
  
})

#Dendrogram
max_height <- attr(dend_list, "height") 
circos.track(ylim = c(0, max_height), bg.border = NA, track.height = 0.6, 
             panel.fun = function(x, y) {
               dend = dend_list
               circos.dendrogram(dend, max_height = max_height)
             })
circos.clear()

library(ComplexHeatmap)
#adding legend key
lgd_links <- Legend(at=c(-5,0,5), col_fun = col_fun, grid_height = unit(2.5, "cm"), 
                    legend_width = unit(6, "cm"), 
                    labels_gp = gpar(fontsize = 20), title_gp = gpar(font=2, fontsize = 25),
                    title_position = "topleft", title = "logFC", direction = "horizontal")
draw(lgd_links, x = unit(0.1, "npc"), y = unit(300, "mm"), 
     just = c("right", "bottom"))
circos.clear()
dev.off()

ggsave(plot = performance_venn_asv, filename = 'dendro.pdf', width = 10, height = 10)

#---------------------------------Differential Abundance over Urbanization Gradient--------------------------------

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

#------------------------------Candelaria-----------------------------------------

write.csv(x = otu_table(phyraresub_candelariaconcolor), file = "phyraresub_candelariaconcolor_heatmap.csv")

#after you download, change the name of the column from code to temperature zone_tree_replicate

otu_table_candelaria=read.csv("phyraresub_candelariaconcolor_heatmap.csv", row.names=1, header = 1, check.names = F)
head(otu_table_candelaria)

otu_table_candelaria <- otu_table_candelaria[rowSums(otu_table_candelaria) > 20 & rowSums(otu_table_candelaria > 0) >= 3, ] #minimal abundance 20 ditemukan di 5 sampel, rule of thumb nya ambil ASV yang banyak ditemukan di banyak sampel untuk mempermudah perhitungan dan visualisasi


otu_table_candelaria=as.matrix(otu_table_candelaria)

otu_table_candelaria <- make_relative(otu_table_candelaria) 

otu_table_candelaria <- na.omit(otu_table_candelaria)

otu_table_candelaria <- t(otu_table_candelaria)

write.csv(x = (otu_table_candelaria), file = "subset_significant_candelaria.csv")

#Subset Significant OTU Over Time
df_successionall_subset_candelaria <- read.csv("subset_significant_candelaria.csv", header=1, check.names = F)

for (i in 2:ncol(df_successionall_subset_candelaria)) {
  res.aov.succession.all_subset_candelaria <- aov(df_successionall_subset_candelaria[,i] ~ df_successionall_subset_candelaria$Group, data = df_successionall_subset_candelaria)
  bb <- data.frame(summary(res.aov.succession.all_subset_candelaria)[[1]])
  if (bb$Pr..F.[1] < 0.05) {
    print(colnames(df_successionall_subset_candelaria)[i])
  }
}

df_successionall_significant_candelaria <- df_successionall_subset_candelaria %>% select(Group,OTU5,OTU56,OTU88,OTU98,OTU143,OTU145,OTU195,OTU215,OTU243,OTU299,OTU359)
write.csv(x = (df_successionall_significant_candelaria), file = "df_successionall_significant_candelaria.csv")

#Make "df_successionall_significant.csv" by transposing 

library(OTUtable)

df_successionall_filter_candelaria <- read.csv("df_successionall_significant_candelaria.csv", row.names=1, header=1, check.names = F)
#df_successionall_filter_candelaria  <- filter_taxa(df_successionall_filter_candelaria, abundance = 1, persistence=12.34)
df_successionall_filter_candelaria =as.matrix(df_successionall_filter_candelaria)
df_successionall_filter_candelaria =t(df_successionall_filter_candelaria )

df_successionall_original <- read.csv("phyraresub_candelariaconcolor_heatmap.csv", row.names=1, header=1, check.names = F)
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

df_filter_candelaria=t(df_successionall_filter_candelaria)
write.csv(x = (df_filter_candelaria), file = "df_filter_candelaria.csv")

# Heatmap candelaria----

heatmap_candelaria_data <- read.csv("df_heatmap_candelaria.csv",
                                    header=T, row.names=1)
Heatmap_candelaria <- heatmap_candelaria_data[,1:11] # data
anno_one_candelaria <- heatmap_candelaria_data[,12] # annotation
anno_two_candelaria <- heatmap_candelaria_data[,13] # annotation
anno_one_candelaria <- factor(anno_one_candelaria, levels = c("temperaturezone1","temperaturezone2","temperaturezone3"))
annotation_col = data.frame(
  temperature_zone = factor(anno_one_candelaria, levels = c("temperaturezone1","temperaturezone2","temperaturezone3")),
  uhi_level = anno_two_candelaria)
heatmap_candelaria_plot <- pheatmap::pheatmap(t(Heatmap_candelaria), 
                                              cluster_row=T,     
                                              cluster_col=F,         
                                              fontsize_col=5,        
                                              fontsize_row=5,         
                                              show_rownames = T,      
                                              show_colnames = F,     
                                              cellwidth=5,
                                              cellheight=5,
                                              color = colorRampPalette(rev(c("#4169E1","#00BFFF","#F8F8FF")))(100),
                                              border_color='grey60',
                                              scale = 'none',
                                              main = '', 
                                              fontsize = 5, 
                                              annotation_col = annotation_col,
                                              gaps_row = c(3,6,9),
                                              treeheight_row = 100
)

heatmap_candelaria_plot

pdf('Heatmap_all.pdf', width = 35, height = 150)
grid::grid.newpage() #kalo kita pengen plot tapi dibikin kosong dulu
grid::grid.draw(heatmap_deseq$gtable) #gtable itu maksudnya struktur data image nya, memvisualisasikan data table menjadi image 
dev.off() #untuk menutup yang dipanggil

#------------------------------Physcia-----------------------------------------

write.csv(x = otu_table(phyraresub_physciaadscendens), file = "phyraresub_physciaadscendens.csv")

#after you download, change the name of the column from code to temperature zone_tree_replicate

otu_table_physcia=read.csv("phyraresub_physciaadscendens_kraken.csv", row.names=1, header = 1, check.names = F)
head(otu_table_physcia)

otu_table_physcia <- otu_table_physcia[rowSums(otu_table_physcia) > 20 & rowSums(otu_table_physcia > 0) >= 3, ] #minimal abundance 20 ditemukan di 5 sampel, rule of thumb nya ambil ASV yang banyak ditemukan di banyak sampel untuk mempermudah perhitungan dan visualisasi


otu_table_physcia=as.matrix(otu_table_physcia)

otu_table_physcia <- make_relative(otu_table_physcia) 

otu_table_physcia <- na.omit(otu_table_physcia)

otu_table_physcia <- t(otu_table_physcia)

write.csv(x = (otu_table_physcia), file = "subset_significant_physcia.csv")

#Subset Significant OTU Over Time
df_successionall_subset_physcia <- read.csv("subset_significant_physcia.csv", header=1, check.names = F)

for (i in 2:ncol(df_successionall_subset_physcia)) {
  res.aov.succession.all_subset_physcia <- aov(df_successionall_subset_physcia[,i] ~ df_successionall_subset_physcia$Group, data = df_successionall_subset_physcia)
  bb <- data.frame(summary(res.aov.succession.all_subset_physcia)[[1]])
  if (bb$Pr..F.[1] < 0.01) {
    print(colnames(df_successionall_subset_physcia)[i])
  }
}

df_successionall_significant_physcia <- df_successionall_subset_physcia %>% select(Group,OTU249, OTU251)
write.csv(x = (df_successionall_significant_physcia), file = "df_successionall_significant_physcia.csv")

#Make "df_successionall_significant.csv" by transposing 

library(OTUtable)

df_successionall_filter_physcia <- read.csv("df_successionall_significant_physcia.csv", row.names=1, header=1, check.names = F)
#df_successionall_filter_physcia  <- filter_taxa(df_successionall_filter_physcia, abundance = 1, persistence=12.34)
df_successionall_filter_physcia =as.matrix(df_successionall_filter_physcia)
df_successionall_filter_physcia =t(df_successionall_filter_physcia )

df_successionall_original <- read.csv("phyraresub_physciaadscendens_heatmap.csv", row.names=1, header=1, check.names = F)
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

df_filter_physcia=t(df_successionall_filter_physcia)
write.csv(x = (df_filter_physcia), file = "df_filter_physcia.csv")

# Heatmap physcia----

heatmap_physcia_data <- read.csv("df_heatmap_physcia.csv",
                                 header=T, row.names=1)
Heatmap_physcia <- heatmap_physcia_data[,1:2] # data
anno_one_physcia <- heatmap_physcia_data[,3] # annotation
anno_two_physcia <- heatmap_physcia_data[,4] # annotation
anno_one_physcia <- factor(anno_one_physcia, levels = c("temperaturezone1","temperaturezone2","temperaturezone3"))
annotation_col_physcia = data.frame(
  temperature_zone = factor(anno_one_physcia, levels = c("temperaturezone1","temperaturezone2","temperaturezone3")),
  uhi_level = anno_two_physcia)
heatmap_physcia_plot <- pheatmap::pheatmap(t(Heatmap_physcia), 
                                           cluster_row=T,     
                                           cluster_col=F,         
                                           fontsize_col=5,        
                                           fontsize_row=5,         
                                           show_rownames = T,      
                                           show_colnames = F,     
                                           cellwidth=5,
                                           cellheight=5,
                                           color = colorRampPalette(rev(c("#4169E1","#00BFFF","#F8F8FF")))(100),
                                           border_color='grey60',
                                           scale = 'none',
                                           main = '', 
                                           fontsize = 5, 
                                           annotation_col = annotation_col_physcia,
                                           gaps_row = c(3,6,9),
                                           treeheight_row = 100
)

#pdf('Heatmap_all.pdf', width = 35, height = 150)
#grid::grid.newpage() #kalo kita pengen plot tapi dibikin kosong dulu
#grid::grid.draw(heatmap_deseq$gtable) #gtable itu maksudnya struktur data image nya, memvisualisasikan data table menjadi image 
dev.off() #untuk menutup yang dipanggil

#------------------------------Xanthoria-----------------------------------------

write.csv(x = otu_table(phyraresub_xanthoriaparietina), file = "phyraresub_xanthoriaparietina_heatmap.csv")

#after you download, change the name of the column from code to temperature zone_tree_replicate

otu_table_xanthoria=read.csv("phyraresub_xanthoriaparietina_kraken.csv", row.names=1, header = 1, check.names = F)
head(otu_table_xanthoria)

otu_table_xanthoria <- otu_table_xanthoria[rowSums(otu_table_xanthoria) > 20 & rowSums(otu_table_xanthoria > 0) >= 3, ] #minimal abundance 20 ditemukan di 5 sampel, rule of thumb nya ambil ASV yang banyak ditemukan di banyak sampel untuk mempermudah perhitungan dan visualisasi


otu_table_xanthoria=as.matrix(otu_table_xanthoria)

otu_table_xanthoria <- make_relative(otu_table_xanthoria) 

otu_table_xanthoria <- na.omit(otu_table_xanthoria)

otu_table_xanthoria <- t(otu_table_xanthoria)

write.csv(x = (otu_table_xanthoria), file = "subset_significant_xanthoria.csv")

#Subset Significant OTU Over Time
df_successionall_subset_xanthoria <- read.csv("subset_significant_xanthoria.csv", header=1, check.names = F)

for (i in 2:ncol(df_successionall_subset_xanthoria)) {
  res.aov.succession.all_subset_xanthoria <- aov(df_successionall_subset_xanthoria[,i] ~ df_successionall_subset_xanthoria$Group, data = df_successionall_subset_xanthoria)
  bb <- data.frame(summary(res.aov.succession.all_subset_xanthoria)[[1]])
  if (bb$Pr..F.[1] < 0.01) {
    print(colnames(df_successionall_subset_xanthoria)[i])
  }
}

df_successionall_significant_xanthoria <- df_successionall_subset_xanthoria %>% select(Group,OTU2,OTU269,OTU312,OTU314,OTU315)
write.csv(x = (df_successionall_significant_xanthoria), file = "df_successionall_significant_xanthoria.csv")

#Make "df_successionall_significant.csv" by transposing 

library(OTUtable)

df_successionall_filter_xanthoria <- read.csv("df_successionall_significant_xanthoria.csv", row.names=1, header=1, check.names = F)
#df_successionall_filter_xanthoria  <- filter_taxa(df_successionall_filter_xanthoria, abundance = 1, persistence=12.34)
df_successionall_filter_xanthoria =as.matrix(df_successionall_filter_xanthoria)
df_successionall_filter_xanthoria =t(df_successionall_filter_xanthoria )

df_successionall_original <- read.csv("phyraresub_xanthoriaparietina_heatmap.csv", row.names=1, header=1, check.names = F)
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
#Procrustes Sum of Squares (m12 squared):        0.3125 
#Correlation in a symmetric Procrustes rotation: 0.8292 
#Significance:  1e-04 

#Permutation: free
#Number of permutations: 9999

plot(vare.proc.all)
plot(vare.proc.all, kind=2)
residuals(vare.proc.all)


#Put the taxonomy Phylum on excel file
#Create another excel file from filter table and replace the OTUCode to OTUID

df_filter_xanthoria=t(df_successionall_filter_xanthoria)
write.csv(x = (df_filter_xanthoria), file = "df_filter_xanthoria.csv")

# Heatmap xanthoria----

heatmap_xanthoria_data <- read.csv("df_xanthoria_heatmap.csv",
                                   header=T, row.names=1)
Heatmap_xanthoria <- heatmap_xanthoria_data[,1:5] # data
anno_one_xanthoria <- heatmap_xanthoria_data[,6] # annotation
anno_two_xanthoria <- heatmap_xanthoria_data[,7] # annotation
anno_one_xanthoria <- factor(anno_one_xanthoria, levels = c("temperaturezone1","temperaturezone2","temperaturezone3"))
annotation_col_xanthoria = data.frame(
  temperature_zone = factor(anno_one_xanthoria, levels = c("temperaturezone1","temperaturezone2","temperaturezone3")),
  uhi_level = anno_two_xanthoria)
heatmap_xanthoria_plot <- pheatmap::pheatmap(t(Heatmap_xanthoria), 
                                             cluster_row=T,     
                                             cluster_col=F,         
                                             fontsize_col=5,        
                                             fontsize_row=5,         
                                             show_rownames = T,      
                                             show_colnames = F,     
                                             cellwidth=5,
                                             cellheight=5,
                                             color = colorRampPalette(rev(c("#4169E1","#00BFFF","#F8F8FF")))(100),
                                             border_color='grey60',
                                             scale = 'none',
                                             main = '', 
                                             fontsize = 5, 
                                             annotation_col = annotation_col_xanthoria,
                                             gaps_row = c(3,6,9),
                                             treeheight_row = 100
)

pdf('Heatmap_all.pdf', width = 35, height = 150)
grid::grid.newpage() #kalo kita pengen plot tapi dibikin kosong dulu
grid::grid.draw(heatmap_deseq$gtable) #gtable itu maksudnya struktur data image nya, memvisualisasikan data table menjadi image 
dev.off() #untuk menutup yang dipanggil

plot_heatmap <- ggarrange(heatmap_candelaria_plot, heatmap_physcia_plot, heatmap_xanthoria_plot,heatmap_candelaria_plot, heatmap_physcia_plot, heatmap_xanthoria_plot,
                          align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                          common.legend = FALSE, legend = "right", ncol = 3, nrow = 2)

plot_heatmap

#Taxa Bar Plot ------------------------
data.phyla.all <- phy_rare %>%
  tax_glom(taxrank = "Species") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>% 
  filter(Abundance > 0.01) %>%
  arrange(Phylum) 
data.phyla.all

data.phyla.all$Class

#data.phyla.all$Treatment1 <- factor(data.phyla.all$Treatment1, levels=c("Soil_Rice_24D", "Soil_Rice", "Soil_24D", "Soil"))

#data.phyla.all[data.phyla.all == "Unidentified"] <- "Zunidentified"

phylum_colors <- c("Ascomycota" = "pink", "Basidiobolomycota" = "blue", "Armatimonadetes" = "beige",
                   "Basidiomycota" = "gold", "Chlamydiae" = "grey", "Glomeromycota" = "red",
                   "Mortierellomycota" = "purple", "Mucoromycota" = "cyan", "Gemmatimonadetes" = "purple",
                   "Fungi_phy_Incertae_sedis" = "lavender", "OD1" = "salmon", "Planctomycetes" = "brown", 
                   "Zoopagomycota" = "darkgreen", "TM7"= "skyblue", "Verrucomicrobia" = "violet",
                   "Unassigned" = "black")

class_colors <- c("Lecanoromycetes" = "pink", "Candelariomycetes" = "blue", "Eurotiomycetes" = "beige",
                  "Agaricomycetes" = "gold", "Arthoniomycetes" = "grey", "Glomeromycota" = "red",
                  "Sordariomycetes" = "purple", "Mucoromycota" = "cyan", "Dothideomycetes" = "purple",
                  "Rhizophydiomycetes" = "lavender", "Leotiomycetes" = "salmon", "Tremellomycetes" = "brown", 
                  "Microbotryomycetes" = "darkgreen", "Fungi_cls_Incertae_sedis"= "skyblue", "Basidiomycota_cls_Incertae_sedis" = "violet",
                  "Unassigned" = "black")

write.csv(x = (data.phyla.all), file = "data.phyla.all.csv")
data.phyla.all=read.csv("data.phyla.all.csv", row.names=1, header = T)

m=ggplot(data.phyla.all, aes(x = lichen_species, y = Abundance, fill = Class)) +
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

m_treatment=ggplot(data.phyla.all, aes(x = temperature_zone, y = Abundance, fill = Class)) +
  facet_grid(.~lichen_species, scales="free_x") +
  geom_bar(stat = "identity",position="fill") + # to equal 1
  scale_fill_manual(values = class_colors) +
  theme(text = element_text(size=15),axis.text.y = element_text(hjust=1)) +
  theme(strip.text = element_text(face="bold", size=18)) +
  scale_x_discrete(drop=TRUE) +
  theme(axis.title.x = element_text(size=12)) + # Remove x axis title
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.5)) + 
  theme(legend.position = "right")

m_treatment

plot_taxa <- ggarrange(m_treatment, m_treatment, m_treatment, m_treatment,align = "v", font.label = list(size = 14, face = "bold", family="serif"), 
                       common.legend = FALSE, legend = "right", ncol = 2, nrow = 2)

plot_taxa

############################## CO-OCCURENCE NETWORK #####################################

#--------- Input data for Candelaria concolor ----
candelariaconcolor <- subset_samples(phy_rare, lichen_species == 'Candelaria concolor')
candelariaconcolor

candelariaconcolor_df <- as.data.frame(otu_table(candelariaconcolor))
candelariaconcolor_df_filt <- candelariaconcolor_df[rowSums(candelariaconcolor_df) > 20 & rowSums(candelariaconcolor_df > 0) >= 3, ] #minimal abundance 20 ditemukan di 3 sampel, rule of thumb nya ambil ASV yang banyak ditemukan di banyak sampel untuk mempermudah perhitungan dan visualisasi
candelariaconcolor_df_filt <- candelariaconcolor_df_filt %>% rownames_to_column('#OTU ID')
write.table(candelariaconcolor_df_filt, 'candelariaconcolor.tsv', quote = F, sep = '\t', col.names = T, row.names = F) 

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table candelariaconcolor.tsv --correlation candelariaconcolor_cor.tsv --covariance candelariaconcolor_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table candelariaconcolor.tsv --number 1000 --prefix bootstrap_counts/candelariaconcolor -t 4
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
fs_Rcandelaria_filt <- fs_Rcandelaria_combine  %>% filter(abs(cor) > 0.8 & pval < 0.01)

# Export to Cytoscape ----
fs_Rcandelaria_filt %>%
  write.table('./candelariaconcolor_net2Cys.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df <- otu_table(candelariaconcolor) %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df[,-1])

otu_annot_df_candelaria <- otu_df %>%
  filter(OTU %in% unique(c(fs_Rcandelaria_filt $otu_1, fs_Rcandelaria_filt $otu2)))

otu_annot_df_candelaria$mean_abun <- rowMeans(otu_annot_df_candelaria[, -1])
otu_annot_df_mean_candelaria <- otu_annot_df_candelaria[,c(1,ncol(otu_annot_df_candelaria))]

tax_df <- tax_table(phy_rare) %>% 
  data.frame() %>%
  rownames_to_column('OTU')

tax_annot_df_candelaria <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rcandelaria_filt $otu_1, fs_Rcandelaria_filt $otu2)))

tax_annot_df_candelaria %>%
  left_join(otu_annot_df_mean_candelaria) %>%
  write.table('candelariaconcolor_net2Cys_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#--------- Input data for Physcia adscendens ----
physciaadscendens <- subset_samples(phy_rare, lichen_species == 'Physcia adscendens')
physciaadscendens

physciaadscendens_df <- as.data.frame(otu_table(physciaadscendens))
physciaadscendens_df_filt <- physciaadscendens_df[rowSums(physciaadscendens_df) > 20 & rowSums(physciaadscendens_df > 0) >= 3, ] #minimal abundance 20 ditemukan di 3 sampel, rule of thumb nya ambil ASV yang banyak ditemukan di banyak sampel untuk mempermudah perhitungan dan visualisasi
physciaadscendens_df_filt <- physciaadscendens_df_filt %>% rownames_to_column('#OTU ID')
write.table(physciaadscendens_df_filt, 'physciaadscendens.tsv', quote = F, sep = '\t', col.names = T, row.names = F) 

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table physciaadscendens.tsv --correlation physciaadscendens_cor.tsv --covariance physciaadscendens_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table physciaadscendens.tsv --number 1000 --prefix bootstrap_counts/physciaadscendens -t 4
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
fs_Rphyscia_filt <- fs_Rphyscia_combine  %>% filter(abs(cor) > 0.8 & pval < 0.01)

# Export to Cytoscape ----
fs_Rphyscia_filt %>%
  write.table('./physciaadscendens_net2Cys.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_physcia <- otu_table(physciaadscendens) %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_physcia[,-1])

otu_annot_df_physcia <- otu_df_physcia %>%
  filter(OTU %in% unique(c(fs_Rphyscia_filt $otu_1, fs_Rphyscia_filt $otu2)))

otu_annot_df_physcia$mean_abun <- rowMeans(otu_annot_df_physcia[, -1])
otu_annot_df_mean_physcia <- otu_annot_df_physcia[,c(1,ncol(otu_annot_df_physcia))]

tax_df <- tax_table(phy_rare) %>% 
  data.frame() %>%
  rownames_to_column('OTU')

tax_annot_df_physcia <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rphyscia_filt $otu_1, fs_Rphyscia_filt $otu2)))

tax_annot_df_physcia %>%
  left_join(otu_annot_df_mean_physcia) %>%
  write.table('physciaadscendens_net2Cys_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#--------- Input data for Xanthoria parietina ----
xanthoriaparietina <- subset_samples(phy_rare, lichen_species == 'Xanthoria parietina')
xanthoriaparietina

xanthoriaparietina_df <- as.data.frame(otu_table(xanthoriaparietina))
xanthoriaparietina_df_filt <- xanthoriaparietina_df[rowSums(xanthoriaparietina_df) > 20 & rowSums(xanthoriaparietina_df > 0) >= 3, ] #minimal abundance 20 ditemukan di 3 sampel, rule of thumb nya ambil ASV yang banyak ditemukan di banyak sampel untuk mempermudah perhitungan dan visualisasi
xanthoriaparietina_df_filt <- xanthoriaparietina_df_filt %>% rownames_to_column('#OTU ID')
write.table(xanthoriaparietina_df_filt, 'xanthoriaparietina.tsv', quote = F, sep = '\t', col.names = T, row.names = F) 

#Next Fastpar in HPC
# conda activate network_analyses

# fastspar --otu_table xanthoriaparietina.tsv --correlation xanthoriaparietina_cor.tsv --covariance xanthoriaparietina_cov.tsv -t 4
# mkdir bootstrap_counts
# fastspar_bootstrap --otu_table xanthoriaparietina.tsv --number 1000 --prefix bootstrap_counts/xanthoriaparietina -t 4
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
fs_Rxanthoria_filt <- fs_Rxanthoria_combine  %>% filter(abs(cor) > 0.8 & pval < 0.01)

# Export to Cytoscape ----
fs_Rxanthoria_filt %>%
  write.table('./xanthoriaparietina_net2Cys.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

#Setelah ini buka cytoscape

## Create annotation file ----
otu_df_xanthoria <- otu_table(xanthoriaparietina) %>% 
  data.frame(check.names = FALSE) %>%
  apply(2,function(x)x/sum(x)) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column('OTU')

colSums(otu_df_xanthoria[,-1])

otu_annot_df_xanthoria <- otu_df_xanthoria %>%
  filter(OTU %in% unique(c(fs_Rxanthoria_filt $otu_1, fs_Rxanthoria_filt $otu2)))

otu_annot_df_xanthoria$mean_abun <- rowMeans(otu_annot_df_xanthoria[, -1])
otu_annot_df_mean_xanthoria <- otu_annot_df_xanthoria[,c(1,ncol(otu_annot_df_xanthoria))]

tax_df <- tax_table(phy_rare) %>% 
  data.frame() %>%
  rownames_to_column('OTU')

tax_annot_df_xanthoria <- tax_df %>%
  filter(OTU %in% unique(c(fs_Rxanthoria_filt $otu_1, fs_Rxanthoria_filt $otu2)))

tax_annot_df_xanthoria %>%
  left_join(otu_annot_df_mean_xanthoria) %>%
  write.table('xanthoriaparietina_net2Cys_annot_abundance.tsv', 
              col.names = T, row.names = F, quote = F, sep = '\t')

bactotraits=read.csv("bactotraits.csv", sep=";", header = T)
write.table(bactotraits, 'bactotraits.tsv', quote = F, sep = '\t', col.names = T, row.names = F) 

####### Pie chart lichen fungi ##########

fungifunction_colors <- c("animal parasite"="#DBFFD6",
                          "litter saprotroph"="#ECD4FF",
                          "mycoparasite"="orange",
                          "soil saprotroph"="#ACE7FF",
                          "plant pathogen"="cyan",
                          "unknown" = "beige",
                          "wood saprotroph" = "#ffa2c5",
                          "lichen parasite"= "#3498db",
                          "ectomycorrhizal" = "#33ffbe",
                          "epiphyte" = "#ff3333",
                          "foliar_endophyte" = "#842704",
                          "lichenized" = "#f4fb00",
                          "sooty mold" = "#a700fb",
                          "foliar endophyte" = "#E7E8E8",
                          "unspecified saprotroph" = "#dab821")

piechart_fungi_lichen=read.csv("deseq_fungi_lichen_function.csv", row.names=1, header = T, sep=";")

piechart_cp_fungilichen <- subset(piechart_fungi_lichen, code == 'CP')
piechart_pc_fungilichen <- subset(piechart_fungi_lichen, code == 'PC')
piechart_cx_fungilichen <- subset(piechart_fungi_lichen, code == 'CX')
piechart_xc_fungilichen <- subset(piechart_fungi_lichen, code == 'XC')
piechart_xp_fungilichen <- subset(piechart_fungi_lichen, code == 'XP')
piechart_px_fungilichen <- subset(piechart_fungi_lichen, code == 'PX')

preference_counts_cp <- table(piechart_cp_fungilichen$guild)
preference_counts_pc <- table(piechart_pc_fungilichen$guild)
preference_counts_cx <- table(piechart_cx_fungilichen$guild)
preference_counts_xc <- table(piechart_xc_fungilichen$guild)
preference_counts_xp <- table(piechart_xp_fungilichen$guild)
preference_counts_px <- table(piechart_px_fungilichen$guild)

print (preference_counts_px)

piechart_all <- ggplot(piechart_fungi_lichen, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(label=abundance), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  facet_grid(.~code, scales="free_x") +
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
grpiechart_cp_fungilichen <- ggplot(piechart_cp_fungilichen, aes(x="", fill=guild, y=abundance))+
  geom_bar(stat="identity", colour="white", linewidth=1, width=1)+
  geom_text(aes(label=abundance), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  scale_fill_manual(values = fungifunction_colors) +
  guides(fill = guide_legend(title = "Guild")) +
  #scale_y_continuous(breaks = piechart_group1_fungilichen$abundance, labels = piechart_group1_fungilichen$guild) +
  theme_bw()+
  theme(text = element_text(size=12),legend.position = "right",
        legend.text = element_text(size=12),
        axis.title=element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 0.5, vjust=0.5),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grpiechart_cp_fungilichen

grpiechart_pc_fungilichen <- ggplot(piechart_pc_fungilichen, aes(x="", fill=guild, y=abundance))+
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

grpiechart_pc_fungilichen

grpiechart_cx_fungilichen <- ggplot(piechart_cx_fungilichen, aes(x="", fill=guild, y=abundance))+
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

grpiechart_cx_fungilichen

grpiechart_xc_fungilichen <- ggplot(piechart_xc_fungilichen, aes(x="", fill=guild, y=abundance))+
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

grpiechart_xc_fungilichen

grpiechart_xp_fungilichen <- ggplot(piechart_xp_fungilichen, aes(x="", fill=guild, y=abundance))+
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

grpiechart_xp_fungilichen

grpiechart_px_fungilichen <- ggplot(piechart_px_fungilichen, aes(x="", fill=guild, y=abundance))+
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

grpiechart_px_fungilichen

plot_piechart <- ggarrange(grpiechart_cp_fungilichen, grpiechart_pc_fungilichen, grpiechart_cx_fungilichen, grpiechart_xc_fungilichen, grpiechart_xp_fungilichen, grpiechart_px_fungilichen,
                           align = "v", font.label = list(size = 12, face = "bold", family="serif"), 
                           common.legend = TRUE, legend = "right", ncol = 3, nrow = 2)

plot_piechart

