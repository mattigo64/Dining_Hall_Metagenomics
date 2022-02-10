library("phyloseq")
library("ggplot2")
library("vegan")
library("DESeq2")
library(readxl)
install.packages("rJava", configure.args="--disable-jri")
#library(venneuler)
library(lme4)
install.packages('gplots')
library(gplots)


#setwd("~/Genomics/Figures.7-7-20")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")



OTU <- read_excel("~/Dropbox/Clear Labs Rutgers University collaboration/Clear Labs Data Extraction - Censored.xlsx", 
                   sheet = "OTU_Test")

TAX <- read_excel("~/Dropbox/Clear Labs Rutgers University collaboration/Clear Labs Data Extraction - Censored.xlsx", 
                  sheet = "Tax_Test")

OTU1=as.matrix(OTU)
rownames(OTU1) <- paste0("OTU", 1:nrow(OTU1))

TAX1 = as.matrix(TAX)
rownames(TAX1) <- rownames(OTU1)

library("phyloseq")
OTU1 = otu_table(OTU1, taxa_are_rows = TRUE)
TAX1 = tax_table(TAX1)

physeq = phyloseq(OTU1, TAX1)

install.packages("ape", dependencies = T)
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))


Sample <- read_excel("~/Dropbox/Clear Labs Rutgers University collaboration/Clear Labs Data Extraction - Censored.xlsx", 
                     sheet = "Samples_Test")

Sample=as.matrix(Sample)
rownames(Sample) = sample_names(physeq)
Sample=as.data.frame(Sample)
sampledata = sample_data(Sample)
sampledata

ps.1 = merge_phyloseq(physeq, sampledata,random_tree)
ps.species = tax_glom(ps.1, taxrank="Genus", NArm=FALSE)



####Everything
jpeg("barplot.jpeg", units="in", width=20, height=20, res=1000)
plot_bar(ps.species, fill = "Order")+theme(legend.text = element_text(size=10),legend.key.size = unit(0.2,"line"),
                                           axis.text.x = element_text(size=10))
dev.off()

#By Category
jpeg("barplot_category.jpeg", units="in", width=20, height=20, res=1000)
plot_bar(ps.species, fill = "Genus") + facet_wrap(~Category, scales= "free_x", nrow=1)+theme(legend.text = element_text(size=3),legend.key.size = unit(0.1,"line"),                                                                                             axis.text.x = element_text(size=5))
dev.off()

#By Category and Holding
jpeg("barplot_category_holding.jpeg", units="in", width=9, height=9, res=1000)
plot_bar(ps.species, fill = "Genus",x='Category') + facet_wrap(~Holding, scales= "free_x", nrow=1)
dev.off()

#Top 10 OTU Plot
top10OTU.names = names(sort(taxa_sums(ps.species), TRUE)[1:20])
top10OTU = prune_taxa(top10OTU.names, ps.species)

jpeg("barplot_top10.jpeg", units="in", width=40, height=40, res=1000)
plot_bar(top10OTU, fill = "Genus") +theme(legend.text = element_text(size=20),legend.key.size = unit(1,"line"),
                                                                                      axis.text.x = element_text(size=20))
dev.off()


#Tree Plot
jpeg("treeplot.jpeg", units="in", width=9, height=9, res=1000)
plot_tree(ps.species, label.tips="Genus", ladderize=T, plot.margin=.1, text.size=1.2,size='abundance',color='Category')
?plot_tree
dev.off()

#HeatMap Plot by Genus
jpeg("heatmap.jpeg", units="in", width=20, height=20, res=1000)
plot_heatmap(ps.species, ylab="Genus",taxa.label="Genus", sample.order='Category',taxa.order="Order")+theme(axis.text.x = element_text(size=6), axis.text.y=element_text(size=6))
dev.off()

#Plot Top 10 HeatMaps
jpeg("heatmap_top10.jpeg", units="in", width=20, height=20, res=1000)
plot_heatmap(top10OTU, ylab="Genus",taxa.label="Genus", sample.order='Category',taxa.order="Order")+theme(axis.text.x = element_text(size=9), axis.text.y=element_text(size=21))
dev.off()

plot_heatmap(top10OTU, ylab="Genus",taxa.label="Genus",taxa.order="Order")+theme(axis.text.x = element_text(size=9), axis.text.y=element_text(size=21))

#Alpha Diversity
jpeg("alpha.category.jpeg", units="in", width=9, height=9, res=1000)
plot_richness(ps.species, x="Category", measures=c("Shannon",'Simpson'))+geom_boxplot()
dev.off()

##Alpha.TPC
jpeg("alpha.TPC.jpeg", units="in", width=9, height=9, res=1000)
plot_richness(ps.species, x="TPC.round", measures=c("Shannon",'Simpson','Chao'))+geom_boxplot()
dev.off()

#Beta Diversity
wunifrac_dist = phyloseq::distance(ps.species, method="unifrac", weighted=F)
ordination = ordinate(ps.species, method="PCoA", distance=wunifrac_dist)

jpeg("beta.category.jpeg", units="in", width=9, height=9, res=1000)
plot_ordination(ps.species, ordination, color="Category") + theme(aspect.ratio=1)
dev.off()


jpeg("beta.ingredient.jpeg", units="in", width=9, height=9, res=1000)
plot1 = plot_ordination(ps.species, ordination, color="Ingredient.Amount") + theme(aspect.ratio=1)
plot1 + 
  stat_ellipse(type = "t")
dev.off()

jpeg("beta.TPC.jpeg", units="in", width=9, height=9, res=1000)
plot1 = plot_ordination(ps.species, ordination, color="TPC.round") + theme(aspect.ratio=1)
dev.off()


adonis(wunifrac_dist ~ sample_data(ps.species)$"TPC.round")

rich=estimate_richness(ps.species,measures=c('Simpson','Shannon'))
pairwise.wilcox.test(rich$Simpson, sample_data(ps.species)$"TPC.round",p.adjust.method = 'bonferroni')

##########ANOVA Simpson
library(car)
aov1 = Anova(lm(rich$Shannon~Sample$Category))
aov1 = aov(rich$Shannon~Sample$Category)
TukeyHSD(aov1)
Anova(lm(rich$Shannon~Sample$Holding))
Anova(lm(rich$Shannon~Sample$`Ingredient Amount`))
Anova(lm(rich$Shannon~Sample$`Dining hall`))
Anova(lm(rich$Shannon~Sample$TPC.round))
Anova(lm(rich$Shannon~Sample$Coli.round))
mod1 = lm(rich$Shannon~Sample$Category+
            Sample$TPC.round+
            Sample$`Dining hall`+
            Sample$Holding + 
            Sample$`Ingredient Amount`+
            Sample$Coli.round
            )
summary(mod1)
summary(Anova(mod1))
##############ANOVA Shannon
Anova(lm(rich$Simpson~Sample$Category))
aov1 = aov(rich$Simpson~Sample$Category)
TukeyHSD(aov1)

Anova(lm(rich$Simpson~Sample$Holding))
Anova(lm(rich$Simpson~Sample$`Ingredient Amount`))
Anova(lm(rich$Simpson~Sample$`Dining hall`))
Anova(lm(rich$Simpson~Sample$TPC.round))
aov1 = aov(rich$Shannon~Sample$TPC.round)
TukeyHSD(aov1)

Anova(lm(rich$Simpson~Sample$Coli.round))


############PRODUCE
ps.species.produce = subset_samples(ps.species, Category=="Produce")
ps.species.produce = tax_glom(ps.species.produce, taxrank="Genus", NArm=FALSE)

#sample_variables(ps.species.produce1)



#library(tidyverse)
#library(forcats)

#plot_bar(ps.species.produce1, fill = "Genus")+theme(legend.text = element_text(size=5),legend.key.size = unit(0.2,"line"),
                                                  #  axis.text.x = element_text(size=5)

#plot_bar(ps.species.produce1, fill = "Genus") + facet_wrap(~Dining.hall, scales= "free_x", nrow=1)+theme(legend.text = element_text(size=5),legend.key.size = unit(0.2,"line"),
               
                                                                                                                                                                                               #    axis.text.x = element_text(size=5))
#ps.species.produce.2 = subset_samples(ps.1, Category=="Produce")
#ps.species.produce.2 = tax_glom(ps.species.produce.2, taxrank="Genus", NArm=FALSE)


#top20OTU.names = names(sort(taxa_sums(ps.species), TRUE)[1:20])
#ps.species.produce.2 = prune_taxa(top20OTU.names, ps.species.produce.2)
#top20OTU=subset_samples(ps.species.produce.2, Category=="Produce")



top10OTU.names = names(sort(taxa_sums(ps.species), TRUE)[1:10])
top10OTU = prune_taxa(top10OTU.names, ps.species.produce)
top20OTU.produce=subset_samples(top10OTU, Category=="Produce")

jpeg("produce.20.jpeg", units="in", width=20, height=20, res=1000)
plot_bar(top20OTU.produce, fill = "Genus")+theme(
                                                 axis.text.x = element_text(size=10),
                                                 strip.text.x = element_text(size=20),
                                                 axis.title.y=element_text(size=24),
                                                 axis.title.x=element_text(size=24),
                                                 legend.key.size = unit(3,"line"),
                                                 legend.text = element_text(size=20),
                                                  legend.title=element_text(size=22) )
dev.off()

jpeg("produce.TPC.jpeg", units="in", width=22, height=22, res=1000)
plot_bar(ps.species.produce, fill = "Genus")+ facet_wrap(~TPC.round, scales= "free_x", nrow=1)+theme(legend.position = 'none',
                                                                                                      axis.text.x = element_text(size=14,face='bold'),
                                                                                                     strip.text.x = element_text(size=20),
                                                                                                     axis.title.y=element_text(size=24),
                                                                                                     axis.title.x=element_text(size=24))
dev.off()

jpeg("produce.Coli.jpeg", units="in", width=22, height=22, res=1000)
plot_bar(ps.species.produce, fill = "Genus")+ facet_wrap(~Coli.round, scales= "free_x", nrow=1)+theme(legend.position = 'none',
                                                                                                    axis.text.x = element_text(size=5),
                                                                                                    strip.text.x = element_text(size=20))
dev.off()


jpeg("produce.20.heat.jpeg", units="in", width=22, height=22, res=1000)
plot_heatmap(top20OTU.produce, ylab="Genus",taxa.label="Genus",taxa.order="Order")+theme(axis.text.x = element_text(size=12, face='bold'), axis.text.y=element_text(size=30), axis.title.y=element_text(size=24),axis.title.x=element_text(size=24))
dev.off()

plot_heatmap(top20OTU.produce, ylab="Genus",taxa.label="Genus",taxa.order="Order")+theme(axis.text.x = element_text(size=12, face='bold'), axis.text.y=element_text(size=30), axis.title.y=element_text(size=24),axis.title.x=element_text(size=24))


####################Else
ps.species.else = subset_samples(ps.1, Category!="Produce")
#else.table=phyloseq(OTU1, TAX1)
#ps.species.else = merge_phyloseq(else.table, ps.species.else)
ps.species.else = tax_glom(ps.species.else, taxrank="Genus", NArm=FALSE)



#plot_bar(ps.species.else, fill = "Genus")+ facet_wrap(~Category, scales= "free_x", nrow=1)+theme(legend.text = element_text(size=5),legend.key.size = unit(0.2,"line")
#)


top10OTU.names = names(sort(taxa_sums(ps.species), TRUE)[1:10])
top10OTU = prune_taxa(top10OTU.names, ps.species)
top20OTU.else = subset_samples(top10OTU, Category!="Produce")


jpeg("else.20.jpeg", units="in", width=22, height=22, res=1000)
plot_bar(top20OTU.else, fill = "Genus")+ facet_wrap(~Category, scales= "free_x", nrow=1)+theme(axis.text.x = element_text(size=10),
                                                                                               strip.text.x = element_text(size=20),
                                                                                               axis.title.y=element_text(size=24),
                                                                                               axis.title.x=element_text(size=24),
                                                                                               legend.key.size = unit(3,"line"),
                                                                                               legend.text = element_text(size=20),
                                                                                               legend.title=element_text(size=22) )
dev.off()


jpeg("else.TPC.jpeg", units="in", width=22, height=22, res=1000)
plot_bar(ps.species.else, fill = "Genus")+ facet_wrap(~TPC.round, scales= "free_x", nrow=1)+theme(legend.position = 'none',
                                                                                                  axis.text.x = element_text(size=5),
                                                                                                  strip.text.x = element_text(size=20))
dev.off()

plot_bar(top20OTU.else, fill = "Genus")+ facet_wrap(~Coli.round, scales= "free_x", nrow=1)+theme(legend.text = element_text(size=5),legend.key.size = unit(0.2,"line")
)



jpeg("else.20.heat.jpeg", units="in", width=22, height=22, res=1000)
plot_heatmap(top20OTU.else, ylab="Genus",taxa.label="Genus", taxa.order="Order")+theme(axis.text.x = element_text(size=12, face='bold'), axis.text.y=element_text(size=30), axis.title.y=element_text(size=24),axis.title.x=element_text(size=24))
dev.off()

##############OrderResults
ps.order = tax_glom(ps.1, taxrank="Order", NArm=FALSE)


plot_bar(ps.order, fill = "Order")+theme(legend.text = element_text(size=10),legend.key.size = unit(0.2,"line"),
                                           axis.text.x = element_text(size=10))


plot_bar(ps.order, fill = "Order") + facet_wrap(~Category, scales= "free_x", nrow=1)+theme(legend.text = element_text(size=3),legend.key.size = unit(0.1,"line"),                                                                                             axis.text.x = element_text(size=5))


ps.order.produce = subset_samples(ps.order, Category=="Produce")
ps.order.produce = tax_glom(ps.order.produce, taxrank="Order", NArm=FALSE)

plot_bar(ps.order.produce, fill = "Order")+theme(legend.text = element_text(size=10),legend.key.size = unit(0.2,"line"),
                                         axis.text.x = element_text(size=10))


row.names(OTU) = TAX$Genus
row.names(OTU)
OTU = as.numeric(OTU)
hm1 = heatmap.2(OTU)
hm1
############
#######


library("phyloseq")
OTU1 = otu_table(OTU1, taxa_are_rows = TRUE)
TAX = tax_table(TAX1)
OTU

physeq = phyloseq(OTU1, TAX)
physeq

plot_bar(physeq, fill = "Family")


library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

plot_heatmap(physeq)


ps.genus = tax_glom(physeq, taxrank="Genus", NArm=FALSE)
ps.genus =  merge_phyloseq(ps.genus, sampledata)

plot_heatmap(ps.genus, ylab="Genus",taxa.label="Genus", sample.order='Category',taxa.order="Order")

random_tree = rtree(ntaxa(ps.genus), rooted=TRUE, tip.label=taxa_names(ps.genus))


plot_bar(ps.genus, fill = "Genus")


Sample <- read_excel("~/Dropbox/Clear Labs Rutgers University collaboration/Clear Labs Data Extraction - Censored.xlsx", 
                  sheet = "Samples_Test")

Sample=as.matrix(Sample)
rownames(Sample) = sample_names(physeq)
Sample=as.data.frame(Sample)
sampledata = sample_data(Sample)
sampledata

physeq1 = merge_phyloseq(ps.genus, sampledata,random_tree)

plot_bar(physeq1, fill="Genus") + facet_wrap(~Holding, scales= "free_x", nrow=1)
plot_tree(physeq1, label.tips="Genus", ladderize=T, plot.margin=0.3)


random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq1 = merge_phyloseq(physeq, sampledata,random_tree)
ps.species = tax_glom(physeq1, taxrank="Species", NArm=FALSE)


plot_richness(ps.species, x="Holding",measures = c('Shannon'))
#rich = estimate_richness(physeq1)
#pairwise.wilcox.test(rich$Shannon, sample_data(ps.species)$Category)

#Beta Diversity
wunifrac_dist = phyloseq::distance(ps.species, method="unifrac", weighted=F)
ordination = ordinate(ps.species, method="PCoA", distance=wunifrac_dist)
plot1= plot_ordination(ps.species, ordination, color="Ingredient.Amount") + theme(aspect.ratio=1)

adonis(wunifrac_dist ~ sample_data(ps.species)$Holding)
library("ggplot2")
plot1 + 
  stat_ellipse(type = "t") +
  theme_bw()





wunifrac_dist = phyloseq::distance(ps.species, method="unifrac", weighted=F)
ordination = ordinate(ps.species, method="PCoA", distance=wunifrac_dist)

jpeg("beta.ingredient.jpeg", units="in", width=9, height=9, res=1000)
plot_ordination(ps.species, ordination, color="Ingredient.Amount") + theme(aspect.ratio=1)
dev.off()




