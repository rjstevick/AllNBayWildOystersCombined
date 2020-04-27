# 16S data from TNC and PJ projects combined
# Analyzed with QIIME
# RJS 3/24/2019

# Figures V2 to V5.

library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(phyloseq)
library(plyr)
library(vegan)
library(openxlsx)
library(ggpubr)


## phylum ----------------
data<-read.xlsx("TNC_PJ_V6_Phyladata_edit.xlsx")
data_gut<-filter(data, SampleType=="gut" & Treatment=="Control")
datag<-gather(data, taxa, value, "Acidobacteria":"Unknown")
# use only gut and control data at first.
datag<-filter(datag, SampleType=="gut" & Treatment=="Control")

#phylum top 20
datatop20<-read.xlsx("TNC_PJ_V6_Phyladata_edit.xlsx", sheet="top20")
datagtopg<-gather(datatop20, taxa, value, "Proteobacteria":"Others")
datatopg<-filter(datagtopg, SampleType=="gut" & Treatment=="Control")
datatopg$taxa = factor(datatopg$taxa, levels = unique(datatopg$taxa))

# define palette
palette<-c("#8fe1c7",
           "#4e736a",
           "#803435",
           "#48211f",
           "#cdd99e",
           "#2b3826",
           "#d8ad93",
           "#4c652d",
           "#c27a7b",
           "#61a178",
           "#7f572f",
           "#72becb",
           "#a79459")


# Number of reads per sample
abund <-ggplot(datag,aes(SampleName,value,fill=Station))+geom_bar(stat="identity")+
  facet_grid(.~SampleType+Station, scales="free")+
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.length = unit(0.6, 'lines'))+
  labs(y="Number of reads \nper sample",x=NULL)+
  scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc"), 
                    labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
                            "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
                            "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond"))

ggplot(data_gut,aes(SampleName,DermoConc,fill=DermoIndex))+geom_bar(stat="identity")+
  facet_grid(.~SampleType+Station, scales="free")+
  scale_y_log10()+
  labs(y="Dermo Concentration \n(cells/oyster)",x=NULL)+theme_grey()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.length = unit(0.3, 'lines'))+
  scale_fill_manual(values=c("red","yellow","orange","transparent"))


# Percent abundances of top taxa
percent <- ggplot(datatopg,aes(SampleName,value,fill=taxa))+
  geom_col(position="fill")+
  theme(legend.text = element_text(size=12, colour="gray20"),
        legend.position = "none", axis.text.x = element_blank(),
        axis.text.y=element_text(angle=90, hjust=0.5),
        strip.text = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.6, 'lines'))+
  facet_grid(.~TypeStation, scales="free")+
  scale_fill_manual(values=palette)+labs(y="Percent \nabundance",x=NULL,fill=NULL)+
  scale_y_continuous(labels = scales::percent)

# Figure V-3
cowplot::plot_grid(abund,percent,nrow=2,align = "hv", rel_heights=c(20,50))


# heatmap of top 20 phyla
datatopg_norm<-ddply(datatopg,.(SampleName),transform,rescale=sqrt(value))

ggplot(datatopg_norm,aes(taxa,SampleName,fill=rescale))+
  geom_tile()+
  theme(legend.text = element_text(size=12, colour="gray20"),
        axis.text.y = element_blank(),axis.text.x = element_text(angle=40, hjust=1),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0.2, 'lines'))+
  facet_grid(SampleType+Station~., scales="free")+
  labs(y="Percent \nabundance",x=NULL,fill=NULL)+
  scale_fill_gradientn(labels = scales::percent,limits=c(0,1),colours=c("white","#fecc5c","#fd8d3c","#f03b20","#bd0026","darkred"))



#order ---------------------------------------
datao<-read.xlsx("TNC_PJ_V6_Orderdata_edit.xlsx", sheet="norm")
dataotab<-datao[,2:389]
row.names(dataotab)<-datao$SAMPLEINDEX
datao_gut<-filter(datao, SampleType=="gut" & Treatment=="Control")
dataog<-gather(datao, taxa, value, "Acidobacteria;AT-s3-28;uncultured.actinobacterium":"Unknown;Unknown;Unknown")

# diversity --------------
dataog$Simpsons<-diversity(dataotab, index="simpson")
dataog<-filter(dataog, SampleType=="gut" | SampleType=="water" & Treatment=="Control")
dataog_gut<-filter(dataog, SampleType=="gut")
dataog_water<-filter(dataog, SampleType=="water")

ggplot(dataog,aes(x=Station,y=Simpsons, fill=Station))+
  geom_boxplot(aes(fill=factor(Station)))+
  facet_grid(~SampleType, labeller=as_labeller(c(gut="Gut samples (n=5)", water="Water samples (n=2)"))) + 
  labs(x=NULL,y="Simpson's Index of Diversity",fill="Station")+
  scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc"), 
                    labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
                            "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
                            "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond"))+
  stat_compare_means(method = "anova", label.y = 0.2, label.x="3.BIS") +        # Add global annova p-value
  theme_bw()+scale_y_continuous(limits=c(0,1))


ggplot(dataog_gut,aes(x=Station,y=Simpsons))+geom_boxplot(aes(fill=factor(Station)))+
  ggtitle("Simpson's Index of diversity for gut samples")+
  labs(x=NULL,y="Simpson's Index of Diversity",fill="Station")+
  scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc"), 
                    labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
                             "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
                             "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond"))
  
ggplot(dataog_water,aes(x=Station,y=Simpsons))+geom_boxplot(aes(fill=factor(Station)))+
  ggtitle("Simpson's Index of diversity for water samples")+
  scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc"), 
                    labels=c("1. Providence River - Bold Point Park", "2. Greenwich Bay - Goddard State Park",
                             "3. Bissel Cove", "4. Narrow River", "5. Point Judith North - Billington Cove", 
                             "6. Point Judith South - Bluff Hill Cove","7. Ninigret Pond"))


# heatmap -----------------------------------

library(reshape2)
library(plyr)
library(scales)
library(viridis)
theme_set(theme_bw())

# use top 40 orders
datao_gut_top40<-datao_gut[,names(sort(colSums(datao_gut[,2:389]), decreasing = TRUE))]
rownames(datao_gut_top40)<-datao_gut$SampleName
datao_gut_top40_merge<-c(datao_gut_top40[,1:40],datao_gut[,390:443])
datao_gut_top40_merge<-as.data.frame(datao_gut_top40_merge)

dataog_gut_top40<-gather(datao_gut_top40_merge, taxa, value, "Unknown.Unknown.Unknown":"Cyanobacteria.Oxyphotobacteria.Nostocales")
dataog_gut_norm<-ddply(dataog_gut_top40,.(SampleName),transform,rescale=sqrt(value))
dataog_gut_norm$taxa = factor(dataog_gut_norm$taxa, levels = rev(unique(dataog_gut_norm$taxa)))


# Figure V-4
# Plot heatmap
ggplot(dataog_gut_norm, aes(SampleName, taxa)) + 
  geom_tile(aes(fill = rescale),colour = "white") + 
  facet_grid(~SampleType+Station,space="free", scales="free")+
#  scale_fill_gradient(na.value = "grey50",low="white",high="red")+
  scale_fill_gradientn(labels = scales::percent,limits=c(0,1),colours=c("white","#fecc5c","#fd8d3c","#f03b20","#bd0026","darkred"))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "bottom",axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, colour="grey40"),
        legend.text = element_text(size=10, colour="grey40"),
        legend.key.size = unit(2, 'lines'))+
  labs(y="Most abundant taxa (Order level)", x=NULL, fill="Relative\nPercent\nAbundance")

pheatmap::pheatmap(t(datao_gut_top40[,1:40]), scale="row")



## PCA with Environmental Parameters, order level -------------------------

library(gridExtra)
library(ggfortify)

pcadataraw <- read.xlsx("TNC_PJ_V6_Orderdata_edit.xlsx", sheet="norm-pca")

pcadatameta<-pcadataraw[1:43,]
pcadata<-as.data.frame(pcadataraw)
rownames(pcadata)<-c(pcadata[,1])
# pcadata<-pcadata[,2:14] # with mass, width, length

pcadataprune<-pcadata[1:43,2:288] # remove random rows
# add in diversity
pcadataprune$BacDiversity<-diversity(dataotab[1:43,], index="simpson")

pca_result <- prcomp(pcadataprune, scale=TRUE)
pca_result$rotation <- -pca_result$rotation
pca_result$rotation

pca_result$x <- - pca_result$x
autoplot(pca_result)


theme_set(theme_grey())
aplot<-autoplot(pca_result,size=6, data=pcadatameta,fill="Station",shape="Station",loadings=FALSE,loadings.label=FALSE,
                frame=TRUE,frame.colour = "Station",level = 0.95,
                loadings.label.fontface="bold",
                loadings.label.repel=TRUE,loadings.label.size=3.5, 
                loadings.label.colour="black",
                loadings.colour="darkgreen")
aplot + scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc"), 
                          labels=c("1. Providence River - Bold Point Park", "2. Greenwich Bay - Goddard State Park",
                                   "3. Bissel Cove", "4. Narrow River", "5. Point Judith North - Billington Cove", 
                                   "6. Point Judith South - Bluff Hill Cove","7. Ninigret Pond"))+
        scale_color_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc"), 
                          labels=c("1. Providence River - Bold Point Park", "2. Greenwich Bay - Goddard State Park",
                                   "3. Bissel Cove", "4. Narrow River", "5. Point Judith North - Billington Cove", 
                                   "6. Point Judith South - Bluff Hill Cove","7. Ninigret Pond"))+
        scale_shape_manual(values = c(21,22,23,24,25,21, 22),
                     labels=c("1. Providence River - Bold Point Park", "2. Greenwich Bay - Goddard State Park",
                              "3. Bissel Cove", "4. Narrow River", "5. Point Judith North - Billington Cove", 
                              "6. Point Judith South - Bluff Hill Cove","7. Ninigret Pond"))
  





# With Phyla-level bacterial data -----------------
# top 20 phyla abundant in 6 or more samples.
pcadataraw <- read.xlsx("TNC_PJ_V6_Phyladata_edit.xlsx", sheet="norm-strip")
theme_set(theme_grey())

pcadatameta<-pcadataraw[1:43,]
pcadata<-as.data.frame(pcadataraw)

pcadataprune<-pcadata[1:43,c(2:22,33:40)]
# add in diversity
pcadataprune$BacterialDiversity<-diversity(pcadata[1:43,c(2:22)], index="simpson")

pca_result <- prcomp(pcadataprune, scale=TRUE)
pca_result$rotation <- -pca_result$rotation
pca_result$x <- - pca_result$x


# Figure V-5
aplot<-autoplot(pca_result,size=6, data=pcadatameta,fill="Station",shape="Station",
                frame=TRUE,frame.colour = "Station",level = 0.95,
                loadings=TRUE,loadings.label=TRUE,frame.alpha=0.3,
                loadings.label.repel=TRUE,loadings.label.size=4, 
                loadings.label.colour=c("black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black",
                                  "magenta4","red","red","red","red","red","red","red","blue"),
                loadings.colour=c("black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black",
                                  "magenta4","red","red","red","red","red","red","red","blue"))
aplot + scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc"), 
                          labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
                                   "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
                                   "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond"))+
        scale_color_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc"), 
                           labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
                                    "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
                                    "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond"))+
        scale_shape_manual(values = c(21,22,23,24,25,21,22),
                labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
                         "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
                         "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond"))+
        theme(legend.position="bottom", legend.text = element_text(size=12, color="grey20"))
