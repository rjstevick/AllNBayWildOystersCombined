# 16S data from RWU, TNC and PJ projects combined
# Analyzed with QIIME
# RJS 3/24/2019

# Figures V-6 to V-8

library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(phyloseq)
library(plyr)
library(vegan)
library(openxlsx)
library(ggpubr)
library(ggbiplot)



## phylum ----------------
data<-read.xlsx("RWU_TNC_PJ_16S_Phyladata.xlsx")
datag<-gather(data, taxa, value, "Acidobacteria":"Unknown")
# use only gut and control data at first.
#datag<-filter(datag, SampleType=="gut" & Treatment=="Control")

#phylum top 30
datatop30<-read.xlsx("RWU_TNC_PJ_16S_Phyladata.xlsx", sheet="norm-strip")
datagtopg<-gather(datatop30, taxa, value, "Acidobacteria":"Unknown")
#datatopg<-filter(datagtopg, SampleType=="gut" & Treatment=="Control")
datagtopg$taxa = factor(datagtopg$taxa, levels = unique(datagtopg$taxa))



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


datagtopgavg<-datagtopg %>%group_by(taxa,Treatment,BucketDay, OysterStage,SpecSampleType,StationTrial,Treatment) %>%
  dplyr::summarise(meanvalue=mean(value)) %>%
  ungroup()

# Percent abundances of top taxa, averaged per rep.
# Figure V-7.
ggplot(arrange(datagtopgavg,taxa),aes(Treatment,meanvalue,fill=taxa))+
  geom_col(position="fill")+
  theme(legend.text = element_text(size=12, colour="gray20"),
        legend.position = "bottom", axis.text.x = element_blank(),
   #     axis.text.y=element_text(angle=90, hjust=0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.6, 'lines'),
        panel.spacing = unit(0.01, "lines"))+
  facet_grid(.~OysterStage+SpecSampleType+StationTrial+BucketDay, scales="free", space="free")+
  scale_fill_manual(values=c(paletteall))+labs(y="Percent \nabundance",x=NULL,fill=NULL)+
  scale_y_continuous(labels = scales::percent, expand = c(0,0))

ggplot(arrange(datagtopg,taxa),aes(SampleName,value,fill=taxa))+
  geom_col(position="fill")+
  theme(legend.text = element_text(size=12, colour="gray20"),
        legend.position = "bottom", axis.text.x = element_blank(),
        #     axis.text.y=element_text(angle=90, hjust=0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.6, 'lines'),
        panel.spacing = unit(0.01, "lines"))+
  facet_grid(.~OysterStage+SpecSampleType+StationTrial+Treatment, scales="free", space="free")+
  scale_fill_manual(values=c(paletteall))+labs(y="Percent \nabundance",x=NULL,fill=NULL)+
  scale_y_continuous(labels = scales::percent, expand = c(0,0))



# heatmap of top 20 phyla
datatopg_norm<-ddply(datagtopg,.(SampleName),transform,rescale=sqrt(value))
datagtopgavg_norm<-ddply(datagtopgavg,.(Group),transform,rescale=sqrt(meanvalue))

# Figure V-6.
ggplot(datatopg_norm,aes(taxa,SampleName,fill=rescale))+
  geom_tile()+coord_flip()+theme_minimal()+
  theme(legend.text = element_text(size=10, colour="gray20"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),legend.position = "bottom",
        axis.ticks.length = unit(0.2, 'lines'),
        panel.spacing = unit(0.01, "lines"))+
  facet_grid(.~OysterStage+SpecSampleType+StationTrial+Treatment, scales="free", space="free")+
  labs(fill="Percent abundance",x=NULL,y=NULL)+
  viridis::scale_fill_viridis(option="B", labels = scales::percent,limits=c(0,1))

ggplot(datagtopgavg_norm,aes(taxa,Group,fill=rescale))+
  geom_tile()+coord_flip()+theme_minimal()+
  theme(legend.text = element_text(size=12, colour="gray20"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.2, 'lines'))+
  facet_grid(.~OysterStage+SpecSampleType+StationTrial, scales="free", space="free")+
  labs(fill="Percent abundance",x=NULL,y=NULL)+
  scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0))+
  scale_fill_gradientn(labels = scales::percent,limits=c(0,1),colours=c("white","#fecc5c","#fd8d3c","#f03b20","#bd0026","darkred"))



# diversity --------------
datagtopg$Simpsons<-diversity(data[,2:30], index="simpson")

ggplot(datagtopg,aes(x=SpecSampleType,y=Simpsons, fill=StationTrial))+
  geom_boxplot(aes(fill=factor(StationTrial)))+
  facet_grid(~OysterStage, scales="free",space="free") + 
  labs(x=NULL,y="Simpson's Index of Diversity",fill="Station")+
#  scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc"), 
#                    labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
#                             "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
#                            "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond"))+
  stat_compare_means(method = "anova", label.y = 0.2, label.x="3.BIS") +        # Add global annova p-value
  theme_bw()+scale_y_continuous(limits=c(0,1))


ggplot(datagtopg,aes(x=StationTrial,y=Simpsons, fill=factor(Treatment))) +
  facet_grid(~OysterStage+SpecSampleType, scales="free",space="free") + 
  geom_point(position=position_dodge(width=0.8), color="grey60")+
  geom_boxplot(alpha=0.5)+
  labs(x=NULL,y="Simpson's Index of Diversity",fill="Station")+
  #  scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc"), 
  #                    labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
  #                             "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
  #                            "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond"))+
  theme_bw()+scale_y_continuous(limits=c(0,1))





# With Phyla-level bacterial data -----------------
pcadataraw <- read.xlsx("RWU_TNC_PJ_16S_Phyladata.xlsx", sheet="norm-strip")
theme_set(theme_grey())

pcadata<-pcadataraw[,1:29]

pca_result <- prcomp(pcadata, scale=TRUE)
pca_result$rotation <- -pca_result$rotation
pca_result$x <- - pca_result$x

aplot<-autoplot(pca_result,size=6, data=pcadataraw,color="StationTrial",
                fill="StationTrial",shape="SpecSampleType",
                frame=FALSE,frame.colour = "StationTrial",level = 0.95,
                loadings=TRUE,loadings.label=TRUE,frame.alpha=0.3,
                loadings.label.repel=TRUE,loadings.label.size=4,
                loadings.label.colour="black",loadings.colour="black")
  
aplot + scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc","black","grey30","grey80"), 
                          labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
                                   "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
                                   "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond",
                                   "1.RWU - Hatchery Trial 1","2.RWU - Hatchery Trial 2",
                                   "3.RWU - Hatchery Trial 3"))+
        scale_color_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc","black","grey40","grey10"), 
                           labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
                                    "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
                                    "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond",
                                    "1.RWU - Hatchery Trial 1","2.RWU - Hatchery Trial 2",
                                    "3.RWU - Hatchery Trial 3"))+
        scale_shape_manual(values = c(24,22,23,25,21,21,9),
                           limits=c("gut","inner swab","outer swab","LarvalOyster","water",
                                    "RearingWater","TankSwab"),
                           labels=c("Adult oyster gut","Adult inner shell",
                                    "Adult outer shell","Larval oyster","Seawater",
                                    "Hatchery Rearing Water",
                                    "Hatchery Tank Biofilm Swab"))+
        theme(legend.position="right", legend.text = element_text(size=12, color="grey20"))


  
  
  
  
# With Averaged Phyla-level bacterial data -----------------
pcadataraw <- read.xlsx("RWU_TNC_PJ_16S_Phyladata.xlsx", sheet="norm-strip")
theme_set(theme_grey())

pcadatarawg<-gather(pcadataraw, taxa, value, "Acidobacteria":"Unknown")
pcadatarawavg<-pcadatarawg %>% group_by(taxa,Treatment,SampleType,BucketDay, OysterStage,SpecSampleType,StationTrial,Treatment) %>%
  dplyr::summarise(meanvalue=mean(value)) %>%
  ungroup()
pcadatarawavgtab<-spread(pcadatarawavg, taxa,meanvalue)

pcadata<-pcadatarawavgtab[,7:35]

pca_result <- prcomp(pcadata, scale=TRUE)
pca_result$rotation <- -pca_result$rotation
pca_result$x <- - pca_result$x

# Figure V-8A
aplot<-autoplot(pca_result,x = 1,y = 2, size=6, data=pcadatarawavgtab,colour="Treatment",
                fill="StationTrial",shape="SpecSampleType", 
                frame=FALSE,frame.colour = "SpecSampleType",level = 0.95,
                loadings=TRUE,loadings.label=TRUE,frame.alpha=0.3,
                loadings.label.repel=TRUE,loadings.label.size=3.7,
                loadings.label.colour="black",loadings.colour="black")
aplot$layers[[1]]$aes_params$stroke <- 1.2
aplot + scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc",
                                    "#2b0045","#7e358a","#da70d6"), 
                          labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
                                   "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
                                   "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond",
                                   "0.RWU.1 - Hatchery Trial 1","0.RWU.2 - Hatchery Trial 2",
                                   "0.RWU.3 - Hatchery Trial 3"))+
  scale_color_manual(values=c("grey60","black","black"), 
                     labels=c("Control","Nutrient Enrichment","Probiotic Treatment"))+
  scale_shape_manual(values = c(24,22,23,25,21,21,23),
                     limits=c("gut","inner swab","outer swab","LarvalOyster","water",
                              "RearingWater","TankSwab"),
                     labels=c("Adult oyster gut","Adult inner shell",
                              "Adult outer shell","Larval oyster","Seawater",
                              "Hatchery Rearing Water",
                              "Hatchery Tank Biofilm Swab"))+
  theme(legend.position="bottom", legend.direction = "vertical", legend.text = element_text(size=12, color="grey20"))

# Figure V-8B.
bplot<-autoplot(pca_result,x = 3,y = 4, size=6, data=pcadatarawavgtab,colour="Treatment",
                fill="StationTrial",shape="SpecSampleType", 
                frame=FALSE,frame.colour = "SpecSampleType",level = 0.95,
                loadings=TRUE,loadings.label=TRUE,frame.alpha=0.3,
                loadings.label.repel=TRUE,loadings.label.size=3.7,
                loadings.label.colour="black",loadings.colour="black")
bplot$layers[[1]]$aes_params$stroke <- 1.2
bplot + scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc",
                                   "#2b0045","#7e358a","#da70d6"), 
                          labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
                                   "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
                                   "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond",
                                   "0.RWU.1 - Hatchery Trial 1","0.RWU.2 - Hatchery Trial 2",
                                   "0.RWU.3 - Hatchery Trial 3"))+
  scale_color_manual(values=c("grey60","black","black"), 
                     labels=c("Control","Nutrient Enrichment","Probiotic Treatment"))+
  scale_shape_manual(values = c(24,22,23,25,21,21,23),
                     limits=c("gut","inner swab","outer swab","LarvalOyster","water",
                              "RearingWater","TankSwab"),
                     labels=c("Adult oyster gut","Adult inner shell",
                              "Adult outer shell","Larval oyster","Seawater",
                              "Hatchery Rearing Water",
                              "Hatchery Tank Biofilm Swab"))+
  theme(legend.position="bottom", legend.direction = "vertical", legend.text = element_text(size=12, color="grey20"))


#1200x1000
cowplot::plot_grid(aplot,bplot)






# Subset of the samples

pcadatarawavg<-pcadatarawg %>% group_by(taxa,Treatment,BucketDay, OysterStage,SpecSampleType,StationTrial,Treatment) %>%
  dplyr::summarise(meanvalue=mean(value)) %>%
  ungroup() %>%
  filter(SpecSampleType!="outer swab" & SpecSampleType!="inner swab")

pcadatarawavgtab<-spread(pcadatarawavg, taxa,meanvalue)

#pcadata<-pcadatarawavgtab[,c(6,7,9:14,16:33)]
pcadata<-pcadatarawavgtab[,6:34]


pca_result <- prcomp(pcadata, scale=TRUE)
pca_result$rotation <- -pca_result$rotation
pca_result$x <- - pca_result$x

aplot<-autoplot(pca_result,size=6, data=pcadatarawavgtab,colour="Treatment",
                fill="StationTrial",shape="SpecSampleType", 
                frame=FALSE,frame.colour = "StationTrial",level = 0.95,
                loadings=TRUE,loadings.label=TRUE,frame.alpha=0.3,
                loadings.label.repel=TRUE,loadings.label.size=4,
                loadings.label.colour="black",loadings.colour="black")


aplot + scale_fill_manual(values=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc","black","mediumorchid4","plum"), 
                          labels=c("1.PVD - Bold Point Park", "2.GB - Goddard State Park",
                                   "3.BIS - Bissel Cove", "4.NAR - Narrow River", "5.PJN - Billington Cove", 
                                   "6.PJS - Bluff Hill Cove","7.NIN - Ninigret Pond",
                                   "1.RWU - Hatchery Trial 1","2.RWU - Hatchery Trial 2",
                                   "3.RWU - Hatchery Trial 3"))+
  scale_color_manual(values=c("grey60","black","black"), 
                     labels=c("Control","Nutrient Enrichment","Probiotic Treatment"))+
  scale_shape_manual(values = c(21,24,22,23,25),
                     limits=c("water","gut","LarvalOyster",
                              "RearingWater","TankSwab"),
                     labels=c("Adult gut","Seawater","Larval oyster",
                              "Hatchery Rearing Water",
                              "Hatchery Tank Biofilm Swab"))+
  theme(legend.position="right", legend.text = element_text(size=12, color="grey20"))

