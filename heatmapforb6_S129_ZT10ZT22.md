## heatmap plot for the ZT22/ZT10 in B6 and S129

## step 1: calling peaks with homer findepeaks
this part can be done with "ChIAPET peks calling from PET" 

## mkdir a new folder (not nessary) and copy peaks file into this folder for furthur analysis
```
mkdir /workspace/rsrch2/panpanliu/23101-02_06302023_173816/combined_2_times/ChIAPET.Tool.V3/heatmap_ZT10vsZT22
```
## step 2: combined these 4 samples(S1_B6_NC_ZT10, S2_B6_NC_ZT22, S5_S129NC_ZT10 and S6_S129NC_ZT22)
```
#/usr/bin/bash
for f in $(ls /workspace/rsrch2/panpanliu/23101-02_06302023_173816/combined_2_times/ChIAPET.Tool.V3/heatmap_ZT10vsZT22/peak_file/*peakcalling.txt)
do
tag_dir=`basename $f .peakcalling.txt`
cat $f | grep -v '#' | awk '{print $2"\t"$3"\t"$4"\t"$1"_'$tag_dir'"}' >>/workspace/rsrch2/panpanliu/23101-02_06302023_173816/combined_2_times/ChIAPET.Tool.V3/heatmap_ZT10vsZT22/B6S129_ZT10ZT22.combined.peak
done
```
## step 3: sort and merge combined peak file
```
sort -k1,1 -k2,2n B6S129_ZT10ZT22.combined.peak | mergeBed -i stdin -d 10 -o collapse -c 4 >B6S129_ZT10ZT22.merged.peak
```
## step 4: annotated peaks with homer
```
annotatePeaks.pl B6S129_ZT10ZT22.merged.peak mm10 -d ../S1/Histone_peakscalling/S1_selfligand_TagDir/ ../S2/Histone_peakscalling/S2_selfligand_TagDir/ \
../S5/Histone_peakscalling/S5_selfligand_TagDir/ ../S6/Histone_peakscalling/S6_selfligand_TagDir/ -size 1000 > B6S129_ZT10ZT22.mergedpeaks.annotate.txt
```

## step 4: seek for the peaks for B6 ZT10vsZT22 >2 with R
/workspace/rsrch2/panpanliu/23101-02_06302023_173816/combined_2_times/ChIAPET.Tool.V3/heatmap_ZT10vsZT22/scripts/filte_peaks_withcondition.R
``rm(list = ls())
wd="/workspace/rsrch2/panpanliu/23101-02_06302023_173816/combined_2_times/ChIAPET.Tool.V3/heatmap_ZT10vsZT22"
setwd(wd)

peak_ann<-read.delim("B6S129_ZT10ZT22.mergedpeaks.annotate.txt",stringsAsFactors = F)
colnames(peak_ann)[1]="PeakID"
colnames(peak_ann)[20:23]=c("S1_B6NC_ZT10","S2_B6NC_ZT22","S5_S129NC_ZT10","S6_S129NC_ZT22")

dim(peak_ann)
##70784    23
library(dplyr)
##filter peaks without gene.name
peak_ann<-filter(peak_ann,Gene.Name !="")

#filter with lowest tag counts
tag_count=0

peak_ann <- peak_ann %>%
            filter(S1_B6NC_ZT10 > tag_count & 
                   S2_B6NC_ZT22 > tag_count &
                   S5_S129NC_ZT10 > tag_count &
                   S6_S129NC_ZT22 > tag_count)

FC=2
#############combined file#######33
peak_ZT10vsZT22_down_FC2 <- peak_ann %>%
  filter(S2_B6NC_ZT22/S1_B6NC_ZT10 <= 1/FC | S6_S129NC_ZT22/S5_S129NC_ZT10 <= 1/FC)
peak_ZT10vsZT22_add_FC4 <- peak_ann %>%
  filter(S2_B6NC_ZT22/S1_B6NC_ZT10 >=FC | S6_S129NC_ZT22/S5_S129NC_ZT10 >=FC)
write.table(peak_ZT10vsZT22_add_FC2,"B6S129_ZT10vs22_upFC2.txt",sep = "\t",quote = F,na="", row.names = F)
write.table(peak_ZT10vsZT22_down_FC2,"B6S129_ZT10vs22_downFC2.txt",sep = "\t",quote = F,na="", row.names = F)


## seperate into different file#########
B6_add_FC2 <- peak_ann %>%
                   filter(S2_B6NC_ZT22/S1_B6NC_ZT10 >=FC)
row.names(B6_add_FC2)=B6_add_FC2$PeakID

S129_add_FC2 <- peak_ann %>%
  filter(S6_S129NC_ZT22/S5_S129NC_ZT10 >=FC )
row.names(S129_add_FC2)=S129_add_FC2$PeakID

B6_down_FC2 <- peak_ann %>%
  filter(S2_B6NC_ZT22/S1_B6NC_ZT10 <= 1/FC)
row.names(B6_down_FC2)=B6_down_FC2$PeakID

S129_down_FC2 <- peak_ann %>%
  filter(S6_S129NC_ZT22/S5_S129NC_ZT10 <= 1/FC )
row.names(S129_down_FC2)=S129_down_FC2$PeakID

###############up peaks for common, specific peak#####3
common <- intersect(B6_add_FC2$PeakID,S129_add_FC2$PeakID)
length(common) ##200
B6_up_spe <-setdiff(B6_add_FC2$PeakID,S129_add_FC2$PeakID)
length(B6_up_spe) ##1066
S129_up_spe <- setdiff(S129_add_FC2$PeakID,B6_add_FC2$PeakID)
length(S129_up_spe) ## 2622

write.table(B6_add_FC2[common,c(1:4)],"./plots_file/up_Z10vs22_common.peak",sep = "\t",quote = F, row.names = F,col.names = F)
write.table(peak_ann[B6_up_spe,c(1:4)],"./plots_file/up_Z10vs22_B6_spe.peak",sep = "\t",quote = F, row.names = F,col.names = F)
write.table(peak_ann[S129_up_spe,c(1:4)],"./plots_file/up_Z10vs22_S129_spe.peak",sep = "\t",quote = F, row.names = F,col.names = F)


#######################down ZT10vsZT22 in B6 and S129 common, specific
common_down <- intersect(B6_down_FC2$PeakID,S129_down_FC2$PeakID)
length(common_down) ##33
B6_down_spe <-setdiff(B6_down_FC2$PeakID,S129_down_FC2$PeakID)
length(B6_down_spe) ##372
S129_down_spe <- setdiff(S129_down_FC2$PeakID,B6_down_FC2$PeakID)
length(S129_down_spe) ## 2259

write.table(peak_ann[common_down,c(1:4)],"./plots_file/down_Z10vs22_common.peak",sep = "\t",quote = F, row.names = F,col.names = F)
write.table(peak_ann[B6_down_spe,c(1:4)],"./plots_file/down_Z10vs22_B6_spe.peak",sep = "\t",quote = F, row.names = F,col.names = F)
write.table(peak_ann[S129_down_spe,c(1:4)],"./plots_file/down_Z10vs22_S129_spe.peak",sep = "\t",quote = F, row.names = F,col.names = F)`

```


## step 5: plot signal plot file for up/down peaks
```
cd /workspace/rsrch2/panpanliu/23101-02_06302023_173816/combined_2_times/ChIAPET.Tool.V3/heatmap_ZT10vsZT22/plots_file
annotatePeaks.pl ../B6S129_ZT10vs22_upFC2.txt mm10 -size 4000 -hist 10 -d ../../S1/Histone_peakscalling/S1_selfligand_TagDir/ ../../S2/Histone_peakscalling/S2_selfligand_TagDir/ ../../S5/Histone_peakscalling/S5_selfligand_TagDir/ ../../S6/Histone_peakscalling/S6_selfligand_TagDir/ >B6S129_ZT10vs22_FC2.4kb.signal.txt
annotatePeaks.pl ../B6S129_ZT10vs22_downFC2.txt mm10 -size 4000 -hist 10 -d ../../S1/Histone_peakscalling/S1_selfligand_TagDir/ ../../S2/Histone_peakscalling/S2_selfligand_TagDir/ ../../S5/Histone_peakscalling/S5_selfligand_TagDir/ ../../S6/Histone_peakscalling/S6_selfligand_TagDir/ >down_B6S129_ZT10vs22_FC2.4kb.signal.txt
```

and then load R to plot file
```
rm(list = ls())

wd="/workspace/rsrch2/panpanliu/23101-02_06302023_173816/combined_2_times/ChIAPET.Tool.V3/heatmap_ZT10vsZT22/plots_file/"
setwd(wd)

library(reshape2)
library(ggplot2)

data<-read.delim("B6S129_ZT10vs22_up_FC2.4kb.signal.txt",stringsAsFactors = F)

data_plot<-data[,c(1,grep("Coverage",colnames(data)))];head(data_plot)

colnames(data_plot)=c("Distance","B6_NC_ZT10","B6_NC_ZT22","Sv_NC_ZT10","Sv_NC_ZT22");head(data_plot)

data_plot<-melt(data_plot, id = "Distance", variable.name = "Condition",value.name = "Coverage")
data_plot$Condition=as.factor(data_plot$Condition)

ggplot(data_plot, aes(x=Distance, y=Coverage, group=Condition))+
    geom_line(aes(color=Condition))+
    scale_x_continuous(breaks=seq(-2000,2000,1000))+
    scale_y_continuous(breaks=seq(0,3,0.5))+
    scale_color_manual(values=c("grey", "green", "blue","red"))+
  labs(title = "signal Plotfile of ZT10vsZT22 up peaks", x="Distance from peaks",y="signal value") +
   theme(legend.position="upper", legend.direction="horizontal",
        legend.title = element_blank())+
   theme_classic()

ggsave("ZT10vsZT22_up_signal.pdf",width = 4,height = 4)
```

