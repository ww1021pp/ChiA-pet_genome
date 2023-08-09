---
attachments: [Clipboard_2023-05-22-11-18-54.png, Clipboard_2023-05-22-16-59-40.png]
title: ChiA-PET V3 tools install
created: '2023-05-22T14:20:57.210Z'
modified: '2023-07-07T16:16:13.830Z'
---

# ChiA-PET V3 tools install

## simple introduce
Chromatin Interaction Analysis with Paired-End Tag (ChIA-PET) sequencing is a technology to study genome-wide long-range chromatin interactions bound by protein factors. ChIA-PET Tool V3, a software package for automatic processing of ChIA-PET sequence data, including: 0. linker detect (--stop_step 0, only run linker detect)

    1.linker filtering
    2.mapping the paired-end reads to a reference genome
    3.purifying the mapped reads
    4.dividing the reads into different categories
    5.peak calling
    6.interaction calling
    7.visualizing the results


## the software depends on the following softwares:


    JDK>=1.8(https://www.oracle.com/technetwork/java/javase/downloads/index.html) ## conda base
    BWA(http://bio-bwa.sourceforge.net/) ## conda install BWA -c bioconda
    SAMtools(http://samtools.sourceforge.net/) ## conda install -c samtools
    BEDTools(https://bedtools.readthedocs.io/en/latest/) ## conda install BEDTools
    R(https://www.r-project.org/) ## install by apt-install R
    R package grid(install.packages("grid")) ## install.packages("grid")
    R package xtable(install.packages("xtable")) 
    R package RCircos(install.packages("RCircos"))

##  excute ChIA-PET_tool
Download the ChIA-PET Tool V3 package from https://github.com/GuoliangLi-HZAU/ChIA-PET_Tool_V3. Unpack the package using the following command in your selected directory:
$ unzip ChIA-PET_Tool_V3.zip  
 
 cd ChIA-PET_Tool_V3, the use ChIA-PET.jar to excute

 ###As the process was mapped with bwa, so need bulid index with bwa

 ```
 bwa index ref.fa 
 ```


 ## pipeline of ChIA-PET
 ![](@attachment/Clipboard_2023-05-22-16-59-40.png)


 can uss "start_step" to choose corresponding steps
