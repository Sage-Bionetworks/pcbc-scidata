library(synapseClient)
library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(tidyr)

synapseLogin()
qr <- synapseQuery("select * from file where projectId=='syn1773109'", blockSize = 250)
d <- qr$collectAll()
save(d, file="syn1773109.RData")

load("syn1773109.RData")

d2 <- d %>% filter(file.public=='true', file.Cell_Type %in% c('PSC', 'somatic'))

keepCols <- c("file.C4_Cell_Line_ID", "file.dataType", "file.dataSubType", "file.fileType",
              "file.Donor_ID", "file.Cell_Line_Type", "file.Cell_Type", "file.Diffname_short")

d.mRNAC4 <- d2 %>% filter(file.dataType == 'mRNA', file.fileType == 'fastq') %>%
  select(one_of(keepCols)) %>% 
  unique

d.miRNAC4 <- d2 %>% filter(file.dataType == 'miRNA', file.fileType == 'fastq') %>%
  select(one_of(keepCols)) %>% 
  unique


d.MethlC4 <- d2 %>% 
  filter(file.dataType == 'methylation', 
         file.fileType == 'idat', file.Channel == 'Grn') %>%
  select(one_of(keepCols)) %>% 
  unique


d.karyoReportC4 <- d2 %>% 
  filter(file.dataType == 'Karyotype', file.fileType=='report') %>%
  select(one_of(keepCols)) %>% 
  unique

d.TeratomaIHCC4 <- d2 %>% 
  filter(file.dataType == 'Teratoma', file.dataSubType=="IHC", 
         file.fileType == 'image') %>%
  select(one_of(keepCols)) %>% 
  unique

d.TeratomaHEC4 <- d2 %>% 
  filter(file.dataType == 'Teratoma', file.dataSubType=="H and E", 
         file.fileType == 'image') %>%
  select(one_of(keepCols)) %>% 
  unique

d.TeratomaReportC4 <- d2 %>% 
  filter(file.dataType == 'Teratoma', file.fileType=="report") %>% 
  select(one_of(keepCols)) %>% 
  unique

d.CNVReportC4 <- d2 %>% 
  filter(file.dataType == 'CNV', file.fileType=="report") %>% 
  select(one_of(keepCols)) %>% 
  unique

dC4 <- rbind(d.mRNAC4, d.miRNAC4, d.MethlC4, d.TeratomaReportC4, 
             d.TeratomaHEC4, d.TeratomaIHCC4, 
             d.karyoReportC4, d.CNVReportC4)

dC4 %>% 
  group_by(file.dataType, file.dataSubType, file.fileType) %>% 
  summarise(n=n_distinct(file.C4_Cell_Line_ID))

dC4 %>% 
  group_by(file.dataType) %>% 
  summarise(n=n_distinct(file.C4_Cell_Line_ID))

dC4 %>% 
  group_by(file.dataType) %>%
  select(file.C4_Cell_Line_ID, file.Diffname_short) %>% 
  unique() %>%
  summarise(n=n())
