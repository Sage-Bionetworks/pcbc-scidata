library(synapseClient)
library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(tidyr)

synapseLogin()

cellMetadataTbl <- synTableQuery('select * from syn2767694')
cellMetadata <- cellMetadataTbl@values %>% filter(public, Cell_Type %in% c('PSC'))


# qr <- synapseQuery("select * from file where projectId=='syn1773109'", blockSize = 250)
# d <- qr$collectAll()
# save(d, file="syn1773109.RData")
# f <- synStore(File("syn1773109.RData", parentId="syn6174635"), forceVersion=FALSE)

f <- synGet("syn6174636")
load(getFileLocation(f))

# Only getting pluripotent stem cell public files
d2 <- d %>% filter(file.public=='true', file.Cell_Type %in% c('PSC'))

keepCols <- c("file.id", "file.UID", "file.name", "file.biologicalSampleName", "file.C4_Cell_Line_ID", 
              "file.dataType", "file.dataSubType", "file.fileType",
              "file.Donor_ID", "file.Cell_Line_Type", "file.Cell_Type", "file.Diffname_short",
              "file.C4_Karyotype_Result", "file.pass_qc", "file.exclude")

# Get raw data for each dataType
d.mRNAC4 <- d2 %>% filter(file.dataType == 'mRNA', file.fileType == 'fastq') %>%
  select(one_of(keepCols)) %>% 
  unique

d.miRNAC4 <- d2 %>% filter(file.dataType == 'miRNA', file.fileType == 'fastq') %>%
  select(one_of(keepCols)) %>% 
  unique

# Methylation raw files have two files per UID (red and green channel)
d.MethylC4 <- d2 %>% 
  filter(file.dataType == 'methylation', 
         file.fileType == 'idat') %>%
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

dC4 <- rbind(d.mRNAC4, d.miRNAC4, d.MethylC4, d.TeratomaReportC4, 
             d.TeratomaHEC4, d.TeratomaIHCC4, 
             d.karyoReportC4, d.CNVReportC4)

dC4Used <- dC4 %>% filter(#file.C4_Karyotype_Result == "normal", 
                          file.pass_qc == "true", 
                          file.exclude == "false")

isa <- dC4Used %>% 
  filter(file.dataType %in% c("mRNA", "miRNA", "methylation")) %>% 
  mutate(`Data Repository`="Synapse") %>% 
  select(file.UID, file.name, `Data Repository`, file.id, file.dataType) %>% 
  rename(`Assay Name`=file.UID, `Raw Data File`=file.name, 
         `Data Record Accession`=file.id, Method=file.dataType)

dC4Used %>% 
  group_by(file.dataType) %>% 
  summarise(n=n_distinct(file.C4_Cell_Line_ID))

dC4Used %>% 
  filter(file.dataType %in% c("mRNA", "miRNA", "methylation")) %>%
  select(file.C4_Cell_Line_ID, file.Diffname_short, file.biologicalSampleName, file.dataType) %>% 
  dcast(file.biologicalSampleName ~ file.dataType)

dC4Used %>% 
  filter(file.dataType %in% c("mRNA", "miRNA", "methylation")) %>%
  select(file.UID)

## Problematic

dC4Used %>% 
  filter(!(file.C4_Cell_Line_ID %in% c("IPS18"))) %>% 
  filter(file.dataType %in% c("mRNA", "miRNA", "methylation")) %>%
  select(file.C4_Cell_Line_ID, file.Diffname_short, file.biologicalSampleName, file.dataType) %>% 
  count(file.biologicalSampleName) %>% filter(n > 3)

dC4 %>% 
  group_by(file.dataType, file.Cell_Line_Type) %>% 
  summarise(n=n_distinct(file.C4_Cell_Line_ID)) %>% 
  reshape2::dcast(file.dataType ~ file.Cell_Line_Type) 

dC4Used %>% 
  group_by(file.dataType) %>% 
  summarise(n=n_distinct(file.C4_Cell_Line_ID)) %>% 
  reshape2::dcast(file.dataType ~ file.Cell_Line_Type) 

dC4Used %>% 
  group_by(file.Cell_Line_Type) %>% 
  summarise(n=n_distinct(file.C4_Cell_Line_ID)) 

dC4Used %>% 
  group_by(file.dataType) %>%
  select(file.C4_Cell_Line_ID, file.Diffname_short) %>% 
  unique() %>%
  summarise(n=n())
