library(openxlsx)
library(biomaRt)
library(stringr)

#### Integrating new variants from inmegen-----
wg.Varnuevas <- read.table("data/Varnuevas_POS.csv", header = T, sep = ",", stringsAsFactors = F)[,-c(5,6)]
colnames(wg.Varnuevas)[1] <- "POS"
wg.Varnuevas$WGid <- paste(wg.Varnuevas$CHROM, wg.Varnuevas$POS, sep = "_")

wg.Consolidadas <- read.xlsx("data/listas_consolidadas3.xlsx", sheet = 1, colNames = T, skipEmptyRows = T, check.names = T)
wg.noFASTA <- droplevels(wg.Consolidadas[grep("POS",wg.Consolidadas$Locus_Name),])
wg.noFASTA$WGid <- paste(wg.noFASTA$Chromosome,wg.noFASTA$Coordinate, sep = "_")

wg.Cons_and_Var <- merge(wg.noFASTA, wg.Varnuevas[,-c(1,4)], by = "WGid")

wg.CT_Edits <- read.xlsx("data/CT_Edits_Original_Input_File.xlsx", sheet = 1, colNames = T, skipEmptyRows = T, check.names = T)
wg.CT_Edits$WGid <- paste0("chr",wg.CT_Edits$Chromosome,"_",wg.CT_Edits$Coordinate)

wg.threelist <- unique(merge(wg.Cons_and_Var, wg.CT_Edits[,c(5,19,1,11:16)], by = "WGid", all.x = T)[,-c(19)])

wg.threelist$ALT <- gsub(",","/",wg.threelist$ALT)
wg.threelist$Sequence <- paste0("N[",wg.threelist$REF,"/",wg.threelist$ALT,"]N")

#wg.threelist$Locus_Name.x <- wg.threelist$Locus_Name.y ##Genrea duplicados
wg.threelist[,9:14] <- wg.threelist[,20:ncol(wg.threelist)]

wg.nuevas <- unique(wg.threelist[,-c(1,17:ncol(wg.threelist))])
wg.nuevas$Genome_Build_Version.x <- 37
wg.nuevas$Source.x <- "new"
wg.nuevas$Source_Version.x <- 0
wg.nuevas$Sequence_Orientation.x <- "unknown"
wg.nuevas$Plus_Minus.x <- "Plus"
wg.nuevas$Species.x <- "Homo sapiens"
wg.nuevas$Force_Infinium_I <- "FALSE"
wg.nuevas$Locus_Name.x <- NA

##wg.dups <- unique(wg.nuevas[,c(3,5,6,7)]) ##Collumn 3 makes duplicate data
##wg.dups <- wg.nuevas[duplicated(wg.nuevas[,c(5:7)]),]


###MEMORY CHECKPOINT
rm(list=setdiff(ls(), "wg.nuevas"))

### ----


wg.0.Hugo_data <- read.xlsx("data/listas_consolidadas3.xlsx", sheet = 1, colNames = T, skipEmptyRows = T, check.names = T)

colnames(wg.nuevas) <- colnames(wg.0.Hugo_data)
###wg.0.Hugo_data <- wg.nuevas ##Just for checking second assignment correct way to do it is via rbind
wg.0.Hugo_data <- rbind(wg.0.Hugo_data, wg.nuevas)

#### DATA STRUCTURE EVALUATION ----
# wg.0.Hugo_data$File <- as.factor(wg.0.Hugo_data$File)
# wg.0.Hugo_data$Locus_Name_LOCAL <- as.factor(wg.0.Hugo_data$Locus_Name_LOCAL)
# wg.0.Hugo_data$Locus_Name <- as.factor(wg.0.Hugo_data$Locus_Name)
# wg.0.Hugo_data$Target_Type <- as.factor(wg.0.Hugo_data$Target_Type)
# wg.0.Hugo_data$Chromosome <- as.factor(wg.0.Hugo_data$Chromosome)
# wg.0.Hugo_data$Coordinate <- as.factor(wg.0.Hugo_data$Coordinate)
# wg.0.Hugo_data$Genome_Build_Version <- as.factor(wg.0.Hugo_data$Genome_Build_Version)
# wg.0.Hugo_data$Source <- as.factor(wg.0.Hugo_data$Source)
# wg.0.Hugo_data$Source_Version <- as.factor(wg.0.Hugo_data$Source_Version)
# wg.0.Hugo_data$Sequence_Orientation <- as.factor(wg.0.Hugo_data$Sequence_Orientation)
# wg.0.Hugo_data$Plus_Minus <- as.factor(wg.0.Hugo_data$Plus_Minus)
# wg.0.Hugo_data$Species <- as.factor(wg.0.Hugo_data$Species)
# wg.0.Hugo_data$Force_Infinium_I <- as.factor(wg.0.Hugo_data$Force_Infinium_I)
# wg.0.Hugo_data$Confidentiality_status <- as.factor(wg.0.Hugo_data$Confidentiality_status)
# 
# summary(wg.0.Hugo_data)


### Removal of rows without valid fasta sequence

wg.noFASTA <- droplevels(wg.0.Hugo_data[is.na(wg.0.Hugo_data$Sequence),])
wg.noFASTA <- droplevels(rbind(wg.noFASTA, wg.0.Hugo_data[grep(">:",wg.0.Hugo_data$Sequence),] ))


wg.1.work_data <- droplevels( subset(wg.0.Hugo_data, !is.na(Sequence) ))
wg.1.work_data <- droplevels( wg.1.work_data[grep(">:",wg.1.work_data$Sequence, invert = T),] )

wg.1.work_data <- unique.data.frame(wg.1.work_data) ##Removes duplicate lines


###MEMORY CHECKPOINT
rm(wg.noFASTA)

###### ----  
### DATA REFORMATING ----

#### Resolver problema de mala antacion en columna Coordinate
wg.1a.work_data_good_coordinate <- wg.1.work_data[ grep("rs", wg.1.work_data$Coordinate, invert = T), ]
wg.1b.work_data_bad_coordinate <- wg.1.work_data[ grep("rs", wg.1.work_data$Coordinate), ]

##Reshuffling bad collumns
# wg.1c.work_data_bad_coordinate <- wg.1b.work_data_bad_coordinate
# wg.1c.work_data_bad_coordinate$Locus_Name_LOCAL <- NA
# wg.1c.work_data_bad_coordinate$Coordinate <- wg.Final_bad_coordinate$Locus_Name_LOCAL

wg.1b.work_data_bad_coordinate$Coordinate <- wg.1b.work_data_bad_coordinate$Locus_Name_LOCAL
wg.1b.work_data_bad_coordinate$Locus_Name_LOCAL <- NA

##Extracting only coordinate

wg.1b.work_data_bad_coordinate$Coordinate <- as.character(wg.1b.work_data_bad_coordinate$Coordinate)

for (i in 1:nrow(wg.1b.work_data_bad_coordinate) ) {
  wg.1b.work_data_bad_coordinate[i,"Coordinate"] <- unlist( strsplit(wg.1b.work_data_bad_coordinate$Coordinate[i], split = "_") )[2]
}

##Remerging 1a and 1b to work_data

wg.1.work_data <- droplevels( rbind(wg.1a.work_data_good_coordinate, wg.1b.work_data_bad_coordinate) )

###MEMORY CHECKPOINT
rm(list = c("wg.1a.work_data_good_coordinate", "wg.1b.work_data_bad_coordinate", "bad_coordinate_rows", "i"))

#####Convert coordinates with rtracklayer DEPRECATED, DONE BY HAND ON EXCEL FILE GIVEN BY HUGO####Correcting GRCh38 coordinates. Translating them to GRCh37
## wg.1c.GRCh37 <- wg.1.work_data[ wg.1.work_data$Genome_Build_Version == 37,  ]
## wg.1c.GRCh38 <- wg.1.work_data[ wg.1.work_data$Genome_Build_Version == 38,  ]  ### Needs to be converted

##CHECKING BIOMART NEEDS
#removing "seq-" from Locus_Name to leave valid SNPid in form rsXXXXX

wg.1.work_data$Locus_Name <- gsub("seq-","",wg.1.work_data$Locus_Name, ignore.case = T)
wg.1.work_data$Locus_Name <- gsub("exm-","",wg.1.work_data$Locus_Name, ignore.case = T)

##summary(wg.1.work_data)

wg.2a.with_rsid <- wg.1.work_data[ grep("rs",wg.1.work_data$Locus_Name) ,]
wg.2b.without_rsid <- wg.1.work_data[ grep("rs",wg.1.work_data$Locus_Name, invert = T) ,] ##Good to go
##wg.2bx.NO_rsid_No_Chr <- wg.2b.without_rsid[ wg.2b.without_rsid$Chromosome == 0 , ] ##THIS ONES ARE LOST. Ideally, should be cero rows in this df

wg.2aa.rsid_AND_Chr <- wg.2a.with_rsid[ wg.2a.with_rsid$Chromosome != 0 , ]
wg.2ab.rsid_No_Chr <- wg.2a.with_rsid[ wg.2a.with_rsid$Chromosome == 0 , ]  ##This needs to be filled with biomaRt

##DEPRECATED##write.table(wg.3a.with_rsid_No_Cromosome$Locus_Name, "tmp/rsid_without_cromosome", quote = F, col.names = F, row.names = F, fileEncoding = "UTF-8")

### DEV OPT: SAMPLING 2ab to test biomaRt
###wg.2ab.rsid_No_Chr <- droplevels( wg.2ab.rsid_No_Chr[c(1:5),] )


##Fill missing coumns with biomaRt. Using rsID, needs Chromosome, Position, Source ----
hg19Mart = useEnsembl(biomart = "snp", GRCh = 37, dataset = "hsapiens_snp")

for (i in 1:nrow(wg.2ab.rsid_No_Chr)) {
  message(paste0("Retrieving dbSNP ", i,"/", nrow(wg.2ab.rsid_No_Chr)  ))
  wg.dbSlice <- getBM(attributes = c("chr_name","chrom_start", "synonym_source"),
                      
                      filters = "snp_filter",
                      values = wg.2ab.rsid_No_Chr$Locus_Name[i],
                      mart = hg19Mart)
  
  wg.2ab.rsid_No_Chr$Chromosome[i] <- wg.dbSlice$chr_name
  wg.2ab.rsid_No_Chr$Coordinate[i] <- wg.dbSlice$chrom_start
  wg.2ab.rsid_No_Chr$Source_Version[i] <- wg.dbSlice$synonym_source
}

wg.2ab.rsid_No_Chr[is.na(wg.2ab.rsid_No_Chr$Source_Version),"Source_Version"] <-"Archive dbSNP"

#### Juntar todos los df que ya tienen el formato correcto ----

wg.1.work_data <- droplevels( rbind( wg.2b.without_rsid, wg.2aa.rsid_AND_Chr, wg.2ab.rsid_No_Chr ) )
###MEMORY CHECKPOINT
rm(list = c (ls(pattern = "wg.2"),"wg.dbSlice") )

### DEV OPT: SAMPLING 1.work to test allele extractor
###wg.1.work_data <- droplevels( wg.1.work_data[c(1:5),] )
###DEV OPT: TESTING MULTIPLE ALT alleles
##wg.1.work_data$Sequence[5] <- "ATTTCTGGTGAGGTCAACAAGGAGTGCCAAT[A/ALT1/ALT2]CATGCATGCAT"

wg.1.work_data$REF <- "."
wg.1.work_data$ALT <- "."

wg.1.work_data$QUAL <- "."
wg.1.work_data$FILTER <- "PASS"
wg.1.work_data$INFO <- "."

for (i in 1:nrow(wg.1.work_data)) {
  Alleles <- unlist( strsplit(wg.1.work_data$Sequence[i], split = "\\[" ) )[2]
  Alleles <- unlist( strsplit(Alleles, split = "\\]" ) )[1]
  REF <- unlist( strsplit(Alleles, split = "/" ) )[1]
  ALT <- paste( unlist( strsplit(Alleles, split = "/" ) )[-1], collapse = "," )
  wg.1.work_data$REF[i] <- REF
  wg.1.work_data$ALT[i] <- ALT
}
rm(i)
###Find multi-allele variants to Solve: 
#Acceptable SNPs:
#  Can only accept [N/N]
#Cannot accept [N/NNNNNN] or [NNN/NNN] or [N/N/N] or anything other than [N/N].

##THIS SOLVES [N/N/N...] ----
wg.1.work_data$number_of_ALTs <- str_count(wg.1.work_data$ALT, ",")

wg.3a.monoALT <- wg.1.work_data[wg.1.work_data$number_of_ALTs == 0,]
wg.3b.multiALT <- wg.1.work_data[wg.1.work_data$number_of_ALTs > 0,]

wg.3ba.deconcat <- wg.3b.multiALT[0,]

for (i in 1:nrow(wg.3b.multiALT)) {
  ALT_number <- wg.3b.multiALT$number_of_ALTs[i] + 1
    for (n in 1:ALT_number) {
      new_row <- wg.3b.multiALT[i,]
      new_row$ALT <- unlist ( strsplit(wg.3b.multiALT$ALT[i], split = ",") )[n]
      wg.3ba.deconcat <- rbind(wg.3ba.deconcat, new_row)
    }
}

wg.1.work_data <- droplevels( rbind(wg.3a.monoALT, wg.3ba.deconcat) )

#### This solves [NNN/NNN] ----

wg.1.work_data$REF_size <- nchar(wg.1.work_data$REF) 
wg.1.work_data$ALT_size <- nchar(wg.1.work_data$ALT) 

wg.3c.fixing_INDEL_tag <- wg.1.work_data

wg.3c.fixing_INDEL_tag$Target_Type <- ifelse( wg.3c.fixing_INDEL_tag$REF_size == wg.3c.fixing_INDEL_tag$ALT_size , "SNP", "INDEL" )
wg.3ca.SNP_only <- wg.3c.fixing_INDEL_tag[wg.3c.fixing_INDEL_tag$Target_Type == "SNP",]
wg.3cb.INDELs <- wg.3c.fixing_INDEL_tag[wg.3c.fixing_INDEL_tag$Target_Type == "INDEL",]

wg.3cb.INDELs$INDEL_type <- ifelse( wg.3cb.INDELs$REF_size > wg.3cb.INDELs$ALT_size , "DEL", "INS" )
wg.3cba.DEL <- wg.3cb.INDELs[ wg.3cb.INDELs$INDEL_type == "DEL", ]
wg.3cba.INS <- wg.3cb.INDELs[ wg.3cb.INDELs$INDEL_type == "INS", ]

# for (i in 1:nrow(wg.3cba.DEL)) {    ## NO GOOD UNTIL INDELS ARE CORRECTLY CALLED, WITH THE CORRECT POSITION OF THE EVENT
#   INDEL_nucleotides <- unlist ( strsplit( wg.3cba.DEL$ALT[i], split = "" ) )
#   INDEL_nucleotides[1] <- "-"
#   INDEL_nucleotides <- paste(INDEL_nucleotides, collapse = "")
#   wg.3cba.DEL$ALT[i] <- INDEL_nucleotides
# }

# for (i in 1:nrow(wg.3cba.INS)) {    ## NO GOOD UNTIL INDELS ARE CORRECTLY CALLED, WITH THE CORRECT POSITION OF THE EVENT
#   INDEL_nucleotides <- unlist ( strsplit( wg.3cba.INS$REF[i], split = "" ) )
#   INDEL_nucleotides[1] <- "-"
#   INDEL_nucleotides <- paste(INDEL_nucleotides, collapse = "")
#   wg.3cba.DEL$REF[i] <- INDEL_nucleotides
#}

wg.3cba.DEL$ALT <- "-"  ##TEMPORARY
wg.3cba.INS$REF <- "-"  ##TEMPORARY

wg.3cb.INDELs <- rbind( wg.3cba.DEL, wg.3cba.INS )

wg.3cb.INDELs <- wg.3cb.INDELs[,-ncol(wg.3cb.INDELs)]

wg.1.work_data <- rbind(wg.3ca.SNP_only, wg.3cb.INDELs)

###SOLVES "Force Infinium I": Easy fix
##This should be set to FALSE  <<<<<<<<<<<<<<<<<<<<We went with this option
##If set to TRUE, it creates two bead types per SNP

wg.1.work_data$Force_Infinium_I <- "FALSE"

write.table (wg.1.work_data, file = "data/full_data_reformatted.tsv", sep = "\t", quote = F, row.names = F, col.names = T, fileEncoding = "UTF8") ## PROBLEM: saved some lines in data with \t\tFALSE instead of just one tab
### PROBLEM SEEMS TO BE IN Species COLLUMN OF SOME ROWS


#### MEMORY CHECKPOINT
rm(list=ls())

wg.1.work_data <- read.table(file = "data/full_data_reformatted.tsv", header = T, sep = "\t", na.strings = "NA", stringsAsFactors = F)

##Preparar Inf_Sequence_Input_File.csv

wg.Inf_Sequence_Input_File <- wg.1.work_data[ ,
                                              c("Locus_Name", "Target_Type", "Sequence", "Chromosome", "Coordinate",
                                                "Genome_Build_Version", "Source", "Source_Version","Sequence_Orientation",
                                                "Plus_Minus", "Species","Force_Infinium_I")]

wg.Inf_Sequence_Input_File$Species <- "Homo sapiens"
wg.Inf_Sequence_Input_File$Chromosome <- gsub("chr","",wg.Inf_Sequence_Input_File$Chromosome)
wg.Inf_Sequence_Input_File$Chromosome <- gsub("XY","X",wg.Inf_Sequence_Input_File$Chromosome)
wg.Inf_Sequence_Input_File[grep("-", wg.Inf_Sequence_Input_File$Sequence) ,"Target_Type"] <- "INDEL"

write.table (wg.Inf_Sequence_Input_File, file = "data/Inf_Sequence_Input_File.csv", sep = ",", quote = F, row.names = F, col.names = T, fileEncoding = "UTF8")


##MEMORY CHECKPOINT
rm(list = ls())

##PASS TO BATCH FOR faidx retrival of fasta

wg.final <- read.csv("Output_files/WG2_Inf_Sequence_Input_File.csv", header = T)

###SUGIERON REALIZAR UNA ANOTACION FINAL CON biomart para anotar los rsIDs de forma automática y correcta1