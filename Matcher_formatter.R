##library(rtracklayer)
library(openxlsx)
library(biomaRt)

wg.0.Hugo_data <- read.xlsx("data/listas_consolidadas3.xlsx", sheet = 1, colNames = T, skipEmptyRows = T, check.names = T)

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
wg.1.work_data <- droplevels( wg.1.work_data[-grep(">:",wg.1.work_data$Sequence),] )

##wg.1.work_data <- unique.data.frame(wg.1.work_data) ##Removes duplicate lines


###MEMORY CHECKPOINT
rm(wg.noFASTA)

###### ----  
### DATA REFORMATING ----

bad_coordinate_rows <- grep("rs", wg.1.work_data$Coordinate)

#### Resolver problema de mala antacion en columna Coordinate
wg.1a.work_data_good_coordinate <- wg.1.work_data[ -bad_coordinate_rows, ]
wg.1b.work_data_bad_coordinate <- wg.1.work_data[ bad_coordinate_rows, ]

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

##DEPRECATED, DONE BY HAND ON EXCEL FILE GIVEN BY HUGO####Correcting GRCh38 coordinates. Translating them to GRCh37
# wg.1c.GRCh37 <- wg.1.work_data[ wg.1.work_data$Genome_Build_Version == 37,  ]
# wg.1c.GRCh38 <- wg.1.work_data[ wg.1.work_data$Genome_Build_Version == 38,  ]  ### Needs to be converted

###Convert coordinates with rtracklayer

##CHECKING BIOMART NEEDS
#removing "seq-" from Locus_Name to leave valid SNPid in forma rsXXXXX

wg.1.work_data$Locus_Name <- gsub("seq-","",wg.1.work_data$Locus_Name, ignore.case = T)
wg.1.work_data$Locus_Name <- gsub("exm-","",wg.1.work_data$Locus_Name, ignore.case = T)

summary(wg.1.work_data)

wg.2a.with_rsid <- wg.1.work_data[ grep("rs",wg.1.work_data$Locus_Name) ,]
wg.2b.without_rsid <- wg.1.work_data[ -grep("rs",wg.1.work_data$Locus_Name) ,] ##Good to go
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

#### Juntar todos los df que ya tienen el formato correcto ----

wg.1.work_data <- droplevels( rbind( wg.2b.without_rsid, wg.2aa.rsid_AND_Chr, wg.2ab.rsid_No_Chr ) )

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


###SUGIERON REALIZAR UNA ANOTACION FINAL CON biomart para anotar los rsIDs de forma automática y correcta1

### Guardar el df Final en el formato correcto para Samtools faidx. El formato es:
#!!## fileformat=VCFv4.1
#!!# CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

fileConn <- file("data/vcf_for_fasta-from-vcf2_script.vcf")
writeLines(c("## fileformat=VCFv4.1",
             "# CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO"),
           fileConn)
close(fileConn)

write.table(wg.1.work_data[,c("Chromosome","Coordinate","Locus_Name","REF","ALT","QUAL","FILTER","INFO")], file = "data/vcf_for_fasta-from-vcf2_script.vcf", append = T, quote = F, sep = "\t", row.names = F, col.names = F, fileEncoding = "UTF-8")
