# Step0 Before starting your project --------------------------------------

## Remove everything in the working environment, not including loaded libraries.
rm(list = objects( all = TRUE ))

if (!is.null( dev.list() )) dev.off() 

clearhistory <- function() {
  write( "", file = ".blank" )
  loadhistory( ".blank" )
  unlink( ".blank" )
}
clearhistory()

## basecal packages
sysPackages <- (.packages())

## data.frame(..., row.names = NULL, check.rows = FALSE,
##            check.names = TRUE, fix.empty.names = TRUE,
##            stringsAsFactors = default.stringsAsFactors())
options( stringsAsFactors = FALSE )


## Now winInet not supported for use in service, but the default setting of 
## download.file.method is "wininet". 
## If your system support "libcurl", set the downloaded method to libcurl.
if ( capabilities( "libcurl" ) == T ) {
  options( download.file.method = "libcurl" )
}
options()$download.file.method


## Change the library location of the packages
## Even your updated your R, you can still use your packages.
.libPaths( c( "G:/R-packages",
              "C:/Program Files/R/R-3.5.2/library") )
.libPaths()




# Step1 download TCGA dateset ---------------------------------------------

if (!file.exists( './data/TCGA-BRCA.htseq_counts.Rdata' )) {
  gzfile <- "./raw_data/TCGA-BRCA.htseq_counts.tsv.gz"
  download.file("https://gdc.xenahubs.net/download/TCGA-BRCA/Xena_Matrices/TCGA-BRCA.htseq_counts.tsv.gz", 
                destfile = gzfile)
  library(R.utils)
  gunzip(gzfile, remove = F)
  library(data.table)
  raw_data <- fread( "./raw_data/TCGA-BRCA.htseq_counts.tsv",
                     sep = '\t', header = T)
  raw_data <- as.data.frame( raw_data )
  raw_data[1:5, 1:6] 
  rownames( raw_data ) <- raw_data[, 1]
  raw_data <- raw_data[, -1]
  raw_data[1:5, 1:6]
  raw_data <- 2^raw_data - 1
  raw_data <- ceiling( raw_data )
  raw_data[1:5, 1:6]
  pick_row <- apply( raw_data, 1, function(x){
    sum(x == 0) < 10
  })
  raw_data <- raw_data[pick_row, ]
  dim(raw_data  )
  save( raw_data, file = './data/TCGA-BRCA.htseq_counts.Rdata' )
}else{
  load('./data/TCGA-BRCA.htseq_counts.Rdata')
}

## firehose-BRCA-RPPA
if (!file.exists( './data/BRCA.Merge_protein_exp_rppa_Level_3_normalization.Rdata' )) {
  ## http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BRCA/20160128/gdac.broadinstitute.org_BRCA.Merge_protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.Level_3.2016012800.0.0.tar.gz"
  untar_F <- "./raw_data/BRCA.protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.data.txt"
  library(data.table)
  raw_data <- fread( untar_F, sep = '\t', header = T)
  raw_data <- as.data.frame( raw_data )
  raw_data[1:5, 1:3] 
  raw_data <- raw_data[-1, ]
  name <- raw_data[ , 1]
  raw_data <- raw_data[ , -1]
  raw_data <- apply(raw_data, 2, as.numeric)
  raw_data[1:5, 1:3]
  rownames( raw_data ) <- name
  raw_data[is.na(raw_data)] <- 0
  raw_data[1:5, 1:3]
  pick_row <- apply( raw_data, 1, function(x){
    sum(x == 0) < 10
  })
  raw_data <- raw_data[pick_row, ]
  dim(raw_data  )
  save( raw_data, file = './data/BRCA.Merge_protein_exp_rppa_Level_3_normalization.Rdata' )
}else{
  load('./data/BRCA.Merge_protein_exp_rppa_Level_3_normalization.Rdata')
}



# Step2 Grouping by special clinical information --------------------------

if (!file.exists( './raw_data/TCGA-BRCA.GDC_phenotype.tsv.gz' )) {
  gzfile <- "./raw_data/TCGA-BRCA.GDC_phenotype.tsv.gz"
  download.file("https://gdc.xenahubs.net/download/TCGA-BRCA/Xena_Matrices/TCGA-BRCA.GDC_phenotype.tsv.gz", 
                destfile = gzfile)
  phenoData <- read.table( gzfile,
                           header = T,
                           sep = '\t',
                           quote = '' )
  save( phenoData, file = './data/TCGA-BRCA.GDC_phenotype.Rdata' )
}else{
  load('./data/TCGA-BRCA.GDC_phenotype.Rdata')
}


pheno_num <- c()
invisible(
  lapply(1:ncol(phenoData), 
         function(col_num){
           ## Assume that the classification project is between 2 and 4
           if (1 < dim(table(phenoData[,col_num])) & 
               dim(table(phenoData[,col_num])) < 5) {
             pheno_num <<- append(pheno_num, col_num, after = length(pheno_num))
           }
         }
  )
)
View(phenoData[, pheno_num])
names(phenoData[, pheno_num])

## Category 1: Triple negative breast cancer
## Pick columns that contains "receptor_status"
receptot_names <- grep('receptor_status', names(phenoData), value = T)
phenoData[1:2, receptot_names]
triPheno <- phenoData[, receptot_names[1:3]]
## Group samples by "ER", "PR", "HER2"
tnbc_rownum <- apply(triPheno, 1, function(x) sum(x == 'Negative'))
tnbc_sample <- phenoData[tnbc_rownum == 3, 1]
tnbc_sample <- intersect(colnames(raw_data), tnbc_sample)
save(tnbc_sample, file = './data/sample_by_tnbc.Rdata')

## Category 2: race
## Pick columns that contains 'race.demographic'
race_names <- grep('race.demographic', names(phenoData), value = T)
racePheno <- phenoData[, race_names]
table(racePheno)
col_asian <- grep("asian", racePheno)
col_black <- grep("black or african american", racePheno)
col_white <- grep("white", racePheno)
asian_sample <- phenoData[col_asian, 1]
black_sample <- phenoData[col_black, 1]
white_sample <- phenoData[col_white, 1]
tnbc_tumor_sample <- tnbc_sample[substr( tnbc_sample,14,15) < 10]
asian_tnbc_sample <- intersect(tnbc_tumor_sample, asian_sample)
black_tnbc_sample <- intersect(tnbc_tumor_sample, black_sample)
white_tnbc_sample <- intersect(tnbc_tumor_sample, white_sample)
save(asian_tnbc_sample, black_tnbc_sample, white_tnbc_sample, 
     file = './data/tnbc_tumor_sample_by_race.Rdata')

## Category 3: BRCA
if (!file.exists( './raw_data/TCGA-BRCA.mutect2_snv.tsv.gz' )) {
  gzfile <- "./raw_data/TCGA-BRCA.mutect2_snv.tsv.gz"
  download.file("https://gdc.xenahubs.net/download/TCGA-BRCA/Xena_Matrices/TCGA-BRCA.mutect2_snv.tsv.gz", 
                destfile = gzfile)
  library(R.utils)
  gunzip(gzfile, remove = F)
  mutype_file <- read.table( "./raw_data/TCGA-BRCA.mutect2_snv.tsv",
                             header = T,
                             sep = '\t',
                             quote = '' )
  save( mutype_file, file = './data/TCGA-BRCA.mutect2_snv.Rdata' )
}else{
  load('./data/TCGA-BRCA.mutect2_snv.Rdata')
}

## Pick columns that contains 'brca'
BRCA <- mutype_file[mutype_file$gene == 'BRCA1' | mutype_file$gene == 'BRCA2',]
brca_sample <- unique( sort( BRCA$Sample_ID ) )
save(brca_sample, file = './data/sample_by_BRCA.Rdata')

## Pick columns that contains 'tp53'
TP53 <- mutype_file[mutype_file$gene == 'tp53' | mutype_file$gene == 'TP53',]
TP53_sample <- unique( sort( TP53$Sample_ID ) )
tumor_sample <- colnames(raw_data)[substr( colnames(raw_data),14,15) < 10]
TP53_sample <- intersect(tumor_sample, TP53_sample)
noTP53_sample <- setdiff(tumor_sample, TP53_sample)
save(TP53_sample, noTP53_sample, file = './data/sample_by_TP53.Rdata')

## Pick columns that contains 'FABP4'
FABP4 <- mutype_file[mutype_file$gene == 'FABP4' | mutype_file$gene == 'fabp4',]
FABP4_sample <- unique( sort( FABP4$Sample_ID ) )
tumor_sample <- colnames(raw_data)[substr( colnames(raw_data),14,15) < 10]
tumor_sample <- substr(tumor_sample, 1, 16)
FABP4_sample <- intersect(tumor_sample, FABP4_sample)
noFABP4_sample <- setdiff(tumor_sample, FABP4_sample)
FABP4_sample <- c("TCGA-A1-A0SO-01A", "TCGA-A1-A0SQ-01A", 
                  "TCGA-BH-A0H5-01A", "TCGA-BH-A0H6-01A")
save(FABP4_sample, file = './data/sample_by_FABP4.Rdata')



# Step3 Filt sample ------------------------------------------------

load('./data/TCGA-BRCA.htseq_counts.Rdata')

## TCGA barcode: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
PairedSample <- function(sample, type){
  sample <- substring( sample, 1, 12 )
  raw_data_sample <- substring( colnames(raw_data), 1, 12 )
  paired_sample <- names(table(raw_data_sample)[table(raw_data_sample) >= 2])
  tnbc_paired_sample <- intersect(sample, paired_sample)
  tnbc_paired_sample <- colnames(raw_data)[raw_data_sample %in% tnbc_paired_sample]
  ## tnbc_paired_sample <- tnbc_paired_sample[-4]
  AssayData <- raw_data[, tnbc_paired_sample]
  group_list <- factor(
    ifelse( as.numeric( substr( 
      colnames( AssayData ),14,15)) < 10, 'tumor', 'normal'))
  filename <- paste('./data/', type, '_paired_AssayData.Rdata', 
                    sep = "", collapse = NULL)
  save(AssayData, group_list, file = filename)
}

PairedSample(sample = tnbc_sample, type = 'tnbc')
PairedSample(sample = brca_sample, type = 'brca')



race_sample <- c(asian_tnbc_sample, black_tnbc_sample, white_tnbc_sample)
AssayData <- raw_data[, race_sample]
dim(AssayData)
group_list <- c(rep('asian', length(asian_tnbc_sample)),
                rep('black', length(black_tnbc_sample)),
                rep('white', length(white_tnbc_sample)))
save(AssayData, group_list, file = './data/tnbc_tumor_race_AssayData.Rdata')


tp53_sample <- c(TP53_sample, noTP53_sample)
AssayData <- raw_data[, tp53_sample]
dim(AssayData)
group_list <- c(rep('TP53', length(TP53_sample)),
                rep('NO_TP53', length(noTP53_sample)))
save(AssayData, group_list, file = './data/tnbc_tumor_TP53_AssayData.Rdata')



AssayData <- raw_data[, substr(colnames(raw_data), 1, 16) %in% FABP4_sample]
dim(AssayData)
group_list <- c("FABP4","contral","FABP4","contral")
save(AssayData, group_list, file = './data/tnbc_tumor_FABP4_AssayData.Rdata')

# Step5 Remove new loaded packages ----------------------------------------

allPackages <- (.packages())
newPackages <- setdiff( allPackages, sysPackages )
lapply( newPackages,
        function(package) {
          package <- paste('package:', package, sep = "", collapse = NULL)
          detach( package, character.only = TRUE )
        }
)

