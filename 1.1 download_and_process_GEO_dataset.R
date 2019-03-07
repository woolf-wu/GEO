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



# Step1 download GEO dateset by GEOquery ----------------------------------

## Basic data
GSE_ID <- ''
GSE_file <- paste('./data/', GSE_ID, '.Rdata', sep = "", collapse = NULL)

library( "GEOquery" )
if (!file.exists( GSE_file )) {
  ## dir.create("./data/")
  download_GEO <- function(GSE_ID) {
    options( download.file.method.GEOquery = 'libcurl' )
    ## dir.create("./raw_data/")
    gset <- getGEO( GSE_ID, getGPL = T, destdir = "./raw_data/" )
    save( gset, file = GSE_file )
  }
  download_GEO( GSE_ID )
}



# Step2 Data loading ---------------------------------------------------

load( GSE_file )
class( gset )
length( gset )

exprSet <- gset[[1]]
str( exprSet, max.level = 2 )

## assayData <- exprSet@assayData$exprs[1:5,1:6]
assayData <- exprs(exprSet)
dim(assayData)
assayData[,1:4]

## phenoData <- exprSet@phenoData@data
phenoData <- pData(exprSet)
dim(phenoData)
phenoData[1:2,1:6]


# Step3 Grouping by special clinical information --------------------------

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

## choose special pheno
group_list <- phenoData[, "neo-adjuvant chemotherapy:ch1"]
table(group_list)
group_list <- ifelse( group_list == 'resistant',
                     'resistant','sensitive')
names(group_list) <- rownames(phenoData)
newAssayData <- assayData[, names(group_list)]
dim(newAssayData)
newAssayData[1:5, 1:6]

if (!file.exists( "./data/AssayData_by_group.Rdata" )) {
  save( newAssayData, group_list, file = "./data/AssayData_by_group.Rdata" )
}



# Step4 Filt gene ---------------------------------------------------------

load( './data/AssayData_by_group.Rdata' )
load( './data/GSE108565.Rdata' )
exprSet <- gset[[1]]
featureData = exprSet@featureData@data
dim( featureData )
colnames( featureData )
View( featureData )
## b <- apply(featureData, 1, 
##       function(i){
##         featureData[i, 1] == featureData[i, 2]
##       })
## table(b)
featureData$max <- apply(newAssayData, 1, max)
featureData <- featureData[order(featureData$ID,
                                 featureData$max,
                                 decreasing = T), ]
dim( featureData )
featureData <- featureData[!duplicated(featureData$ID), ]
dim( featureData )
AssayData <- newAssayData[featureData$ID, ]
dim(AssayData)
AssayData[1:5, 1:6]

# AssayData = log2(AssayData)
# dim(AssayData)
# AssayData[1:5,1:6]
# group_list <- unname(group_list)

save(AssayData, group_list, file = './data/final_AssayData.Rdata')



# Step5 Remove new loaded packages ----------------------------------------

allPackages <- (.packages())
newPackages <- setdiff( allPackages, sysPackages )
lapply( newPackages,
        function(package) {
          package <- paste('package:', package, sep = "", collapse = NULL)
          detach( package, character.only = TRUE )
        }
)

