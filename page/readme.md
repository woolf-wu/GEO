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



# Step1 Setting CRAN mirror -----------------------------------------------

local({
  options( repos  = "https://mirrors.ustc.edu.cn/CRAN/" )
  options( BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/" )
})



# Step2 List of the used packages  ----------------------------------------

bioPackages <- 
c( 
  "dplyr", "stringi", "purrr", ## ERROR
  "R.utils", "data.table", ## unzip and read table
  "GEOquery", ## download
  "FactoMineR", "factoextra", "ggfortify", ## PCA
  "pheatmap", ## heatmap
  "ggplot2", ## Volcano plot
  "limma", "DESeq2", "edgeR", ## DEG
  "clusterProfiler", "org.Hs.eg.db", ## annotation
  "pathview" ## kegg
)



# Step3 Install the packages ----------------------------------------------

lapply( bioPackages, 
  function(bioPackage) {
    if ( !require( bioPackage, character.only = T ) ) {
      CRANpackages <- available.packages()
      
      ## install packages by CRAN
      if ( bioPackage %in% rownames( CRANpackages ) ) {
        install.packages( bioPackage )
      
      }else{
        ## install packages by bioconductor
          ## R version >= 3.5 ===> BiocManager
        if ( as.character( sessionInfo()$R.version$minor ) >= 3.5 ) {
          if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
          BiocManager::install(bioPackage, update = TRUE, ask = FALSE)
        
        }else{
          ## R version < 3.5 ===> BiocInstaller
          if (!requireNamespace("BiocInstaller", quietly = TRUE))
            source( "https://bioconductor.org/biocLite.R" )
          BiocInstaller::biocLite( bioPackage, ask = FALSE)
        }
      }
    }
  }
)



# Step4 Remove new loaded packages ----------------------------------------

allPackages <- (.packages())
newPackages <- setdiff( allPackages, sysPackages )
lapply( newPackages,
        function(package) {
          package <- paste('package:', package, sep = "", collapse = NULL)
          detach( package, character.only = TRUE )
        }
)
