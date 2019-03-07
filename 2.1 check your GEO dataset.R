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




# Step1 check data of PCA -------------------------------------------------

## dir.create("./fig/")
load( "./data/tnbc_tumor_FABP4_AssayData.Rdata" )

pca_data <- as.data.frame( t(AssayData) )

library("FactoMineR")
library("factoextra")
dat.pca <- PCA(pca_data, graph = FALSE)
get_eig(dat.pca)
fviz_screeplot(dat.pca, addlabels = TRUE, ylim = c(0, 50))
## Contributions of variables to PC1
fviz_contrib(dat.pca, choice = "var", axes = 1, top = 10)
fviz_pca_ind( dat.pca,
              geom.ind = "point",
              col.ind = group_list,
              palette = c("#00AFBB", "#FC4E07"),
              addEllipses = TRUE,
              legend.title = "Groups"
)
ggsave('./fig/all_gene_PCA_FactoMineR.png')

library("ggfortify")
pca_data$group = group_list
png( './fig/all_gene_PCA_ggfortify.png', res = 120 )
autoplot( 
  prcomp( pca_data[, 1:(ncol( pca_data ) - 1)] ), 
  data = pca_data, 
  colour = 'group') + theme_bw()
dev.off()




# Step2 check data of hclust ----------------------------------------------

colnames( AssayData ) = paste( group_list, 1:ncol( AssayData ), sep = '_' )
nodePar <- list( lab.cex = 0.7, pch = c( NA, 15 ), cex = 0.7, col = "red" )
hc = hclust( dist( t( AssayData ) ) )
png('./fig/all_gene_hclust.png', res = 100, height = 600)
plot( as.dendrogram( hc ), nodePar = nodePar, horiz = TRUE )
dev.off()



# Step3 check data of haetmap ---------------------------------------------

library( "pheatmap" )
choose_gene <- names(tail(sort(apply(AssayData,1,sd)),100))
choose_matrix <- t( scale( t( AssayData[ choose_gene, ] ) ) )

choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix < -2] = -2

annotation_col <- data.frame( CellType = factor( group_list ) )
rownames( annotation_col ) <- colnames( AssayData )

pheatmap( fontsize = 6, choose_matrix, annotation_col = annotation_col, 
          show_rownames = T, show_colnames = F,
          annotation_legend = T, cluster_cols = F, 
          filename = "./fig/heatmap_top100_sd.png")



# Step4 Remove new loaded packages ----------------------------------------

allPackages <- (.packages())
newPackages <- setdiff( allPackages, sysPackages )
lapply( newPackages,
        function(package) {
          package <- paste('package:', package, sep = "", collapse = NULL)
          detach( package, character.only = TRUE )
        }
)
