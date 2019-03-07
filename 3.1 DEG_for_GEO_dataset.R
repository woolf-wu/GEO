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




# Step1 DEG by limma ------------------------------------------------------

load( "./data/final_exprSet.Rdata" )

library( "limma" )
{
  design <- model.matrix( ~0 + factor( group_list ) )
  colnames( design ) = levels( factor( group_list ) )
  rownames( design ) = colnames( AssayData )
}
design

contrast.matrix <- makeContrasts( "tumor-normal", levels = design )
contrast.matrix

{
  fit <- lmFit( AssayData, design )
  fit2 <- contrasts.fit( fit, contrast.matrix ) 
  fit2 <- eBayes( fit2 )
  nrDEG <-  topTable( fit2, coef = 1, n = Inf )
  save(nrDEG, file = './data/nrDEG.Rdata')
}
head(nrDEG)



# Step2 Heatmap -----------------------------------------------------------

library( "pheatmap" )
{
  nrDEG_Z = nrDEG[ order( nrDEG$logFC ), ]
  nrDEG_F = nrDEG[ order( -nrDEG$logFC ), ]
  choose_gene = c( rownames( nrDEG_Z )[1:50], rownames( nrDEG_F )[1:50] )
  choose_matrix = AssayData[ choose_gene, ]
  choose_matrix = t( scale( t( choose_matrix ) ) )
  
  choose_matrix[choose_matrix > 2] = 2
  choose_matrix[choose_matrix < -2] = -2
  
  annotation_col = data.frame( CellType = factor( group_list ) )
  rownames( annotation_col ) = colnames( AssayData )
  pheatmap( fontsize = 6, choose_matrix, annotation_col = annotation_col, 
            show_rownames = T, show_colnames = F,
            annotation_legend = T, cluster_cols = T, 
            filename = "./fig/heatmap_top100_logFC.png")
}



# Step3 Volcano plot ------------------------------------------------------

library( "ggplot2" )
logFC_cutoff <- with( nrDEG, mean( abs( logFC ) ) + 2 * sd( abs( logFC ) ) )
logFC_cutoff
## logFC_cutoff = 0.5

{
  nrDEG$change = as.factor( ifelse( 
    nrDEG$P.Value < 0.01 & abs(nrDEG$logFC) > logFC_cutoff,
      ifelse( nrDEG$logFC > logFC_cutoff, 'UP', 'DOWN' ), 'NOT' ) )
  
  save( nrDEG, file = "./data/nrDEG_by_logFC.Rdata" )
  
  this_tile <- paste0( 
    'Cutoff for logFC is ', round( logFC_cutoff, 3 ),
    '\nThe number of up gene is ', nrow(nrDEG[ nrDEG$change == 'UP', ] ),
    '\nThe number of down gene is ', nrow(nrDEG[ nrDEG$change == 'DOWN', ] ) )
  
  volcano = ggplot(data = nrDEG, 
                   aes( x = logFC, y = -log10(P.Value), color = change)) +
            geom_point( alpha = 0.4, size = 1.75 ) +
            theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
            xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
            ggtitle( this_tile ) + 
            theme( plot.title = element_text( size = 15, hjust = 0.5 )) +
            scale_colour_manual( values = c('blue','black','red') )
  print( volcano )
  ggsave( volcano, filename = './fig/volcano.png' )
  dev.off()
}



# Step4 Remove new loaded packages ----------------------------------------

allPackages <- (.packages())
newPackages <- setdiff( allPackages, sysPackages )
lapply( newPackages,
        function(package) {
          package <- paste('package:', package, sep = "", collapse = NULL)
          detach( package, character.only = TRUE )
        }
)
