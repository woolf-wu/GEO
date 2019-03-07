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

## function of KEGG pathway
kegg_plot <- function(type) {
  kk.up <- enrichKEGG(   gene          =  gene_up    ,
                         organism      =  'hsa'      ,
                         universe      =  gene_all   ,
                         pvalueCutoff  =  0.8        ,
                         qvalueCutoff  =  0.8        )
  kk.down <- enrichKEGG( gene          =  gene_down  ,
                         organism      =  'hsa'      ,
                         universe      =  gene_all   ,
                         pvalueCutoff  =  0.8        ,
                         qvalueCutoff  =  0.8        )
  library( "ggplot2" )
  kegg_down_dt <- as.data.frame( kk.down )
  kegg_up_dt <- as.data.frame( kk.up )
  down_kegg <- kegg_down_dt[ kegg_down_dt$pvalue < 0.05, ]
  down_kegg$group <- -1
  up_kegg <- kegg_up_dt[ kegg_up_dt$pvalue < 0.05, ]
  up_kegg$group <- 1
  dat = rbind( up_kegg, down_kegg )
  dat$pvalue = -log10( dat$pvalue )
  dat$pvalue = dat$pvalue * dat$group
  dat = dat[ order( dat$pvalue, decreasing = F ), ]
  g_kegg <- ggplot( dat, 
                    aes(x = reorder( Description, order( pvalue, decreasing = F )), 
                        y = pvalue, fill = group)) + 
    geom_bar( stat = "identity" ) + 
    scale_fill_gradient( low = "blue", high = "red", guide = FALSE ) + 
    scale_x_discrete( name = "Pathway names" ) +
    scale_y_continuous( name = "log10P-value" ) +
    coord_flip() + theme_bw() + 
    theme( plot.title = element_text( hjust = 0.5 ), 
           axis.text.x = element_text(size = 10),
           axis.text.y = element_text(size = 7)) +
    ggtitle( "Pathway Enrichment" ) 
  print( g_kegg )
  filename <- paste('./fig/kegg_up_down_', type, '.png', sep = "", collapse = NULL)
  ggsave( g_kegg, filename = filename )
}

## function of GO pathway
go_plot <- function(type) {
  go_enrich_results <- lapply( g_list, function(gene) {
    lapply( c( 'BP', 'MF', 'CC' ) , function(ont) {
      cat( paste( 'Now process', ont ) )
      ego <- enrichGO( gene          =  gene,
                       universe      =  gene_all,
                       OrgDb         =  org.Hs.eg.db,
                       ont           =  ont ,
                       pAdjustMethod =  "BH",
                       pvalueCutoff  =  0.99,
                       qvalueCutoff  =  0.99,
                       readable      =  TRUE)
      print( head( ego ) )
      return( ego )
    })
  })
  gofilename <- paste('./data/go_enrich_result', type, '.Rdata', 
                      sep = "", collapse = NULL)
  save( go_enrich_results, file = gofilename )
  
  n1 = c( 'gene_up', 'gene_down', 'gene_diff' )
  n2 = c( 'BP', 'MF', 'CC' ) 
  for (i in 1:3) {
    for (j in 1:3) {
      fn = paste0( './fig/dotplot_', n1[i], '_', n2[j], '_', type, '.png' )
      cat( paste0( fn, '\n' ) )
      png( fn, res = 150, width = 1080 )
      print( dotplot( go_enrich_results[[i]][[j]] ) )
      dev.off()
    }
  }
}



# Step1 annotation --------------------------------------------------------

library( "clusterProfiler" )
library( "org.Hs.eg.db" )
keytypes(org.Hs.eg.db)
library("stringr")

load( "./data/DESeq2_nrDEG_by_logFC.Rdata" )
load( "./data/edgeR_nrDEG_by_logFC.Rdata" )
load( "./data/limma_voom_nrDEG_by_logFC.Rdata" )

## tans1: ENSEMBL2ENTREZID
table( nrDEG$change )
rownames( nrDEG ) <- str_sub(rownames( nrDEG ), start = 1, end = 15)
nrDEG$ENSEMBL <- rownames( nrDEG )
df <- bitr( rownames( nrDEG ), fromType = "ENSEMBL", toType = c( "ENTREZID" ), 
            OrgDb = org.Hs.eg.db )
head( df )
nrDEG = merge( nrDEG, df, by = 'ENSEMBL' )
head( nrDEG )

gene_up = nrDEG[ nrDEG$change == 'UP', 'ENTREZID' ] 
gene_down = nrDEG[ nrDEG$change == 'DOWN', 'ENTREZID' ]
gene_diff = c( gene_up, gene_down )
gene_all = as.character(nrDEG[ ,'ENTREZID'] )
g_list = list( gene_up = gene_up, gene_down = gene_down, gene_diff = gene_diff)


# Step2 pathway analysis ------------------------------------------

kegg_plot("DESeq2")
go_plot("DESeq2")

kegg_plot("edgeR")
go_plot("edgeR")

kegg_plot("limma_voom")
go_plot("limma_voom")

library(pathview)
geneList <- nrDEG$logFC
names( geneList ) <- nrDEG$ENTREZID
geneList <- sort( geneList, decreasing = T )
pathview( gene.data = geneList, 
          pathway.id = dat$ID, 
          species = "hsa", 
          limit = list(gene = 5, cpd = 1))



# Step3 Remove new loaded packages ----------------------------------------

allPackages <- (.packages())
newPackages <- setdiff( allPackages, sysPackages )
lapply( newPackages,
        function(package) {
          package <- paste('package:', package, sep = "", collapse = NULL)
          detach( package, character.only = TRUE )
        }
)

