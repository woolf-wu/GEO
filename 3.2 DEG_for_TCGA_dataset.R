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

load( "./data/tnbc_tumor_FABP4_AssayData.Rdata" )


## function for paint fig
draw_heatmap <- function(nrDEG, type){
  library( "pheatmap" )
  nrDEG_Z = nrDEG[ order( nrDEG$logFC ), ]
  nrDEG_F = nrDEG[ order( -nrDEG$logFC ), ]
  choose_gene = c( rownames( nrDEG_Z )[1:50], rownames( nrDEG_F )[1:50] )
  choose_matrix = AssayData[ choose_gene, ]
  choose_matrix = t( scale( t( choose_matrix ) ) )
  
  choose_matrix[choose_matrix > 2] = 2
  choose_matrix[choose_matrix < -2] = -2
  
  annotation_col = data.frame( CellType = factor( group_list ) )
  rownames( annotation_col ) = colnames( AssayData )
  filename <- paste('./fig/', type, '_heatmap_top100_logFC.png',
                    sep = "", collapse = NULL)
  pheatmap( fontsize = 6, choose_matrix, annotation_col = annotation_col, 
            show_rownames = T, show_colnames = F,
            annotation_legend = T, cluster_cols = T, 
            filename = filename)
}

draw_volcano <- function(nrDEG, type){
  library( "ggplot2" )
  logFC_cutoff <- with( nrDEG, mean( abs( logFC ) ) + 2 * sd( abs( logFC ) ) )
  
  nrDEG$change = as.factor( ifelse( 
    nrDEG$P.Value < 0.01 & abs(nrDEG$logFC) > logFC_cutoff,
    ifelse( nrDEG$logFC > logFC_cutoff, 'UP', 'DOWN' ), 'NOT' ) )
  nrDEGfile <- paste('./data/', type, '_nrDEG_by_logFC.Rdata',
                    sep = "", collapse = NULL)
  save( nrDEG, file = nrDEGfile )
  
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
  filename <- paste('./fig/', type, '_volcano_logFC.png',
                    sep = "", collapse = NULL)
  ggsave( volcano, filename = filename )
  dev.off()
}



# Step1 Firstly run DESeq2 ------------------------------------------------

library( DESeq2 )
colData <- data.frame( row.names = colnames(AssayData), 
                       group_list = group_list)
## results of an analysis of differential expression
dds <- DESeqDataSetFromMatrix( countData = AssayData,
                               colData = colData,
                               design = ~ group_list)
ddsFile <- "./data/TCGA-BRCA-mRNA-DESeq2-dds.Rdata"
if ( !file.exists( ddsFile )) {
  ## Differential expression analysis based on the Negative Binomial 
  ## (a.k.a. Gamma-Poisson) distribution
  dds <- DESeq(dds)
  save(dds, file = ddsFile )
}else{
  load( file = ddsFile )
}

resultsNames(dds)
res <- results( dds )
resOrdered <- res[ order(res$padj), ]
head(resOrdered)
nrDEG_DESeq2 <- as.data.frame( resOrdered )
nrDEG_DESeq2 <- na.omit( nrDEG_DESeq2 )
colnames(nrDEG_DESeq2)[2] <- c("logFC") 
colnames(nrDEG_DESeq2)[5] <- c("P.Value") 
draw_heatmap(nrDEG = nrDEG_DESeq2, type = 'DESeq2')
draw_volcano(nrDEG = nrDEG_DESeq2, type = 'DESeq2')
  


# Step2 Then run edgeR ----------------------------------------------------

library(edgeR)
## A list-based S4 class for storing read counts and associated information 
## from digital gene expression or sequencing technologies.
DGElist <- DGEList( counts = AssayData, group = factor(group_list) )
## Counts per Million or Reads per Kilobase per Million
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
## Calculate Normalization Factors to Align Columns of a Count Matrix
DGElist <- calcNormFactors( DGElist )
DGElist$samples

design <- model.matrix( ~0 + factor(group_list) )
rownames(design) <- colnames(DGElist)
colnames(design) <- levels(factor(group_list))

## Estimate Common Dispersion for Negative Binomial GLMs
DGElist <- estimateGLMCommonDisp(DGElist, design)
## Estimate Trended Dispersion for Negative Binomial GLMs
DGElist <- estimateGLMTrendedDisp(DGElist, design)
## Empirical Bayes Tagwise Dispersions for Negative Binomial GLMs
DGElist <- estimateGLMTagwiseDisp(DGElist, design)

## glmFit fits genewise negative binomial glms, all with the same design matrix 
## but possibly different dispersions, offsets and weights
fit <- glmFit(DGElist, design)
## https://www.biostars.org/p/110861/
## glmLRT conducts likelihood ratio tests for one or more coefficients in the 
## linear model.
results <- glmLRT(fit, contrast = c(-1, 1)) 
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
head(nrDEG_edgeR)
colnames(nrDEG_edgeR)[4] <- c("P.Value") 
draw_heatmap(nrDEG = nrDEG_edgeR, type = 'edgeR')
draw_volcano(nrDEG = nrDEG_edgeR, type = 'edgeR')


# Step3 Lastly run voom from limma ----------------------------------------

library(limma)
## A list-based S4 class for storing read counts and associated information 
## from digital gene expression or sequencing technologies.
DGElist <- DGEList( counts = AssayData, group = factor(group_list) )
## Counts per Million or Reads per Kilobase per Million
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
## Calculate Normalization Factors to Align Columns of a Count Matrix
DGElist <- calcNormFactors( DGElist )
DGElist$samples

design <- model.matrix( ~0 + factor(group_list) )
rownames(design) <- colnames(DGElist)
colnames(design) <- levels(factor(group_list))

## Transform RNA-Seq Data Ready for Linear Modelling
v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
## Fit linear model for each gene given a series of arrays
fit <- lmFit(v, design)

## Construct the contrast matrix corresponding to specified contrasts of a set 
## of parameters.
cont.matrix <- makeContrasts(contrasts = c('TP53-NO_TP53'), levels = design)
## Given a linear model fit to microarray data, compute estimated coefficients 
## and standard errors for a given set of contrasts.
fit2 <- contrasts.fit(fit, cont.matrix)
## Empirical Bayes Statistics for Differential Expression
fit2 <- eBayes(fit2)
  
nrDEG_limma_voom = topTable(fit2, coef = 'TP53-NO_TP53', n = Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
draw_heatmap(nrDEG = nrDEG_limma_voom, type = 'limma_voom')
draw_volcano(nrDEG = nrDEG_limma_voom, type = 'limma_voom')


# Step4 Compare three methods ---------------------------------------------

mi <- unique(c(rownames(nrDEG_DESeq2),
               rownames(nrDEG_edgeR),
               rownames(nrDEG_limma_voom)))
lf <- data.frame(DESeq2 = nrDEG_DESeq2[mi, 2],
                 edgeR = nrDEG_edgeR[mi, 1],
                 limma_voom = nrDEG_limma_voom[mi, 1])
cor(na.omit(lf))
mi <- unique(c(rownames(nrDEG_edgeR),
               rownames(nrDEG_limma_voom)))
lf <- data.frame(edgeR = nrDEG_edgeR[mi, 1],
                 limma_voom = nrDEG_limma_voom[mi, 1])
cor(na.omit(lf))

## firehose-BRCA-RPPA
if (!file.exists( './data/BRCA.antibody_annotation.Rdata' )) {
  ## http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BRCA/20160128/gdac.broadinstitute.org_BRCA.RPPA_AnnotateWithGene.Level_3.2016012800.0.0.tar.gz"
  untar_F <- "./raw_data/gdac.broadinstitute.org_BRCA.RPPA_AnnotateWithGene.Level_3.2016012800.0.0/BRCA.antibody_annotation.txt"
  phenoData <- read.table( untar_F,
                           header = T,
                           sep = '\t',
                           quote = '' )
  save( phenoData, file = './data/BRCA.antibody_annotation.Rdata' )
}else{
  load('./data/BRCA.antibody_annotation.Rdata')
}

# Step5 Remove new loaded packages ----------------------------------------

allPackages <- (.packages())
newPackages <- setdiff( allPackages, sysPackages )
lapply( newPackages,
        function(package) {
          package <- paste('package:', package, sep = "", collapse = NULL)
          detach( package, character.only = TRUE )
        }
)
