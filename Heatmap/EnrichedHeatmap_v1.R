## Wrapper for EnrichedHeatmap, producing heatmaps for one target and one or many bigwigs.

## FIRST DRAFT, not finished yet!

require(rtracklayer)
require(GenomicRanges)
require(EnrichedHeatmap)

set.seed(2019)

plot.EnrHtmp <- function(bigwigs,            ## list of bigwigs
                         bigwig.names = NULL, ## a list of length(BigWig) with the names be used for plotting and lengend title
                         target.gr,          ## GRanges with target regions
                         Target.name,        ## a name for the regions to use while plotting
                         window.size,        ## total window size around the centers of target.gr, so 10000 would be 5kb each side
                         pdf.dir = "./",     ## directory to save plots to
                         scale.all = TRUE,   ## scale all plots to min of the minest and max of maxest sample
                         do.log2 = FALSE,    ## whether to log2 transform the bigwig
                         log2.prior = 1,     ## prior count in case of log
                         col.low = "darkblue",
                         col.high = "darkgoldenrod1",
                         main.titles = NULL,  ## NULL or a lsit of length(bigwigs)
                         font.fam = "sans",
                         font.size.title = 15,
                         legend.title = "ATAC-seq [CPM]"
                         ) {
  
  #######################################################################################################################################
  
  ## check for correct formats:
  stop_quietly <- function() {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
  
  if (length(bigwigs) != length(bigwig.names)){
    message("[Error]: Length of BigWig and bigwig.names is not the same!")
    stop_quietly()
  }
  
  if (class(target.gr) != "GRanges") {
    message("[Error]: target.gr is not a GRanges object!")
    stop_quietly()
  }
  
  if (!is.numeric(window.size)) {
    message("[Error]: window.size must be an integer!")
    stop_quietly()
  }
  
  window.size <- ceiling(window.size)
  
  ## check if bigwigs exist:
  noBig <- c()
  for (i in bigwigs){
    if (!file.exists(i)) noBig <- c(noBig, i)
  }; rm(i)
  if (length(noBig) > 0) {
    message("[Error]: Those bigwigs do not exist: ", paste(noBig, collapse = " "))
    stop_quietly()
  }
  
  if (!is.null(main.titles) && (length(main.titles) != length(bigwigs))){
    message("[Error]: main.titles can either be NULL or same length as bigwigs!")
    stop_quietly()
  }
  
  if (is.null(legend.title)){
    message("[Error]: Enter a title for the legend!")
    stop_quietly()
  }
   
  #######################################################################################################################################
  
  ## resize the target ranges according to window.size:
  tmp.target0 <- GenomicRanges::resize(target.gr, fix = "center", width = 1)
  tmp.target1 <- GenomicRanges::resize(target.gr, fix = "center", width = window.size)
  
  ## Now prepare normalizedMatrix for every bigwig:
  
  tmp.nM <<- mclapply(bigwigs, mc.cores=detectCores()-1, function(x){
    
    ## use rtracklayer to load only the relevant part of the bigwig:
    tmp.bigwig <- rtracklayer::import(x, format = "BigWig", selection = BigWigSelection(tmp.target1))
    
    ## normalizedMatrix:
    nM <- normalizeToMatrix(signal = tmp.bigwig, 
                            target = tmp.target0, 
                            background = 0, 
                            keep = c(0, 0.99), ## to exclude any outliers from skewing the plot 
                            target_ratio = 0,
                            mean_mode = "w0", 
                            value_column = "score", 
                            extend = window.size/2)
    return(nM)
    }
  )
  
  ## If scale.all then all plots will have the same range/axis limits.
  ## For this get these limits:
  if (scale.all){
    
    ## for the profile plot => min and max of colMeans
    tmp.r <<- sapply(tmp.nM, function(r) {
      
      tmp.cm <- colMeans(r)
      tmp.min <- min(tmp.cm)
      tmp.max <- max(tmp.cm)
      data.frame(Min=tmp.min, Max=tmp.max)
    })
    
    ## Use 0 and 10% added to tmp.upper as limits:
    tmp.ylims.profile <- c( round(( min(unlist(tmp.r[1,])) - (min(unlist(tmp.r[1,])) * 0.2) ), digits = 3), 
                            round(( max(unlist(tmp.r[2,])) + (max(unlist(tmp.r[2,])) * 0.2) ), digits = 3) 
    )

    ## now for the heatmap ranges
    
    tmp.quant <<- sapply(tmp.nM, function(r) {
                  qu <- quantile(r, c(0, .99))
                  data.frame(lower=qu[1], upper=qu[2])
                  })
    ## the lowest of the 0th quantile (=the min) and the max of the 99th qpercentile are the limits:
    tmp.lower.htmp  <- min(unlist(tmp.quant[1,]))
    tmp.upper.htmp  <- max(unlist(tmp.quant[2,]))
    col_fun = colorRamp2(c(tmp.lower.htmp, tmp.upper.htmp), c(col.low, col.high))
    
  }
    
  ## Now plot everything:
  trashy.can <- mclapply(1:length(tmp.nM), mc.cores = detectCores()-1,function(p){
      
    if (is.null(main.titles)) tmp.title = ""
    if (!is.null(main.titles)) tmp.title <- main.titles[p]
      
    ## heatmap function:
    enrHtmp <- EnrichedHeatmap(
      mat = tmp.nM[[p]], 
                               
      pos_line = FALSE, ## no dashed lines around the start
                               
      border = FALSE,   ## no ugly box around heatmap
                               
      col = col_fun,    ## color gradients
                               
      column_title = tmp.title,## column title 
      column_title_gp = gpar(fontsize = font.size.title, fontfamily = font.fam),
                       
      ## use raster to make high quality but powerpoint-ready plot as pdf
      use_raster = TRUE, raster_quality = 10, raster_device = "png",
                               
      rect_gp = gpar(col = "transparent"), ## no bg color (transparent)
                               
      ## legend:
      heatmap_legend_param = list(
        legend_direction = "horizontal",
        title = legend.title),
                               
      ## profile plot on top
        top_annotation = HeatmapAnnotation(
          enriched = anno_enriched(
            gp = gpar(col = "black", lty = 1, lwd=2),
            ylim = tmp.ylims.profile,
            col="black",
            axis_param = list(
              at = tmp.ylims.profile,
              labels = as.character(tmp.ylims.profile),
              side = "right",
              facing = "outside")))
    )
    
    ## draw it defining padding distances and annotation/heatmap labels
    if (!is.null(bigwig.names)) tmp.name <- bigwig.names[p]
    
    ## if no names are given use the basename of the bigwig as plot/legend name:
    if (is.null(bigwig.names)) tmp.name <- tools::file_path_sans_ext(basename(bigwigs[p]))
    
    GetDate <- function(){ gsub("^20", "", format(Sys.Date(), "%Y%m%d")) }
    
    pdf(paste0(pdf.dir, "/", GetDate(), "_", tmp.name, "_on_", Target.name, ".pdf"), width = 2, height = 6)
      
    draw(enrHtmp,
         heatmap_legend_side = "bottom", 
         annotation_legend_side = "bottom",
         padding = unit(c(3, 3, 3, 3), "mm")
    ); dev.off()
      
    }
  ) ## end of plotting mclapply
    
} ## end of if(scale.all == TRUE)

  
gr <- makeGRangesFromDataFrame(fread("~/IMTB/Project_Leukemia/ATACseq_May2019/Lists/Clusters/190804_cluster2.bed"), seqnames.field = "V1", start.field = "V2", end.field = "V3")

plot.EnrHtmp(bigwigs = c("~/IMTB/C1_mean_TMM.bigwig",
                         "~/IMTB/C2_TMM.bigwig",
                         "~/IMTB/C3_TMM.bigwig",
                         "~/IMTB/C4_TMM.bigwig"),
             bigwig.names = c("Celltype1", "Celltype2", "Celltype3", "Celltype4"), 
             target.gr = gr, Target.name = "cluster2", window.size = 5000, scale.all = TRUE, main.titles = c("Celltype1", "Celltype2", "Celltype3", "Celltype4"), pdf.dir = "./")
 
  
