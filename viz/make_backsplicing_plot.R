# XD: modified based on leafcutter/leafcutter/R/make_cluster_plot.R

#' Make cluster level plot
#'
#' @import ggplot2
#' @import foreach
#' @importFrom gridExtra grid.arrange
#' @import intervals
#' @export

library(tidyverse);
library(intervals);
library(foreach)
library(gridExtra)
library(Hmisc);

make_backsplicing_plot <- function(
  cluster_to_plot, 
  main_title = NA, 
  exons_table = NULL, 
  meta = NULL, 
  counts = NULL, 
  introns = NULL,
  snp_pos=NA){
  
  if( is.null(cluster_to_plot)){
    print("no cluster selected!")
  }
  #for testing!
  #cluster_to_plot <- "clu_8845"
  #cluster_to_plot <- "clu_57214"; main_title = c("FAM47E-STBD1", "clu_57214"); meta = meta; exons_table = exons_table; counts = counts;introns = introns;  snp_pos=NA
  #cluster_to_plot <- "clu_1234567"; main_title = c("test", "clu_1234567"); meta = meta; exons_table = exons; counts = counts;introns = circRNA;  snp_pos=NA
  #cluster_to_plot <- "chr2_231940224_231951895"; main_title = c("chr2_231940224_231951895","kgp9431896_A:G"); meta = meta; exons_table = exons; counts = counts;introns = circRNAs;  snp_pos=NA
  
  meta$group=as.factor(meta$group)
  group_names=levels(meta$group)
  
  stopifnot(cluster_to_plot %in% rownames(counts))
  # create variables for later
  y <- t(counts[cluster_to_plot, ,drop = FALSE])
  #x <- numeric(nrow(y))+1
  x <- meta$group
  length_transform <- function(g){ log(g+1) }
  introns_to_plot <- introns[ introns$clusterID == cluster_to_plot, ,drop = FALSE] 
  summary_func <- colSums
  legend_title <- "Mean counts"
  alphabet <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z")
  
  
  junction_colour <- "red"
  cryptic_colour <- "pink"

  # convert colnames(y) into intron meta data
  intron_meta = do.call(rbind, strsplit(colnames(y), "_"))
  colnames(intron_meta) = c("chr", "start", "end")
  intron_meta = as.data.frame(intron_meta, stringsAsFactors = F)
  intron_meta$start = as.numeric(intron_meta$start)
  intron_meta$end = as.numeric(intron_meta$end)
  
  intron_meta$verdict <- introns_to_plot$verdict[match(paste(intron_meta$start,intron_meta$end), paste(introns_to_plot$start, introns_to_plot$end) ) ]

  ## XD: comment out
  # # give alphabetical rank based on dPSI
  # intron_meta$dPSI <- introns_to_plot$deltapsi[match(paste(intron_meta$start,intron_meta$end), paste(introns_to_plot$start, introns_to_plot$end) ) ]
  # ranks <- alphabet[1:nrow(intron_meta)]
  # absPSI <- intron_meta$dPSI[ order( abs(intron_meta$dPSI), decreasing=TRUE) ]
  # intron_meta$rank <- ranks[match(intron_meta$dPSI, absPSI)]
  
  # make sure intron_meta has "chr" in front of chromosome name so it plays nice with the exon table
  if( all( ! grepl("chr", intron_meta$chr))){
    intron_meta$chr <- paste0("chr", as.character(intron_meta$chr))
  }

  #print(intron_meta)

  new_theme_empty <- theme_bw(base_size = 15 )
  new_theme_empty$panel.background = element_rect(fill="white", colour = "white")
  new_theme_empty$line <- element_blank()
  new_theme_empty$rect <- element_blank()
  new_theme_empty$strip.text <- element_blank()
  new_theme_empty$axis.text <- element_blank()
  
  groups=sort(unique(x))
  
  max_log=.5*ceiling(2*log10( 1+max( unlist( foreach (tis=groups) %do% { intron_meta$counts=summary_func(y[ tis==x,,drop=F]) } ) ) ))
  
  breaks=if (max_log <= 2.5) seq(0,max_log,by=0.5) else seq(0,ceiling(max_log),by=1)
  limits=c(0.0,max_log)
  
  intron_meta$id=as.factor(1:nrow(intron_meta)) # number each junction
  temp=intron_meta[,c("id","start","end")]
  m=reshape2::melt(temp, id.vars = "id") # melt to get list of all coordinates
  
  s=unique(m$value) # get unique start and end values
  if (!is.na(snp_pos)) s=c(s,snp_pos)
  s=sort(s) 
  d=s[2:length(s)]-s[1:length(s)-1] # get the difference between pairs of coordinates
  trans_d <- length_transform(d) # e.g. log(d+1), sqrt(d), d ^ .33 # apply a trasnformation function - doesn't work on negative numbers!
  coords <- c(0,cumsum(trans_d))
  names(coords)=s
  
  snp_coord=coords[as.character(snp_pos)]
  
  total_length=sum(trans_d) # ==max(coords)
  my_xlim=c(-.2*total_length,1.2*total_length)

  ############## 
  # PLOT SETTINGS
  ###############
  
  min_height=0
  max_height=0
  curv <- 0.45
  min_exon_length <- 0.5
  maxratio=0
  minratio=1.0
  yFactor = 0.65   # originally set as 0.65
  yConstant = -0.25 # originally set as 0.5
  labelTextSize=3.5 # orignally set as 5
  curveMax = 10
  curveExponent = 2
  yOffset = 0
  centreLineWidth = 3   # horizontal white line to clean up the edges of the curves
  len = 500
  
  
  mainPalette <- c(junction_colour, cryptic_colour)
  
  # XD: below funciton is good for multiple columns (e.g. intron cluster), not for single column
  # sweep is dividing each entry in each row by the sum of all entries in that row and then apply is finding the mean value of each column
  #summary_func=function(a) apply( sweep(a,1,rowSums(a),"/"),2, function(g) mean(g, na.rm=T) ) 
  
  # XD: just take the mean expression
  summary_func=function(a) apply(a,2,mean)

  last_group=groups[length(groups)]
  
  plots <- list()
  for( fancyVar in 1:length(groups) ){

    intron_meta$counts=summary_func(y[ groups[fancyVar]==x,,drop=F])
    intron_meta$prop=intron_meta$counts#/sum(intron_meta$counts)  # this has been changed

    group_sample_size=sum(groups[fancyVar]==x)
    #print(intron_meta$counts)

    allEdges=do.call(rbind,foreach (i=1:nrow(intron_meta)) %do% {
      if (intron_meta$counts[i]==0) return(NULL)
      start=coords[ as.character(intron_meta$start[i]) ]
      end=coords[ as.character(intron_meta$end[i]) ]
      l=end-start
      h=(1+sqrt(l)) * ( (i %% 2)*2-1 )
      min_height=min(min_height,h)
      max_height=max(max_height,h)
      edge = data.frame(Hmisc::bezier(x=c(start, start-0.3*total_length, (start + end)/2, end+0.3*total_length, end), y=c(0, h*2/3, h, h*2/3, 0),evaluation = len))
      edge$Sequence <- (1+intron_meta$counts[i]) * sin( seq(0,pi,length.out=len) ) # For size and colour weighting in plot
      edge$log10counts=(1+intron_meta$counts[i])
      edge$verdict <- intron_meta$verdict[i]
      edge$Group <- i
      edge
    })
    
    # XD: unused?
    # MAXcounts <- max(c(max(with(allEdgesP,log10counts)),max(with(allEdges,log10counts))))

    # add 1 or -1 to the function, in case there is only Positive or Negative introns
    #YLIMP <- 1.25 * max( allEdgesP$ytext, .5)  
    #YLIMN <- 1.25 * min( allEdges$ytext, -.5)

    #print(paste("loop:",fancyVar, length(allEdgesP), length(allEdges), group_sample_size))
    
    g <- ggplot()
      
    g <- g + geom_path(data=allEdges, aes(x = x, y = y, group = Group, colour=verdict, size = Sequence, alpha=.9)) +
      new_theme_empty +
      # make the y axis label the group
      ylab(paste0(groups[fancyVar],"\n(n=",group_sample_size,")")) +
      xlab("") +
      xlim(my_xlim) +

      # horizontal line - smooth out the ends of the curves
      geom_hline(yintercept=0, size = centreLineWidth, colour = "white") +
      geom_hline(yintercept=0,alpha=.9, size=1) 
      
      # label the junctions
    if(length(allEdges)>0) g <- g +
      geom_label(data=data.frame(x=mean(range(allEdges$x)),y=max(allEdges$y),label=unique(allEdges$log10counts)),
                 aes(x=x,y=y,label=label),
                 size= labelTextSize, label.size = NA, parse=TRUE,
                 fill = "white", colour = "black", label.r = unit(0.3,"lines"),
                 label.padding = unit(0.3,"lines") )

      g <- g +
      #ylim(YLIMN,YLIMP) +
      scale_size_continuous(limits=c(0,10),guide='none')
    # is this used for anything? color is currently set to clu which doesn't change for each junction

    if (!is.na(snp_coord)) {
      df=data.frame(x=snp_coord,xend=snp_coord,y=0,yend=max_height*1.1)
      g <- g + geom_segment(data=df,aes(x=x,y=y,xend=xend,yend=yend)) #+geom_vline(xintercept=snp_coord)
    }
      
    plots[[fancyVar]] <- g  # return the plot
  }

  # ADDING EXTRA STUFF TO THE PLOTS
  
  df <- data.frame(x=coords, xend=total_length*(s-min(s))/(max(s)-min(s)), y=0, yend=min_height)
  # Control segment between diagram to exon
  # if (! is.null( exons_table) ){
  # plots[[length(plots)]] <- plots[[length(plots)]] + geom_segment(data=df, aes(x=x,y=y,xend=xend,yend=yend),alpha=0.2, lty=2)
  # }
  
  # ADDING EXON ANNOTATION
  
  if (!is.null(exons_table)) {
    
    # find the exons
    exons_chr <- exons_table[ exons_table$chr==intron_meta$chr[1], ] # subset by chr
    stopifnot( nrow(exons_chr) > 0)
    exons_here=unique(exons_chr[ exons_chr$start <= max(s) & min(s) <= exons_chr$end, ])
    # XD: this is to control which exons to display

    ## XD: code below are optional
    # # if exons found - remove exons that don't actually start or end with a junction
    # # and any repeated exons or any exons larger than 500bp
    if( nrow(exons_here) > 0 ){
      exons_here <-  unique(
        exons_here[ ( exons_here$start %in% intron_meta$start |
                        exons_here$end %in% intron_meta$end ) &
                      ( exons_here$end - exons_here$start <= 500 |
                        exons_here$start == min(intron_meta$start) |
                        exons_here$end == max(intron_meta$end) ), ]
      )
    }
    
    # if any exons survive the cull
    if ( nrow( exons_here) > 0) {
      exons_here$gene_name=factor(exons_here$gene_name)
      n_genes <- seq(1, length( levels(exons_here$gene_name) ) )
      gene_name_df <- data.frame( x= 0.2*total_length,#(n_genes * total_length) / (max(n_genes) + 1), 
                                  y=YLIMN, 
                                  label=rev(levels(exons_here$gene_name))
                                )
      # count the occurences of each gene's exons and order by it
      gene_name_df$N <- table(exons_here$gene_name)[ gene_name_df$label ]
      gene_name_df <- gene_name_df[ order(gene_name_df$N, decreasing = TRUE), ]
      
      # for gene name colouring (black by default, other colours for multiple genes)
      cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
      # sort out colour pallette - most common gene goes first!
      cbbPalette <- cbbPalette[ 1:length(gene_name_df$label) ]
      names(cbbPalette) <- gene_name_df$label
      # add junction colours 
      mainPalette <- c(cbbPalette, junction_colour, cryptic_colour)
      names(mainPalette)[ (length(mainPalette)-1):length(mainPalette) ] <- c("circRNA","ciRNA")

      # fit exons within the cluster scale
      # s is the sorted unique values of start and end
      invert_mapping=function(pos){
        if (pos %in% s) coords[as.character(pos)] else
          if (pos < min(s)) my_xlim[1] else
            if (pos > max(s)) my_xlim[2] else {
              w=which( pos < s[2:length(s)] & pos > s[1:(length(s)-1)] )
              stopifnot(length(w)==1)
              coords[w] + (coords[w+1]-coords[w])*(pos - s[w])/(s[w+1]-s[w])
            }
      }
      

      exon_df <- data.frame( x=sapply(exons_here$start,invert_mapping), 
                               xend=sapply(exons_here$end,invert_mapping), 
                               y=0, 
                               yend=0,
                               label = exons_here$gene_name)
        
      # alter exon sizes to conform to a minimum exon length
      # XD: extend the end of exon if it's not part of the edge
      exon_df[ (exon_df$xend - exon_df$x) < min_exon_length & exon_df$x %in% c(allEdgesP$start,allEdges$start), ]$xend <- exon_df[ (exon_df$xend - exon_df$x) < min_exon_length & exon_df$x %in% c(allEdgesP$start,allEdges$start), ]$x + min_exon_length
      exon_df[ (exon_df$xend - exon_df$x) < min_exon_length & exon_df$xend %in% c(allEdgesP$end,allEdges$end), ]$x <- exon_df[ (exon_df$xend - exon_df$x) < min_exon_length & exon_df$xend %in% c(allEdgesP$end,allEdges$end), ]$xend - min_exon_length
      
      # remove exons that are now duplicates - overlapping exons that extend outside of the plotting space and are squished to same coordinates
      exon_df <- exon_df[ !duplicated( paste(exon_df$x, exon_df$xend) ),]
      
      ## STRANDING
      # pick principal gene - the one that contributes the most exons to the cluster will be at the top of the df
      principal_gene <- gene_name_df$label[1]
      gene_strand <- unique(exons_here[exons_here$gene_name == principal_gene,]$strand)
      if( length(gene_strand)  > 1 & !is.na(gene_strand) ){
        gene_strand <- NULL
      }else{
        # assign strand arrows based on lengths of intron
        exon_intervals <- Intervals( matrix(data = c(exon_df$x, exon_df$xend), ncol = 2) )
        #exon_intervals <- intervals::interval_union( exon_intervals )
        intron_intervals <- intervals::interval_complement( exon_intervals )
        intron_intervals <- intron_intervals[ 2:(nrow(intron_intervals)-1),]
        # strand arrows should be placed between exons
        strand_df <- as.data.frame(intron_intervals)
        # remove strand arrows where the introns are too small
        strand_df <- strand_df[ (strand_df$V2 - strand_df$V1) > 0.025 * total_length,]
        strand_df$midpoint <- strand_df$V1 + (strand_df$V2 - strand_df$V1) / 2
        
        group <- c()
        for(i in 1:nrow(strand_df)){
          group <- c(group, rep(i,2))
        }
        # make strand_df
        if( gene_strand == "+"){
          strand_df <- data.frame( x = c(rbind(strand_df$V1, strand_df$midpoint)),
                                   group = group, 
                                   y = 0);
          strand_pos = "last";
        }
        if( gene_strand == "-"){
          strand_df <- data.frame( x = c(rbind(strand_df$midpoint, strand_df$V2)),
                                   group = group,
                                   y = 0);
          strand_pos = "first";
        }
        
        # add strand arrows to plot
        for (i in 1:length(plots) ){
          plots[[i]] <- plots[[i]] +
            geom_line( data = strand_df, 
                       aes( x = x, y = y, group = group ), colour = "black", size=1, 
                       arrow = arrow(ends = strand_pos, type = "open", angle = 30, length = unit(0.1, units = "inches" ))) 
        }
      }

      # add exons to plots
      for (i in 1:length(plots) ){ 
        plots[[i]] <- plots[[i]] +
          geom_segment( data=exon_df, aes(x=x,y=y,xend=xend,yend=yend, colour = label), alpha=1, size=6) +
          geom_segment( data = exon_df, aes(x = x, xend = x+0.01, y = y, yend = yend), colour = "white", size = 6, alpha = 1) +
          geom_segment( data = exon_df, aes(x = xend-0.01, xend=xend, y = y, yend = yend), colour = "white", size = 6, alpha = 1)
      }
    }
    
  }

  # TITLE
  
  for (i in 1:length(plots) ){
    plots[[i]] = plots[[i]] + scale_colour_manual("", values = mainPalette ) + scale_alpha(guide="none",range = c(0.1, 1)) + 
guides(colour='none') +  # hide color palette
      ggtitle(paste(main_title, collapse = " at ")) + theme(plot.title = element_text(colour="white", size = 15))
    
    if(i==1) plots[[i]] = plots[[i]] + theme(plot.title = element_text(colour="black"))
    if(i==length(plots)) plots[[i]] = plots[[i]] + guides(colour = guide_legend(override.aes = list(size=.8))) + 
        theme(legend.box = "horizontal", legend.direction = "horizontal",
              legend.justification=c(1,0), legend.position=c(1,-.4)) 

  }
  
  # ARRANGE PLOTS
  
  do.call( gridExtra::grid.arrange, c(plots, list(ncol=1)))

}