#' dend_meta_plot
#' 
#' 
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns. Requires cluster_id which needs to be sequential in order of the dendrogram.
#' @param dend Dendrogram object. 
#' @param section_wedges Default is NULL. Use annotation to separate  Can be used to generate lines between sections to divide leaves of dendrogram, e.g. separating subclass, class. 
#' @param bar_variables base name of variables to be represented as bargraphs below dendrogram. Annotation variables need to be represented as \_id, \_label, \_color in anno.
#' @param nr_hist plotting of cluster size below dend. Default is TRUE.
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' @example_data:
#'  
#' load("data/dend_meta_example/dend.rda")
#' load("data/dend_meta_example/anno.df.rda")
#' 
#' 
#' @usage dend_meta_plot <- plot_dend_meta(anno=anno.df, dend=dend, section_wedges="class", bar_variables=c("region", "gender", "batch"), nr_hist=TRUE)
#' 
#'  
#'    

## to remove later ##
#load("//allen/programs/celltypes/workgroups/mct-t200/Manuscript/2019_Yao/common/dend.rda")
#load("//allen/programs/celltypes/workgroups/mct-t200/Manuscript/2019_Yao/common/anno.df.rda")

#load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_TH_20190513/dend.rda")
#load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_TH_20190513/anno.df.rda")



plot_dend_meta <- function(anno,
                           dend,
                           section_wedges=NULL,
                           bar_variables=NULL,
                           nr_hist=TRUE,
                           panel_width=0.2
                           ) {
  
  
  # required libraries
  #library(dendextend)
  library(dplyr)
  library(ggplot2)
  
  rect.offset= c()
  panel_pad <- 0.02
  
  
  # cluster_id order doesn't match dendrogram order
  n_clusters <- max(anno$cluster_id)
  
  # convert to ggdend
  dend_gg <- dendextend::as.ggdend(dend)
  # extract segments for separate plotting later
  dend_seg <- dend_gg$segments
  
  
  if(!is.null(section_wedges)){
    var_wedge <- paste0(section_wedges,"_id")
    sections <- anno %>%
    group_by_at(var_wedge) %>%
    summarise(x = min(cluster_id - 0.5),
              xend = n_clusters + 0.5)
  
    sections = sections %>% arrange(x)
    tmp = sections$x
    section.length = c(tmp[-1], max(sections$xend)) - tmp
    select = which(section.length >= 3)
    select = sort(union(select, select+1))
    sections = sections[select,]
    
    wedge_lines <- data.frame(x = unique(c(sections$x, sections$xend)),
                              y = 0,
                              yend = -1.8) %>%
                  mutate(xend = x)
  }
  
  # set padding of initial graph below dend
    offset= panel_pad
    rect.offset= c(rect.offset, offset)
    panel_width = panel_width 

    
  ############
  # generating stacked bar graphs for selected metadata
  ############ 
    
  if(!is.null(bar_variables)) {
    var_rects <- list()
    

    
    grouping_id <- paste0(bar_variables[1], "_id")
    grouping_label <- paste0(bar_variables[1], "_label")
    grouping_color <- paste0(bar_variables[1], "_color")
    cluster_id <- "cluster_id"
    # platform rectangles
    rect <- anno %>%
      select(cluster_id, cluster_label, cluster_color, grouping_id,  grouping_label, grouping_color) %>%        
      group_by_(cluster_id, grouping_id,  grouping_label) %>%
      mutate(ly_n = n()) %>%
      ungroup() %>%
      group_by(cluster_id) %>%
      arrange_(grouping_id) %>%
      mutate(cluster_n = n(),
             ly_frac = ly_n/cluster_n) %>%
      unique() %>%
      arrange_(grouping_id) %>%
      mutate(ly_cum_frac = cumsum(ly_frac)) %>%
      ungroup() %>%
      arrange_(cluster_id, grouping_id) %>%
      group_by(cluster_id) %>%
      mutate(xmin = cluster_id - 0.5,
             xmax = cluster_id + 0.5,
             ymax = -offset - lag(ly_cum_frac, default = 0)* panel_width,
             ymin = -offset - ly_cum_frac * panel_width)
    
    var_rects[[1]] <- as.data.frame(rect)  
    
    if(length(bar_variables) >1){
      for(i in 2:length(bar_variables)) {
      
        # add panel padding for each bargraph to be plotted
        offset = offset + panel_width + panel_pad
        rect.offset= c(rect.offset, offset)
        panel_width =0.2
        
        grouping_id <- paste0(bar_variables[i], "_id")
        grouping_label <- paste0(bar_variables[i], "_label")
        grouping_color <- paste0(bar_variables[i], "_color")
        cluster_id <- "cluster_id"
        # platform rectangles
        rect <- anno %>%
              select(cluster_id, cluster_label, cluster_color, grouping_id,  grouping_label, grouping_color) %>%        
              group_by_(cluster_id, grouping_id,  grouping_label) %>%
              mutate(ly_n = n()) %>%
              ungroup() %>%
              group_by(cluster_id) %>%
              arrange_(grouping_id) %>%
              mutate(cluster_n = n(),
                     ly_frac = ly_n/cluster_n) %>%
              unique() %>%
              arrange_(grouping_id) %>%
              mutate(ly_cum_frac = cumsum(ly_frac)) %>%
              ungroup() %>%
              arrange_(cluster_id, grouping_id) %>%
              group_by(cluster_id) %>%
              mutate(xmin = cluster_id - 0.5,
                     xmax = cluster_id + 0.5,
                     ymax = -offset - lag(ly_cum_frac, default = 0)* panel_width,
                     ymin = -offset - ly_cum_frac * panel_width)
          
        var_rects[[i]] <- as.data.frame(rect)  
      }
    }
    
  offset = offset + panel_width + panel_pad
  rect.offset= c(rect.offset, offset)
  panel_width =1.5*panel_width
  }
  
    
    
  ############
  # generating cell nr. histogram data
  ############ 
    
  if(nr_hist==TRUE){
    n_rects <- anno  %>%
      group_by(cluster_id, cluster_color, cluster_label) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(adj_n = log10(n)) %>%
      mutate(f = adj_n/3) %>%
      mutate(xmin = cluster_id - 0.5,
             xmax = cluster_id + 0.5,
             ymin = -offset - f  * panel_width,
             ymax = -offset)
    
    
    ## still to fix guide labels 
    n_guides <- data.frame(y = seq(-offset - panel_width, -offset, by = 1/5 * panel_width),
                           x = 0.5,
                           xend = n_clusters + 1,
                           label = seq(5, 0, by = -1)) %>%
              mutate(yend = y)
  
  offset = offset + panel_width + panel_pad*2
  }
  

  
  dend_leaves <- dend_gg$labels %>%
    mutate(cluster_label = label,
           y = - offset, cex=0.3)
  
  ############
  # plotting of the dendrogram 
  ############
  
  
  # dend segments
  flat_plot <- ggplot() +
    geom_segment(data = dend_seg,
                 aes(x = x,
                     xend = xend,
                     y = y,
                     yend = yend,
                     size = lwd,
                     color = col),
                 lineend = "square") 
  
  
  # different bar plots with metadata
  if(!is.null(bar_variables)) {
    lapply(var_rects, function(df) {
      
      flat_plot <<- flat_plot + geom_rect(data = df,
                                         aes(xmin = xmin,
                                             xmax = xmax,
                                             ymin = ymin,
                                             ymax = ymax,
                                             fill = df[,6]) )
      
      } )
    }
  
  if(nr_hist==TRUE){
    flat_plot <- flat_plot +
      geom_rect(data = n_rects,
                aes(xmin = xmin,
                  xmax = xmax,
                  ymin = ymin,
                  ymax = ymax,
                  fill = cluster_color)) 
    }
    
  # N Cells labels
  flat_plot <- flat_plot +  
    geom_text(data = n_guides,
              aes(x = 0,
                  y = y,
                  label = label),
              size = 2,
              hjust = 1) +
    geom_segment(data = n_guides,
                 aes(x=x,
                     xend = xend,
                     y=y,
                     yend = yend),
                 linetype="dashed") +
    # Leaf Labels
    geom_text(data = dend_leaves,
              aes(x = x,
                  y = y,
                  label = label,
                  color = col),
              angle = 90,
              hjust = 1,
              vjust = 0.3,
              size = 2) 
  
  if(!is.null(section_wedges)){ 
    #Vertical class separators
    flat_plot <- flat_plot +
      geom_segment(data = wedge_lines,
                   aes(x = x,
                       xend = xend,
                       y = y,
                       yend = yend)) 
    
  }  
  
  flat_plot <- flat_plot +
        scale_size(range = c(0.5, 1)) +
        scale_color_identity() +
        scale_fill_identity() +
        scale_y_continuous(expand = c(0.25,0)) +
        scale_x_continuous(limits = c(-1,n_clusters + 1)) +
        theme_void()
  
  
  
  
  return(flat_plot)
  
  
}
