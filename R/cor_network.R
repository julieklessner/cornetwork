#' Correlation Network
#'
#' This function allows you to build a network based on correlations
#' @param data (required) A phyloseq object.
#' @param threshold.cor Threshold value for the correlation (default: 0.8).
#' @param threshold.count Threshold value for abundance reads (default: 10).
#' @param show.label Show vertex names as labels (default: TRUE).
#' @param scale.abund Scale the node size according to average abundance (default: FALSE).
#' @param node.size Size of nodes if scale.abund=FALSE (default: 5) 
#' @param scale.cor Scale the edge color according to strength of correlation (default: FALSE).
#' @param edge.size Size of egdes if scale.cor=FALSE (default: 0.5)
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (default: "Phylum")
#' @param tax.class Converts a specific phyla to class level instead, e.g. "p__Proteobacteria" (default: "none").
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: "best").
#' @param highlight.OTU Highlight OTUs using see-throgh nodes, e.g. c("MiDAS_37", "MiDAS_73") (default: none).
#' @param highlight.color Color of highlights, e.g. c("purple", "green") (default: none).
#' @param highlight.size Size of highlights (default: 10).
#' @param highlight.label Highlight OTUs using node labels, e.g. c("MiDAS_37", "MiDAS_73") (default:none).
#' @param add.tax.info Show taxonomic information (default: FALSE).
#' @param show.edge.label Show correlation values as edge labels (default: FALSE).
#' @keywords network correlation ggnetwork
#' @import network
#' @import ggnetwork
#' @import sna
#' @export

cor_network <- function(data, threshold.cor= 0.8, threshold.count=10, show.label=TRUE, scale.abund=FALSE, node.size=5, scale.cor=FALSE, edge.size = 0.5, highlight.OTU=NULL, highlight.color=NULL, highlight.size=10, highlight.label=NULL, tax.aggregate="Phylum", tax.add=NULL, tax.class=NULL, tax.empty="best", add.tax.info=FALSE, show.edge.label=FALSE){
  
  data0 <- list(abund = as.data.frame(otu_table(data)@.Data),
                tax = data.frame(tax_table(data)@.Data, OTU =   rownames(tax_table(data))),
                sample = suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))
  
  # Clean and rename taxonomy -------------------------------------------------------
  
  data1 <- amp_rename(data = data0,
                      tax.class=tax.class, 
                      tax.empty=tax.empty, 
                      tax.level = tax.aggregate)
  
  # Divide data to seperate data frames --------------------------------------------
  
  abund <- data1[["abund"]]
  tax <- data1[["tax"]]
  sample <- data1[["sample"]]
  
  # Subsetting by threshold.count-----------------------------------------------------  
  
  abund[abund <= threshold.count] <- NA
  
  # Producing correlation matrix-----------------------------------------------------
  
  cormat <- cor(t(abund), method = "pearson")
  
  # Revoming self-correlation-------------------------------------------------------  
  
  diag(cormat) <- NA  
  
  # Removing "doubles"---------------------------------------------------------------
  
  cormat[lower.tri(cormat)] <- NA
  
  # Converting to dataframe for easier subsetting------------------------------------
  
  cor2 <- melt(cormat, value.name="Rval")
  cor3 <- cor2[!is.na(cor2$Rval), ]
  
  # Subsetting by threshold.cor------------------------------------------------------
  
  cor4 <- subset(cor3, abs(cor3$Rval) >= threshold.cor)
  
  #Building network ----------------------------------------------------------------
  
  net <- network(cor4[ ,1:2 ], directed=FALSE)
  
  #Assigning Rval to edges
  
  set.edge.attribute(net, "Rval", cor4$Rval)
  
  #Assigning abundance values to nodes---------------------------------------------
  
  names = ggnetwork(net)$vertex.names %>% unique() %>% as.character()
  names2 = ggnetwork(net)$vertex.names %>% as.data.frame()
  
  abund_cor <- abund[names, ] 
  abund_mean <- rowMeans(abund_cor) 
  abund_mean2 <- abund_mean %>% as.data.frame()
  
  net %v% "Avg_Abundance" <- abund_mean
  
  #Assigning taxonomic information to nodes
  
  tax <- data.frame(tax, Display = tax[,tax.aggregate])
  rownames(tax) <- tax$OTU
  
  tax_info = tax[names, ] 
  tax_info_display = tax_info[,"Display"] %>% as.character()
  
  set.vertex.attribute(net, "taxo_info", tax_info_display)
  
  #Visualizing using ggnetwork-------------------------------------------------------
  
  
  p <- ggplot(net, arrow.gap = 0, aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_nodes()+
    theme_blank()
  
  #Options for edges---------------------------------------------------------------- 
  if (scale.cor == T){
    p <- p + geom_edges(alpha = 0.6, size=2, aes(color = Rval)) #+ scale_color_gradient2()
  } else {
    p <- p + geom_edges(alpha = 0.6, size=edge.size)
  }
  
  #Options for nodes----------------------------------------------------------------
  if (scale.abund == T){
    p <- p + geom_nodes(aes(size= Avg_Abundance))
  } else {
    p <- p + geom_nodes(size = node.size)
  }
  
  if(add.tax.info == T & scale.cor == T){
    print("Mapping a color to both a vertex attribute and an edge attribute will violate the grammer of graphics")
  } else if (add.tax.info == T & scale.abund == T) {
    p <- p + geom_nodes(aes(size= Avg_Abundance, col=taxo_info))
  } else if (add.tax.info == T & scale.abund == F) {
    p <- p + geom_nodes(aes(size= node.size, col=taxo_info))
  }
  
  #Options for labels----------------------------------------------------------------
  if(show.edge.label == T){
    p <- p + geom_edgelabel(aes(label = signif(Rval, digits = 2)))
  }
  
  if (show.label == T){
    p <- p + geom_nodetext(aes(label=vertex.names), nudge_x = 0.05) 
  }
  
  #Options for highlights------------------------------------------------------------
  
  if(!is.null(highlight.OTU)){
    p <- p + geom_nodes(data = function(x) { x[ x$vertex.names == highlight.OTU, ]}, col=highlight.color, size= 10, alpha = 0.2)
  }
  
  if(!is.null(highlight.label)){
    p <- p + geom_nodelabel_repel(aes(label = vertex.names),
                                  box.padding = unit(1, "lines"),
                                  data = function(x) { x[ x$vertex.names == highlight.label, ]})
  }  
p
  }