#' Correlation Network
#'
#' This function allows you to build a network based on correlations
#' @param data (required) A phyloseq object.
#' @param threshold.cor Threshold value for the correlation (default: 0.8).
#' @param threshold.abund Threshold value for relative abundance (default: 0.01).
#' @param threshold.obs Threshold value for the average relative abundance per observation (default: 0.1).
#' @param threshold.pval Threshold value for the p-value (default: 0.05).
#' @param show.label Show vertex names as labels (default: FALSE).
#' @param scale.abund Scale the node size according to average abundance (default: FALSE).
#' @param node.size Size of nodes if scale.abund=FALSE (default: 5)
#' @param scale.cor Scale the edge color according to strength of correlation (default: FALSE).
#' @param edge.size Size of egdes if scale.cor=FALSE (default: 0.5)
#' @param tax.display The taxonomic level that the data should be displayed (default: "Phylum").
#' @param tax.class Converts a specific phyla to class level instead, e.g. "p__Proteobacteria" (default: "none").
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: "best").
#' @param highlight.OTU Highlight OTUs using see-throgh nodes, e.g. c("MiDAS_37", "MiDAS_73") (default: none).
#' @param highlight.color Color of highlights, e.g. c("purple", "green") (default: none).
#' @param highlight.size Size of highlights (default: 10).
#' @param highlight.label Highlight OTUs using node labels, e.g. c("MiDAS_37", "MiDAS_73") (default:none).
#' @param add.tax.info Show taxonomic information (default: FALSE).
#' @param show.edge.label Show correlation values as edge labels (default: FALSE).
#' @param correlation Choose the method for calculation of correlation from Pearson, Spearman, SparCC, CClasso or Rebacca (default: "Pearson").
#' @param show.top Number of correlations returned if return.output is set to top (default: 25).
#' @param return.output Either "top" or "total" or "plot" (default: "top").
#' @keywords network correlation ggnetwork
#' @import network
#' @import ggnetwork
#' @import sna
#' @import viridis
#' @import glmnet
#' @import MCMCpack
#' @import Hmisc
#' @import pi0
#' @export

amp_network <- function(data, threshold.cor= 0.8, threshold.abund=0.01, threshold.obs=0.1, threshold.pval=0.05, show.label=FALSE, scale.abund=FALSE, node.size=5, scale.cor=FALSE, edge.size = 0.5, highlight.OTU=NULL, highlight.color=NULL, highlight.size=10, highlight.label=NULL, tax.display="Phylum", tax.add=NULL, tax.class=NULL, tax.empty="best", add.tax.info=FALSE, show.edge.label=FALSE, correlation = "Pearson", return.output = "top", show.top = 25){

  data0 <- list(abund = as.data.frame(otu_table(data)@.Data),
                tax = data.frame(tax_table(data)@.Data, OTU =   rownames(tax_table(data))),
                sample = suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))

  # Clean and rename taxonomy -------------------------------------------------------

  data1 <- amp_rename(data = data0,
                      tax.class=tax.class,
                      tax.empty=tax.empty,
                      tax.level = tax.display)

  # Divide data to seperate data frames --------------------------------------------

  abund <- data1[["abund"]]
  tax <- data1[["tax"]]
  sample <- data1[["sample"]]

  # Normalising abundance
  abund_n <- apply(as.matrix(abund), 2, function(x) x / sum(x) * 100) %>%  as.data.frame()

  # Subsetting by threshold.obs------------------------------------------------
  abund_n1 <- abund_n[rowMeans(abund_n) >= threshold.obs,]

  # Subsetting by threshold.abund
  abund_n1[abund_n1 <= threshold.abund] <- 0

  # Subsetting abundance count
  otus <- rownames(abund_n1)
  abunds <- abund[otus, ]

  ##################################
  # Producing correlation matrix   #
  ##################################


  if(correlation == "Pearson"){
    
    # Calculating correlations
    resp <- rcorr(as.matrix(t(abunds)), type = "pearson")
    
    # Extracting correlation matrix
    cormat1 <- resp$r %>% as.matrix()
    
    # Extracting correlation matrix
    pmat <- resp$P %>% as.matrix()
    cormat <- ifelse(pmat <= threshold.pval, cormat1, NA)
  }

  if(correlation == "Spearman"){
    
    # Calculating correlations
    resp <- rcorr(as.matrix(t(abunds)), type = "spearman")
    
    # Extracting correlation matrix
    cormat1 <- resp$r %>% as.matrix()
    
    # Extracting p-values
    pmat <- resp$P %>% as.matrix()
    cormat <- ifelse(pmat <= threshold.pval, cormat1, NA)
  }

  if(correlation == "SparCC"){
    
    # Calculating correlations
    res_cc <- sparcc(t(abunds)) # From sparcc.R
    
    # Extracting correlation matrix
    cormat1 <- res_cc$Cor
    rownames(cormat1) <- rownames(abunds)
    colnames(cormat1) <- rownames(abunds)
    
    # Calculating and extracting p-values
    mat <- cormat1
    mat <- as.matrix(mat)
    n <- ncol(mat)
    pmat <- matrix(NA, n, n)
    diag(pmat) <- 0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        tmp <- t.test(mat[,i], mat[,j], conf.level = 0.95, paired=T)
        pmat[i,j] <- pmat[j,i] <- tmp$p.value
      }}
    
    rownames(pmat) <- rownames(abunds)
    colnames(pmat) <- rownames(abunds)
    
    cormat <- ifelse(pmat <= threshold.pval, cormat1, NA)
  }

  if(correlation == "CClasso"){
    
    # Calculating correlations
    res_cc <- cclasso(t(abunds), counts = TRUE, n_boot = 20) # From CClasso.R
    
    # Extracting correlation matrix
    cormat1 <- res_cc$cor_w %>% as.matrix()
    rownames(cormat1) <- rownames(abunds)
    colnames(cormat1) <- rownames(abunds)
    
    # Extracting p-values
    pmat <- res_cc$p_vals %>% as.matrix
    rownames(pmat) <- rownames(abunds)
    colnames(pmat) <- rownames(abunds)
  
    cormat <- ifelse(pmat <= threshold.pval, cormat1, NA)
  }

  if(correlation == "Rebacca"){
    
    # Calculating correlations
    x.rslt <- rebacca(as.matrix(abunds), nbootstrap=20, N.cores=1) # From REBACCA.R
    
    #Producing p-values
    tau = stability_cutoff(x.rslt$Stability, x.rslt$q, B=50, FWER=threshold.pval) # From REBACCA.R
    x.adj = sscore2adjmatrix(x.rslt$Stability, tau) # From REBACCA.R
    
    #Extracting correlation matrix
    x.est = rebacca_adjm2corr(abunds, x.adj) # From REBACCA.R
    cormat <- x.est$corr %>% as.matrix() # From REBACCA.R
    rownames(cormat) <- rownames(abunds)
    colnames(cormat) <- rownames(abunds)
  }

  # Revoming self-correlation-------------------------------------------------------

  diag(cormat) <- NA

  # Removing "doubles"---------------------------------------------------------------

  cormat[lower.tri(cormat)] <- NA

  # Converting to dataframe for easier subsetting------------------------------------

  cor2 <- melt(cormat, value.name="Cor_val")
  cor3 <- cor2[!is.na(cor2$Cor_val), ]

  # Subsetting by threshold.cor------------------------------------------------------

  cor4 <- subset(cor3, abs(cor3$Cor_val) >= threshold.cor)

  ####################
  # Building network #
  ####################

  net <- network(cor4[ ,1:2 ], directed=FALSE)

  # Assigning correlation values to edges ---------------------------------------------

  set.edge.attribute(net, "Cor_val", cor4$Cor_val)

  ##################################################
  # Assigning normalised abundance values to nodes #
  ##################################################

  names = ggnetwork(net)$vertex.names %>% unique() %>% as.character()
  names2 = ggnetwork(net)$vertex.names %>% as.data.frame()


  ## Matching abundances to correlated abundance------------------------------------
  abund_cor <- abund_n[names, ]

  ## Calculating mean abundance -------------------------------------------------------

  abund_mean <- rowMeans(abund_cor)
  #abund_mean2 <- abund_mean %>% as.data.frame()

  ## Asigning abundance values to vertices --------------------------------------------
  net %v% "Avg_Abundance" <- abund_mean

  ############################################
  # Assigning taxonomic information to nodes #
  ############################################

  tax <- data.frame(tax, Display = tax[,tax.display])
  rownames(tax) <- tax$OTU

  tax_info = tax[names, ]
  tax_info_display = tax_info[,"Display"] %>% as.character()
  tax_info_display2 = tax_info[,"Display"] %>% as.data.frame()

  set.vertex.attribute(net, "taxo_info", tax_info_display)

  ################################
  # Visualizing using ggnetwork  #
  ################################

  p <- ggplot(net, arrow.gap = 0, aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_nodes()+
    theme_blank()

  # Options for edges----------------------------------------------------------------
  if (scale.cor == T){
    p <- p + geom_edges(alpha = 0.6, size=2, aes(color = Cor_val)) +
      scale_color_viridis("Correlation value") #+ scale_color_gradient2()
  } else {
    p <- p + geom_edges(alpha = 0.6, size=edge.size)
  }

  # Options for highlights------------------------------------------------------------

  if(!is.null(highlight.OTU)){
    p <- p + geom_nodes(data = function(x) { x[ x$vertex.names %in% highlight.OTU, ]},
                        col=highlight.color,
                        size= highlight.size,
                        alpha = 0.2)
  }

  if(!is.null(highlight.label)){
    p <- p + geom_nodelabel_repel(aes(label = vertex.names),
                                  box.padding = unit(1, "lines"),
                                  data = function(x) { x[ x$vertex.names %in% highlight.label, ]})
  }

  # Options for nodes----------------------------------------------------------------
  if (scale.abund == T){
    p <- p + geom_nodes(aes(size= Avg_Abundance)) +
      scale_size(name="Mean % Read Abundance", range=c(node.size, (node.size+4)))
  } else {
    p <- p + geom_nodes(size = node.size)
  }

  if(add.tax.info == T & scale.cor == T){
    print("Mapping a color to both a vertex attribute and an edge attribute will violate the grammer of graphics")
  } else if (add.tax.info == T & scale.abund == T) {
    p <- p + geom_nodes(aes(size= Avg_Abundance, col=taxo_info))
  } else if (add.tax.info == T & scale.abund == F) {
    p <- p + geom_nodes(size= node.size, aes(col=taxo_info))
  }

  # Options for labels----------------------------------------------------------------
  if(show.edge.label == T){
    p <- p + geom_edgelabel(aes(label = signif(Rval, digits = 2)))
  }

  if (show.label == T){
    p <- p + geom_nodetext(aes(label=vertex.names), nudge_x = 0.05)
  }


  ##########################################
  # Matching taxonomy to correlation pairs #
  ##########################################

  #Retrieving taxonomic information for Var2
  Var2 <- cor4[, "Var2", drop=FALSE]
  Tax_Var2 <- tax[rownames(tax) %in% cor4$Var2,]
  Tax_Var2 <- Tax_Var2[ , c("OTU", "Display")]
  cor5 <- merge(cor4, Tax_Var2, by.x=2, by.y=1)

  #Retrieving taxonomic information for Var1
  Var1 <- cor4[, "Var1", drop=FALSE]
  Tax_Var1 <- tax[rownames(tax) %in% cor4$Var1,]
  Tax_Var1 <- Tax_Var1[ , c("OTU", "Display")]
  cor5 <- merge(cor5, Tax_Var1, by.x=2, by.y=1)

  ###################
  # Defining Output #
  ###################

  #cor4$Rel_abund <- abund_n
  Result_list_total <- data.frame(
    Cor1 = cor5$Var1,
    Tax1 = cor5$Display.y,
    Cor2 = cor5$Var2,
    Tax2 = cor5$Display.x,
    Cor_val = cor5$Cor_val)

  if(return.output == "total"){
    outlist <- list(network = p, Cor_result = Result_list_total)
    return(outlist)
  }

  if(return.output == "top"){
    Result_list_top <- arrange(Result_list_total, desc(abs(Cor_val))) %>% head(n=show.top)
    outlist <- list(network = p, Cor_result = Result_list_top)
    return(outlist)
  }

if(return.output == "plot"){
  return(p)
}


  }
