#' Correlation Network
#'
#' This function allows you to build a network based on correlations
#' @param data (required) A phyloseq object.
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (default: "Phylum")
#' @param tax.add Additional taxonomic levels to display for each entry, e.g. "Phylum" (default: "none").
#' @param tax.class Converts a specific phyla to class level instead, e.g. "p__Proteobacteria" (default: "none").
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: "best").
#' @param threshold.cor Threshold value for the correlation (default: 0.8).
#' @param threshold.count Threshold value for abundance reads (default: 10). 
#' @keywords network correlation
#' @import network
#' @import ggnetwork
#' @import sna
#' @export

cor_network <- function(data, tax.aggregate="Phylum", tax.add=NULL, tax.class=NULL, tax.empty="best", threshold.cor= 0.8, threshold.count=10){
  
  # Extracting data from phyloseq------------------------------------------------------
  
  data0 <- list(abund = as.data.frame(otu_table(data)@.Data),
                tax = data.frame(tax_table(data)@.Data, OTU =   rownames(tax_table(data))),
                sample = suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))
  
  # Clean and rename taxonomy ---------------------------------------------------------
  
  data1 <- amp_rename(data = data0,
                      tax.class=tax.class, 
                      tax.empty=tax.empty, 
                      tax.level = tax.aggregate)
  
  # Divide data to seperate data frames -----------------------------------------------
  
  abund <- data1[["abund"]]
  tax <- data1[["tax"]]
  sample <- data1[["sample"]]
  
  # Subsetting by threshold.count-----------------------------------------------------  
  
  abund[abund <= threshold.count] <- NA
  
  # Producing correlation matrix-------------------------------------------------------
  
  cormat <- cor(t(abund), method = "pearson")
  
  # Revoming self-correlation---------------------------------------------------------  
  
  diag(cormat) <- NA  
  
  # Removing "doubles"---------------------------------------------------------------
  
  cormat[lower.tri(cormat)] <- NA
  
  # Converting to dataframe for easier subsetting-------------------------------------
  
  cor2 <- melt(cormat, value.name="Rval")
  cor3 <- cor2[!is.na(cor2$Rval), ]
  
  # Subsetting by threshold.cor------------------------------------------------------
  
  cor4 <- subset(cor3, abs(cor3$Rval) >= threshold.cor)
  
  # Building network -----------------------------------------------------------------
  
  net <- network(cor4[ ,1:2 ], directed=TRUE)
  
  #Assigning Rval to edges
  
  set.edge.attribute(net, "Rval", cor4$Rval)
  
  #Assigning abundance values to nodes------------------------------------------------
  
  values <- cor4$Var1
  abund_cor <- abund[values, ] 
  abund_mean <- rowMeans(abund_cor) 
  abund_mean2 <- abund_mean %>% as.data.frame()
  
  net %v% "Avg_Abundance" = abund_mean
  
  #Visualizing using ggnetwork--------------------------------------------------------
  ggplot(net, arrow.gap = 0, aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_edges(alpha = 0.6, size=2, aes(color = Rval)) + 
    geom_nodes(aes(size= Avg_Abundance)) +
    geom_nodetext(aes(label=vertex.names), nudge_x = 0.05) +
    theme_blank()
}