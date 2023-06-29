# =========================================
# =========================================
#################################   
#plot the network of clusters - use posterior_probability value to plot IBD network
# Have look at Fig4-IBD code for details how to generate input files 
##################################################



# calculate the posterior probability of IBD sharing
samples238_final_posterior_probability <- getIBDposterior(ped.genotypes=samples238_final_genotypes,
                                                          parameters = samples238_final_parameters, 
                                                          number.cores =4, 
                                                          error = 0.001)


# Calculate IBD metrics

IBD_metrics <- list()
IBD_metrics[["fraction_IBD"]] <- (samples238_final_posterior_probability[,5:ncol(samples238_final_posterior_probability)]*2) %>% 
  colMeans() %>% as.data.frame() %>% 
  tibble::rownames_to_column() %>% 
  tidyr::extract(rowname, c("iid1", "iid2"),
                 regex="[A-Za-z0-9_()//.-]+/([A-Za-z0-9_()//.-]+)/[A-Za-z0-9_()-]+/([A-Za-z0-9_()-]+)") %>%
  mutate(iida=pmin(iid1, iid2), iidb=pmax(iid1, iid2)) %>%
  dplyr::select(-iid1, -iid2)
colnames(IBD_metrics[["fraction_IBD"]]) <- c("fraction_IBD", "iida", "iidb")
IBD_metrics[["fraction_IBD"]]$fraction_IBD <- IBD_metrics[["fraction_IBD"]]$fraction_IBD*100

#write.csv(IBD_metrics, "IBD_metrics_final.csv")

# ======================== metadata Network analysis-grouped provinces into bigger region to minimize sample size variation =============================

metadata <- samples238_final_genotypes[["pedigree"]][, c("iid", "fid")] %>% 
  tidyr::extract(fid, c("Region", "Pop"), convert=TRUE, remove=FALSE,
                 regex="([A-Za-z-_]+)([0-9]+)") %>%
  mutate(Region=gsub("-|_", "", Region) %>% stringr::str_to_title(),
         Region=ifelse(Region=="Regionone", "Luapula-Northern",Region),
         Region=ifelse(Region=="Regiontwo", "Eastern-Muchinga",Region),
         Region=ifelse(Region=="Regionthree", "Copperbelt-NorthWestern",Region),
         Region=ifelse(Region=="Regionfour", "Western",Region),
         Pop=ifelse(Pop=="101", "Copperbelt",Pop),
         Pop=ifelse(Pop=="102", "Eastern",Pop),
         Pop=ifelse(Pop=="103", "Luapula",Pop),
         Pop=ifelse(Pop=="104", "Muchinga",Pop),
         Pop=ifelse(Pop=="105", "North-Western",Pop),
         Pop=ifelse(Pop=="106", "Northern",Pop),
         Pop=ifelse(Pop=="107", "Western",Pop),
         Province=ifelse(Region=="Luapula-Northern","Regionone", "Others")) %>%
  dplyr::rename(Sample=iid, Study=fid) 

rownames(metadata) <- metadata$Sample

population_levels <- c("Study", "Region", "Pop", "Province")

pop_colours <- pop_shapes <- list()

pop_colours[["Pop"]] <- data.frame(row.names=c("Copperbelt", "Eastern", "Luapula", "Muchinga", "North-Western", "Northern","Western"),
                                   Color=c("red","orange", "tomato4","mediumvioletred","darkgreen", "darkblue", "dodgerblue"))

pop_colours[["Study"]] <- read.delim("zambia_province_colours_yearformated_sorted.txt", col.names=c("Pop", "Color"), header=FALSE, stringsAsFactors=FALSE) %>% tibble::column_to_rownames(var="Pop")


pop_shapes[["Region"]] <- data.frame(row.names=c("Luapula-Northern", "Eastern-Muchinga", "Copperbelt-NorthWestern", "Western"),
                                     Shape=c(16, 16, 16, 16))


# ======================== generate IBD networks =================

# generate layout and igraph object for pairs of isolates with IBD sharing over
# an arbitrary threshold
set_pairwise_threshold <- function(threshold, IBD_metrics, metric, metadata, population_levels, province=NULL) {
  pair_net <- list()
  if (is.null(province)) {
    ibd_sum <- IBD_metrics[[metric]]
  } else {
    keep_samples <- (metadata %>% subset(Province==province))$Sample
    ibd_sum <- IBD_metrics[[metric]] %>% 
      subset(iida %in% keep_samples & iidb %in% keep_samples)
  }
  ibd_sum$length <- ibd_sum[, metric]
  pairs_ibd <- ibd_sum %>% subset(length>=threshold) %>%
    dplyr::select(iida, iidb, length)
  
  # represent all isolates on network, irrespective of IBD sharing
  nodes <- metadata[as.character(unique(c(ibd_sum$iida, ibd_sum$iidb))), 
                    c("Sample", population_levels)] %>%
    subset(!is.na(Sample))
  
  # generate graph with weighted edges
  pair_net[["graph"]] <- igraph::graph.data.frame(pairs_ibd, vertices=nodes, directed=F)
  
  #zambia.network <- ggnetwork::ggnetwork(pair_net[["graph"]])
  #zambia.network.i.network <- zambia.network[!duplicated(zambia.network[,c(1,2,4)]),]
  #pair_net[["layout"]] <- as.matrix(zambia.network.i.network[,1:2], ncol=2)
  
  # isolates with no IBD sharing relegated to the periphery
  pair_net[["layout"]] <- layout_with_fr(pair_net[["graph"]], weights = NULL)
  colnames(pair_net[["layout"]]) <- c("x", "y")
  
  E(pair_net[["graph"]])$weight <- pairs_ibd$length
  return(pair_net)
}

# given a predefined graph/layout for pairs of isolates with IBD sharing above
# a particular threshold, generate IBD network by a particular population level
plot_pairwise_network <- function(pair_net, pop_colours, pop_shapes, title="") {
  my_nodes <- pair_net[["layout"]] %>% as.data.frame %>% 
    mutate(Sample=as_ids(V(pair_net[["graph"]])),
           Pop=as.character(vertex.attributes(pair_net[["graph"]])[["Pop"]]),
           Region=vertex.attributes(pair_net[["graph"]])[["Region"]])
  
  my_edges <- igraph::as_data_frame(pair_net[["graph"]]) %>%
    merge(my_nodes, by.x="from", by.y="Sample") %>%
    merge(my_nodes, by.x="to", by.y="Sample") %>%
    transmute(x=x.x, xend=x.y, y=y.x, yend=y.y, weight=weight)
  
  keep_Regions <- my_nodes$Region %>% unique %>% as.character %>% sort
  keep_Pops <- my_nodes$Pop %>% unique %>% as.character %>% sort
  
  my_net_test <- ggplot() +
    geom_point(data=my_nodes, aes(x=x, y=y, color=Pop, pch=Region), size=3) +
    geom_segment(data=my_edges, aes(x=x, xend=xend, y=y, yend=yend),
                 colour="black", alpha=0.5, size=0.5) +
    scale_shape_manual(values=(pop_shapes[["Region"]][keep_Regions,"Shape"])) +
    scale_color_manual(values=as.character(pop_colours[["Pop"]][keep_Pops,"Color"])) +
    ggtitle(title) +
    #scale_alpha_continuous(range = c(0, 1), limits=c(0, 100)) +
    theme_void() + theme(plot.title = element_text(face="bold", hjust=0.5))
  
  return(my_net_test)
}

# EXAMPLE USAGE

network_cutoffs <- list(ibd_1= 1, ibd_5= 5, ibd_10= 10, ibd_25=25, ibd_35=35,
                        ibd_45=45, ibd_50=50, ibd_75=75, ibd_90=90)

poster_networks <- lapply(network_cutoffs, function(x) {
  set_pairwise_threshold(x, IBD_metrics, "fraction_IBD", metadata, population_levels)})

poster_network_plots <- lapply(rev(names(network_cutoffs)), function(x) {
  plot_pairwise_network(poster_networks[[x]], pop_colours,
                        pop_shapes, title=paste0("Isolates sharing IBD% >", network_cutoffs[[x]]))})

View(poster_network_plots) 

# save the figures .pdf format