
packages <- c("ggplot2", "cowplot", "dplyr", "svglite", "ape","phytools","ggtree")
install.packages(setdiff(packages, rownames(installed.packages()))) 

library(ape)
library(phytools)
library(dplyr)
library(ggplot2)
library(ggtree)
library(cowplot)

packageVersion("ggtree")
packageVersion("ape")
packageVersion("phytools")


args <-  commandArgs(trailingOnly = TRUE)
results_dir <- args[1]
visualizations_dir <- args[2]
image_format <- args[3]
tip_label_size <- args[4]

setwd(results_dir)

## Troubleshooting ##
#setwd('/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset4/test2/mexb/internal_data') #debugging
#visualizations_dir <- '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset4/test2/mexb/visualizations' #debugging
#image_format <- 'png' # debugging
#tip_label_size <- 3 # debugging
## ## ## ## ## ## ## ##

tree_in <- ape::read.tree('representative_seed_sequences.nwk')
tree_in <- phytools::midpoint.root(tree_in)

tree_out <- ggtree(tree_in)+
  geom_tiplab(aes(label=as.character(label)),size=as.numeric(tip_label_size), align=TRUE, hjust=-0.05)+
  geom_tippoint(size=1)+
  xlim(0,2)

cowplot::save_plot(paste(visualizations_dir,'/representative_seed_phylogeny','.',image_format,sep='') ,tree_out,base_height=8,base_aspect_ratio = 2)

