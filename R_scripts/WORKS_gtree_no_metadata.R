library(dplyr)
library(ape)
library(ggtree)
library(treeio)
library(ggplot2)
library(stringr)

# 1. Load and Process Tree Support
tree <- read.tree(file.choose())
tree$node.label[tree$node.label == ""] <- NA
support <- str_split(tree$node.label, "/", simplify = TRUE)

# Extract SH-aLRT and UFBoot from labels
SHaLRT <- as.numeric(support[,1])
UFBoot <- as.numeric(support[,2])

# Create a data frame for plotting support values
df <- data.frame(
  node = (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode),
  SHaLRT = SHaLRT,
  UFBoot = UFBoot
)

# 2. Rooting
# Note: Ensure "outgroup_name" exactly matches a label in your tree
#outgroup_name <- "EU582663.1/38-665"
#rooted_tree <- root(tree, outgroup = outgroup_name, resolve.root = TRUE)

# 3. Plotting IQ-TREE Results
#p <- ggtree(rooted_tree) %<+% df
p <- ggtree(tree) %<+% df


final_p <- p +
  geom_text2(aes(
    subset = !isTip & UFBoot >= 70,
    label = paste0(SHaLRT, "/", UFBoot)
  ), hjust = 1.1, vjust = -1, size = 2) +
  geom_tiplab(size = 3) + # Uses original tree labels
  hexpand(0.2) +
  geom_treescale(x = 0, y = 0, width = 0.01, offset = 0.02, fontsize = 2)
final_p
# Save Plot
setwd("C:/Users/nasia/Lab/trees")
ggsave("COX1_BEAST.jpg", plot = final_p, width = 6, height = 5)

############################### Mr Bayes Section ###############################

# Simple load and plot for MrBayes consensus trees
mb_file <- file.choose()
mb <- read.mrbayes(mb_file)

p_mb <- ggtree(mb) +
  geom_text2(aes(label = round(as.numeric(prob), 2)), hjust = -0.3, size = 2) + 
  geom_tiplab(align = TRUE, size = 2.5) + 
  hexpand(0.2)

ggsave("partition_merge.pdf", plot = p_mb, width = 10, height = 20)