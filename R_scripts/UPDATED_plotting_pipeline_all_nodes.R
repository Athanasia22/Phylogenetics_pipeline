# ---------------------------------------------------------
# 1. SETUP
# ---------------------------------------------------------
rm(list = ls()) 
library(ape); library(ggtree); library(ggplot2); library(dplyr); library(readr); library(stringr)


# --- USER TOGGLE ---
SHOW_NUMBERS <- FALSE  # TRUE = text | FALSE = Shape Dots (Black/Gray/White)
GENE_NAME    <- "COX1"

# ---------------------------------------------------------
# 2. LOAD FILES
# ---------------------------------------------------------
print("Select your .treefile...")
tree <- read.tree(file.choose())

######################DON'T FORGET TO CHANGE THE OUTGROUP#####################
rooted_tree <- root(tree, outgroup = "Hap_8_n1", resolve.root = TRUE) #as shown in the fasta/treefile

print("Select your Python-generated Metadata TSV...")
metadata <- read_tsv(file.choose(), show_col_types = FALSE)

# ---------------------------------------------------------
# 3. SUPPORT LOGIC (Data Extraction)
# ---------------------------------------------------------
# Clean empty labels and split SH-aLRT/UFboot
rooted_tree$node.label[rooted_tree$node.label == ""] <- NA
supp <- str_split(rooted_tree$node.label, "/", simplify = TRUE)

node_data <- data.frame(
  node = (Ntip(rooted_tree) + 1):(Ntip(rooted_tree) + rooted_tree$Nnode),
  label_text = rooted_tree$node.label,
  SH = as.numeric(supp[,1]),
  UF = as.numeric(supp[,2])
) %>%
  mutate(Support_Type = case_when(
    SH >= 75 & UF >= 75 ~ "Congruent (SH/UF >= 75)",
    UF >= 75            ~ "High UFboot only",
    SH >= 75            ~ "High SH-aLRT only",
    TRUE                ~ "Low"
  ))



# ---------------------------------------------------------
# 4. LABEL BUILDING (Standard R Plotmath)
# ---------------------------------------------------------
metadata <- metadata %>%
  mutate(new_label = paste0('bold("', acc_versioned, '") ~ italic("', species_name, '") ~ "|" ~ "', ncbi_country, ' [n=', Total_n, ']"'))

# Re-match to tree tips
current_tips <- rooted_tree$tip.label
clean_tips <- gsub("_n[0-9]+", "", current_tips) 
rooted_tree$tip.label <- metadata$new_label[match(clean_tips, metadata$Final_Hap_ID)]

# ---------------------------------------------------------
# 5. THE PLOT (Standard Fonts + Red Greece)
# ---------------------------------------------------------
p <- ggtree(rooted_tree) %<+% node_data +
  # Draw the tip labels
  geom_tiplab(aes(label = label, 
                  # This line turns Greece Red and others Black
                  color = grepl("Greece", label)), 
              parse = TRUE, 
              size = 3.2, 
              offset = 0.001) +
  # Define the colors: TRUE (Greece) = Red, FALSE (Others) = Black
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none")

if(SHOW_NUMBERS) {
 # ---------------------------------------------------------
# OPTION A: Show Support Numbers (Strict Filter)
# ---------------------------------------------------------
p <- p + geom_text2(
    aes(
      # ONLY show labels for nodes where both values are >= 75
      subset = (Support_Type == "Congruent (SH/UF >= 75)"), 
      label = label_text
    ), 
    size = 2.5,       # Slightly larger for readability
    vjust = -0.7,     # Position slightly above the branch
    hjust = 1.3,      # Position slightly to the left of the node split
    fontface = "plain" # Makes the numbers pop in the PDF
  

)} else {
  # OPTION B: Show Professional Shapes
  p <- p + geom_nodepoint(data = td_filter(!isTip & Support_Type != "Low"),
                          aes(shape = Support_Type, fill = Support_Type), size = 3) +
    scale_shape_manual(values = c("Congruent (SH/UF >= 75)" = 21, 
                                  "High UFboot only" = 24, 
                                  "High SH-aLRT only" = 22), 
                       name = "Node Support") +
    scale_fill_manual(values = c("Congruent (SH/UF >= 75)" = "black", 
                                 "High UFboot only" = "gray", 
                                 "High SH-aLRT only" = "white"), 
                      name = "Node Support")
}

# ---------------------------------------------------------
# 6. EXPORT & THEME
# ---------------------------------------------------------
final_plot <- p + 
  hexpand(0.35) + 
  geom_treescale(x = 0, y = -1, offset = 0.4, fontsize = 3, linesize = 0.6, width = 0.01) +
  theme(legend.position = c(0.85, 0.15), #x,y down right is 1,0 #ITS: 0.15, 0.35, COI: 0.15, 0.2, 18S: 0.85, 0.15
        legend.title = element_text(face="bold"),
        legend.background = element_rect(color="black", fill="white"),
        legend.margin = margin(6, 6, 6, 6))
print(final_plot)

#########CHANGE WD###########################
while (!is.null(dev.list())) dev.off()

setwd("C:/Users/nasia/Lab/trees/18S")   #18S--> W=9.5, H=5, COI--> W=11, H=18, ITS--> W=13, H=10
# Save with standard proportions
ggsave(paste0(GENE_NAME, "_Phylogeny_Final.jpg"), 
       plot = final_plot, width = 9.5, height = 5, units = "in", dpi = 600)


# --- SAVE AS EPS (Vector for Scientific Reports) ---
ggsave(paste0(GENE_NAME, "_Phylogeny_Final.eps"), 
       plot = final_plot, 
       device = cairo_ps,        # Use cairo_ps for better font handling
       width = 9.5, height = 5, units = "in")

# --- SAVE AS TIFF (600 DPI High-Res) ---
ggsave(paste0(GENE_NAME, "_Phylogeny_Final.tiff"), 
       plot = final_plot, 
       width = 9.5, height = 5, units = "in",
       #limitsize = FALSE, 
       dpi = 600,
       compression = "lzw")    