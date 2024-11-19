library(tidyverse)
library(here)
#BiocManager::install("topGO")
library(topGO)

clean_gene_annotations <- function(input_file) {
    # Read the file
    data <- read.table(input_file, sep="\t", quote="", fill=TRUE, stringsAsFactors=FALSE)
    #DOES NOT WORK: data <- read_table(input_file, col_names = FALSE)
    
    # Extract the column containing gene information (based on our format, it's column 13)
    gene_info <- data[, 12]
    
    # Function to extract information from each entry
    extract_info <- function(entry) {
        # Extract Gene ID
        gene_id <- gsub(".*ID=(.*?);.*", "\\1", entry)
        
        # Extract gene name and organism
        name_org <- gsub(".*Name=(.*?);.*", "\\1", entry)
        # Extract organism name from square brackets
        organism <- gsub(".*\\[(.*?)\\].*", "\\1", name_org)
        # Get the UniProt ID (everything before the :)
        uniprot_id <- gsub("^(.*?):.*", "\\1", organism)
        # Clean up organism name to remove any UniProt IDs
        organism <- gsub("^.*?:", "", organism)
        
        # Extract gene name (everything before the square bracket)
        gene_name <- gsub("\\s*\\[.*\\].*$", "", name_org)
        
        return(data.frame(ID = gene_id, gene_name = gene_name, organism= organism, full_name= name_org, uniprot_id = uniprot_id, original_entry = entry))
    }
    
    # Apply the extraction function to each row
    result <- map(gene_info, extract_info, .progress = TRUE) 
    
    # Convert to data frame with proper column names
    result_df <- bind_rows(result)

    # Remove duplicates (multiple SNPs that fall in the same gene)
    result_df <- distinct(result_df)

    # Replace rows where ID is not present with NA
    result_df$organism <- ifelse(grepl("ID=", result_df$organism), NA, result_df$organism)
    result_df$uniprot_id <- ifelse(grepl("ID=", result_df$uniprot_id), NA, result_df$uniprot_id)
    
    return(result_df)
}

# Read in genes 
# CHECK WHY THERE ARE PARSING FAILURES
genes <- clean_gene_annotations(here("analysis", "gea", "outputs", "gea_genes.bed"))
genes_org <- 
  genes %>% 
  drop_na(organism) %>% 
  mutate(ID = gsub("GNX-", "", ID))  %>% 
  # remove version from uniprot_id
  mutate(uniprot_id = gsub("\\..*", "", uniprot_id)) 

write_csv(genes_org, here("analysis", "gea", "outputs", "genes_list.csv"))
genes_org <- read_csv(here("analysis", "gea", "outputs", "genes_list.csv"))
unique(genes_org$Gene_Name)
unique(genes_org$Organism)

#https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html
#install.packages("gprofiler2")
library(gprofiler2)

genes_org %>% count(organism) %>% arrange(desc(n)) %>% head(20)

run_go <- function(org_key){  
    search_key <- substring(org_key, 2)

    ids <- 
        genes_org %>% 
        filter(grepl(search_key, organism)) %>%
        pull(uniprot_id)

    go_result <-
        gost(
            query = ids, 
            organism = org_key, ordered_query = FALSE, 
            multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
            measure_underrepresentation = FALSE, evcodes = TRUE, 
            user_threshold = 0.05, correction_method = "g_SCS", 
            domain_scope = "annotated", custom_bg = NULL, 
            numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE
        )

    return(go_result)
}

# Column descriptions for gprofiler output:
# 
# query: The name or identifier of the query. This indicates which query the results are associated with.
# significant: A logical value (TRUE or FALSE) indicating whether the term is significantly enriched.
# p_value: The p-value for the enrichment of the term. It indicates the statistical significance of the enrichment.
# term_size: The total number of genes associated with the GO term in the background set.
# query_size: The total number of genes in the query set.
# intersection_size: The number of genes in the query set that are also associated with the GO term.
# *precision:* The precision of the enrichment, calculated as the ratio of intersection_size to query_size. Precision measures the proportion of relevant genes (i.e., genes associated with the GO term) in the query set.
# *recall:* The recall of the enrichment, calculated as the ratio of intersection_size to term_size. Recall measures the proportion of genes associated with the GO term that are present in the reference set.
# term_id: The identifier of the GO term or other ontology term.
# source: The source of the term, such as "GO:BP" for Gene Ontology Biological Process, "GO:MF" for Gene Ontology Molecular Function, or "GO:CC" for Gene Ontology Cellular Component.
# term_name: The name of the GO term or other ontology term.
# effective_domain_size: The effective domain size, which represents the number of genes in the background set that are considered for the enrichment analysis.
# source_order: The order of the source in the results, which can be used for sorting or prioritizing results.
# parents: The parent terms of the GO term or other ontology term. This provides hierarchical context for the term.
# highlighted: A logical value (TRUE or FALSE) indicating whether the term is highlighted in the results. This can be used for visualization or emphasis in the output.
# Run GO analysis for human genes
go_hsapien <- run_go("hsapiens")
write_csv(go_hsapien$result, here("analysis", "gea", "outputs", "hsapiens.csv"))
head(go_hsapien$result, 3)
#gostplot(go_hsapien, capped = FALSE, interactive = TRUE)
p <- gostplot(go_hsapien, capped = FALSE, interactive = FALSE)
pdf(here("analysis", "gea", "plots", "hsapiens.pdf"))
print(p)
dev.off()

bp_hsapien <- 
    go_hsapien$result %>% 
    filter(source == "GO:BP") %>% 
    filter(p_value < 0.05) %>%
    arrange(p_value) %>%
    as_tibble() %>%
    dplyr::select(p_value, precision, recall, term_name)

genes_hsapien <- bp_hsapien$intersection
nrow(bp_hsapien)
nrow(bp_hsapien %>% filter(grepl("regulation", term_name)))
bp_hsapien %>% arrange(p_value) %>% head(10)
bp_hsapien %>% arrange(desc(recall)) %>% head(10)

# Run GO analysis for chickens
go_ggallus <- run_go("ggallus")
write_csv(go_ggallus$result, here("analysis", "gea", "outputs", "ggallus.csv"))
head(go_ggallus$result, 3)
#gostplot(go_hsapien, capped = FALSE, interactive = TRUE)
p <- gostplot(go_ggallus, capped = FALSE, interactive = FALSE)
pdf(here("analysis", "gea", "plots", "ggallus.pdf"))
print(p)
dev.off()
bp_ggallus <- 
    go_ggallus$result %>% 
    filter(source == "GO:BP") %>% 
    filter(p_value < 0.05) %>%
    arrange(p_value) %>%
    as_tibble() %>%
    dplyr::select(p_value, precision, recall, term_name)
nrow(bp_ggallus)
bp_ggallus %>% arrange(p_value) %>% head(10)
bp_ggallus %>% arrange(desc(recall)) %>% head(10)

go_mmusculus <- run_go("mmusculus")
write_csv(go_mmusculus$result, here("analysis", "gea", "outputs", "mmusculus.csv"))
bp_mmusculus <- 
    go_mmusculus$result %>% 
    filter(source == "GO:BP") %>% 
    filter(p_value < 0.05) %>%
    arrange(p_value) %>%
    as_tibble() %>%
    dplyr::select(p_value, precision, recall, term_name)
nrow(bp_mmusculus)
bp_mmusculus %>% arrange(p_value) %>% head(10)
bp_mmusculus %>% arrange(desc(recall)) %>% head(10)


bp_hsapien %>% arrange(p_value) %>% head(10)
bp_ggallus %>% arrange(p_value) %>% head(10)

bp_ggallus %>% arrange(desc(recall)) %>% head(10)
bp_hsapien %>% arrange(desc(recall)) %>% head(10)


# PLOTS
plot_bp <- function(bp){
  gg_df <- 
    bp %>%
    arrange(recall) %>%
    mutate(factor = factor(term_name, levels = unique(term_name))) 
  
  ggplot(gg_df, aes(y = -log(p_value), x = term_name, col = precision)) +
    geom_point() +
    ggrepel::geom_text_repel(data = head(arrange(gg_df, p_value), 10), aes(label = term_name)) +
    theme_classic() +
    labs(x = "", y = "-log(p)") +
    scale_color_viridis_c(option = "plasma", direction = -1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(20, 20, 20, 20) ) 
}


pdf(here("analysis", "gea", "plots", "hsapiens.pdf"),width = 12, height = 10)
plt1 <- plot_bp(bp_hsapien)
plt2 <- plot_bp(bp_ggallus)
plt3 <- plot_bp(bp_mmusculus)
plot_grid(plt1, plt2, plt3, nrow = 3, labels = c("Human", "Chicken", "Mouse"))
# Precision: Precision measures the proportion of relevant genes (i.e., genes associated with the GO term) in the query set
dev.off()


# Proteins of interest

# Heat shock proteins
genes_org %>% filter(grepl("HSP", full_name) | grepl("heat shock", full_name))
genes_org %>% filter(grepl(" heat", full_name))


search <- function(x, y = NULL){
  if (is.null(y)) y <- genes_org 
  y %>% filter(grepl(x, full_name)) %>% pull(full_name)
}
#Humans
#HSP70 (Heat Shock Protein 70): This protein helps protect cells from stress-induced damage, including heat stress. It plays a crucial role in protein folding and repair.
search("HSP70")
search("HSP")
search("Heat shock") # many of these!
length(search("Heat shock"))
search("heat shock")
#UCP1 (Uncoupling Protein 1): Found in brown adipose tissue, UCP1 is involved in thermogenesis, helping humans generate heat in cold environments.
search("UCP1")
search("UCP")
search("Uncoupling")
search("uncoupling")
#MC1R (Melanocortin 1 Receptor): This receptor is involved in pigmentation and can affect skin's response to UV radiation, which is a crucial adaptation in different climates.
search("MC1R")
search("Melanocortin")
search("melanocortin")
#EPAS1 (Endothelial PAS Domain Protein 1): This gene, also known as HIF-2α, is involved in the response to hypoxia and is associated with high-altitude adaptation.
search("EPAS1") 
search("PAS")# a hit!
search("Endothelial")

#Mice
#FGF21 (Fibroblast Growth Factor 21): This protein helps regulate metabolism and energy expenditure, playing a role in cold adaptation by promoting brown fat activity.
search("FGF21")
search("FGF")
search("Fibroblast growth factor 21")
search("Fibroblast growth factor")
#PRDM16 (PR Domain Containing 16): Important in the development of brown adipose tissue, which is crucial for thermoregulation in cold environments.
search("PRDM16")
search("PRDM")
#HIF1A (Hypoxia-Inducible Factor 1 Alpha): Similar to humans, this protein helps mice adapt to low oxygen levels by regulating genes involved in oxygen homeostasis.
search("HIF1A")
search("HIF")
search("Hypoxia")
# ADRB3 (Beta-3 Adrenergic Receptor): This receptor is involved in the regulation of lipolysis and thermogenesis in brown fat, aiding in cold adaptation.
search("ADRB3")
search("Adrenergic")
search("ADR")

#Chickens
#HSP70 (Heat Shock Protein 70): Like in humans, HSP70 in chickens helps protect cells from heat stress.
search("HSP70")
#UCP (Uncoupling Protein): Chickens also have uncoupling proteins that play a role in thermogenesis, helping them cope with temperature fluctuations.
search("UCP")
#TRPV4 (Transient Receptor Potential Vanilloid 4): This protein is involved in temperature sensation and regulation.
#TRP in general is associated with thermosensation
search("TRPV4")
search("TRPV")
search("Transient receptor potential")
search("TRP")
#ADRB2 (Beta-2 Adrenergic Receptor): This receptor plays a role in the regulation of metabolic processes and is involved in the adaptation to cold stress.
search("ADRB2")

#Additional 
#Aquaporins (AQPs): These proteins are involved in water transport across cell membranes and play a role in osmoregulation, which can be crucial in adapting to varying humidity levels.
search("Aquaporin")
#Leptin: This hormone, involved in regulating energy balance, can influence how organisms respond to different climatic conditions, especially in terms of energy storage and expenditure.
search("Leptin")
search("leptin")


# Look at genes associated with biological regulation in humans
br_genes <- 
  bp_hsapien %>% 
  filter(grepl("biological regulation", term_name))%>%
  pull(intersection) %>%
  # split by comma
  strsplit(",") %>%
  unlist()
br_genes <- genes_org %>% filter(uniprot_id %in% br_genes)
br_genes$gene_name
search("Heat shock", bp_genes)
search("PAS", bp_genes)

# Look at p-values
rda_genes <-
  # file from intersect_genes.R
  read_csv(here("analysis", "gea", "outputs", "gene_snp.csv")) %>%
  dplyr::rename(original_entry = full_name)  %>%
  left_join(genes_org) %>%
  # remove genes without IDs
  drop_na(gene_name)# %>%
  # filter to only include humans
 # filter(grepl("Homo sapien", organism)) 

# Get top 5 genes based on p-value
arrange(rda_genes, p.adj) %>% head(10) %>% dplyr::select(gene_name, organism, uniprot_id)

# Manhattan plots
# Prepare the data for the Manhattan plot
gg_df <- 
  read_csv(here("analysis", "gea", "outputs", "gene_snp.csv"))  %>%
  mutate(logp = -log10(p.adj)) %>%
  mutate(scaffold = case_when(grepl("Scaffold", scaffold) ~ "Anonymous scaffolds", TRUE ~ scaffold)) %>%
  mutate(scaffold = factor(scaffold, levels = c(paste0("chr", 1:11), "Anonymous scaffolds"))) %>%
  arrange(scaffold, start) %>%
  group_by(scaffold) %>%
  mutate(chr_len = max(start)) %>%
  ungroup() %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  mutate(BPcum = start + tot) %>%
  filter(grepl("chr", scaffold))

# Get top 5 genes based on p-value
top <- 
  head(arrange(rda_genes, p.adj), 5) %>%
  mutate(logp = -log10(p.adj)) %>%
  mutate(scaffold = case_when(grepl("Scaffold", scaffold) ~ "Anonymous scaffolds", TRUE ~ scaffold)) %>%
  mutate(scaffold = factor(scaffold, levels = c(paste0("chr", 1:11), "Anonymous scaffolds"))) %>%
  left_join(dplyr::select(gg_df, scaffold, start, BPcum), by = c("scaffold", "start")) 
print(dplyr::select(top, gene_name, organism))

# Axis labels
axisdf <- 
  gg_df %>%
  group_by(scaffold) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

png(here("analysis", "gea", "plots", "manhattan.png"), width = 25, height = 3, units = "in", res = 300)
# Manhattan plot
ggplot(gg_df, aes(x = BPcum, y = logp, color = as.factor(scaffold))) +
  geom_point(alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_x_continuous(label = axisdf$scaffold, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, NA)) +
  # add labels for top 5 genes
  geom_text(data = top, aes(x = BPcum, y = logp, label = gene_name), vjust = -0.5, hjust = 0.5, size = 3, col = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Chromosome", y = "-log10(p-value)")
dev.off()
