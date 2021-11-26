library(tidyverse)
library(ggrepel)
library(patchwork)
library(phyloseq)
library(ggtree)
library(plotly)


theme_set(theme_minimal() + theme(
  legend.background = element_rect(fill=scales::alpha("white",0.3), color="grey"),
  legend.justification = c(1,1)))



# Shared functions ---------------------------------------------------------------

format_taxon_name <- function(taxon_name) {
  taxon_name %>% 
    str_remove("[a-z]__") %>% 
    str_replace("_L", " (Lachnospiraceae)") %>% 
    str_replace("_R", " (Ruminococcaceae)") #this stuff is just for histogram_ and density_data
}
format_biomarker_label <- function(biomarker_label) {
  case_when(grepl("KLRG1", biomarker_label)  ~ "Naive vs. KLRG1 high CD8+",
            grepl("CD57", biomarker_label) ~ "CD57+ vs. - NK",
            grepl("CRP", biomarker_label) ~ "C-reactive protein", 
            grepl("DNA", biomarker_label) ~ "Hallmark DNA Repair",
            tolower(biomarker_label) == "hba1c" ~ "HbA1c",
            grepl("ffector CD8", biomarker_label) ~ "Naive vs. Effector CD8+",
            grepl("Th1", biomarker_label) ~ "Th1 vs. Naive CD4+",
            grepl("Transcriptomic", biomarker_label) ~ "Transcriptomic Age",
            grepl("Treg", biomarker_label) ~ "Treg vs. Naive CD4+",
            TRUE ~ biomarker_label)
}
make_biomarker_groups <- function(biomarker_label) {
  case_when(biomarker_label %in% DBS_BIOMARKERS ~ "DBS", 
            grepl("Hallmark", biomarker_label) ~ "Hallmark", 
            TRUE ~ "Immunological")
}
filter_biomarker_groups <- function(data, group) {
  if(group == "All") data else filter(data, biomarker_group==group)
}
drop_taxa_labels_except <- function(data, taxa_names) {
  if("All" %in% taxa_names) {
    return(data)
  } else {
    data %>%
      mutate(#taxon_is_labelled = taxon_is_labelled & taxon_label %in% taxa_names,
             taxon_label = ifelse(taxon_label %in% taxa_names, taxon_label, NA))
  }
}
is_taxon_labelled <- function(log_fold_change, neg_log_sval, lfc_cutoff, sval_cutoff) {
  abs(log_fold_change) > lfc_cutoff & neg_log_sval > sval_cutoff
}
get_sorted_list_of_labelled_taxa <- function(data, sort_by="appearance") {
  taxa <- data %>% 
    count(taxon_label) %>% 
    drop_na() %>% 
    arrange(-n) %>% 
    pull(taxon_label)
  
  if(sort_by == "alphabetic") sort(taxa) else taxa
}

collapse_biomarker_labels <- function(data, collapse=TRUE, other_label = "Others") {
  if(collapse) {
    data %>% 
      group_by(biomarker_label) %>%
      mutate(any_taxa_labeled = sum(taxon_is_labelled) > 0) %>%
      ungroup %>%
      mutate(biomarker_label = ifelse(any_taxa_labeled, biomarker_label, other_label))
  } else {
    data
  }
}
label_by_lfc_sval <- function(data, lfc_cutoff, sval_cutoff) {
  data %>% mutate(taxon_is_labelled = is_taxon_labelled(log_fold_change = logfc_HDL_reversed, 
                                                        neg_log_sval = neg_log_sval,
                                                        lfc_cutoff = lfc_cutoff,
                                                        sval_cutoff = sval_cutoff),
                  taxon_label = ifelse(taxon_is_labelled, format_taxon_name(Taxon), NA))
}
make_key_variables <- function(raw_data) {
  out_data <- raw_data %>% 
    mutate(neg_log_sval = -log10(`S-value`),
           biomarker_label = format_biomarker_label(Biomarker), 
           logfc_HDL_reversed = ifelse(biomarker_label=="HDL", -lfc, lfc),
           biomarker_source = ifelse(biomarker_label %in% DBS_BIOMARKERS, "dbs", "gene"))
  
  out_data
}

add_formatted_ci <- function(data, dig=2) {
  fmt <- function(x) formatC(x, digits=dig, format="f")
  data %>% mutate(`90% HPDI` = paste0(fmt(`Lower 90% HPDI`), ", ", fmt(`Upper 90% HPDI`)))
}

# Data and constants ------------------------------------------------------

table_s1 <- readxl::read_excel(path = "eTable1.xlsx", sheet = 1, skip = 1) %>%
  filter(`taxon level` %in% c("Family", "Genus","Species","Phylum")) #I don't use any OTU results here

raw_manhattan_data <- read_csv("manhattan_data.csv")
tree <- dget(file = "tree.txt")
taxonomy <- read_csv("taxonomy.csv") %>% rename(id = X1)

glommed_trees <- dget(file = "glommed trees.txt")
propionate <- readxl::read_excel("Reichardt 2014 propionate.xlsx") %>%
  separate(name, into = c("Genus", "Species"), sep = " " ) %>%
  mutate(propio_species = 1)
taxonomy_more_info <- readxl::read_excel("taxonomy_more_info.xlsx") %>% 
  mutate(across(everything(), format_taxon_name)) %>%
  group_by(Family) %>%
  mutate(Butyrate = case_when(any(Butyrate == "Yes") ~ ifelse(rank == "Family", "Yes", Butyrate), 
                           Butyrate == "NA" ~ "No",
                           TRUE ~ Butyrate),
         across(everything(), ~ifelse(.=="NA", NA, .))) %>%
  left_join(propionate %>% mutate(Species = NA, Genus = NA) %>% rename(propio_family = propio_species)) %>%
  left_join(propionate %>% mutate(Species = NA) %>% rename(propio_genus = propio_species)) %>%
  left_join(propionate) %>%
  mutate(across(starts_with("propio"), ~ifelse(is.na(.), 0, .)),
         Propionate = ifelse(rowSums(across(starts_with("propio"))) > 0, "Yes", "No")) %>%
  select(-propio_family:-propio_species)



glommed_taxonomies <- read_csv("full_taxonomy.csv") %>%
  mutate(across(everything(), format_taxon_name)) %>%
  left_join(taxonomy_more_info %>% select(tax_id, Butyrate, Propionate, `More info`))


full_raw_data <- table_s1[!duplicated(table_s1$`log2 Fold Change`), ]  %>%
  rename(lfc=`log2 Fold Change`) %>%
  inner_join(raw_manhattan_data) %>%
  mutate(Taxon = str_remove(Taxon, "[a-z]__"))

histogram_data <- read_csv("histograms.csv")[,-1] %>% mutate(tax_id = format_taxon_name(tax_id))
histogram_nonzero_data <- read_csv("histograms_nonzero.csv")[,-1] %>% mutate(tax_id = format_taxon_name(tax_id))
density_data <- read_csv("densities.csv")[,-1] %>% mutate(tax_id = format_taxon_name(tax_id))
density_nonzero_data <- read_csv("densities_nonzero.csv")[,-1] %>% mutate(tax_id = format_taxon_name(tax_id))



DBS_BIOMARKERS <- c("HDL","LDL","Glucose","HbA1c", "C-reactive protein")
COLOR_PALETTE <- RColorBrewer::brewer.pal(12, "Paired")[-11]
TAXON_CHOICES <- table_s1 %>% 
  filter(abs(`log2 Fold Change`) > 0.8 & -log10(`S-value`) > 1.5) %>%
  mutate(Taxon = format_taxon_name(Taxon)) %>%
  group_by(Taxon) %>%
  count() %>% 
  arrange(-n) %>%
  pull(Taxon)
TREE_SQUISHER <- "
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBB#
"  
# Volcano plot stuff ------------------------------------------------------



make_volcano_plot <- function(data, 
                              biomrkr_source = "dbs", 
                              legend_position = c(0.6, 1)){
  
  data %>%
    filter(biomarker_source == biomrkr_source) %>%
    ggplot(aes(x=logfc_HDL_reversed, y=neg_log_sval, color=biomarker_label)) +
    geom_point(aes(size=`Taxon prevalence`, alpha=taxon_is_labelled, fill=biomarker_label), shape=21) +
    geom_label_repel(aes(label=taxon_label), alpha=0.8, force=10, 
                     show.legend = FALSE, min.segment.length = 0)  +
    scale_fill_manual(values=COLOR_PALETTE, name="Biomarker") +
    scale_color_manual(values=COLOR_PALETTE, guide="none") +
    scale_size_continuous(name="Taxon prevalence", breaks=c(0,0.1,0.3,0.6), guide="none") +
    scale_alpha_manual(values=c(1/4, 1), guide="none") + 
    theme(strip.text = element_blank(), 
          legend.position = legend_position) +
    labs(x="log2-Fold Change", y="-log10(s-value)")+
    guides(fill = guide_legend(override.aes = list(size=4)))
  
}



# manhattan plot stuff ----------------------------------------------------


 
#functions
get_rare_phyla <- function(taxonomy, minimum_n = 5) {
  taxonomy %>% 
    filter(id %in% raw_manhattan_data$id) %>%
    count(Phylum) %>% 
    filter(n < minimum_n) %>% 
    pull(Phylum) %>%
    format_taxon_name() %>% 
    c("SR1","TM7")
}
collapse_rare_phyla <- function(phyla, taxonomy, minimum_n = 5) {
  phyla %>% 
    as.character() %>% 
    ifelse(. %in% get_rare_phyla(taxonomy, minimum_n), "", .) 
}
format_and_collapse_tree_phyla <- function(tree, taxonomy) {
  attr(tree, "Phylum") <- attr(tree, "Phylum") %>%
    format_taxon_name() %>% 
    collapse_rare_phyla(taxonomy) 
  tree
}
add_phylum_to_manh_data <- function(manhattan_data, taxonomy_data) {
  taxonomy_data %>% 
    mutate(Phylum = format_taxon_name(Phylum)) %>%
    select(taxon_label = tax_id, Phylum) %>%
    right_join(manhattan_data)
}
make_tree_plot <- function(tree) {
  ggtree(tree, aes(color=Phylum), ladderize = FALSE) +
    theme_tree()
}
add_labels_to_tree_plot <- function(tree_plot, 
                                    tree_distance_from_xaxis = 1/2,
                                    label_distance = 200) {
  tree_plot_label_data <- tree_plot$data %>% 
    group_by(Phylum) %>%
    summarise(y_mean = mean(y)-label_distance)  %>%
    mutate(y_mean = y_mean - (Phylum == "Actinobacteria")*200  + (Phylum == "Fusobacteria")*100)
  tree_plot + geom_text(data=tree_plot_label_data, 
                        mapping=aes(x=-tree_distance_from_xaxis, y=y_mean, color=Phylum, label=Phylum), 
                        inherit.aes = FALSE,
                        show.legend = FALSE,
                        angle=15, vjust=-1) +
    theme(legend.position = "none") + 
    coord_flip()
}
order_manhattan_data <- function(data) {
  data %>% 
    arrange(order) %>%
    select(-order) %>%
    group_by(id) %>%
    nest() %>%
    ungroup() %>%
    mutate(order = row_number()) %>%
    unnest(everything()) %>%
    ungroup() 
}
make_manhattan_plot <- function(data) {
  
  data %>% 
    ggplot(aes(x=order, y=neg_log_sval, fill=biomarker_label, group=biomarker_label)) + 
    geom_point(aes(alpha=taxon_is_labelled), shape=21, color="grey", size=4) + 
    geom_label_repel(aes(label=taxon_label, color=biomarker_label), fill="white", show.legend = FALSE) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.position = c(0,1), legend.justification = c(0,1),
          panel.grid.major.x = element_blank()) +
    labs(x=NULL, y="-log10(S-value)") +
    scale_fill_manual(values = COLOR_PALETTE) +
    scale_color_manual(values = COLOR_PALETTE) +
    scale_alpha_manual(values=c(0.1, 1), guide="none") + 
    guides(fill=guide_legend(title="Biomarker")) 
}



# Density / histogram stuff -----------------------------------------------

make_density_plot <- function(density_data, taxon) {
  density_data %>% 
    filter(tax_id == taxon) %>%
    mutate(y = y/max(y)) %>%
    ggplot(aes(x=x)) + 
    geom_ribbon(aes(ymin=y, ymax=1), fill="white") + 
    geom_ribbon(aes(ymin=0, ymax=y), fill="grey70", color="black")  +
    scale_x_continuous(labels = scales::percent_format(scale = 100)) +
    theme(axis.text.y = element_blank()) +
    labs(y="Density", x=paste("Relative abundance of", taxon, "in each sample"))
}

make_histogram_plot <- function(histogram_data, taxon) {
  histogram_data %>%
    filter(tax_id == taxon) %>%
    rename(right_break = breaks) %>%
    mutate(left_break =right_break - min(right_break),
           widths = right_break - left_break) %>%
    ggplot() + 
    geom_rect(aes(xmin=left_break, xmax=right_break, ymin=0, ymax=counts),
              fill="grey60", color="grey30") +
    labs(y="Number of samples (participants)", x=paste("Count of", taxon, "in each sample"))
}



# Taxonomy information ----------------------------------------------------


get_taxonomy_table <- function(glommed_taxonomies, taxon) {
  glommed_taxonomies %>%
    filter(tax_id == taxon) %>%
    select(-tax_id) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("rowname") %>%
    mutate(rowname = paste0(rowname, ": ")) %>%
    set_names(c("",""))
}


# Pie chart ---------------------------------------------------------------


make_pie_data <- function(data, taxon) {
  data %>%
    select(taxon_label, prevalence = `Taxon prevalence`) %>%
    filter(taxon_label == taxon) %>%
    unique() %>%
    {bind_rows(list(., tibble(taxon_label = "Others", prevalence = 1 - .$prevalence)))} %>%
    mutate(taxon_label = fct_inorder(taxon_label))
}

make_pie_chart <- function(pie_data) {
  ggplot(pie_data, aes(x="", y=prevalence, fill=taxon_label)) +
    geom_bar(stat="identity", width=1, color="grey50", show.legend = FALSE) +
    coord_polar("y", start=0) +
    theme_void() +
    theme(rect = element_rect(fill = "transparent", color="transparent")) +
    labs(title = paste0("Prevalence = ", round(pie_data$prevalence[1] * 100,1), "%"),
         y = NULL, x = NULL) +
    scale_fill_manual(values = c("darkgreen", "transparent"))
}

# Tests -------------------------------------------------------------------


#actual program
filteredResultsData <- 
full_raw_data %>% 
  make_key_variables() %>%
  label_by_lfc_sval(0.8, 2) %>%
  collapse_biomarker_labels() %>%
  drop_taxa_labels_except("All")

filteredResultsData %>%
    make_volcano_plot() 

filteredResultsData %>%
  order_manhattan_data() %>%
  add_phylum_to_manh_data(glommed_taxonomies) %>%
  add_formatted_ci() %>%
  make_manhattan_plot() 

filteredResultsData %>%
  make_pie_data("Acinetobacter") %>%
  make_pie_chart()

 
