#Load the packages/libraries required 
library("qiime2R")
library("microbiome")
library("vegan")
library("ggplot2")
library("RColorBrewer")
library("gridExtra")
library("cowplot")
library("phyloseq")
library("dplyr")
library("tidyr")
library("ggsignif")
library("ape")
require("GUniFrac")
require("ade4")
library("reshape2")
library("stringr")
library("grid")

####Make Phyloseq Object with Qiime Outputs ####

#Import data
pseq <- qiime2R::qza_to_phyloseq(
  features = "table-with-phylum.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.sklearn.qza",
  metadata = "metadata.tsv"
)

#Summary of phyloseq object
microbiome::summarize_phyloseq(pseq)

#Number of reads per sample
sample_depths <- microbiome::readcount(pseq)
hist(sample_depths, main = "Histogram of read depths")

#Extract ASV table
phyloseq::otu_table(pseq)

#make vector to keep track
num_asvs_vec <- c(nrow(phyloseq::otu_table(pseq)))
names(num_asvs_vec)[1] <- "abundance"
num_asvs_vec

#### Relative Abundance ####

#Transform abundance table to a relative abundance table
pseq_relabund <- microbiome::transform(pseq, "compositional")

#Summarise and check sample counts which should each amount to 1
microbiome::summarize_phyloseq(pseq_relabund)
microbiome::readcount(pseq_relabund)


###Phylum phyloseq###
phylum_pseq <- microbiome::aggregate_taxa(pseq_relabund, "Phylum", verbose = FALSE)

head(phyloseq::otu_table(phylum_pseq))
paste0("Number of phyla: ", nrow(phyloseq::otu_table(phylum_pseq)))
#Summarise
microbiome::summarize_phyloseq(phylum_pseq)
microbiome::readcount(phylum_pseq)

# Generate a colour palette
palette_phy <- colorRampPalette(brewer.pal(12, "Paired"))(51)

# Subset by seaweed species
phylum_sub_dulse <- subset_samples(phylum_pseq, species == "Palmaria palmata")
phylum_sub_fucus <- subset_samples(phylum_pseq, species == "Fucus serratus")

#Plot P. palamaria Graph
phylum_dulse <-microbiome::plot_composition(phylum_sub_dulse, 
                                            otu.sort = NULL,
                                            group_by = "rep") +
  scale_fill_manual("Phylum", values = palette_phy) +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6")) +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, lineheight = 0.5, margin = margin(t = 0.5, b = 0.5)),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) + 
  xlab("") +
  ylab("Relative abundance")

#Plot F. serratus Graph
phylum_fucus <-microbiome::plot_composition(phylum_sub_fucus, 
                                            otu.sort = NULL,
                                            group_by = "rep") +
  scale_fill_manual("Phylum", values = palette_phy) +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6")) +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, lineheight = 0.5, margin = margin(t = 0.5, b = 0.5)),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) + 
  xlab("") +
  ylab("Relative abundance")

#Make plot for legend
phylum_fucus_legend <- plot_composition(phylum_sub_fucus, 
                                        otu.sort = NULL,
                                        group_by = "rep") +
  scale_fill_manual("Phylum", values = palette_phy) +
  guides(fill = guide_legend(title = "",
                             override.aes = list(size = 4),
                             ncol = 4)) +
  theme(legend.box = "vertical", 
        legend.key.size = unit(0.2, "cm"))

# Extract the legend
legend1 <- get_legend(phylum_fucus_legend)

# Assuming you have the grid of plots and the saved legend
grid_plot1 <- grid.arrange(phylum_dulse, phylum_fucus, ncol = 1)

#lengend underneath
combined_plot1 <- plot_grid(grid_plot1,
                            legend1,
                            ncol = 1,
                            rel_widths = c(0.6, 0.2),
                            rel_heights = c(0.8, 0.2))
combined_plot1

###Phyloseq Family###
family_pseq <- microbiome::aggregate_taxa(pseq_relabund, "Family", verbose = FALSE)
head(phyloseq::otu_table(family_pseq))
paste0("Number of family: ", nrow(phyloseq::otu_table(family_pseq)))
microbiome::summarize_phyloseq(family_pseq)
microbiome::readcount(family_pseq)

#Aggregate rare families
family_rareaggregate_pseq <- microbiome::aggregate_rare(
  pseq_relabund, level = "Family",
  detection = 0.025, prevalence = 5/100
)

# Subset by Seaweed Species
family_sub_dulse <- subset_samples(family_rareaggregate_pseq, species == "Palmaria palmata")
family_sub_fucus <- subset_samples(family_rareaggregate_pseq, species == "Fucus serratus")

set.seed(NULL)  # Optional: Ensures randomness each time
palette_fam <- sample(colorRampPalette(brewer.pal(12, "Paired"))(51))

#Plot P. palmaria Graph
family_dulse <-microbiome::plot_composition(family_sub_dulse, 
                                            otu.sort = NULL,
                                            group_by = "rep") +
  scale_fill_manual("Family", values = palette_fam) +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6")) +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, lineheight = 0.5, margin = margin(t = 0.5, b = 0.5)),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) + 
  xlab("") +
  ylab("Relative abundance")

#Plot F. serratus Graph
family_fucus <-microbiome::plot_composition(family_sub_fucus, 
                                            otu.sort = NULL,
                                            group_by = "rep") +
  scale_fill_manual("Family", values = palette_fam) +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6")) +  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, lineheight = 0.5, margin = margin(t = 0.5, b = 0.5)),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) + 
  xlab("") +
  ylab("Relative abundance")

# Plot graph with legend present
family_fucus_legend <- plot_composition(family_sub_fucus, 
                                 otu.sort = NULL,
                                 group_by = "rep") +
  scale_fill_manual("Family", values = palette_fam) +
  guides(fill = guide_legend(title = "",
                             override.aes = list(size = 4),
                             ncol = 4)) +
  theme(legend.box = "vertical",
        legend.key.size = unit(0.2, "cm"))

# Extract the legend
legend <- get_legend(family_fucus_legend)
#Plot species together
grid_plot <- grid.arrange(family_dulse, family_fucus, ncol = 1)
#Add Legend
combined_plot <- plot_grid(grid_plot,
                           legend,
                           ncol = 1,
                           rel_widths = c(0.6, 0.2),
                           rel_heights = c(0.8, 0.2))
combined_plot

#Combine the Family & Phylum Plots
taxa_together <- plot_grid(combined_plot1, combined_plot, ncol = 2)
taxa_together

png("Figure2", width = 3500, height = 3500, res = 300)
print(taxa_together)  # Print the plot to the file

dev.off() 
####Rarefaction Curves ####

#Extract ASV table as data frame
asv_abund_df <- as.data.frame(t(phyloseq::otu_table(pseq)))

# Minimum sample count to be rarefied
raremax <- min(microbiome::readcount(pseq))

# Define colors for the plot
col <- c("darkred", "forestgreen", "hotpink", "blue", "orange", "black")
set.seed(3)
grp <- factor(sample(seq_along(col), nrow(asv_abund_df), replace = TRUE))
cols <- col[grp]

# Perform rarefaction
out <- rarecurve(asv_abund_df, step=100, sample=raremax, label=FALSE)

# Extract Nmax and Smax for the plot ranges
Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
Smax <- sapply(out, max)

# Set up an empty plot
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Read Depth",
     ylab = "ASVs", type = "n")

# Add vertical line for rarefaction threshold
abline(v = raremax)

# Add rarefaction lines for each sample
for (i in seq_along(out)) {
  N <- attr(out[[i]], "Subsample")
  lines(N, out[[i]], col = cols[i])
}

#Rarefy to minimum depth
pseq_rarefy <- phyloseq::rarefy_even_depth(
  pseq, sample.size = min(microbiome::readcount(pseq)),
  rngseed = 1000
)
#Summarise and check sample counts
microbiome::summarize_phyloseq(pseq_rarefy)
microbiome::readcount(pseq_rarefy)

#Add relative abundance ASV count
num_asvs_vec["rarefied"] <- nrow(phyloseq::otu_table(pseq_rarefy))
num_asvs_vec


####Alpha Diversity####

### GENUS ###
#Produce data frame of all alpha diversity values
alpha_df <- phyloseq::estimate_richness(physeq = pseq_rarefy)
metadata <- sample_data(pseq_rarefy)
alpha_df1 <- cbind(metadata, alpha_df)

#Shannon's
hist(alpha_df1$Shannon)
shapiro.test(alpha_df1$Shannon)
kruskal.test(Shannon ~ rep, data=alpha_df1)
pairwise.wilcox.test(alpha_df1$Shannon, alpha_df1$rep, p.adjust.method = "fdr")

summary_genus <- alpha_df1 %>%
  group_by(rep) %>%
  summarise(
    mean_Shannon = mean(Shannon, na.rm = TRUE),   
    sd_Shannon = sd(Shannon, na.rm = TRUE)          
  )

####Phylum###
pseq_rarefy_phy <- microbiome::aggregate_taxa(pseq_rarefy, "Phylum", verbose = FALSE)
#Produce data frame of all alpha diversity values
alpha_df_phy <- phyloseq::estimate_richness(physeq = pseq_rarefy_phy, measures = c("Shannon", "Simpson", "Chao1", "Observed"))
head(alpha_df_phy)
metadata_phy <- sample_data(pseq_rarefy_phy)
alpha_df1_phy <- cbind(metadata_phy, alpha_df_phy)

#Shannon's
hist(alpha_df1_phy$Shannon)
shapiro.test(alpha_df1_phy$Shannon)
kruskal.test(Shannon ~ rep, data=alpha_df1_phy)
pairwise.wilcox.test(alpha_df1_phy$Shannon, alpha_df1_phy$rep, p.adjust.method = "fdr")

summary_phylum <- alpha_df1_phy %>%
  group_by(rep) %>%
  summarise(
    mean_Shannon = mean(Shannon, na.rm = TRUE),   
    sd_Shannon = sd(Shannon, na.rm = TRUE)          
  )

###Class ###
pseq_rarefy_cla <- microbiome::aggregate_taxa(pseq_rarefy, "Class", verbose = FALSE)
#Produce data frame of all alpha diversity values
alpha_df_cla <- phyloseq::estimate_richness(physeq = pseq_rarefy_cla, measures = c("Shannon", "Simpson", "Chao1", "Observed"))
head(alpha_df_cla)
metadata_cla <- sample_data(pseq_rarefy_cla)
alpha_df1_cla <- cbind(metadata_cla, alpha_df_cla)

#Shannon's
hist(alpha_df1_cla$Shannon)
shapiro.test(alpha_df1_cla$Shannon)
kruskal.test(Shannon ~ rep, data=alpha_df1_cla)
pairwise.wilcox.test(alpha_df1_cla$Shannon, alpha_df1_cla$rep, p.adjust.method = "fdr")

summary_class <- alpha_df1_cla %>%
  group_by(rep) %>%
  summarise(
    mean_Shannon = mean(Shannon, na.rm = TRUE),   
    sd_Shannon = sd(Shannon, na.rm = TRUE)          
  )

###Order ###
pseq_rarefy_ord <- microbiome::aggregate_taxa(pseq_rarefy, "Order", verbose = FALSE)
#Produce data frame of all alpha diversity values
alpha_df_ord <- phyloseq::estimate_richness(physeq = pseq_rarefy_ord, measures = c("Shannon", "Simpson", "Chao1", "Observed"))
head(alpha_df_ord)
metadata_ord <- sample_data(pseq_rarefy_ord)
alpha_df1_ord <- cbind(metadata_ord, alpha_df_ord)

#Shannon's
hist(alpha_df1_ord$Shannon)
shapiro.test(alpha_df1_ord$Shannon)
kruskal.test(Shannon ~ rep, data=alpha_df1_ord)
pairwise.wilcox.test(alpha_df1_ord$Shannon, alpha_df1_ord$rep, p.adjust.method = "fdr")

summary_order <- alpha_df1_ord %>%
  group_by(rep) %>%
  summarise(
    mean_Shannon = mean(Shannon, na.rm = TRUE),   
    sd_Shannon = sd(Shannon, na.rm = TRUE)          
  )

###Family ###
pseq_rarefy_fam <- microbiome::aggregate_taxa(pseq_rarefy, "Family", verbose = FALSE)
#Produce data frame of all alpha diversity values
alpha_df_fam <- phyloseq::estimate_richness(physeq = pseq_rarefy_fam, measures = c("Shannon", "Simpson", "Chao1", "Observed"))
head(alpha_df_fam)
metadata_fam <- sample_data(pseq_rarefy_fam)
alpha_df1_fam <- cbind(metadata_fam, alpha_df_fam)

#Shannon's
hist(alpha_df1_fam$Shannon)
shapiro.test(alpha_df1_fam$Shannon)
kruskal.test(Shannon ~ rep, data=alpha_df1_fam)
pairwise.wilcox.test(alpha_df1_fam$Shannon, alpha_df1_fam$rep, p.adjust.method = "fdr")

summary_family <- alpha_df1_fam %>%
  group_by(rep) %>%
  summarise(
    mean_Shannon = mean(Shannon, na.rm = TRUE),   
    sd_Shannon = sd(Shannon, na.rm = TRUE)          
  )

### Create Plot ###
#Merge Shannon Scores into 1 dataframe
# Convert row names to a column in each dataframe
alpha_df1$SampleID <- row.names(alpha_df1)
alpha_df1_cla$SampleID <- row.names(alpha_df1_cla)
alpha_df1_fam$SampleID <- row.names(alpha_df1_fam)
alpha_df1_ord$SampleID <- row.names(alpha_df1_ord)
alpha_df1_phy$SampleID <- row.names(alpha_df1_phy)

combined_df <- alpha_df1[, c("SampleID", "rep", "Shannon")]
colnames(combined_df)[3] <- "Genus"

combined_df <- merge(combined_df, alpha_df1_cla[, c("SampleID", "Shannon")], by = "SampleID", all = TRUE)
colnames(combined_df)[ncol(combined_df)] <- "Class"

combined_df <- merge(combined_df, alpha_df1_fam[, c("SampleID", "Shannon")], by = "SampleID", all = TRUE)
colnames(combined_df)[ncol(combined_df)] <- "Family"

combined_df <- merge(combined_df, alpha_df1_ord[, c("SampleID", "Shannon")], by = "SampleID", all = TRUE)
colnames(combined_df)[ncol(combined_df)] <- "Order"

combined_df <- merge(combined_df, alpha_df1_phy[, c("SampleID", "Shannon")], by = "SampleID", all = TRUE)
colnames(combined_df)[ncol(combined_df)] <- "Phylum"

# Convert the data to long format
shannon_long <- combined_df %>%
  pivot_longer(cols = c(Genus, Family, Order, Class, Phylum),
               names_to = "Metric",
               values_to = "Value")

# Reorder the Metric factor levels
shannon_long$Metric <- factor(shannon_long$Metric, 
                              levels = c("Phylum", "Class", "Order", "Family", "Genus"))
# Reorder the rep factor levels
shannon_long$rep <- factor(shannon_long$rep, 
                           levels = c("Fucus: Uncultured", "Fucus: SSM", "Fucus: MB", "Palmaria: Uncultured", "Palmaria: SSM", "Palmaria: MB"))

# Define the list of comparisons
comparisons <- list(
  c("Fucus: Uncultured", "Palmaria: Uncultured"),
  c("Fucus: Uncultured", "Fucus: SSM"),
  c("Fucus: Uncultured", "Fucus: MB"),
  c("Fucus: MB", "Fucus: SSM"),
  c("Palmaria: Uncultured", "Palmaria: SSM"),
  c("Palmaria: MB", "Palmaria: Uncultured"),
  c("Palmaria: MB", "Palmaria: SSM")
)


fig1 <-ggplot(shannon_long, aes(x = rep, y = Value, fill = rep)) +
  geom_boxplot() +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +  # Set nrow to 1
  theme_bw() +
  scale_fill_manual(values = c("Fucus: Uncultured" = "#669966", 
                               "Fucus: SSM" = "#DC9596",
                               "Fucus: MB" = "#999999",
                               "Palmaria: SSM" = "#660033",
                               "Palmaria: MB" = "#666666",
                               "Palmaria: Uncultured" = "#003300")) +
  labs(x = "Treatments", y = "Shannon Diversity") +
  geom_signif(comparisons = comparisons, 
              annotations = c("ns", "**","**","**","**","**","*"),
              y_position = c(1.5,0.7,0.8,0.4,0.7,0.8,0.4),
              textsize = 4, data = shannon_long %>% filter(Metric == "Phylum"), tip_length = 0) +
  geom_signif(comparisons = comparisons, 
              annotations = c("ns", "**","**","**","**","**","*"),
              y_position = c(1.8,0.9,1.1,0.4,0.9,1.1,0.4),
              textsize = 4, data = shannon_long %>% filter(Metric == "Class"), tip_length = 0) +
  geom_signif(comparisons = comparisons, 
              annotations = c("ns", "**","**","**","**","**","**"),
              y_position = c(2.8,1.5,1.8,0.5,1.5,1.8,0.5),
              textsize = 4, data = shannon_long %>% filter(Metric == "Order"), tip_length = 0) +
  geom_signif(comparisons = comparisons, 
              annotations = c("ns", "**","**","ns","**","**","ns"),
              y_position = c(3.1,2,2.2,1.65,2,2.2,1.65),
              textsize = 4, data = shannon_long %>% filter(Metric == "Family"), tip_length = 0) +
  geom_signif(comparisons = comparisons, 
              annotations = c("ns", "**","**","ns","**","**","ns"),
              y_position = c(5.5,3.3,3.7,2.9,3.3,3.7,2.9),
              textsize = 4, data = shannon_long %>% filter(Metric == "Genus"), tip_length = 0) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "lightgrey"),
    legend.title = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
  )
png("Figure1.png", width = 4000, height = 2000, res = 300)
print(fig1)  # Print the plot to the file
dev.off() 

####BETA Diversity ####

#Prune tree after rarefication
ps_tree = phy_tree(pseq_rarefy)
sprintf("Is tree binary: %s", is.binary(ps_tree))
phy_tree(pseq_rarefy) = multi2di(ps_tree)
sprintf("Is tree binary: %s", is.binary(phy_tree(pseq_rarefy)))
dist = phyloseq::distance(pseq_rarefy, method = "unifrac", weighted = F)
sprintf("`dist` class: %s", class(dist))

ord.mds.wunifrac <-phyloseq::ordinate(pseq_rarefy, method = "MDS", distance = "wunifrac")

wUF.ordu = ordinate(pseq_rarefy, method="NMDS", distance="unifrac", weighted=TRUE)

#plot shannons diversity as point size
shannon <- (alpha_df$Shannon)

NMDS_Plot <- phyloseq::plot_ordination(pseq_rarefy, ord.mds.wunifrac,
                                       color = "type", shape = "species") +
  geom_point(aes(size = alpha_df$Shannon), alpha = 0.5) +  # Mapping size to "Shannon" column
  theme_minimal() +
  scale_size_continuous(name = "Shannon's Diversity", range = c(3.5, 7.5)) +  # Adjust range for visibility
  scale_color_manual(name = "Treatment", values = c("MB" = "#DC9596", "Uncultured" = "#1F271B", "SSM" = "#52489C")) +  # Custom colors
  guides(size = FALSE) +
  theme_q2r() +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)  # Set larger size for legend titles
  ) +
  labs(x = "Axis 1 (63.7%)", y = "Axis 2 (18.9%)", shape = "Species") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_shape_discrete(labels = c(expression(italic("Fucus serratus")), 
                                  expression(italic("Palmaria palmata"))))

###PERMANOVA###

#Produce distance matrix
wunifrac_dist_mat <- phyloseq::distance(pseq_rarefy, method="wunifrac")
#check homogeneity
anova(betadisper(wunifrac_dist_mat, metadata$rep))

#Extract metadata data frame from phyloseq
metadf <- data.frame(phyloseq::sample_data(pseq_rarefy))
bray.bd = betadisper(wunifrac_dist_mat, metadata$rep)
permutest(bray.bd, pairwise=TRUE, permutations=1000)

#PERMANOVA/ADONIS of rep
wunifrac_adonis <- vegan::adonis2(wunifrac_dist_mat ~ rep, data = metadf, by = "margin", permutations = 9999)
wunifrac_adonis

#Pairwise comparison loop
#Get combinations of rep values
cbn <- combn(x = unique(metadf$rep), m = 2)
#Create empty final data frame with 4 columns
pairwise_permanova_df <- as.data.frame(matrix(data = NA, nrow = ncol(cbn), ncol = 4))
colnames(pairwise_permanova_df) <- c("1","2","p","p.adj")

#Loop through the combinations
for(i in 1:ncol(cbn)){
  pseq_rarefy_subset <- phyloseq::subset_samples(pseq_rarefy, rep %in% cbn[,i])
  metadf_subset <- data.frame(phyloseq::sample_data(pseq_rarefy_subset))
  wunifrac_dist_mat <- phyloseq::distance(pseq_rarefy_subset, method="wunifrac")
  wunifrac_pairwise_adonis <- vegan::adonis2(wunifrac_dist_mat ~ rep, 
                                             data = metadf_subset, by = "margin")
  pairwise_permanova_df[i,1:2] <- cbn[,i]
  pairwise_permanova_df[i,3] <- wunifrac_pairwise_adonis[1,"Pr(>F)"]
}

#Add adjusted P-values
pairwise_permanova_df$p.adj <- p.adjust(pairwise_permanova_df$p, method = "BH")

###Plot as Bar Chart###

otu.tab <- t(otu_table(pseq_rarefy))
treefile <- phy_tree(pseq_rarefy)

# calculate the UniFracs & create dataframe
unifracs <- GUniFrac(otu.tab, treefile, alpha = c(0, 0.5, 1))$unifracs
d5 <- unifracs[, , "d_0.5"] 

# convert to long form
df <- melt(d5, id.vars = "row_names", variable.name = "col_names", value.name = "score")

#Create columns for each variables group comparison
group_vector1 <- character(nrow(df))
group_vector2 <- character(nrow(df))

# Loop through each row of the data frame
for (i in 1:nrow(df)) {
  # Check conditions and assign group name for Var1
  group_vector1[i] <- ifelse(grepl("Dulse", df$Var1[i]) & grepl("Agar", df$Var1[i]), "Palmaria SSM",
                             ifelse(grepl("Dulse", df$Var1[i]) & grepl("MB", df$Var1[i]), "Palmaria MB",
                                    ifelse(grepl("Dulse", df$Var1[i]) & grepl("Microbiome", df$Var1[i]), "Palmaria Uncultured",
                                           ifelse(grepl("Fucus", df$Var1[i]) & grepl("Agar", df$Var1[i]), "Fucus SSM",
                                                  ifelse(grepl("Fucus", df$Var1[i]) & grepl("MB", df$Var1[i]), "Fucus MB",
                                                         ifelse(grepl("Fucus", df$Var1[i]) & grepl("Microbiome", df$Var1[i]), "Fucus Uncultured", NA)))))) 
  # Check conditions and assign group name for Var2
  group_vector2[i] <- ifelse(grepl("Dulse", df$Var2[i]) & grepl("Agar", df$Var2[i]), "Palmaria SSM",
                             ifelse(grepl("Dulse", df$Var2[i]) & grepl("MB", df$Var2[i]), "Palmaria MB",
                                    ifelse(grepl("Dulse", df$Var2[i]) & grepl("Microbiome", df$Var2[i]), "Palmaria Uncultured",
                                           ifelse(grepl("Fucus", df$Var2[i]) & grepl("Agar", df$Var2[i]), "Fucus SSM",
                                                  ifelse(grepl("Fucus", df$Var2[i]) & grepl("MB", df$Var2[i]), "Fucus MB",
                                                         ifelse(grepl("Fucus", df$Var2[i]) & grepl("Microbiome", df$Var2[i]), "Fucus Uncultured", NA)))))) 

}

# Add the group vectors as new columns to the data frame
df$group1 <- group_vector1
df$group2 <- group_vector2

# Dulse Barplot
subset_df_dulse <- subset(df, group1 == "Palmaria Uncultured" &
                      (group2 %in% c("Palmaria SSM", "Palmaria MB", "Palmaria Uncultured", "Fucus Uncultured")))

# Select columns 3 to 5
subset_df_dulse <- subset_df_dulse[, 3:5]
subset_df_dulse <- subset(subset_df_dulse, score != 0.0000000)
#Make groups factors
subset_df_dulse$group1 <- as.factor(subset_df_dulse$group1)
subset_df_dulse$group2 <- as.factor(subset_df_dulse$group2)
summary(subset_df_dulse)

subset_df_dulse$group2 <- factor(subset_df_dulse$group2, levels = c("Palmaria Uncultured", "Palmaria SSM", "Palmaria MB", "Fucus Uncultured"))

barplot_dulse <- subset_df_dulse %>%
  ggplot(aes(x = group2, y = score, fill = group2)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 0.9)) +
  facet_grid(~group1, labeller = labeller(group1 = c(
    "Palmaria Uncultured" = "Compared to \n Uncultured P. palmaria",
    "Fucus" = "Compared to \n Uncultured F. serratus"      
  ))) +
  ylab("Weighted Unifrac Distance") +
  theme_bw() +
  scale_fill_manual(values = c("Palmaria Uncultured" = "#660033", 
                               "Palmaria SSM" = "#1F271B", 
                               "Palmaria MB" = "#DC9596", 
                               "Fucus Uncultured" = "#669966")) +
  scale_x_discrete(labels = function(x) str_wrap(
    case_when(
      x == "Palmaria Uncultured" ~ "Uncultured P. palmaria",
      x == "Palmaria SSM" ~ "SSM",
      x == "Palmaria MB" ~ "MB",
      x == "Fucus Uncultured" ~ " Uncultured F. serratus"
    ), width = 15)) +  
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, face = "bold.italic"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 10),  
        axis.title.y = element_text(size = 12),  
        strip.text = element_text(size = 12, face = "bold.italic"),  
        plot.title = element_text(size = 16),
        strip.background = element_rect(fill = "lightgrey"))

# Fucus Barplot

subset_df_fucus <- subset(df, group1 == "Fucus Uncultured" &
                            (group2 %in% c("Fucus SSM", "Fucus MB", "Palmaria Uncultured", "Fucus Uncultured")))

# Select columns 3 to 5
subset_df_fucus <- subset_df_fucus[, 3:5]
subset_df_fucus <- subset(subset_df_fucus, score != 0.0000000)
#Make groups factors
subset_df_fucus$group1 <- as.factor(subset_df_fucus$group1)
subset_df_fucus$group2 <- as.factor(subset_df_fucus$group2)
summary(subset_df_fucus)
subset_df_fucus$group2 <- factor(subset_df_fucus$group2, levels = c("Fucus Uncultured", "Fucus SSM", "Fucus MB", "Palmaria Uncultured"))

barplot_fucus <-subset_df_fucus %>%
  ggplot(aes(x = group2, y = score, fill = group2)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 0.9)) +
  facet_grid(~group1, labeller = labeller(group1 = c(
    "Palmaria Uncultured" = "Compared to \n Uncultured P. palmaria",  
    "Fucus Uncultured" = "Compared to \n Uncultured F. serratus"      
  ))) +
  ylab("Weighted Unifrac Distance") +
  theme_bw() +
  scale_fill_manual(values = c("Fucus Uncultured" = "#669966", 
                               "Fucus SSM" = "#1F271B", 
                               "Fucus MB" = "#DC9596", 
                               "Palmaria Uncultured" = "#660033")) +
  scale_x_discrete(labels = function(x) str_wrap(
    case_when(
      x == "Fucus Uncultured" ~ "F. serratus Uncultured",
      x == "Fucus SSM" ~ "SSM",
      x == "Fucus MB" ~ "MB",
      x == "Palmaria Uncultured" ~ "P. palmaria Uncultured"
    ), width = 15)) +  
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, face = "bold.italic"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 10),  
        axis.title.y = element_text(size = 12),  
        strip.text = element_text(size = 12),  
        plot.title = element_text(size = 16),
        strip.background = element_rect(fill = "lightgrey"))+ theme(strip.text = element_text(face = "bold.italic"))

barplots <-grid.arrange(barplot_dulse, barplot_fucus, ncol=2)

# Create plot labels
label_a <- textGrob("a)", x = unit(0, "npc"), y = unit(1, "npc"), 
                    just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
label_b <- textGrob("b)", x = unit(0, "npc"), y = unit(1, "npc"), 
                    just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))

# Add labels to the plots
NMDS_Plot_labelled <- arrangeGrob(NMDS_Plot, top = label_a)
barplots_labelled <- arrangeGrob(barplots, top = label_b)

# Arrange plots without immediately plotting
Fig3 <- arrangeGrob(NMDS_Plot_labelled, barplots_labelled)

# Save the arranged plot as a PNG
png("Figure3.png", width = 3000, height = 3500, res = 300)
grid.draw(Fig3)  # Draw the grob to the file
dev.off()
####ANCOMBC####
ancom<- read.delim("ANCOMBC.txt")
# Create a column to distinguish between enriched and depleted
ancom <- ancom %>%
  mutate(status = ifelse(log_fold_change > 0, "enriched", "depleted"))

# Plot the data
figS3 <-ggplot(ancom, aes(x = log_fold_change, y = reorder(Bacteria, log_fold_change), 
               fill = status)) +
  geom_bar(stat = "identity") +
  geom_errorbarh(aes(xmin = log_fold_change - se, xmax = log_fold_change + se), 
                 height = 0.2) +
  scale_fill_manual(values = c("enriched" = "#669966", "depleted" = "#660033"),
                    labels = c("depleted" = expression("Depleted in " * italic("P. palmata")), 
                               "enriched" = expression("Enriched in " * italic("F. serratus")))) +
  labs(x = "Log Fold Change (LFC)", y = "", fill = NULL) +  # Set fill to NULL to
#remove legend title
  theme_minimal()
figS3
png("FigureS3.png", width = 2500, height = 1500, res = 300)
grid.draw(figS3)  # Draw the grob to the file
dev.off()
