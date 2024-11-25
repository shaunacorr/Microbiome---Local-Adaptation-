#Load the packages/libraries required 
library("qiime2R")
library("vegan")
library("ggplot2")
library("phyloseq")
library("microbiome")
library("dplyr")
library("tidyr")
library("viridis")
library("car")
library("ggsignif")
library("multcomp")


####Import data####
pseq_with_blanks <- qiime2R::qza_to_phyloseq(
  features = "table-with-phylum.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.sklearn.qza",
  metadata = "Metadata.txt"
)

#Filter out blank samples
sample_data(pseq_with_blanks)
pseq <- subset_samples(pseq_with_blanks, Rep != "Blank")

####Rarefaction Curves####

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
# plot
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Read Depth",
     ylab = "ASVs", type = "n")
abline(v = raremax)
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

####Alpha Diversity#####

#Produce data frame of all alpha diversity values
alpha_df <- phyloseq::estimate_richness(physeq = pseq_rarefy)
head(alpha_df)

metadata <- sample_data(pseq_rarefy)

alpha_df1 <- cbind(metadata, alpha_df)

#Get averages and SDs for all treatments
alpha_summary <- alpha_df1 %>%
  group_by(Rep, Microbiome, Agar) %>%
  summarise(
    Mean_Shannon = mean(Shannon),
    SD_Shannon = sd(Shannon),
    Mean_Chao1 = mean(Chao1),
    SD_Chao1 = sd(Chao1),
    Mean_Observed = mean(Observed),
    SD_Observed = sd(Observed)
  )

#####Diversity Graph - Fig. 5####

# Reverse alphabetic order for Microbiome
alpha_summary$Microbiome <- factor(alpha_summary$Microbiome, levels = rev(sort(unique(alpha_summary$Microbiome))))
# Alphabetical order for Agar
alpha_summary$Agar <- factor(alpha_summary$Agar, levels = sort(unique(alpha_summary$Agar)))

#Scale by sympatric environment
# Create the reference values
palmaria_ref <- alpha_summary[alpha_summary$Microbiome == "Palmaria" & alpha_summary$Agar == "Palmaria", "Mean_Shannon"]
calli_ref <- alpha_summary[alpha_summary$Microbiome == "Calliblepharis" & alpha_summary$Agar == "Calliblepharis", "Mean_Shannon"]
cera_ref <- alpha_summary[alpha_summary$Microbiome == "Ceramium" & alpha_summary$Agar == "Ceramium", "Mean_Shannon"]
lome_ref <- alpha_summary[alpha_summary$Microbiome == "Lomentaria" & alpha_summary$Agar == "Lomentaria", "Mean_Shannon"]
osmu_ref <- alpha_summary[alpha_summary$Microbiome == "Osmundea" & alpha_summary$Agar == "Osmundea", "Mean_Shannon"]
chon_ref <- alpha_summary[alpha_summary$Microbiome == "Chondrus" & alpha_summary$Agar == "Chondrus", "Mean_Shannon"]

# Create a new ref_column using ifelse for each Microbiome
alpha_summary$ref_column <- NA

alpha_summary$ref_column <- ifelse(alpha_summary$Microbiome == "Palmaria", palmaria_ref,
                                     ifelse(alpha_summary$Microbiome == "Calliblepharis", calli_ref,
                                            ifelse(alpha_summary$Microbiome == "Ceramium", cera_ref,
                                                   ifelse(alpha_summary$Microbiome == "Lomentaria", lome_ref,
                                                          ifelse(alpha_summary$Microbiome == "Osmundea", osmu_ref,
                                                                 ifelse(alpha_summary$Microbiome == "Chondrus", chon_ref, NA))))))

#Create scaling column

alpha_summary$Mean_Shannon <- as.numeric(as.character(alpha_summary$Mean_Shannon))
alpha_summary$ref_column <- as.numeric(as.character(alpha_summary$ref_column))
alpha_summary$shannon_sympatric<-alpha_summary$Mean_Shannon/alpha_summary$ref_column
alpha_summary$SD_Shannon <- round(alpha_summary$SD_Shannon, 2)
alpha_summary$sig <- c("","","***","","","","","","**","","","","","","","","","","","","**","","","","***","","***","","","*","","","***","","","")
alpha_summary$sig <- as.character(alpha_summary$sig)

#Shannon's Diversity Plot

ggplot(alpha_summary, aes(x = Agar, y = Microbiome)) +
  geom_point(aes(size = shannon_sympatric, fill = Mean_Shannon), shape = 21, color = "NA") +
  geom_text(aes(label = SD_Shannon), vjust = 3.5, size = 3, hjust = 0.5) +
  geom_text(aes(label = sig), vjust = -1, size = 5, hjust = 0.5) +
  scale_fill_viridis_c(option = "magma") +
  scale_size(range = c(3, 15), guide = "none") +
  labs(title = NULL,
       x = "Seaweed Derived Agar",
       y = "Innoculated Microbiome",
       fill = "Shannon's Diversity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(vjust = 0.5, face = "italic"),
    axis.text.y = element_text(face = "italic"),
    axis.title.x = element_text(vjust = -7, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(position = "top")


#####Summary of Host Environment vs. Microbiome####

#Get averages and SDs for all agar environments
alpha_summary_Agar <- alpha_df1 %>%
  group_by(Agar) %>%
  summarise(
    Mean_Shannon = mean(Shannon),
    SD_Shannon = sd(Shannon),
    Mean_Chao1 = mean(Chao1),
    SD_Chao1 = sd(Chao1),
    Mean_Observed = mean(Observed),
    SD_Observed = sd(Observed)
  )

shapiro.test(alpha_df1$Shannon)
kruskal.test(Shannon ~ Agar, data=alpha_df1)
pairwise.wilcox.test(alpha_df1$Shannon, alpha_df1$Agar, p.adjust.method="fdr")


# Define the list of comparisons
comparisons_agar_sha <- list(
  c("Calliblepharis", "Ceramium"),
  c("Calliblepharis", "Chondrus"),
  c("Calliblepharis", "Palmaria"),
  c("Osmundea", "Ceramium"),
  c("Osmundea", "Chondrus"),
  c("Osmundea", "Lomentaria"),
  c("Osmundea", "Palmaria")
)

#Plot
ggplot(alpha_df1, aes(x = Agar, y = Shannon, fill = Agar)) +  
  geom_boxplot(color = "black", alpha = 0.7) +  
  geom_jitter(aes(color = Agar), width = 0.2, size = 2, alpha = 0.8) +  
  scale_fill_viridis_d(option = "D") +  
  scale_color_viridis_d(option = "D") +  
  labs(title = "", 
       x = "Seaweed Derived Agars", 
       y = "Shannon's Diveristy") +
  geom_signif(comparisons = comparisons_agar_sha, 
              annotations = c("*","*","**","**","**","**","**"),
              y_position = c(1.4,1.6,2.6,1.1,1.2,1.3,1.35),
              textsize = 4, data = alpha_df1,tip_length = 0) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),  
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic")
  )


#### Compare by Sympatric Environment #####

#subset by microbiome
df_os <- subset(alpha_df1, Microbiome == "Osmundea")
df_pal <- subset(alpha_df1, Microbiome == "Palmaria")
df_cal <- subset(alpha_df1, Microbiome == "Calliblepharis")
df_cho <- subset(alpha_df1, Microbiome == "Chondrus")
df_cer <- subset(alpha_df1, Microbiome == "Ceramium")
df_lom <- subset(alpha_df1, Microbiome == "Lomentaria")

###Calliblepharis###
cal_shan <- lm(Shannon ~ Agar, data = df_cal)
summary(cal_shan)
cal_shan_1 <- glht(cal_shan, linfct = mcp(Agar = "Dunnett"))
summary(cal_shan_1)

# Test normality of residuals
shapiro.test(residuals(cal_shan))
# Q-Q plot
qqnorm(residuals(cal_shan))
qqline(residuals(cal_shan), col = "blue")

#Test for homogeneity of variances
leveneTest(Shannon ~ Agar, data = df_cal)
plot(fitted(cal_shan), residuals(cal_shan))
abline(h = 0, col = "red")


####Ceramium###
df_cer$Agar <- relevel(df_cer$Agar, ref = "Ceramium")
cer_shan <- lm(Shannon ~ Agar, data = df_cer)
summary(cer_shan)
cer_shan_1 <- glht(cer_shan, linfct = mcp(Agar = "Dunnett"))
summary(cer_shan_1)

#Test for normality of residuals
shapiro.test(residuals(cer_shan))
qqnorm(residuals(cer_shan))
qqline(residuals(cer_shan), col = "blue")

#Test for homogeneity of variances
leveneTest(Shannon ~ Agar, data = df_cer)
plot(fitted(cer_shan), residuals(cer_shan))
abline(h = 0, col = "red")


###Chondrus###
df_cho$Agar <- relevel(df_cho$Agar, ref = "Chondrus")
cho_shan <- lm(Shannon ~ Agar, data = df_cho)
summary(cho_shan)
cho_shan_1 <- glht(cho_shan, linfct = mcp(Agar = "Dunnett"))
summary(cho_shan_1)

# Test for normality of residuals
shapiro.test(residuals(cho_shan))
qqnorm(residuals(cho_shan))
qqline(residuals(cho_shan), col = "blue")
#Test for homogeneity of variances
leveneTest(Shannon ~ Agar, data = df_cho)
plot(fitted(cho_shan), residuals(cho_shan))
abline(h = 0, col = "red")

kruskal.test(Shannon ~ Agar, data = df_cho) 
pairwise.wilcox.test(df_cho$Shannon, df_cho$Agar, p.adjust.method="fdr")

###Lomentaria###
df_lom$Agar <- relevel(df_lom$Agar, ref = "Lomentaria")
lo_shan <- lm(Shannon ~ Agar, data = df_lom)
summary(lo_shan)
lo_shan_1 <- glht(lo_shan, linfct = mcp(Agar = "Dunnett"))
summary(lo_shan_1)

#Test for normality of residuals
shapiro.test(residuals(lo_shan))
qqnorm(residuals(lo_shan))
qqline(residuals(lo_shan), col = "blue")
#Test for homogeneity of variances
leveneTest(Shannon ~ Agar, data = df_lom)
plot(fitted(lo_shan), residuals(lo_shan))
abline(h = 0, col = "red")

kruskal.test(Shannon ~ Agar, data = df_lom) 
pairwise.wilcox.test(df_lom$Shannon, df_lom$Agar, p.adjust.method="fdr")

#Osmundea
df_os$Agar <- relevel(df_os$Agar, ref = "Osmundea")
os_shan <- lm(Shannon ~ Agar, data = df_os)
summary(os_shan)
os_shan_1 <- glht(os_shan, linfct = mcp(Agar = "Dunnett"))
summary(os_shan_1)

# Test for normality of residuals
shapiro.test(residuals(os_shan))
qqnorm(residuals(os_shan))
qqline(residuals(os_shan), col = "blue")
# Test for homogeneity of variances
leveneTest(Shannon ~ Agar, data = df_os)
plot(fitted(os_shan), residuals(os_shan))
abline(h = 0, col = "red")

kruskal.test(Shannon ~ Agar, data = df_os) 

###Palmaria###
df_pal$Agar <- relevel(df_pal$Agar, ref = "Palmaria")
pal_shan <- lm(Shannon ~ Agar, data = df_pal)
summary(pal_shan)
pal_shan_1 <- glht(pal_shan, linfct = mcp(Agar = "Dunnett"))
summary(pal_shan_1)

# Test for normality of residuals
shapiro.test(residuals(pal_shan))
qqnorm(residuals(pal_shan))
qqline(residuals(pal_shan), col = "blue")
# Test for homogeneity of variances
leveneTest(Shannon ~ Agar, data = df_pal)
plot(fitted(pal_shan), residuals(pal_shan))
abline(h = 0, col = "red")

##### Plot Together ##########

# Custom labels
custom_labels <- c(
  Calliblepharis = "Calliblepharis\nMicrobiome", 
  Ceramium = "Ceramium\nMicrobiome", 
  Chondrus = "Chondrus\nMicrobiome", 
  Lomentaria = "Lomentaria\nMicrobiome", 
  Palmaria = "Palmaria\nMicrobiome", 
  Osmundea = "Osmundea\nMicrobiome"
)

comparisons_micro_cal <- list(
  c("Calliblepharis", "Ceramium"),
  c("Calliblepharis", "Chondrus"),
  c("Calliblepharis", "Lomentaria"),
  c("Calliblepharis", "Palmaria"),
  c("Calliblepharis", "Osmundea")
)

comparisons_micro_lom <- list(
   c("Lomentaria", "Osmundea")
)
comparisons_micro_pal <- list(
  c("Palmaria", "Osmundea")
)


ggplot(alpha_df1, aes(x = Agar, y = Shannon, fill = Agar)) +
  geom_boxplot() +
  facet_wrap(~Microbiome, scales = "free_y", nrow = 1,labeller = as_labeller(custom_labels)) +
  theme_bw() +
  scale_fill_viridis_d(option = "D") +  
  scale_color_viridis_d(option = "D") + 
  labs(x = "Seaweed Derived Agar", y = "Shannon's Diversity") +
  geom_signif(comparisons = comparisons_micro_cal, 
              annotations = c("**", "**","***","***","***"),
              y_position = c(2.05,1.9,1.8,1.6,1.7),
              textsize = 4, data = alpha_df1 %>% filter(Microbiome == "Calliblepharis"), tip_length = 0) +
  geom_signif(comparisons = comparisons_micro_lom, 
              annotations = c("***"),
              y_position = c(1.75),
              textsize = 4, data = alpha_df1 %>% filter(Microbiome == "Lomentaria"), tip_length = 0) +
  geom_signif(comparisons = comparisons_micro_pal, 
              annotations = c("*"),
              y_position = c(1.75),
              textsize = 4, data = alpha_df1 %>% filter(Microbiome == "Palmaria"), tip_length = 0) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
    axis.ticks.x = element_line(),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "lightgrey"),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
  )

