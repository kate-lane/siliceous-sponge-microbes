library(ggplot2)
library(vegan)
library(dplyr)
library(pairwiseAdonis)
library(tidyr)
library(car)

# setwd('')
data <- read.csv("Glass_Sponges_Archaeal_Bacterial_Phylum_Community_Composition_COVERAGE.csv")

# Clean Data
community_data <- data[, !(names(data) %in% c("Dive", "Location", "Sample", "Smithsonian", "Letter", "Genus", "Family", "Group"))]
community_data[is.na(community_data)] <- 0   #Replace NaN values with 0
data$Location <- as.factor(data$Location)
data$Group <- as.factor(data$Group)

## Individual PERMANOVAs: Location and Sponge Family

# Sponge Family Taxa 
adonis2(community_data~data$Group,data=data,permutations = 999,method="bray")
# Location of Specimen
adonis2(community_data~data$Location,data=data,permutations = 999,method="bray")

## Two-way crossed PERMANOVA

group_loc <- data.frame(Location = data$Location, Group = data$Group) # Extract Location and Group fields 

# Due to low degrees of freedom this collapses the model/ does not show Location, Group, and Location:Group
adonis2(community_data ~ Location * Group, data = group_loc, permutations = 999, method = "bray")

# Run with by margin to accommodate low degrees of freedom and low number of samples
adonis2(community_data ~ Location * Group, data = group_loc, permutations = 999, method = "bray", by = "margin")

# Interaction term between Location and Group
adonis2(community_data ~ Location:Group, data = group_loc, permutations = 999, method = "bray")



## Pairwise post hoc test

pairwise.adonis(community_data,data$Group) #Because Group is significant, while Location is not

#test pairwise combinations of group and Location
pairwise.adonis(community_data, factors = interaction(data$Group, data$Location))



## Simper analysis

simper(community_data,data$Group,permutations=999,parellel=1)
# Check for homogeneity of variance using bray-curtis dissimilarity matrix:
bc_dist <- vegdist(community_data, method = "bray")

dispersion_group <- betadisper(bc_dist, data$Group)
permutest(dispersion_group)

dispersion_location <- betadisper(bc_dist, data$Location)
permutest(dispersion_location)

####################################################


## Simpsons Diversity Index (1-D)

# Calculate Simpsons Index for each sample 
simpson_index <- diversity(community_data, index = "simpson")

diversity_summary <- data.frame(Sample = data$Sample, 
                                 Simpson = simpson_index, 
                                 Group = data$Group,
                                 Location = data$Location)  

# Fit to linear model
diversity_model <- lm(Simpson ~ Group * Location, data = diversity_summary)
summary(diversity_model)

# Plot normality of residuals - distribution appears normal
qqnorm(residuals(diversity_model))
qqline(residuals(diversity_model), col = "red")

# Levene's test for homogeneity of variance
# Assumption of homogeneity of variances is met; p is greater than .05 so proceed with anova
leveneTest(Simpson ~ Group * Location, data = diversity_summary)

# Simpson Index data is univariate so using aov() instead of adonis because it isn't community data 
anova_test<-aov(Simpson ~ Group + Location + Group*Location, data = diversity_summary)
summary(anova_test)

# Post-hoc test 
TukeyHSD(anova_test)


####################################################


## PERMANOVAs on sponge families collected at multiple sites
# one last thing to check for: Is variance within each family driven by site?

# pull sponge families that were recovered at multiple locations
multi_site_families <- data %>%
  group_by(Group) %>%
  summarize(n_sites = n_distinct(Location)) %>%
  filter(n_sites > 1) %>%
  pull(Group)

# One-way PERMANOVA for Location within each 'multi-site' family

for (family in multi_site_families) {
  cat("PERMANOVA for", family)
  subset_meta <- data[data$Group == family, ] # Subset data
  subset_comm <- community_data[data$Group == family, ]
  permanova_results <- adonis2(subset_comm ~ Location, data = subset_meta, permutations = 999, method = "bray")
  print(permanova_results)
}


####################################################


# NMDS

nmds_result <- metaMDS(community_data, distance = "bray", autotransform = FALSE) 
nmds_df <- data.frame(NMDS1 = nmds_result$points[, 1], NMDS2 = nmds_result$points[, 2], 
                      Sample = data$Sample, Group = data$Group, Location = data$Location, Letter = data$Letter)

NMDS_Figure <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Group, shape = Location, label = Letter)) + 
    geom_point(aes(size = Location), show.legend = TRUE) +  
    geom_text(nudge_y = 0.025, size = 3)  +  # adjust position of letter labels next to point
    labs(x = "NMDS1", y = "NMDS2") +
    scale_color_manual(values = c("Farreidae" = "#E69F00", "Hyalonematidae" = "#56B4E9", 
                                  "Cladorhizidae" = "#009E73", "Euplectellidae" = "#D55E00", 
                                  "Phloeodictyidae" = "#CC79A7", "Euretidae" = "#005A8E")) +  
    scale_shape_manual(values = c("Desecheo Ridge" = 15, "Esperanza Ridge" = 17, 
                                  "Whiting Seamount" = 18, "Ridge SW of Vieques" = 19)) +  
    scale_size_manual(values = c("Desecheo Ridge" = 4, "Esperanza Ridge" = 4, 
                                 "Whiting Seamount" = 5, "Ridge SW of Vieques" = 4)) +  # make the diamond a little larger
    theme_bw() +
    theme(legend.position = "right", 
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),   
          axis.text.y = element_text(size = 10)) + 
    guides(color = guide_legend(title = "Sponge Taxonomic Family", override.aes = list(size = 3), order = 1))

ggsave("NMDS_Figure.pdf", plot = NMDS_Figure, width = 9, height = 6, units = "in")
#Figure S1

# Run NMDS without Phloeodictyidae outlier

data_noP <- data[data$Group != "Phloeodictyidae",]
community_data_noP <- community_data[community_data$Group !=  "Phloeodictyidae",]
nmds_result <- metaMDS(community_data_noP, distance = "bray", autotransform = FALSE) 
nmds_df <- data.frame(NMDS1 = nmds_result$points[, 1], NMDS2 = nmds_result$points[, 2], 
                      Sample = data_noP$Sample, Group = data_noP$Group, Location = data_noP$Location, Letter = data_noP$Letter)

NMDS_Figure_noP <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Group, shape = Location, label = Letter)) + 
    geom_point(aes(size = Location), show.legend = TRUE) +  
    geom_text(nudge_y = 0.025, size = 3)  +  # adjust position of letter labels next to point
    labs(x = "NMDS1", y = "NMDS2") +
    scale_color_manual(values = c("Farreidae" = "#E69F00", "Hyalonematidae" = "#56B4E9", 
                                  "Cladorhizidae" = "#009E73", "Euplectellidae" = "#D55E00", 
                                  "Euretidae" = "#005A8E")) +  
    scale_shape_manual(values = c("Desecheo Ridge" = 15, "Esperanza Ridge" = 17, 
                                  "Whiting Seamount" = 18, "Ridge SW of Vieques" = 19)) +  
    scale_size_manual(values = c("Desecheo Ridge" = 4, "Esperanza Ridge" = 4, 
                                 "Whiting Seamount" = 5, "Ridge SW of Vieques" = 4)) +  # make the diamond a little larger
    theme_bw() +
    theme(legend.position = "right", 
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10),   
          axis.text.y = element_text(size = 10)) + 
    guides(color = guide_legend(title = "Sponge Taxonomic Family", override.aes = list(size = 3), order = 1))

ggsave("NMDS_Figure_noP.pdf", plot = NMDS_Figure_noP, width = 9, height = 6, units = "in")
#Figure 4


## Simpsons Diversity Figure

diversity_plot <- ggplot(diversity_summary, aes(x = Group, y = Simpson, fill = Group)) + 
    geom_boxplot(alpha = 0.3) + 
    geom_point(size = 3, aes(color = Group), position = position_dodge(width = 0.75), na.rm = TRUE) +  
    labs(
        x = "Sponge Taxonomic Family",
        y = "Simpson's Diversity Index") +
    theme_bw() +
    scale_fill_manual(values = c("Farreidae" = "#E69F00", 
                                 "Hyalonematidae" = "#56B4E9", 
                                 "Cladorhizidae" = "#009E73", 
                                 "Euplectellidae" = "#D55E00", 
                                 "Phloeodictyidae" = "#CC79A7", 
                                 "Euretidae" = "#005A8E")) +  
    scale_color_manual(values = c("Farreidae" = "#E69F00", 
                                  "Hyalonematidae" = "#56B4E9", 
                                  "Cladorhizidae" = "#009E73", 
                                  "Euplectellidae" = "#D55E00", 
                                  "Phloeodictyidae" = "#CC79A7", 
                                  "Euretidae" = "#005A8E")) +  
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 13),  
          axis.title.y = element_text(size = 13, margin = margin(r = 10)))  

ggsave("Simpson_Diversity_BoxPlot.pdf", plot = diversity_plot, width = 8, height = 6, units = "in")
#Figure 5
