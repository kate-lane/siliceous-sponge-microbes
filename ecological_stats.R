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

# // Which groups differ? 
pairwise.adonis(community_data,data$Group) #Because Group is significant, while Location is not

#test pairwise combinations of group and Location
pairwise.adonis(community_data, factors = interaction(data$Group, data$Location))



## Simper analysis

# // Which microbes drive the differences between sponge species?
simper(community_data,data$Group,permutations=999,parellel=1)

# Also checked for homogeneity of variance using bray-curtis dissimilarity matrix:
bc_dist <- vegdist(community_data, method = "bray")

dispersion_group <- betadisper(bc_dist, data$Group)
permutest(dispersion_group)

dispersion_location <- betadisper(bc_dist, data$Location)
permutest(dispersion_location)
#P values are both greater than 0, so dispersion is homogenous


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
# // Assumption of homogeneity of variances is met; p is greater than .05 so proceed with anova
leveneTest(Simpson ~ Group * Location, data = diversity_summary)

# Simpson Index data is univariate so using aov() instead of adonis because it isn't community data 
anova_test<-aov(Simpson ~ Group + Location + Group*Location, data = diversity_summary)
summary(anova_test)

# Post-hoc test 
# // Which groups differ? 
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
  
  # Subset metadata and community data
  subset_meta <- data[data$Group == family, ]
  subset_comm <- community_data[data$Group == family, ]
  
  # Run PERMANOVA
  permanova_results <- adonis2(subset_comm ~ Location, data = subset_meta, permutations = 999, method = "bray")
  print(permanova_results)

}
