# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Accompanying code for the paper "Dramatic plant-pollinator network change 
# across more than a century in the subarctic" 
# authored by L. Zoller, J. Bennett, and T. Knight
#
# Script authored by L. Zoller
# Please refer to the paper for detailed methodological explanations 
# and relevant citations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Load packages and data ----------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load required packages
library(data.table)
library(dplyr)
library(bipartite)
library(bootstrapnet)
library(ggplot2)
library(car)
library(tidyverse)
library(vegan)
library(grid)
library(gridGraphics)
library(RColorBrewer)
library(cowplot)


{ 
  # Load pre-prepared data file of plant-pollinator interactions
  kittila <- read.csv("data/InteractionData_Zoller_et_al.csv")               
  
  # Prepare data: Concatenate variables
  # Create variable: "plantBinomen" (Binomial name of the plant species without the authorship part)
  kittila <- kittila %>%
    unite("plantBinomen", plantGenericName, plantSpecificEpithet, sep = " ", remove = FALSE)
  # Create variable: "plantPerEra"
  kittila <- kittila %>%
    unite("plantPerEra", plantBinomen, era, sep = " ", remove = FALSE)
  # Create variable: "floralFormPerEra"
  kittila <- kittila %>%
    unite("floralFormPerEra", floralForm, era, sep = " ", remove = FALSE)
  # Create variable: "plantPerYear"
  kittila <- kittila %>%
    unite("plantPerYear", plantBinomen, year, sep = " ", remove = FALSE)
  # Create variable: "floralFormPerYear"
  kittila <- kittila %>%
    unite("floralFormPerYear", floralForm, year, sep = " ", remove = FALSE)
  # Create variable: "interaction"
  kittila <- kittila %>%
    unite("interaction", animalFinestCommonTaxon, plantBinomen, sep = " : ", remove = FALSE)
  
  setDT(kittila)
}

set.seed(66)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Changes in pollinator community composition across eras -------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# * Proportional composition of pollinator taxonomic groups across eras ---- 

# Sort names according to taxonomy and define colors for plotting
namesTaxa <- c("moths", "butterflies", "other wasps","bumblebees", 
               "solitary bees", "hoverflies", "other flies", "muscoid flies" )
Colours <- c("#68957B","#003D45","#9B0E1D","#FFC309",
             "#F7961D", "#9EC0E3", "#007FC4", "#153F87")


{ 
  # Sum up number of observations (indv) of pollinator functional groups in each era
  dt_counts <- kittila[era %in% c("past", "present"), 
                       .(indv = sum(animalQuantityConservativeEstimate, na.rm = FALSE)),
                       by = c("animalFunctionalGroup", "era")]
  dt_counts$animalFunctionalGroup <- factor(dt_counts$animalFunctionalGroup, levels = namesTaxa)
  
  # Visualize results
  pA <- ggplot(dt_counts, aes(x = era, y = indv, fill = animalFunctionalGroup)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = Colours) +               # customize colors
    labs(fill = "Taxonomic group") +                    # customize legend title
    theme_classic() +                                   # customize theme
    theme(axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
          axis.text.y = element_text(face = "bold", size = 12, color = "black",
                                     angle = 90, hjust = 0.5), 
          axis.line = element_line(colour = 'black', size = 0.25),
          legend.title = element_blank(),
          legend.text  = element_text(size = 12, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    coord_flip()
}


# * Dissimilarity of pollinator communities across eras ----

{
  # Create a community matrix of plants and their pollinators in each era
  dt_counts <- kittila[era %in% c("past","present"),
                       .(indv = sum(animalQuantityConservativeEstimate)),
                       by = c("plantPerEra", "animalFinestCommonTaxon")]
  com_mat <- dcast(dt_counts, plantPerEra ~ animalFinestCommonTaxon,
                   fun.aggregate = length,
                   value.var = "indv", fill = 0)
  
  # For the calculation of NMDS values, the community matrix can only contain abundance information, 
  # No other variables, so the first row from the community matrix needs to be removed
  com_mat <- as.matrix(com_mat[,-1])  
  
  # Create sorted data set with variables for plotting
  dt_var <- sort(kittila[era %in% c("past","present"),
                         .(indv = sum(animalQuantityConservativeEstimate)),
                         by = c("plantPerEra","era", "floralForm", "plantBinomen")],
                 decreasing = FALSE)
}

# Calculate NMDS scores
nmds <- metaMDS(com_mat, distance ="bray")

# Extract data for plotting and add variables relevant for plotting to data frame 
{ 
  data.scores = as.data.frame(scores(nmds))
  data.scores$era = dt_var$era
  data.scores$floralForm = dt_var$floralForm
  data.scores$plantBinomen = dt_var$plantBinomen
}

# Visualize results
{ 
  pB = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(size = 4, stroke = 0.8,  alpha = 0.7, 
               aes(shape = floralForm, color = era)) +
    theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.y = element_text(colour = "black", size = 12, face = "bold", 
                                     angle = 90, hjust = 0.5),  
          legend.text = element_text(size = 12, face = "bold", colour = "black"), 
          legend.position = "right", 
          axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = 0.25),
          #panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          legend.key = element_blank()) + 
    scale_shape_manual(values = c(0:8)) + 
    scale_color_manual(values = c("#a0a0a0", "#25AAE5")) +
    labs(x = "NMDS1", shape = "Floral form", color = "era", y = "NMDS2")
}

# Save figures
#  {
#   CombinedPlot <- plot_grid(pA, pB,
#                             align = "vh",
#                             axis = "tblr",
#                             nrow = 2)
#   CombinedPlot
#  
#   ggsave(path = "figures/",                     # specify path
#          filename = paste("Fig_2AB", ".pdf"),
#          device = "pdf",
#          width = 24, height = 16, units = "cm",
#          colormodel = 'cmyk', bg = "white")
#   }

# Similarity analyses (ANOSIM)
summary(anosim(com_mat, dt_var$era,
               distance = "bray", permutations = 9999))


# * Dissimilarity of pollinator communities across years and plant floral forms within each era ----

# Past:
# Create a community matrix of plants and their pollinators in each year of the past
{
  dt_counts <- kittila[era %in% "past",
                       .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                       by = c("plantPerYear", "animalFinestCommonTaxon")]
  com_mat <- dcast(dt_counts, plantPerYear ~ animalFinestCommonTaxon,
                   value.var = "indv", fill = 0)
  com_mat <- (com_mat[,-1])
}

# Create sorted dataset with relevant variables
dt_var <- sort(df <- kittila[era %in% "past",
                            .(indv = sum(animalQuantityConservativeEstimate)),
                            by = c("plantPerYear", "plantBinomen", "floralForm", "era", "year")],
              decreasing = FALSE)

# Similarity analyses (ANOSIM)
summary(anosim(com_mat, dt_var$year, distance = "bray", permutations = 9999))
summary(anosim(com_mat, dt_var$floralForm, distance = "bray", permutations = 9999))

# Present:
# Create a community matrix of plants and their pollinators in each year of the present
{
  dt_counts <- kittila[era %in% "present",
                       .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                       by = c("plantPerYear", "animalFinestCommonTaxon")]
  com_mat <- dcast(dt_counts, plantPerYear ~ animalFinestCommonTaxon,
                   value.var = "indv", fill = 0)
  com_mat <- (com_mat[,-1])
}

# Create sorted dataset with relevant variables
dt_var <-sort(df <- kittila[era %in% "present",
                            .(indv = sum(animalQuantityConservativeEstimate)),
                            by = c("plantPerYear", "plantBinomen", "floralForm", "era", "year")],
              decreasing = FALSE)

# Similarity analyses (ANOSIM)
summary(anosim(com_mat, dt_var$year, distance = "bray", permutations = 9999))
summary(anosim(com_mat, dt_var$floralForm, distance = "bray", permutations = 9999))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Interaction networks ------------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# * Past interaction network ----

{ 
  # Prepare data: For each era, create a community matrix of plants and their pollinators  
  present <- kittila[era %in% "present",
                     .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                     by = c("plantBinomen", "animalOrder", "animalFunctionalGroup", 
                            "animalFamily", "animalGenericName", "animalFinestCommonTaxon")]
  com_matPresent <- as.data.frame(dcast(present, plantBinomen ~ animalFinestCommonTaxon,
                                        value.var = "indv", fill = 0))
  rownames(com_matPresent) <- com_matPresent$plantBinomen
  com_matPresent <- as.matrix(com_matPresent[,-1])
  
  
  past <- kittila[era %in% "past",
                  .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                  by = c("plantBinomen", "animalOrder", "animalFunctionalGroup", 
                         "animalFamily", "animalGenericName", "animalFinestCommonTaxon")]
  com_matPast<- as.data.frame(dcast(past, plantBinomen ~ animalFinestCommonTaxon,
                                    value.var = "indv", fill = 0))
  rownames(com_matPast) <- com_matPast$plantBinomen
  com_matPast <- as.matrix(com_matPast[,-1])
}

# Prepare data for plotting: Define sequence of plotting
{ 
  PlantList <<- rownames(com_matPast)
  PlantOrder <<- PlantList[order(PlantList)] 
  
  # Sort pollinators by functional group
  kittilaPast <- kittila[order(animalFunctionalGroup) & era == "past",]
  kittilaGroupedTab <- kittilaPast %>% group_by(animalFunctionalGroup) %>% arrange(animalFinestCommonTaxon)
  setDT(kittilaGroupedTab)
  PolliOrderTab <- kittilaGroupedTab[, unique(animalFinestCommonTaxon), keyby=animalFunctionalGroup]  
  
  PolliList <<- colnames(com_matPast)
  PolliTab <- merge(PolliOrderTab, as.data.frame(PolliList), by.x = "V1", by.y = "PolliList")
  Polliorder <- PolliTab[order(animalFunctionalGroup),]$V1
  
  countplants <<- length(PlantList)
  
  # Define colors for each pollinator functional group
  kittilaCol <- rbind(PolliTab[PolliTab$animalFunctionalGroup == "hoverflies",],
                      PolliTab[PolliTab$animalFunctionalGroup == "other flies",],
                      PolliTab[PolliTab$animalFunctionalGroup == "muscoid flies",],
                      PolliTab[PolliTab$animalFunctionalGroup == "solitary bees",],
                      PolliTab[PolliTab$animalFunctionalGroup == "bumblebees",],
                      PolliTab[PolliTab$animalFunctionalGroup == "other wasps",],
                      PolliTab[PolliTab$animalFunctionalGroup == "moths",],
                      PolliTab[PolliTab$animalFunctionalGroup == "butterflies",])
  
  table(kittilaCol$animalFunctionalGroup)
  
  # Adding a column specifying the colors for plotting
  kittilaCol <- kittilaCol %>%
    mutate(colGroup = case_when(
      endsWith(animalFunctionalGroup, "bumblebees") ~ "#FFC100",
      endsWith(animalFunctionalGroup, "butterflies") ~ "#053225",
      endsWith(animalFunctionalGroup, "hoverflies") ~ "#A2D6F9",
      endsWith(animalFunctionalGroup, "moths") ~ "#60A561",
      endsWith(animalFunctionalGroup, "muscoid flies") ~ "#072AC8",
      endsWith(animalFunctionalGroup, "other flies") ~ "#1E96FC",
      endsWith(animalFunctionalGroup, "other wasps") ~ "#96031A",
      endsWith(animalFunctionalGroup, "solitary bees") ~ "#FF8200"
    ))
  
  # Generate color vector for pollinators
  colPoll  <- kittilaCol$colGroup
  
  # Generate one color for plants
  colPlant <- rep("darkolivegreen3", length(PlantList))
  
  # Generate sequence of species
  seqList <- list(seq.higher = kittilaCol$V1, seq.lower=PlantOrder)
}

# Visualize network
{ 
  networkPast <- grid.grabExpr(grid.echo(function() plotweb(sortweb(com_matPast, sort.order = "seq", sequence = seqList), 
                                                            method = "normal",
                                                            labsize = 0.7, 
                                                            col.high = colPoll,
                                                            col.low = colPlant, 
                                                            bor.col.interaction = NA,
                                                            y.width.low = 0.025,
                                                            y.width.high = 0.025,
                                                            high.y = 1.25,           
                                                            low.y = 0.7,              
                                                            text.rot = 90)))
  
  label_grobs <- getGrob(gTree = networkPast,
                         gPath = "graphics-plot-1-text-.*",
                         grep = TRUE,
                         global = TRUE)
  
  # Get species labels
  sp <- rep(NA_character_, length(label_grobs))
  for (i in 1:length(label_grobs)){
    sp[i] <- label_grobs[[i]][["label"]]
  }
  grid.newpage(); grid.draw(networkPast)
  
}

### * Present interaction network (rarefied) ----

{
  # Prepare data: Rarefy present interaction matrix to contain same number of interactions as past interaction matrix
  # Get sample size of each row from past community matrix
  obsPast <- .rowSums(com_matPast, (dim(com_matPast))[[1]],
                      (dim(com_matPast))[[2]], 
                      na.rm = FALSE)
  
  obsPres <- .rowSums(com_matPresent, (dim(com_matPresent))[[1]],
                      (dim(com_matPresent))[[2]], 
                      na.rm = FALSE)
  
  # Select the smallest observations
  obs_min <- data.frame(cbind(obsPast, obsPres))
  obs_min$min <- with(obs_min, pmin(obsPast, obsPres))
  
  # Randomly rarefy community matrix
  com_matPresent_rarefied <- rrarefy(com_matPresent, sample = obs_min$min)
  com_matPresent_rarefied <- as.data.frame.matrix(com_matPresent_rarefied)
  com_matPresent_rarefied <- as.matrix(com_matPresent_rarefied)
  
  # Check if pollinator species are absent in the rarefied matrix
  apply(com_matPresent_rarefied, 2, sum)
  
  # Remove absent pollinators
  com_matPresent_rarefied = com_matPresent_rarefied[,colSums(com_matPresent_rarefied != 0) != 0]
}

{
  # Define sequence of plotting
  PlantList <<- rownames(com_matPresent_rarefied)
  Plantorder <<- PlantList[order(PlantList)] 
  
  # Order pollinators by group for plot
  kittilaPresent <- kittila[order(animalFunctionalGroup) & era == "present",]
  kittilaGroupedTab <- kittilaPresent %>% group_by(animalFunctionalGroup) %>% arrange(animalFinestCommonTaxon)
  setDT(kittilaGroupedTab)
  PolliOrderTab <- kittilaGroupedTab[, unique(animalFinestCommonTaxon), keyby = animalFunctionalGroup]  
  
  PolliList <- colnames(com_matPresent_rarefied)
  PolliTab <- merge(PolliOrderTab, as.data.frame(PolliList), by.x = "V1", by.y = "PolliList", all=F)
  PolliOrder <- PolliTab[order(animalFunctionalGroup),]$V1
  
  countPlants <<- length(PlantList)
  
  # Define Color for each pollinator group
  kittilaCol <- rbind(PolliTab[PolliTab$animalFunctionalGroup == "hoverflies",],
                      PolliTab[PolliTab$animalFunctionalGroup == "other flies",],
                      PolliTab[PolliTab$animalFunctionalGroup == "muscoid flies",],
                      PolliTab[PolliTab$animalFunctionalGroup == "solitary bees",],
                      PolliTab[PolliTab$animalFunctionalGroup == "bumblebees",],
                      PolliTab[PolliTab$animalFunctionalGroup == "other wasps",],
                      PolliTab[PolliTab$animalFunctionalGroup == "moths",],
                      PolliTab[PolliTab$animalFunctionalGroup == "butterflies",])
  
  table(kittilaCol$animalFunctionalGroup)
  
  # Add a column specifying the colors for plotting
  kittilaCol <- kittilaCol %>% 
    mutate(col_group = case_when(
      endsWith(animalFunctionalGroup, "bumblebees") ~ "#FFC100",
      endsWith(animalFunctionalGroup, "butterflies") ~ "#053225",
      endsWith(animalFunctionalGroup, "hoverflies") ~ "#A2D6F9",
      endsWith(animalFunctionalGroup, "moths") ~ "#60A561",
      endsWith(animalFunctionalGroup, "muscoid flies") ~ "#072AC8",
      endsWith(animalFunctionalGroup, "other flies") ~ "#1E96FC",
      endsWith(animalFunctionalGroup, "other wasps") ~ "#96031A",
      endsWith(animalFunctionalGroup, "solitary bees") ~ "#FF8200"))
  
  # Generate color vector for pollinators
  colPoll <- kittilaCol$col_group
  
  # Generate 1 color for plants
  colPlant <- rep("darkolivegreen3", length(PlantList))
  
  # Generate sequence of species
  seqList <- list(seq.higher=kittilaCol$V1, seq.lower=Plantorder)
}

# Plot network
{
  networkPresent <- grid.grabExpr(grid.echo(function() plotweb(
    sortweb(com_matPresent_rarefied, 
            sort.order = "seq", sequence = seqList),
    method = "normal",
    labsize = 0.7, 
    col.high = colPoll,# color pollinators
    col.low = colPlant, # color plants
    bor.col.interaction = NA,
    y.width.low = 0.025,
    y.width.high = 0.025,
    high.y = 1.25, # position of high boxes
    low.y = 0.7, # position of low boxes
    text.rot = 90)))
  
  label_grobs <- getGrob(gTree = networkPresent,
                         gPath = "graphics-plot-1-text-.*",
                         grep = TRUE,
                         global = TRUE)
  
  # Get species labels
  sp <- rep(NA_character_, length(label_grobs))
  for (i in 1:length(label_grobs)){
    sp[i] <- label_grobs[[i]][["label"]]
  }
  
  networkPresent <- editGrob(grob = networkPresent,
                             gPath = "graphics-plot-1-text-2[8-9]|3[0-6]", 
                             gp = gpar(font = 3L),
                             grep = TRUE,
                             global = TRUE)
  
  grid.newpage(); grid.draw(networkPresent)
}

# Save figures
# {
#  CombinedNetworks <- plot_grid(networkPresent, networkPast,
#                                align = "vh",
#                                axis = "tblr",
#                                nrow = 2)
#  CombinedNetworks
# 
#  ggsave(path = "Figures/",
#         filename = paste("Fig_3AB", ".pdf"),
#         device = "pdf",
#         width = 25, height = 17, units = "cm",
#         colormodel = 'cmyk', bg = "white")
#   }


# * Interactions unique and common across eras ----

# Data table containing all interactions from past and present
interactionsCount <- kittila[era %in% c("past", "present"),
                             .(num = sum(animalQuantityConservativeEstimate)),
                             by = c("era","interaction")]

length(unique(interactionsCount$interaction))   # total number of unique interactions

# Interactions from the past
interactionsPast <- interactionsCount[interactionsCount$era == "past",]
length(unique(interactionsPast$interaction))    # number of unique interactions in the past

interactionsPresent <- interactionsCount[ interactionsCount$era == "present",]
length(unique(interactionsPresent$interaction)) # number of unique interactions in the present

# Select interactions from the past that still persist in the present
interactionsPersist <- interactionsPast$interaction %in% interactionsPresent$interaction
table(interactionsPersist)
interactionsCommon <- interactionsPast[interactionsPersist == "TRUE",] 

# Percentage of interactions present in both eras
length(interactionsCommon$interaction)/length(unique(interactionsCount$interaction)) 

# Percentage of past interactions that persist into the present
length(interactionsCommon$interaction)/length(unique(interactionsPast$interaction)) 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Network structure and specialization - Full network -----------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                   

# * Bootstrapped network level specialization (H2´) ----

metric = list(Past = com_matPast, Present = com_matPresent) %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower",   # column name for plants
                    col_higher = "higher", # column name for insects
                    level = "both",
                    index = "H2",          # select network metric
                    start = 20,
                    step = 10,
                    n_boot = 1000,         # number of bootstraps
                    n_cpu = 4)             # number of CPU-s to use

# Visualize metric
p <- gg_networklevel(metric) # link to ggplot

# pdf(paste("Figures/", "Fig_4_H2", ".pdf", sep = ""),  # specify path
#     width = 5, height = 4, colormodel = "cmyk")

p$H2 + theme(panel.background = element_blank(),
             axis.text.y = element_text(colour = "black", size = 12, face = "bold",
                                        angle = 90, hjust = 0.5),
             axis.text.x = element_text(colour = "black", face = "bold", size = 12),
             axis.line = element_line(colour = 'black', size = 0.25),
             axis.title = element_text(face = "bold", size = 14, colour = "black"),
             legend.key = element_blank(),
             legend.position = "right",
             legend.title = element_blank(),
             legend.text = element_text(size = 12,face ="bold", colour ="black")) + 
  xlab("sample size") + 
  ylab("H2´") +
  scale_color_manual(values = c("#a0a0a0", "#25AAE5"))

# dev.off()


# * Regression of change in relative abundance and specialization index (d´) ----

{
  # Create community matrix including both past and present observations
  all <- kittila[era %in%  c("past", "present"),
                 .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                 by = c("plantBinomen", "animalFinestCommonTaxon")]
  com_matAll <- as.data.frame(dcast(all, plantBinomen ~ animalFinestCommonTaxon,
                                    value.var = "indv", fill = 0))
  rownames(com_matAll) <- com_matAll$plantBinomen
  com_matAll <- as.matrix(com_matAll[,-1])
  
  # Create data frame containing variables for community matrix
  dt_All <- kittila[era %in% c("past","present"),
                    .(indv = sum(animalQuantityConservativeEstimate)),
                    by = c("animalOrder", "animalFunctionalGroup",
                           "animalFamily", "animalFinestCommonTaxon")]
  
  # Order animal taxa alphabetically
  dt_All <- dt_All[order(dt_All$animalFinestCommonTaxon),]
}

# Calculate species-level specialization (d´)
dtSpec <- data.frame(specieslevel(com_matAll, index = "d", level = "higher"))

{
  # Add variables to data frame:
  dtSpec$animalOrder <- dt_All$animalOrder
  dtSpec$animalFamily <- dt_All$animalFamily
  dtSpec$animalFunctionalGroup <- dt_All$animalFunctionalGroup
  dtSpec$animalFinestCommonTaxon <-dt_All$animalFinestCommonTaxon
  
  # Count number of observations for each taxon in past and present:
  dtEra <- kittila[era %in% c("past","present"),
                   .(indv = sum(animalQuantityConservativeEstimate)),
                   by = c("era", "animalOrder","animalFunctionalGroup",
                          "animalFamily", "animalFinestCommonTaxon" )]
  
  dtEra <- dtEra[order(dtEra$animalFinestCommonTaxon),]  # order pollinator taxa alphabetically
  
  No_obsPast <- data.frame(animalFinestCommonTaxon = dtEra$animalFinestCommonTaxon[dtEra$era =="past"],
                           No_obsPast = dtEra$indv[dtEra$era =="past"])
  No_obsPresent <- data.frame(animalFinestCommonTaxon = dtEra$animalFinestCommonTaxon[dtEra$era =="present"], 
                              No_obsPresent = dtEra$indv[dtEra$era =="present"])
  Obs <- merge(No_obsPast, No_obsPresent, by = "animalFinestCommonTaxon", all = TRUE)
  Obs[is.na(Obs)] <- 0
  
  # Add observations from past and present to dtSpec data frame:
  dtSpec$No_obsPast <- Obs$No_obsPast  # number of observations in past
  dtSpec$No_obsPresent <- Obs$No_obsPresent  # number of observations in present
  dtSpec$no_obsPooled <- dtSpec$No_obsPast + dtSpec$No_obsPresent  # total number of observations
  dtSpec$relAbundPast <- dtSpec$No_obsPast/sum(dtSpec$No_obsPast)  # relative abundance in past
  dtSpec$relAbundPresent <- dtSpec$No_obsPresent/sum(dtSpec$No_obsPresent)  # relative abundance in present
  dtSpec$changeRelAbund <- dtSpec$relAbundPresent - dtSpec$relAbundPast  # change in relative abundance
  setDT(dtSpec)
}


# * Histogram showing the frequency of observations of pollinators ----

ggplot(dtSpec, aes(x = no_obsPooled)) + 
  geom_histogram(stat = "bin", 
                 color = "black",
                 breaks = seq(0, max(dtSpec$no_obsPooled) + 5, 5),
                 fill = "grey60") +
  scale_x_continuous(limits = c(0, 350), 
                     breaks = seq(0, 350, by = 5), 
                     labels = c(   0, rep("", 9),  50, rep("", 9), 
                                   100, rep("", 9), 150, rep("", 9),
                                   200, rep("", 9), 250, rep("", 9),
                                   300, rep("", 9), 350)) + 
  xlab("Number of observations") +
  ylab("Frequency") +
  theme_classic() +
  theme(axis.text = element_text(face= "bold", size = 12, color = "black"),
        legend.text  = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")) 

# ggsave(path = "Figures/",                     # specify path
#        filename = paste("Fig_S7", ".pdf"),
#        device = "pdf",
#        width = 12, height = 12, units = "cm",
#        colormodel = 'cmyk', bg = "white")


# Linear regression:
dtSpec <- dtSpec[dtSpec$no_obsPooled > 10]       # select only species with >= 10 observations total

m <- lm((dtSpec$d) ~ (dtSpec$changeRelAbund))
summary(m)

# Visualize regression:
# customize ggplot for linear regression:
ggplotRegression <- function (fit) {
  ggplot(fit$model, aes_string(x = names(fit$model)[1], y = names(fit$model)[2])) + 
    geom_point(col = "grey60") +
    stat_smooth(method = "lm", col = "black") +
    theme(plot.title = element_text(size = 12),
          axis.text.y = element_text(colour = "black", size = 12, face = "bold",
                                     angle = 90, hjust = 0.5), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank())
}


p <- ggplotRegression(m) +
  xlab("d`") +
  ylab("Δ rel. abundance") + 
  geom_hline(yintercept = 0, 
               colour = "#808080", 
               lty = 5)

# ggsave(path = "Figures/",
#       filename = paste("Fig_4_Regression.pdf"),
#       width = 9,
#       height = 6.1,
#       units = "cm",
#       dpi = 300,
#       device = "pdf")


# * Dissimilarity of interaction networks between eras ----

betaDiv <- betalinkr(webs2array(list(past = com_matPast,
                                     present = com_matPresent)), 
                     index = "bray", # use Bray-Curtis dissimilarity index 
                     binary = FALSE, # use abundance info
                     partitioning = "commondenom",
                     partition.st = FALSE,
                     partition.rr = FALSE)
betaDiv



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Network structure and specialization - Taxonomic subsets ------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# assign "Diptera", Hymenoptera" or "Lepidoptera" to TaxonSelect for all flies, bees & wasps or moths & butterflies subsets,
# assign "Syrphidae" to TaxonSelect and replace animalOrder with animalFamily for hoverfly subset
TaxonSelect <- "Diptera" 
Taxon <-  kittila[animalOrder %in% TaxonSelect]

{ 
  # Prepare community matrices for single taxon
  TaxonPast <- Taxon[era %in% "past",
                     .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                     by = c("plantBinomen", "animalFinestCommonTaxon")]
  com_matTaxPast <- as.data.frame(dcast(TaxonPast, plantBinomen ~ animalFinestCommonTaxon,
                                        value.var = "indv", fill = 0))
  rownames(com_matTaxPast) <- com_matTaxPast$plantBinomen
  com_matTaxPast <- as.matrix(com_matTaxPast[,-1])
  
  TaxonPresent <- Taxon[era %in% "present",
                        .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                        by = c("plantBinomen", "animalFinestCommonTaxon")]
  com_matTaxPresent <- as.data.frame(dcast(TaxonPresent, plantBinomen ~ animalFinestCommonTaxon,
                                           value.var = "indv", fill = 0))
  rownames(com_matTaxPresent) <- com_matTaxPresent$plantBinomen
  com_matTaxPresent <- as.matrix(com_matTaxPresent[,-1])
}


# * Bootstrapped network level specialization (H2´) ----

metricTaxon = list(Past = com_matTaxPast, Present = com_matTaxPresent) %>%
  lapply(web_matrix_to_df) %>%
  boot_networklevel(col_lower = "lower",   # column name for plants
                    col_higher = "higher", # column name for insects
                    level = "both",
                    index = "H2",          # select network metric
                    start = 20,
                    step = 10,
                    n_boot = 1000,         # number of bootstraps
                    n_cpu = 4)             # number of CPU-s to use


{
 # Visualize bootstrapped network metric
 p <- gg_networklevel(metricTaxon)

 # pdf(paste("Figures/", "Fig_5_H2_", TaxonSelect, ".pdf", sep = ""), 
 #     width = 5, height = 4, colormodel = 'cmyk')

 p$H2 + theme(panel.background = element_blank(),
              axis.text.y = element_text(colour = "black", size = 12, face = "bold",
                                         angle = 90, hjust = 0.5),
              axis.text.x = element_text(colour = "black", face = "bold", size = 12),
              axis.line = element_line(colour = 'black', size = 0.25),
              axis.title = element_text(face = "bold", size = 14, colour = "black"),
              legend.key = element_blank(),
              legend.position = "right",
              legend.title = element_blank(),
              legend.text = element_text(size = 12,face ="bold", colour ="black")) + 
        xlab("sample size") + 
        ylab("H2´") +
        scale_color_manual(values = c("#a0a0a0", "#25AAE5"))

# dev.off()
}

# * Regressions of change in relative abundance and specialization index d´  ----

# assign "Diptera", Hymenoptera" or "Lepidoptera" to TaxonSelect for all flies, bees & wasps or moths & butterflies subsets,
# assign "Syrphidae" to TaxonSelect and replace animalOrder with animalFamily for hoverfly subset
TaxonSelect <- "Diptera" 
Taxon <-  kittila[animalOrder %in% TaxonSelect]

{ 
  # Prepare community matrices for single taxon
  dtMat <- Taxon[era %in% c("past", "present"),
                 .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                 by = c("plantBinomen", "animalFinestCommonTaxon")]
  com_matTax <- as.data.frame(dcast(dtMat, plantBinomen ~ animalFinestCommonTaxon,
                                    value.var = "indv", fill = 0))
  rownames(com_matTax) <- com_matTax$plantBinomen
  com_matTax <- as.matrix(com_matTax[,-1])
}

{
  # Create data frame containing variables for the community matrix
  dtTax <- Taxon[era %in% c("past", "present"),
                 .(indv = sum(animalQuantityConservativeEstimate)),
                 by = c("animalOrder", "animalFunctionalGroup" , 
                        "animalFamily", "animalFinestCommonTaxon" )]
  
  # Order taxa alphabetically
  dtTax <- dtTax[order(dtTax$animalFinestCommonTaxon),]
  
  # Calculate species-level specialization (d)
  dtSpec <- data.frame(specieslevel(com_matTax, index = "d", level = "higher"))
  
  # Add variables to data frame:
  dtSpec$animalOrder <- dtTax$animalOrder
  dtSpec$animalFamily <- dtTax$animalFamily
  dtSpec$animalFunctionalGroup <- dtTax$animalFunctionalGroup
  dtSpec$animalFinestCommonTaxon <-dtTax$animalFinestCommonTaxon
  
  # Count number of observations in past and present:
  dtEra <- Taxon[era %in% c("past", "present"),
                 .(indv = sum(animalQuantityConservativeEstimate)),
                 by = c("era", "animalOrder","animalFunctionalGroup",
                        "animalFamily", "animalFinestCommonTaxon" )]
  
  dtEra <- dtEra[order(dtEra$animalFinestCommonTaxon),]  # order pollinator taxa alphabetically
  
  No_obsPast <- data.frame(animalFinestCommonTaxon = dtEra$animalFinestCommonTaxon[dtEra$era =="past"],
                           No_obsPast = dtEra$indv[dtEra$era =="past"])
  No_obsPresent <- data.frame(animalFinestCommonTaxon = dtEra$animalFinestCommonTaxon[dtEra$era =="present"], 
                              No_obsPresent = dtEra$indv[dtEra$era =="present"])
  Obs <- merge(No_obsPast, No_obsPresent, by = "animalFinestCommonTaxon", all = TRUE)
  Obs[is.na(Obs)] <- 0
  
  # Add observations from past and present to d data frame:
  dtSpec$No_obsPast <- Obs$No_obsPast                                        # number of observations in past
  dtSpec$No_obsPresent <- Obs$No_obsPresent                                  # number of observations in present
  dtSpec$no_obsPooled <- Obs$No_obsPast + Obs$No_obsPresent                  # total number of observations
  dtSpec$relAbundPast <- Obs$No_obsPast/sum(Obs$No_obsPast)                  # relative abundance in past
  dtSpec$relAbundPresent <- Obs$No_obsPresent/sum(Obs$No_obsPresent)         # relative abundance in present
  dtSpec$changeRelAbund <- dtSpec$relAbundPresent - dtSpec$relAbundPast      # change in relative abundance
  setDT(dtSpec)
}


# Linear regression:
dtSpec <- dtSpec[dtSpec$no_obsPooled > 10] # select only species with > 10 observations total

m <- lm((dtSpec$d) ~ (dtSpec$changeRelAbund))
summary(m)

# Visualize regression:
p <- ggplotRegression(m) +
  xlab("d`") +
  ylab("Δ rel. abundance") +
  geom_hline(yintercept = 0, 
             colour = "#808080", 
             lty = 5) 

# ggsave(path = "Figures/",
#        filename = paste("Fig_5_Regression_", TaxonSelect, ".pdf", sep=""),
#        width = 9,
#        height = 5.8,
#        units = "cm",
#        dpi = 300,
#        device = "pdf")


# * Dissimilarity of interaction networks between eras ----

# assign "Diptera", Hymenoptera" or "Lepidoptera" to TaxonSelect for all flies, bees & wasps or moths & butterflies subsets,
# assign "Syrphidae" to TaxonSelect and replace animalOrder with animalFamily for hoverfly subset
TaxonSelect <- "Diptera"
Taxon <-  kittila[animalOrder %in% TaxonSelect] 

# Prepare community matrices for single taxon in the past and present
{
  TaxonPast <- Taxon[era %in% "past",
                     .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                     by = c("plantBinomen", "animalFinestCommonTaxon")]
  com_matTaxPast <- as.data.frame(dcast(TaxonPast, plantBinomen ~ animalFinestCommonTaxon,
                                        value.var = "indv", fill = 0))
  rownames(com_matTaxPast) <- com_matTaxPast$plantBinomen
  com_matTaxPast <- as.matrix(com_matTaxPast[,-1])
  
  TaxonPresent <- Taxon[era %in% "present",
                        .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                        by = c("plantBinomen", "animalFinestCommonTaxon")]
  com_matTaxPresent <- as.data.frame(dcast(TaxonPresent, plantBinomen ~ animalFinestCommonTaxon,
                                           value.var = "indv", fill = 0))
  rownames(com_matTaxPresent) <- com_matTaxPresent$plantBinomen
  com_matTaxPresent <- as.matrix(com_matTaxPresent[,-1])
}

betaDiv <- betalinkr(webs2array(list(past = com_matTaxPast,
                                     present = com_matTaxPresent)), 
                     index = "bray", # use Bray-Curtis dissimilarity index 
                     binary = FALSE, # use abundance info
                     partitioning = "commondenom",
                     partition.st = FALSE,
                     partition.rr = FALSE)
betaDiv



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Species-level metric - Pollination Service Index (PSI) --------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

{
  # Create community matrices for past and present with plant species as rows and visitor families as columns:
  FunGroupPast <- kittila[era %in% "past",
                          .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                          by = c("plantBinomen", "animalFunctionalGroup")]
 com_matFunGoupPast <- as.data.frame(dcast(FunGroupPast, plantBinomen ~ animalFunctionalGroup,
                                           value.var = "indv", fill = 0))
 rownames(com_matFunGoupPast) <- com_matFunGoupPast$plantBinomen
 com_matFunGoupPast <- as.matrix(com_matFunGoupPast[,-1])

 FunGroupPresent <- kittila[era %in% "present",
                            .(indv = sum(animalQuantityConservativeEstimate, na.rm = TRUE)),
                            by = c("plantBinomen", "animalFunctionalGroup")]
 com_matFunGoupPresent <- as.data.frame(dcast(FunGroupPresent, plantBinomen ~ animalFunctionalGroup,
                                              value.var = "indv", fill = 0))
 rownames(com_matFunGoupPresent) <- com_matFunGoupPresent$plantBinomen
 com_matFunGoupPresent <- as.matrix(com_matFunGoupPresent[,-1])
}


# Bootstrap metric
metricSL <- list()
p <- list()

metricSL[[i]] = list(Past = com_matFunGoupPast, Present = com_matFunGoupPresent) %>%
  lapply(web_matrix_to_df) %>%
  boot_specieslevel(col_lower = "lower",   # column name for plants
                    col_higher = "higher", # column name for insects
                    level = "higher",
                    index = "PSI",         # Select metric to bootstrap
                    start = 20,
                    step = round(dim(com_matFunGoupPast)[1] * 0.2),
                    n_boot = 1000,         # number of bootstraps
                    n_cpu = 3)             # number of CPU-s to use

# Visualize bootstrapped PSI
p[[i]] <- gg_specieslevel_web_by_web(metricSL[[i]],
                                     sp_higher = c("bumblebees", "muscoid flies", "moths", "hoverflies")) 

cols <- c("moths" = "#60A561", "bumblebees" = "#FFC100", 
          "hoverflies" = "#A2D6F9", "muscoid flies" = "#072AC8")

# Visualize metric
# Past
# pdf(paste("Figures/Figure_6_PSI_past.pdf"), 
#     width = 5, height = 4, colormodel = 'cmyk')

p[[i]]$Past.higher_level + 
  theme(panel.background = element_blank(), 
        axis.text.y = element_text(colour = "black", size = 12, face = "bold",
                                   angle = 90, hjust = 0.5), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        axis.line = element_line(colour = 'black', size = 0.25),
        axis.title = element_text(face = "bold", size = 14, colour = "black"),
        legend.key = element_blank(), 
        legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size = 12,face ="bold", colour ="black")) + 
  xlab("sample size")+
  ylab("PSI") +
  scale_color_manual(values = cols)

# dev.off()

# Present
# pdf(paste("Figures/Figure_6_PSI_present.pdf"), 
#     width = 5, height = 4, colormodel = 'cmyk')

p[[i]]$Present.higher_level + 
  theme(panel.background = element_blank(), 
        axis.text.y = element_text(colour = "black", size = 12, face = "bold",
                                   angle = 90, hjust = 0.5), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        axis.line = element_line(colour = 'black', size = 0.25),
        axis.title = element_text(face = "bold", size = 14, colour = "black"),
        legend.key = element_blank(), 
        legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size = 12,face ="bold", colour ="black")) + 
  xlab("sample size")+
  ylab("PSI") +
  scale_color_manual(values = cols)

# dev.off()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Sensitivity analysis ------------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Replace all "animalQuantityConservativeEstimate" in the above analyses with "animalQuantityGenerousEstimate"


