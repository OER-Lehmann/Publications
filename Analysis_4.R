### ANALYSIS 4 OF LEHMANN ET AL. (2018) ########################################
### See the associated README.txt file for the description of this analysis. ###

##### SET UP. ------------------------------------------------------------------

# Clean the workspace.

# rm(list = ls())

# This script is based in version 0.2 of Claddis and the data structures and
# functions it provides. THIS SCRIPT WILL NOT WORK FOR VERSIONS > 0.2. In order
# to install version 0.2, use the following code (commented for safety).

# library("devtools")
# 
# install_github("graemetlloyd/Claddis",
#                ref = "ffa42d63d1f28eed387949b0aa247e9f2f49d6ee")

### Read libraries and functions.

# Library.

library(Claddis)

# Functions.

source(file = "scripts/disp_functions.R")

### Definition of vectors and simulation parameters.

# Get the name of folder where the matrices are.

folder <- "Group_matrices"

# Get the names of the .nex files.

files <- dir(path = paste(getwd(), "/", folder, sep = ""), pattern = "\\.nex")

# Get the clean names of the matrices.

good_names <- gsub(pattern = "\\.nex", x = files, replacement = "")
good_names <- gsub(pattern = "_", x = good_names, replacement = " ")

# Read the morphological matrices, creating of list with all of them.

matrices <- lapply(paste0(folder, "/", files), ReadMorphNexus)

# Set the number of replications.

reps <- 200

### FUNCTION DEFINITION. ------------------------------------------------------

### This is the helper function that constructs the random groups of taxa and
### drops the necessary amount of taxa in order to make the matrix divisible by
### the desired number of groups.

groupSimulation <- function (mtx, n_groups) {
  
  m <- mtx$matrix
  
  nt <- nrow(m)
  
  if ((nt %% n_groups) > 0) {
    
    m <- m[-sample(x = seq_len(nt), size = (nt %% n_groups), replace = FALSE), ]
    
    nt <- nrow(m)
    
  }
  
  taxa <- sample(x = rownames(m), size = nt, replace = FALSE)
  
  size <- nt / n_groups
  
  groups <- split(x = taxa, f = ceiling(seq_along(taxa) / size))
  
  names(groups) <- paste0("g", seq_len(n_groups))
  
  mtx$matrix <- m
  
  out <- list(matrix = mtx, groups = groups)
  
  return(out)
  
}

### This is just the function of Claddis but with no printing when there was no
### trimming of the distance matrix. It just crammed the console and I could not
### clearly see the replication and the elapsed time.

silent_TrimMorphDistMatrix <- function (dist.matrix, tree = NULL) {
  
  if (class(dist.matrix) == "dist") 
    dist.matrix <- as.matrix(dist.matrix)
  if (!is.matrix(dist.matrix)) 
    stop("ERROR: Input must be a distance matrix (i.e., either an object of class \"dist\" or a square matrix).")
  if (is.null(tree)) {
    if (length(which(is.na(dist.matrix))) == 0) {
      removed.taxa <- NULL
      out <- list(dist.matrix, tree, removed.taxa)
      names(out) <- c("dist.matrix", "tree", "removed.taxa")
      return(out)
    }
    else {
      removes <- vector(mode = "character")
      while (length(which(is.na(dist.matrix))) > 0) {
        na.lengths <- vector(mode = "numeric")
        for (i in 1:length(dist.matrix[, 1])) na.lengths[i] <- length(which(is.na(dist.matrix[i, 
                                                                                              ])))
        taxon.to.delete <- rownames(dist.matrix)[which(na.lengths == 
                                                         max(na.lengths))[1]]
        delete.row <- which(rownames(dist.matrix) == 
                              taxon.to.delete)
        dist.matrix <- dist.matrix[-delete.row, -delete.row]
        removes <- c(removes, taxon.to.delete)
      }
      out <- list(dist.matrix, tree, removes)
      names(out) <- c("dist.matrix", "tree", "removed.taxa")
      return(out)
    }
  }
  else {
    if (length(which(is.na(dist.matrix))) == 0) {
      print("There are no gaps in the distance matrix")
      removed.taxa <- NULL
      out <- list(dist.matrix, tree, removed.taxa)
      names(out) <- c("dist.matrix", "tree", "removed.taxa")
      return(out)
    }
    else {
      node.nos <- as.character((Ntip(tree) + 1):(Ntip(tree) + 
                                                   Nnode(tree)))
      for (i in match(node.nos, rownames(dist.matrix))) colnames(dist.matrix)[i] <- rownames(dist.matrix)[i] <- paste(sort(tree$tip.label[FindDescendants(rownames(dist.matrix)[i], 
                                                                                                                                                          tree)]), collapse = "%%")
      removes <- vector(mode = "character")
      temp.dist.matrix <- dist.matrix
      while (length(which(is.na(dist.matrix))) > 0) {
        na.lengths <- vector(mode = "numeric")
        for (i in 1:length(dist.matrix[, 1])) na.lengths[i] <- length(which(is.na(dist.matrix[i, 
                                                                                              ])))
        taxon.to.delete <- rownames(dist.matrix)[which(na.lengths == 
                                                         max(na.lengths))[1]]
        delete.row <- which(rownames(dist.matrix) == 
                              taxon.to.delete)
        dist.matrix <- dist.matrix[-delete.row, -delete.row]
        removes <- c(removes, taxon.to.delete)
      }
      dist.matrix <- temp.dist.matrix
      tips.to.remove <- removes[sort(match(tree$tip.label, 
                                           removes))]
      nodes.to.remove <- removes[grep("%%", removes)]
      tip.name.nodes.to.remove <- vector(mode = "numeric")
      for (i in tips.to.remove) {
        originating.node <- tree$edge[match(which(tree$tip.label == 
                                                    i), tree$edge[, 2]), 1]
        originating.node.name <- paste(sort(tree$tip.label[FindDescendants(originating.node, 
                                                                           tree)]), collapse = "%%")
        tip.name.nodes.to.remove <- unique(c(tip.name.nodes.to.remove, 
                                             originating.node.name))
      }
      if (length(setdiff(nodes.to.remove, tip.name.nodes.to.remove)) > 
          0) {
        nodes.left.to.remove <- setdiff(nodes.to.remove, 
                                        tip.name.nodes.to.remove)
        for (i in nodes.left.to.remove) {
          originating.node <- FindAncestor(strsplit(i, 
                                                    "%%")[[1]], tree)
          descendant.nodes <- tree$edge[which(tree$edge[, 
                                                        1] == originating.node), 2]
          if (length(which(descendant.nodes <= Ntip(tree))) > 
              0) {
            taxon.to.exclude <- tree$tip.label[min(descendant.nodes)]
          }
          else {
            if (length(FindDescendants(descendant.nodes[1], 
                                       tree)) >= length(FindDescendants(descendant.nodes[2], 
                                                                        tree))) {
              taxon.to.exclude <- tree$tip.label[FindDescendants(descendant.nodes[2], 
                                                                 tree)]
            }
            else {
              taxon.to.exclude <- tree$tip.label[FindDescendants(descendant.nodes[1], 
                                                                 tree)]
            }
          }
          removes <- c(removes, taxon.to.exclude)
          tips.to.remove <- c(tips.to.remove, taxon.to.exclude)
        }
      }
      dist.matrix <- dist.matrix[-match(removes, rownames(dist.matrix)), 
                                 -match(removes, rownames(dist.matrix))]
      for (i in 1:length(dist.matrix[, 1])) {
        nms <- strsplit(rownames(dist.matrix)[i], "%%")[[1]]
        rownames(dist.matrix)[i] <- colnames(dist.matrix)[i] <- paste(sort(nms[which(is.na(match(nms, 
                                                                                                 tips.to.remove)))]), collapse = "%%")
      }
      redundant.rows <- c(which(duplicated(rownames(dist.matrix))), 
                          which(rownames(dist.matrix) == ""))
      if (length(redundant.rows) > 0) {
        dist.matrix <- dist.matrix[-redundant.rows, -redundant.rows]
      }
      tree <- drop.tip(tree, tips.to.remove)
      node.names <- rownames(dist.matrix)[grep("%%", rownames(dist.matrix))]
      for (j in 1:length(node.names)) names(node.names)[j] <- FindAncestor(strsplit(node.names[j], 
                                                                                    "%%")[[1]], tree)
      colnames(dist.matrix)[match(node.names, colnames(dist.matrix))] <- rownames(dist.matrix)[match(node.names, 
                                                                                                     rownames(dist.matrix))] <- names(node.names)
      out <- list(dist.matrix, tree, removes)
      names(out) <- c("dist.matrix", "tree", "removed.taxa")
      return(out)
    }
  }
}

### This is the core function for this script. It encompasses a whole simulation
### incluiding all the treatments. It is not fully optimized, and it is quite
### redundant, but it works.

rankListing <- function(grouped_mtx) {
  
  ### Construct the first groups.
  
  m <- grouped_mtx$matrix
  
  groups <- grouped_mtx$groups
  
  ### Get the base results
  
  df <- data.frame(group = paste0("g", 1:5),
                   GED_MPD = rep(0, 5),
                   GED_WMPD = rep(0, 5),
                   GED_SoV = rep(0, 5),
                   GED_SoR = rep(0, 5),
                   MORD_MPD = rep(0, 5),
                   MORD_WMPD = rep(0, 5),
                   MORD_SoV = rep(0, 5),
                   MORD_SoR = rep(0, 5))
  
  for (i in seq_along(groups)) {
    
    df[i, 2:9] <- unlist(stdPipe(m, g = groups[[i]]), use.names = FALSE)
    
  }
  
  # Rank the groups according to their 
  
  rank_list <- list(GED_MPD = groups[sort(df$GED_MPD, index.return = TRUE)$ix],
                    GED_WMPD = groups[sort(df$GED_WMPD, index.return = TRUE)$ix],
                    GED_SoV = groups[sort(df$GED_SoV, index.return = TRUE)$ix],
                    GED_SoR = groups[sort(df$GED_SoR, index.return = TRUE)$ix],
                    MORD_MPD = groups[sort(df$MORD_MPD, index.return = TRUE)$ix],
                    MORD_WMPD = groups[sort(df$MORD_WMPD, index.return = TRUE)$ix],
                    MORD_SoV = groups[sort(df$MORD_SoV, index.return = TRUE)$ix],
                    MORD_SoR = groups[sort(df$MORD_SoR, index.return = TRUE)$ix])
  
  # Uncorrelate the groups.
  
  for (i in seq_along(rank_list)) {
    
    rank_list[[i]] <- rank_list[[i]][c(2, 5, 3, 1, 4)]
    
  }
  
  return(rank_list)
  
}

### .-----

rankPipe <- function(grouped_mtx, treatment, rank_list) {
  
  out <- list(GED_MPD = NA,
              GED_WMPD = NA,
              GED_SoV = NA,
              GED_SoR = NA,
              MORD_MPD = NA,
              MORD_WMPD = NA,
              MORD_SoV = NA,
              MORD_SoR = NA)
  
  m <- missify(grouped_mtx$matrix, treatment, rank_list[[1]])
  
  dms <- MorphDistMatrixFast(m, "GED")
  
  out$GED_MPD <- c(mpd(dms[rank_list[[1]][[1]],
                           rank_list[[1]][[1]]]),
                   mpd(dms[rank_list[[1]][[2]],
                           rank_list[[1]][[2]]]),
                   mpd(dms[rank_list[[1]][[3]],
                           rank_list[[1]][[3]]]),
                   mpd(dms[rank_list[[1]][[4]],
                           rank_list[[1]][[4]]]),
                   mpd(dms[rank_list[[1]][[5]],
                           rank_list[[1]][[5]]]))
  
  m <- missify(grouped_mtx$matrix, treatment, rank_list[[2]])
  
  dms <- MorphDistMatrixFast(m, c("GED", "Comp"))
  
  out$GED_WMPD <- c(wmpd(dms[[1]][rank_list[[2]][[1]],
                                  rank_list[[2]][[1]]],
                         dms[[2]][rank_list[[2]][[1]],
                                  rank_list[[2]][[1]]]),
                    wmpd(dms[[1]][rank_list[[2]][[2]],
                                  rank_list[[2]][[2]]],
                         dms[[2]][rank_list[[2]][[2]],
                                  rank_list[[2]][[2]]]),
                    wmpd(dms[[1]][rank_list[[2]][[3]],
                                  rank_list[[2]][[3]]],
                         dms[[2]][rank_list[[2]][[3]],
                                  rank_list[[2]][[3]]]),
                    wmpd(dms[[1]][rank_list[[2]][[4]],
                                  rank_list[[2]][[4]]],
                         dms[[2]][rank_list[[2]][[4]],
                                  rank_list[[2]][[4]]]),
                    wmpd(dms[[1]][rank_list[[2]][[5]],
                                  rank_list[[2]][[5]]],
                         dms[[2]][rank_list[[2]][[5]],
                                  rank_list[[2]][[5]]]))
  
  m <- missify(grouped_mtx$matrix, treatment, rank_list[[5]])
  
  dms <- MorphDistMatrixFast(m, "Max")
  
  out$MORD_MPD <- c(mpd(dms[rank_list[[5]][[1]], rank_list[[5]][[1]]]),
                    mpd(dms[rank_list[[5]][[2]], rank_list[[5]][[2]]]),
                    mpd(dms[rank_list[[5]][[3]], rank_list[[5]][[3]]]),
                    mpd(dms[rank_list[[5]][[4]], rank_list[[5]][[4]]]),
                    mpd(dms[rank_list[[5]][[5]], rank_list[[5]][[5]]]))
  
  m <- missify(grouped_mtx$matrix, treatment, rank_list[[6]])
  
  dms <- MorphDistMatrixFast(m, c("Max", "Comp"))
  
  out$MORD_WMPD <- c(wmpd(dms[[1]][rank_list[[6]][[1]], rank_list[[6]][[1]]],
                          dms[[2]][rank_list[[6]][[1]], rank_list[[6]][[1]]]),
                     wmpd(dms[[1]][rank_list[[6]][[2]], rank_list[[6]][[2]]],
                          dms[[2]][rank_list[[6]][[2]], rank_list[[6]][[2]]]),
                     wmpd(dms[[1]][rank_list[[6]][[3]], rank_list[[6]][[3]]],
                          dms[[2]][rank_list[[6]][[3]], rank_list[[6]][[3]]]),
                     wmpd(dms[[1]][rank_list[[6]][[4]], rank_list[[6]][[4]]],
                          dms[[2]][rank_list[[6]][[4]], rank_list[[6]][[4]]]),
                     wmpd(dms[[1]][rank_list[[6]][[5]], rank_list[[6]][[5]]],
                          dms[[2]][rank_list[[6]][[5]], rank_list[[6]][[5]]]))
  
  ### Post-ordination measures.
  
  m <- missify(grouped_mtx$matrix, treatment, rank_list[[3]])
  
  dms <- MorphDistMatrixFast(m, "GED")
  
  gp <- L.PCoA(D = as.dist(dms),
               correction = "lingoes",
               do.stress = FALSE,
               do.R2likeratio = FALSE)$results$vectors[, 1:2]
  
  out$GED_SoV <- c(sum_of_variance(gp[rank_list[[3]][[1]], ], axes = "all"),
                   sum_of_variance(gp[rank_list[[3]][[2]], ], axes = "all"),
                   sum_of_variance(gp[rank_list[[3]][[3]], ], axes = "all"),
                   sum_of_variance(gp[rank_list[[3]][[4]], ], axes = "all"),
                   sum_of_variance(gp[rank_list[[3]][[5]], ], axes = "all"))
  
  m <- missify(grouped_mtx$matrix, treatment, rank_list[[4]])
  
  dms <- MorphDistMatrixFast(m, "GED")
  
  gp <- L.PCoA(D = as.dist(dms),
               correction = "lingoes",
               do.stress = FALSE,
               do.R2likeratio = FALSE)$results$vectors[, 1:2]
  
  out$GED_SoR <- c(sum_of_ranges(gp[rank_list[[4]][[1]], ], axes = "all"),
                   sum_of_ranges(gp[rank_list[[4]][[2]], ], axes = "all"),
                   sum_of_ranges(gp[rank_list[[4]][[3]], ], axes = "all"),
                   sum_of_ranges(gp[rank_list[[4]][[4]], ], axes = "all"),
                   sum_of_ranges(gp[rank_list[[4]][[5]], ], axes = "all"))
  
  m <- missify(grouped_mtx$matrix, treatment, rank_list[[7]])
  
  dms <- silent_TrimMorphDistMatrix(MorphDistMatrixFast(m, "Max"))[[1]]
  
  gp <- L.PCoA(D = as.dist(dms),
               correction = "lingoes",
               do.stress = FALSE,
               do.R2likeratio = FALSE)$results$vectors[, 1:2]
  
  tax <- rownames(dms)
  
  out$MORD_SoV <- c(sum_of_variance(gp[intersect(rank_list[[7]][[1]], tax), ], axes = "all"),
                    sum_of_variance(gp[intersect(rank_list[[7]][[2]], tax), ], axes = "all"),
                    sum_of_variance(gp[intersect(rank_list[[7]][[3]], tax), ], axes = "all"),
                    sum_of_variance(gp[intersect(rank_list[[7]][[4]], tax), ], axes = "all"),
                    sum_of_variance(gp[intersect(rank_list[[7]][[5]], tax), ], axes = "all"))
  
  m <- missify(grouped_mtx$matrix, treatment, rank_list[[8]])
  
  dms <- silent_TrimMorphDistMatrix(MorphDistMatrixFast(m, "Max"))[[1]]
  
  tax <- rownames(dms)
  
  gp <- L.PCoA(D = as.dist(dms),
               correction = "lingoes",
               do.stress = FALSE,
               do.R2likeratio = FALSE)$results$vectors[, 1:2]
  
  out$MORD_SoR <- c(sum_of_ranges(gp[intersect(rank_list[[8]][[1]], tax), ], axes = "all"),
                    sum_of_ranges(gp[intersect(rank_list[[8]][[2]], tax), ], axes = "all"),
                    sum_of_ranges(gp[intersect(rank_list[[8]][[3]], tax), ], axes = "all"),
                    sum_of_ranges(gp[intersect(rank_list[[8]][[4]], tax), ], axes = "all"),
                    sum_of_ranges(gp[intersect(rank_list[[8]][[5]], tax), ], axes = "all"))
  
  return(out)
  
}

treatments <- matrix(data = c(0.15, 0.125, 0.1, 0.075, 0.05,
                              0.2, 0.15, 0.1, 0.05, 0,
                              0.25, 0.225, 0.2, 0.175, 0.15,
                              0.3, 0.25, 0.2, 0.15, 0.1,
                              0.35, 0.325, 0.3, 0.275, 0.25,
                              0.4, 0.35, 0.3, 0.25, 0.2,
                              0.5, 0.4, 0.3, 0.2, 0.1,
                              rep(0.1, 5),
                              rep(0.2, 5),
                              rep(0.3, 5),
                              rep(0.4, 5),
                              rep(0.5, 5),
                              rep(0.6, 5)),
                     nrow = 13,
                     ncol = 5,
                     byrow = TRUE,
                     dimnames = list(Treatment = c("dot1_low",
                                                   "dot1_high",
                                                   "dot2_low",
                                                   "dot2_high",
                                                   "dot3_low",
                                                   "dot3_med",
                                                   "dot3_high",
                                                   paste0("dot", 1:6, "_homo")),
                                     Group = paste0("Group ", 1:5)))

df_list <- vector("list", length(files) * reps)

index <- 1

time_0 <- Sys.time()

for (i in seq_along(files)) {
  
  x <- ReadMorphNexus(paste0(folder, "/",files[i]))
  
  gx <- groupSimulation(x, 5)
  
  ranks <- rankListing(gx)
  
  for (j in seq_len(reps)) {
    
    output <- list()
    
    for (k in seq_len(nrow(treatments))) {
      
      output <- list(output, rankPipe(gx, treatments[k, ], ranks))
      
      time_out <- as.numeric(round(difftime(Sys.time(),
                                            time_0,
                                            units = "min"),
                                   digits = 2))
      
      message("Matrix ", i, "/", length(files), ". Rep ", j, "/", reps,
              ". Treatment ", k, "/", nrow(treatments), ". ",
              time_out, " minutes elapsed.")
      
    }
    
    res <- unlist(output, use.names = FALSE)
    
    treats <- rownames(treatments)
    
    treats <- factor(x = treats,
                     levels = c("dot1_homo",
                                "dot1_low",
                                "dot1_high",
                                "dot2_homo",
                                "dot2_low",
                                "dot2_high",
                                "dot3_homo",
                                "dot3_low",
                                "dot3_med",
                                "dot3_high",
                                "dot4_homo",
                                "dot5_homo",
                                "dot6_homo"),
                     labels = c("0.1 - Homogeneous",
                                "0.1 - Low dispersion",
                                "0.1 - High dispersion",
                                "0.2 - Homogeneous",
                                "0.2 - Low dispersion",
                                "0.2 - High dispersion",
                                "0.3 - Homogeneous",
                                "0.3 - Low dispersion",
                                "0.3 - Medium dispersion",
                                "0.3 - High dispersion",
                                "0.4 - Homogeneous",
                                "0.5 - Homogeneous",
                                "0.6 - Homogeneous"))
    
    df <- data.frame(matrix = good_names[[i]],
                     rep = rep(paste0("rep", j), 520),
                     treatment = rep(treats, each = 40),
                     measure = rep(rep(c("MPD", "WMPD", "SoV", "SoR"), each = 5), 26),
                     distance = rep(rep(c(rep("GED", 4), rep("MORD", 4)), each = 5), 13),
                     group = rep(c("2", "5", "3", "1", "4"), 104),
                     values = res)
    
    df_list[[index]] <- df
    
    index <- index + 1
    
  }
  
}

save(file = "", df_list)

##### FIGURE. -----

# setwd("")
# 
# load("raw_groups_SetI.R")
# 
# groupsDataframe_I <- data.table::rbindlist(df_list)
# 
# load("raw_groups_SetII.R")
# 
# groupsDataframe_II <- data.table::rbindlist(df_list)
# 
# load("raw_groups_SetIII.R")
# 
# groupsDataframe_III <- data.table::rbindlist(df_list)
# 
# load("raw_groups_SetIV.R")
# 
# groupsDataframe_IV <- data.table::rbindlist(df_list)
# 
# groupsDataframe_all <- rbind.data.frame(groupsDataframe_I,
#                                         groupsDataframe_II,
#                                         groupsDataframe_III,
#                                         groupsDataframe_IV)
# 
# save(groupsDataframe_all, file = "groupsDataframe_all.R")

### SET-UP. --------------------------------------------------------------------

# Clean the workspace (commented for safety).

# rm(list = ls())

# Load the data.

load("Results/groupsDataframe_all.R")

# Load the required libraries.

library(dplyr)
library(ggplot2)

# Get the names of the matrices and reorder them alphabetically.

matricesNames <- sort(levels(groupsDataframe_all$matrix))

# Relabel the treatments.

groupsDataframe_all$treatment <- factor(x = groupsDataframe_all$treatment,
                                        
                                        levels = c("0.1 - Homogeneous",
                                                   "0.1 - Low dispersion",
                                                   "0.1 - High dispersion",
                                                   "0.2 - Homogeneous",
                                                   "0.2 - Low dispersion",
                                                   "0.2 - High dispersion",
                                                   "0.3 - Homogeneous",
                                                   "0.3 - Low dispersion",
                                                   "0.3 - Medium dispersion",
                                                   "0.3 - High dispersion",
                                                   "0.4 - Homogeneous",
                                                   "0.5 - Homogeneous",
                                                   "0.6 - Homogeneous"),
                                        
                                        labels = c("10% - Homogeneous",
                                                   "10% - Low dispersion",
                                                   "10% - High dispersion",
                                                   "20% - Homogeneous",
                                                   "20% - Low dispersion",
                                                   "20% - High dispersion",
                                                   "30% - Homogeneous",
                                                   "30% - Low dispersion",
                                                   "30% - Medium dispersion",
                                                   "30% - High dispersion",
                                                   "40% - Homogeneous",
                                                   "50% - Homogeneous",
                                                   "60% - Homogeneous"))

# Relabel the disparity measures.

groupsDataframe_all$measure <- factor(x = groupsDataframe_all$measure,
                                      
                                      levels = c("MPD",
                                                 "WMPD",
                                                 "SoV",
                                                 "SoR"),
                                      
                                      labels = c("MPD",
                                                 "WMPD",
                                                 "Sum of variances",
                                                 "Sum of ranges"))

# Create the data frame with the Spearman's correlations.

correlationDataframe <- groupsDataframe_all %>% 
  
  group_by(matrix, rep, treatment, measure, distance) %>% 
  
  summarize(correlation = cor(values, 1:5, method = "spearman"))

# Create a data frame for the figure rectangles.

rectDataframe <- data.frame(xmin = rep(-1, 13),
                            xmax = rep(1, 13),
                            ymin = seq(0.5, 12.5, 1),
                            ymax = seq(1.5, 13.5, 1))

# Create the bootstrap function.

meanBootstrap <- function (x, reps = 999, conf_level = 0.95, seed = NULL) {
  
  set.seed(seed)
  
  sample_mean <- mean(x)
  
  resamples <- replicate(n = reps,
                         expr = sample(x, length(x), TRUE))
  
  delta_stars <- apply(X = resamples,
                       MARGIN = 2,
                       FUN = function (x) {mean(x) - sample_mean})
  
  res <- sort(sample_mean - quantile(delta_stars,
                                     c((1 - conf_level) / 2,
                                       1 - ((1 - conf_level) / 2)),
                                     names = FALSE))
  
  set.seed(NULL)
  
  return(list(mean = sample_mean, interval = res))
  
}

# Create a variable with random numbers to act as seeds for the bootsrap.

correlationDataframe <- correlationDataframe %>% 
  mutate(randomSeed = floor(runif(n = n(), min = 1, max = 200000)))

# Summarize the correlation information..

correlationDataframe <- correlationDataframe %>%
  
  ungroup() %>% 
  
  group_by(matrix, treatment, distance, measure) %>% 
  
  summarize(meanCorrelation = mean(correlation),
            ciLower = meanBootstrap(correlation,
                                    seed = randomSeed[1])$interval[1],
            ciUpper = meanBootstrap(correlation,
                                    seed = randomSeed[1])$interval[2]) %>% 
  
  ungroup()

# Create a list for the plots.

gList <- vector(mode = "list", length = length(matricesNames))

# Fill the list with the plots.

for (i in seq_along(matricesNames)) {
  
  plotDataframe <- correlationDataframe %>% 
    filter(matrix == matricesNames[i])
  
  gList[[i]] <- ggplot(plotDataframe) +
    
    geom_errorbarh(aes(x = meanCorrelation, 
                       y = treatment,
                       colour = distance,
                       xmin = ciLower,
                       xmax = ciUpper), 
                   size = 0.55,
                   height = 0.35,
                   position = position_nudge(x = 0,
                                             y = c(0.25, -0.25))) +
    
    geom_vline(aes(xintercept = 0), lty = 2) +
    
    geom_point(aes(x = meanCorrelation,
                   y = treatment,
                   colour = distance),
               position = position_nudge(x = 0,
                                         y = c(0.25, -0.25)),
               cex = 1.9) +
    
    ggtitle(label = matricesNames[i]) +
    
    xlab("Correlation coefficient") +
    
    ylab(NULL) +
    
    xlim(c(-1, 1)) +
    
    geom_rect(data = rectDataframe,
              aes(xmin = xmin,
                  xmax = xmax,
                  ymin = ymin,
                  ymax = ymax),
              fill = NA,
              colour = "black") +
    
    facet_wrap(~measure) +
    
    theme_bw() +
    
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank())
  
}

pdf(file = "Groups supplementary figures.pdf",
    width = 10,
    height = 7,
    paper = "a4r")

for (i in seq_along(gList)) {
  
  plot(gList[[i]])
  
  message(i, " / ", length(gList))
  
}

dev.off()
