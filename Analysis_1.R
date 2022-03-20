### ANALYSIS 1 OF LEHMANN ET AL. (2018) ########################################
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

# Load required libraries.

library("Claddis")
library("data.table")

# Load functions.

source("scripts/disp_functions.R")

##### DATA IMPORT. -------------------------------------------------------------

# Set the folder where the matrices are.

folder <- "Matrices"

# Get the names of the .nex files and how many there are.

files <- dir(path = paste(getwd(), "/", folder, sep = ""), pattern = "\\.nex")
n_files <- length(files)

# Remove the .nex for the filenames.

good_names <- gsub(pattern = "\\.nex", replacement = "", x = files)
good_names <- gsub(pattern = "_", replacement = " ", x = good_names)

# Read all the matrices.

matrices <- lapply(X = paste0(folder, "/", files),
                   FUN = ReadMorphNexus)

##### List for the distances with 3 axes and for all the axes. -----------------

dataFrameList <- vector(mode = "list", length = length(matrices))

### Main for loop. -------------------------------------------------------------

timeStart <- Sys.time()

for (i in seq_along(matrices)) {
  
  ### General information about the matrix.
  
  # Number of taxa.
  
  nTax <- nrow(matrices[[i]]$matrix)
  
  # Number of characters.
  
  nChar <- ncol(matrices[[i]]$matrix)
  
  # Proportion of missing entries.
  
  generalMD <- propmd(matrices[[i]]$matrix)
  
  # Proportion of missing entries per taxon.
  
  perTaxonMD <- apply(matrices[[i]]$matrix, 1, propmd)
  
  ### Calculate the distance matrices and do the PCoA.
  
  # Calculate the distance matrices with the GED and the MORD.
  
  distMatrices <- MorphDistMatrixFast(matrices[[i]], c("GED", "Max"))
  
  # Do the PCoA for the GED matrix.
  
  gedPCoA <- L.PCoA(D = as.dist(distMatrices[[1]]),
                    correction = "lingoes",
                    do.stress = FALSE,
                    do.R2likeratio = FALSE)
  
  # Trim the MORD matrix.
  
  distMatrices[[2]] <- silent_TrimMorphDistMatrix(distMatrices[[2]])[[1]]
  
  # Do the PCoA for the MORD matrix.
  
  mordPCoA <- L.PCoA(D = as.dist(distMatrices[[2]]),
                     correction = "lingoes",
                     do.stress = FALSE,
                     do.R2likeratio = FALSE)
  
  # Extract the vectors from the PCoAs.
  
  gedVectors <- gedPCoA$results$vectors
  
  mordVectors <- mordPCoA$results$vectors
  
  # Extract the explained variance for the first 3 axes.
  
  ged3DVar <- gedPCoA$quality$cumvar[3]
  
  mord3DVar <- mordPCoA$quality$cumvar[3]
  
  ### Distance to centroid calculation.
  
  # For the GED PCoA.
  
  gedDistances3D <- distanceToCentroid(gedVectors, axes = 3)
  
  gedDistancesAll <- distanceToCentroid(gedVectors)
  
  # For the MORD PCoA
  
  mordDistances3D <- distanceToCentroid(mordVectors, axes = 3)
  
  mordDistancesAll <- distanceToCentroid(mordVectors, axes = "all")
  
  ### Create the data frame for the current matrix.
  
  # Get the number of taxa for the GED.
  
  gedTax <- nrow(distMatrices[[1]])
  
  # Get the number of taxa for the MORD.
  
  mordTax <- nrow(distMatrices[[2]])
  
  # Get the names of the taxa not trimmed.
  
  mordKeep <- names(mordDistances3D)
  
  # Create the data frame.
  
  df <- data.frame(matrix = rep(x = good_names[i],
                                times = 2 * (gedTax + mordTax)),
                   nTax = rep(x = nTax,
                              times = 2 * (gedTax + mordTax)),
                   nChar = rep(x = nChar,
                               times = 2 * (gedTax + mordTax)),
                   distance = c(rep("GED", 2 * gedTax),
                                rep("MORD", 2 * mordTax)),
                   axes = c(rep("3 axes", gedTax),
                            rep("All axes", gedTax),
                            rep("3 axes", mordTax),
                            rep("All axes", mordTax)),
                   taxon = c(names(gedDistances3D),
                             names(gedDistancesAll),
                             names(mordDistances3D),
                             names(mordDistancesAll)),
                   expvar = c(rep(ged3DVar, gedTax),
                              rep(1, gedTax),
                              rep(mord3DVar, mordTax),
                              rep(1, mordTax)),
                   propMD = c(perTaxonMD,
                              perTaxonMD,
                              perTaxonMD[mordKeep],
                              perTaxonMD[mordKeep]),
                   value = c(gedDistances3D,
                             gedDistancesAll,
                             mordDistances3D,
                             mordDistancesAll))
  
  dataFrameList[[i]] <- df
  
  timeNow <- round(difftime(Sys.time(), timeStart, units = "min"), digits = 2)
  
  message("Matrix ", i, "/", length(matrices), " completed. ",
          timeNow, " minutes elapsed.")
  
}

### Data frame assembly. -------------------------------------------------------

message("Assembling data frames...")

distances_df <- rbindlist(dataFrameList)

save(distances_df, file = "Distances_dataframe.R")

message("Data frame assembling complete.")

##### CORRELATION COEFICIENT V. 0.05 TO 1 EXPLAINED VARIANCE.-------------------

# Clean the workspace (commented for safety).

# rm(list = ls())

# Load requiered libraries.

library("Claddis")
library("data.table")

# Source auxiliary functions.

source("scripts/disp_functions.R")

# Get the filenames of the matrices.

files <- dir(path = "matrices/", pattern = "\\.nex")

# Read all the matrices.

matricesList <- lapply(X = paste0("matrices/", files),
                       FUN = ReadMorphNexus)

# Clean the matrices names.

goodNames <- gsub(pattern = "\\.nex",
                  replacement = "",
                  x = files)

goodNames <- gsub(pattern = "_",
                  replacement = " ",
                  x = goodNames)

### MAIN FUNCTION (FOR ONE MATRIX). --------------------------------------------

mainFunction <- function (morph.matrix, matrix.name) {
  
  ### Set-up the data frame.
  
  res <- data.frame(matrix = rep(x = matrix.name,
                                 times = 40),
                    
                    distance = rep(x = c("GED", "MORD"),
                                   each = 20),
                    
                    explainedVarianceStep = rep(x = seq(5, 100, 5),
                                                times = 2),
                    
                    explainedVarianceReal = rep(x = 0,
                                                times = 40),
                    
                    correlationCoefficient = rep(x = 0,
                                                 times = 40),
                    
                    significance = rep(x = NA,
                                       times = 40))
  
  ### Calculate the distance matrices and the PCoAs.
  
  # Calculate the distance matrices. [[1]] = GED; [[2]] = MORD.
  
  distanceMatrices <- MorphDistMatrixFast(morph.matrix = morph.matrix,
                                          distance = c("GED", "Max"))
  
  # Trim the MORD matrix.
  
  distanceMatrices[[2]] <- silent_TrimMorphDistMatrix(distanceMatrices[[2]])[[1]]
  
  # Trim the GED matrix accordingly.
  
  toKeep <- rownames(distanceMatrices[[2]])
  
  distanceMatrices[[1]] <- distanceMatrices[[1]][toKeep, toKeep]
  
  # Calculate the PCoAs.
  
  pcoaGED <- L.PCoA(D = as.dist(distanceMatrices[[1]]),
                    correction = "lingoes",
                    do.stress = FALSE,
                    do.R2likeratio = FALSE)
  
  pcoaMORD <- L.PCoA(D = as.dist(distanceMatrices[[2]]),
                     correction = "lingoes",
                     do.stress = FALSE,
                     do.R2likeratio = FALSE)
  
  ### Get the number of PCOs for each level of variance.
  
  # Function that obtains the number of PCOs that accumulate that variance.
  
  getPCOs <- function (PCoA.obj, value) {
    
    res <- sum(PCoA.obj$quality$cumvar < value) + 1
    
    if (res > length(PCoA.obj$quality$cumvar)) {
      
      return(length(PCoA.obj$quality$cumvar))
      
    } else {
      
      return(res)
      
    }
    
  }
  
  # Get the indices of the requiered levels of variance (GED).
  
  indicesGED <- sapply(X = seq(0.05, 1, 0.05),
                       FUN = getPCOs,
                       PCoA.obj = pcoaGED,
                       USE.NAMES = FALSE,
                       simplify = TRUE)
  
  # Get the indices of the requiered levels of variance (MORD).
  
  indicesMORD <- sapply(X = seq(0.05, 1, 0.05),
                        FUN = getPCOs,
                        PCoA.obj = pcoaMORD,
                        USE.NAMES = FALSE,
                        simplify = TRUE)
  
  # Get the real value of accumulated variance. Add to the results as %.
  
  res$explainedVarianceReal[1:20] <- pcoaGED$quality$cumvar[indicesGED] * 100
  
  res$explainedVarianceReal[21:40] <- pcoaMORD$quality$cumvar[indicesMORD] * 100
  
  ### Correlation calculation.
  
  # Get the proportion of missing entries for the kept taxa.
  
  propMissingData <- apply(X = morph.matrix$matrix[toKeep, ],
                           MARGIN = 1,
                           FUN = propmd)
  
  # Get the distance to the centroid for the taxa as a list (GED).
  
  distancesGED <- lapply(X = indicesGED,
                         FUN = distanceToCentroid,
                         vectors = pcoaGED$results$vectors)
  
  # Get the distance to the centroid for the taxa as a list (MORD).
  
  distancesMORD <- lapply(X = indicesMORD,
                          FUN = distanceToCentroid,
                          vectors = pcoaMORD$results$vectors)
  
  # Calculate the correlations (GED).
  
  correlationsGED <- lapply(X = distancesGED,
                            FUN = cor.test,
                            y = propMissingData,
                            method = "spearman")
  
  # Calculate the correlations (MORD).
  
  correlationsMORD <- lapply(X = distancesMORD,
                             FUN = cor.test,
                             y = propMissingData,
                             method = "spearman")
  
  ### Extract the results.
  
  # Extract the correlation coefficients. (4th element in the cor.test object)
  
  res$correlationCoefficient[1:20] <- sapply(X = correlationsGED,
                                             FUN = `[[`,
                                             4)
  
  res$correlationCoefficient[21:40] <- sapply(X = correlationsMORD,
                                              FUN = `[[`,
                                              4)
  
  # Extract the significance. (3rd element in the cor.test object)
  
  res$significance[1:20] <- ifelse(test = sapply(X = correlationsGED,
                                                 FUN = `[[`,
                                                 3) < 0.05,
                                   yes = TRUE,
                                   no = FALSE)
  
  res$significance[21:40] <- ifelse(test = sapply(X = correlationsMORD,
                                                  FUN = `[[`,
                                                  3) < 0.05,
                                    yes = TRUE,
                                    no = FALSE)
  
  ### Return.
  
  return(res)
  
}

### Apply mainFunction() to all the matrices. ----------------------------------

# Create a list for the results.

resList <- vector(mode = "list",
                  length = length(goodNames))

# Get the total number of matrices.

totalMatrices <- length(matricesList)

# Get the start time.

startTime <- Sys.time()

## Loop.

for (i in seq_along(matricesList)) {
  
  # Use mainFunction().
  
  resList[[i]] <- mainFunction(morph.matrix = matricesList[[i]],
                               matrix.name = goodNames[i])
  
  # Get the elapsed time.
  
  elapsedTime <- round(difftime(Sys.time(), startTime, units = "secs"), digits = 2)
  
  # Report the progress of the loop.
  
  message(i, " / ", totalMatrices, ". ", elapsedTime, " seconds elapsed.")
  
}

# Consolidate the results in a single data frame.

cumvarDataframe <- rbindlist(resList)

# Save the data frame (commented for safety).

# save(file = "Results/cumvarDataframe.R", cumvarDataframe)

##### FIGURES 1, 2, AND 3. -----------------------------------------------------


# rm(list = ls())

# Load libraries.

library("Claddis")
library("reshape2")
library("ggplot2")
library("dplyr")

# Source disparity functions.

source("scripts/disp_functions.R")

# Load the data frame.

load("Results/Distances_dataframe.R")

# Create a vector with the names of the matrices.

mtx <- levels(distances_df$matrix)

### CALCULATE THE CORRELATIONS. ------------------------------------------------

# Create a data frame with the correlations.

correlationsDataframe <- distances_df %>% 
  group_by(matrix, distance, axes) %>% 
  summarize(correlation = cor(propMD, value, method = "spearman"),
            pValue = cor.test(propMD, value, method = "spearman")$p.value)

correlationsDataframe <- correlationsDataframe %>% 
  mutate(significant = ifelse(pValue < 0.05,
                              TRUE,
                              FALSE),
         signifSign = ifelse(pValue < 0.05,
                             ifelse(correlation < 0,
                                    "sig-neg",
                                    "sig-pos"),
                             "non-sig"))

# Create a data frame with the proportion of significant results and their
# sign.

frequenciesDataframe <- correlationsDataframe %>%
  group_by(distance, axes) %>% 
  summarize(sig_neg = sum(signifSign == "sig-neg") / n(),
            non_sig = sum(signifSign == "non-sig") / n(),
            sig_pos = sum(signifSign == "sig-pos") / n())

frequenciesDataframe <- melt(frequenciesDataframe)

frequenciesDataframe$value <- round(frequenciesDataframe$value, 3) * 100

frequenciesDataframe

# Create a data frame with the explained variance for the cases with 3 axes.

expvar <- distances_df %>% 
  filter(axes == "3 axes") %>% 
  group_by(distance) %>% 
  summarize(meanExpVar = mean(expvar),
            sdExpVar = sd(expvar))

### PLOTS .---------------------------------------------------------------------

### Figure 1. ------------------------------------------------------------------

plotDataframe <- distances_df %>% 
  filter(matrix %in% c("Lamsdell (2015)",
                       "Toljagic and Butler (2013)",
                       "Baron et al. (2017)")) %>% 
  group_by(matrix, distance, axes) %>% 
  mutate(scaledValue = value / max(value))

svg(filename = "Figure 1 raw.svg", width = 12, height = 12)

ggplot(plotDataframe) +
  
  geom_point(aes(x = propMD,
                 y = scaledValue,
                 fill = distance),
             pch = 21,
             cex = 4) +
  
  xlim(c(0, 1)) +
  
  scale_fill_manual(name = "Distance measure",
                    labels = c("GED", "MORD"),
                    values = c("red2",
                               "skyblue")) +
  
  scale_colour_manual(name = "Distance measure",
                      labels = c("GED", "MORD"),
                      values = c("red2",
                                 "skyblue")) +
  
  geom_smooth(aes(x = propMD,
                  y = scaledValue,
                  colour = distance),
              method = "lm",
              se = FALSE) +
  
  facet_wrap(~matrix + axes,
             nrow = 3,
             ncol = 2) +
  
  theme_bw()

dev.off()

### Figure 2. ------------------------------------------------------------------

svg(filename = "Figure 2 raw.svg", width = 12, height = 6)

ggplot(frequenciesDataframe) + 
  geom_col(aes(x = distance,
               y = value,
               fill = variable),
           colour = "black",
           position = "dodge") +
  scale_fill_manual(name = "Correlation coefficient", 
                    labels = c("Significant & negative",
                               "Not significant",
                               "Significant & positive"),
                    values = c("black",
                               "grey95",
                               "darkgrey")) +
  facet_wrap(~axes) +
  theme_bw()

dev.off()

### Figure 3. ------------------------------------------------------------------

# Load data frames.

load("Results/cumvarDataframe.R")

# Keep what we want.

plotDataframe <- cumvarDataframe %>% 
  
  filter(matrix %in% c("Baron et al. (2017)",
                       "Toljagic and Butler (2013)",
                       "Lamsdell (2015)"))

# Clean the devices.

dev.off()

# Plot.

svg(filename = "Figure 3 raw.svg", width = 10, height = 3.5)

par(mfrow = c(1, 3))

# TOLJAGIC AND BUTLER (2013).

plot(x = NA,
     y = NA,
     xlim = c(0, 100),
     ylim = c(-1, 1),
     xlab = NA,
     ylab = NA,
     col = "white")

abline(h = 0, lty = 2)

df <- plotDataframe %>%
  filter(matrix == "Toljagic and Butler (2013)",
         distance == "GED")

points(x = df$explainedVarianceReal,
       y = df$correlationCoefficient,
       pch = ifelse(df$significance, 16, 10),
       cex = 1.5,
       col = "red2")

df <- plotDataframe %>%
  filter(matrix == "Toljagic and Butler (2013)",
         distance == "MORD")

points(x = df$explainedVarianceReal,
       y = df$correlationCoefficient,
       pch = ifelse(df$significance, 16, 10),
       cex = 1.5,
       col = "skyblue")

# LAMSDELL (2015).

plot(x = NA, y = NA,
     xlim = c(0, 100),
     ylim = c(-1, 1),
     xlab = NA,
     ylab = NA,
     col = "white")

abline(h = 0, lty = 2)

df <- plotDataframe %>%
  filter(matrix == "Lamsdell (2015)",
         distance == "GED")

points(x = df$explainedVarianceReal,
       y = df$correlationCoefficient,
       pch = ifelse(df$significance, 16, 10),
       cex = 1.5,
       col = "red2")

df <- plotDataframe %>%
  filter(matrix == "Lamsdell (2015)",
         distance == "MORD")

points(x = df$explainedVarianceReal,
       y = df$correlationCoefficient,
       pch = ifelse(df$significance, 16, 10),
       cex = 1.5,
       col = "skyblue")

# BARON ET AL. (2017).

plot(x = NA,
     y = NA,
     xlim = c(0, 100),
     ylim = c(-1, 1),
     xlab = NA,
     ylab = NA,
     col = "white")

abline(h = 0, lty = 2)

df <- plotDataframe %>% 
  filter(matrix == "Baron et al. (2017)",
         distance == "GED")

points(x = df$explainedVarianceReal,
       y = df$correlationCoefficient,
       pch = ifelse(df$significance, 16, 10),
       cex = 1.5,
       col = "red2")

df <- plotDataframe %>%
  filter(matrix == "Baron et al. (2017)",
         distance == "MORD")

points(x = df$explainedVarianceReal,
       y = df$correlationCoefficient,
       pch = ifelse(df$significance, 16, 10),
       cex = 1.5,
       col = "skyblue")

dev.off()
