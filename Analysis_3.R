### ANALYSIS 3 OF LEHMANN ET AL. (2018) ########################################
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

# Load the required libraries.

library(Claddis)
library(crayon)

# Source the requiered functions.

source("scripts/disp_functions.R")

# Get the matrices files.

files <- dir(path = "Complete_matrices/", pattern = "\\.nex")

# Read the matrices.

matricesList <- lapply(X = paste0("Complete_matrices/", files),
                       FUN = ReadMorphNexus)

# Clean the names of the matrices.

goodNames <- gsub(pattern = "\\.nex",
                  replacement = "",
                  x = files)

goodNames <- gsub(pattern = "_",
                  replacement = " ",
                  x = goodNames)

### mainFunction(). -------------------------------------------------------------

mainFunction <- function (morph.matrix, matrix.name, missing.data, n.reps) {
  
  ### Create the data frame that will hold the results. ########################
  
  # Get the number of missing data levels.
  
  nMD <- length(missing.data)
  
  # Create the data frame.
  
  res <- data.frame(matrix = rep(matrix.name,
                                 times = 8 + (8 * nMD * n.reps)),
                    
                    missingData = factor(c(rep(0,
                                               times = 8),
                                           rep(missing.data,
                                               each = 8 * n.reps))),
                    
                    replication = factor(c(rep("rep1", times = 8),
                                           rep(rep(paste0("rep",
                                                          seq_len(n.reps)),
                                                   each = 8),
                                               times = nMD))),
                    
                    distance = factor(rep(rep(c("GED", "MORD"),
                                              each = 4),
                                          times = 1 + (nMD * n.reps))),
                    
                    measure = factor(rep(c("MPD", "WMPD", "SoV", "SoR"),
                                         times = 2 + (2 * n.reps * nMD))),
                    
                    value = rep(0,
                                times = 8 + (8 * nMD * n.reps)))
  
  # Create a list for the results.
  
  resList <- vector("list", length = 1 + (nMD * n.reps))
  
  # Set an index for its filling (starts with 2 to skip the base-case).
  
  index <- 2
  
  ### Base-case resolution. ####################################################
  
  # Create a vector for the results.
  
  values <- rep(0, 8)
  
  # Calculate the distance matrices and the comparable character matrix.
  
  distanceMatrices <- MorphDistMatrixFast(morph.matrix = morph.matrix,
                                          distance = c("GED", "Max", "Comp"))
  
  ## Pre-ordination measures.
  
  # Calculate the MPD (GED).
  
  values[1] <- mpd(D = distanceMatrices[[1]])
  
  # Calculate the WMPD (GED).
  
  values[2] <- wmpd(D = distanceMatrices[[1]],
                    comp_chars = distanceMatrices[[3]])
  
  # Calculate the MPD (MORD).
  
  values[5] <- mpd(D = distanceMatrices[[2]])
  
  # Calculate the WMPD (MORD).
  
  values[6] <- wmpd(D = distanceMatrices[[2]],
                    comp_chars = distanceMatrices[[3]])
  
  ## Post-ordination measures.
  
  # Trim the MORD matrix.
  
  distanceMatrices[[2]] <- silent_TrimMorphDistMatrix(distanceMatrices[[2]])[[1]]
  
  # Get the taxa to keep.
  
  toKeep <- rownames(distanceMatrices[[2]])
  
  # Trim the GED matrix.
  
  distanceMatrices[[1]] <- distanceMatrices[[1]][toKeep, toKeep]
  
  # Calculate the PCoAs.
  
  pcoaGED <- L.PCoA(D = as.dist(distanceMatrices[[1]]),
                    correction = "lingoes",
                    do.stress = FALSE,
                    do.R2likeratio = FALSE)$results$vectors
  
  pcoaMORD <- L.PCoA(D = as.dist(distanceMatrices[[2]]),
                     correction = "lingoes",
                     do.stress = FALSE,
                     do.R2likeratio = FALSE)$results$vectors
  
  # Calculate the sum of variances (GED).
  
  values[3] <- sum_of_variance(vectors = pcoaGED, axes = "all")
  
  # Calculate the sum of ranges (GED).
  
  values[4] <- sum_of_ranges(vectors = pcoaGED, axes = "all")
  
  # Calculate the sum of variances (MORD).
  
  values[7] <- sum_of_variance(vectors = pcoaMORD, axes = "all")
  
  # Calculate the sum of ranges (MORD).
  
  values[8] <- sum_of_ranges(vectors = pcoaMORD, axes = "all")
  
  ## Add the results to the list.
  
  resList[[1]] <- values
  
  message("From function: Base-case complete.")
  
  ### Missify and replications. ################################################
  
  for (i in seq_along(missing.data)) {
    
    for (j in seq_len(n.reps)) {
      
      # Create a vector for the values.
      
      values <- rep(0, 8)
      
      # Missify the matrix.
      
      missifiedMatrix <- missify(morph.matrix, missing.data[i])
      
      # Calculate the distance matrices and the comparable character matrix.
      
      distanceMatrices <- MorphDistMatrixFast(morph.matrix = missifiedMatrix,
                                              distance = c("GED", "Max", "Comp"))
      
      ## Pre-ordination measures.
      
      # Calculate the MPD (GED).
      
      values[1] <- mpd(D = distanceMatrices[[1]])
      
      # Calculate the WMPD (GED).
      
      values[2] <- wmpd(D = distanceMatrices[[1]],
                        comp_chars = distanceMatrices[[3]])
      
      # Calculate the MPD (MORD).
      
      values[5] <- mpd(D = distanceMatrices[[2]])
      
      # Calculate the WMPD (MORD).
      
      values[6] <- wmpd(D = distanceMatrices[[2]],
                        comp_chars = distanceMatrices[[3]])
      
      ## Post-ordination measures.
      
      # Trim the MORD matrix.
      
      distanceMatrices[[2]] <- silent_TrimMorphDistMatrix(distanceMatrices[[2]])[[1]]
      
      # Get the taxa to keep.
      
      toKeep <- rownames(distanceMatrices[[2]])
      
      # Trim the GED matrix.
      
      distanceMatrices[[1]] <- distanceMatrices[[1]][toKeep, toKeep]
      
      # Calculate the PCoAs.
      
      pcoaGED <- L.PCoA(D = as.dist(distanceMatrices[[1]]),
                        correction = "lingoes",
                        do.stress = FALSE,
                        do.R2likeratio = FALSE)$results$vectors
      
      pcoaMORD <- L.PCoA(D = as.dist(distanceMatrices[[2]]),
                         correction = "lingoes",
                         do.stress = FALSE,
                         do.R2likeratio = FALSE)$results$vectors
      
      # Calculate the sum of variances (GED).
      
      values[3] <- sum_of_variance(vectors = pcoaGED, axes = "all")
      
      # Calculate the sum of ranges (GED).
      
      values[4] <- sum_of_ranges(vectors = pcoaGED, axes = "all")
      
      # Calculate the sum of variances (MORD).
      
      values[7] <- sum_of_variance(vectors = pcoaMORD, axes = "all")
      
      # Calculate the sum of ranges (MORD).
      
      values[8] <- sum_of_ranges(vectors = pcoaMORD, axes = "all")
      
      ## Add the results to the list.
      
      # Add to the list.
      
      resList[[index]] <- values
      
      message("From function: step ", index - 1, " of ", length(resList) - 1,
              " completed.")
      
      # Increase the index.
      
      index <- index + 1
      
    }
    
  }
  
  ### Complete the data frame. #################################################
  
  # Complete the data frame.
  
  res$value <- unlist(x = resList, use.names = FALSE)
  
  # Return the data frame.
  
  return(res)
  
}

### APPLY THE MAIN FUNCTION TO ALL THE MATRICES. ###############################

resultsList <- vector(mode = "list",
                      length = length(matricesList))

startTime <- Sys.time()

for (i in seq_along(matricesList)) {
  
  resultsList[[i]] <- mainFunction(morph.matrix = matricesList[[i]],
                                   matrix.name = goodNames[i],
                                   missing.data = seq(0.1, 0.6, 0.1),
                                   n.reps = 200)
  
  elapsedTime <- round(difftime(Sys.time(), startTime, units = "mins"),
                       digits = 2)
  
  message(green(i, " / ", length(matricesList), ". ",
                elapsedTime, " minutes elapsed."))
  
}

save(file = "homoDispMeasuresData.R", resultsList)

##### FIGURES 5 and SI 9. ------------------------------------------------------

# Clear the workspace (commented for safety).

# rm(list = ls())

# Load required libraries.

library("dplyr")
library("ggplot2")
library("data.table")

# Load data.

load("Results/homoDispMeasuresData.R")

# Create the data frame.

homoDispMeasures <- data.table::rbindlist(resultsList)

# Get the names of the matrices.

matricesNames <- levels(homoDispMeasures$matrix)

# Order the disparity measures.

homoDispMeasures$measure <- factor(x = homoDispMeasures$measure,
                                   levels = c("MPD", "WMPD", "SoV", "SoR"),
                                   labels = c("MPD",
                                              "WMPD",
                                              "Sum of variances",
                                              "Sum of ranges"))

# Relabel the proportion of missing data as percentages.

homoDispMeasures$missingData <- factor(x = homoDispMeasures$missingData,
                                       levels = levels(homoDispMeasures$missingData),
                                       labels = c("0", as.character(seq(10, 60, 10))))

# Calculate the mean and the confidence interval for each variable.

processedDataframe <- homoDispMeasures %>% 
  
  group_by(matrix, missingData, distance, measure) %>% 
  
  summarize(mean = mean(value),
            
            lci = {ifelse(n() == 1,
                          NA,
                          t.test(value)$conf.int[1])},
            
            uci = {ifelse(n() == 1,
                          NA,
                          t.test(value)$conf.int[2])}) %>% 
  
  ungroup() %>% 
  
  group_by(matrix, distance, measure) %>% 
  
  mutate(scaledMean = mean / max(c(mean, lci, uci), na.rm = TRUE),
         scaledLCI = lci / max(c(mean, lci, uci), na.rm = TRUE),
         scaledUCI = uci / max(c(mean, lci, uci), na.rm = TRUE))

### FIGURE SI 9. ---------------------------------------------------------------

gList <- vector(mode = "list", length = length(matricesNames))

for (i in seq_along(matricesNames)) {
  
  plotDataframe <- processedDataframe %>% 
    filter(matrix == matricesNames[[i]])
  
  gList[[i]] <- ggplot(plotDataframe) +
    
    geom_point(aes(x = missingData,
                   y = scaledMean,
                   colour = distance),
               position = position_dodge(width = 0.25),
               cex = 2) +
    
    scale_colour_manual(name = "Distance",
                        labels = c("GED", "MORD"),
                        values = c("red2", "skyblue")) +
    
    geom_errorbar(aes(x = missingData,
                      ymin = scaledLCI,
                      ymax = scaledUCI,
                      colour = distance),
                  position = position_dodge(width = 0.25),
                  width = 0.35) +
    
    ggtitle(label = matricesNames[[i]]) +
    
    xlab("Percentage of missing data added") +
    
    ylim(c(0, 1)) +
    
    ylab("Scaled disparity value") +
    
    facet_wrap(~measure) +
    
    theme_bw()
  
}

pdf(file = "homoDispMeasures.pdf", width = 11, height = 8, paper = "a4r")

for (i in seq_along(matricesNames)) {
  
  plot(gList[[i]])
  
  message(i)
  
}

dev.off()

### FIGURE 5. ------------------------------------------------------------------

plotDataframe <- processedDataframe %>% 
  filter(matrix == "Perez (2017)")

ggplot(plotDataframe) +
  
  geom_point(aes(x = missingData,
                 y = scaledMean,
                 colour = distance),
             position = position_dodge(width = 0.25),
             cex = 2) +
  
  scale_colour_manual(name = "Distance",
                      labels = c("GED", "MORD"),
                      values = c("red2", "skyblue")) +
  
  geom_errorbar(aes(x = missingData,
                    ymin = scaledLCI,
                    ymax = scaledUCI,
                    colour = distance),
                position = position_dodge(width = 0.25),
                width = 0.35) +
  
  xlab("Percentage of missing data added") +
  
  ylim(c(0, 1)) +
  
  ylab("Scaled disparity value") +
  
  facet_wrap(~measure) +
  
  theme_bw()
