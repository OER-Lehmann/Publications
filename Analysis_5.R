rm(list = ls())

library(Claddis)

source("scripts/disp_functions.R")

folder <- "Ciampaglio_matrices/"

# Get the names of the .nex files.

files <- dir(path = folder, pattern = "\\.nex")

# Get the clean names of the matrices.

good_names <- gsub(pattern = "\\.nex", x = files, replacement = "")
good_names <- gsub(pattern = "_", x = good_names, replacement = " ")

# Read the morphological matrices, creating of list with all of them.

matrices <- lapply(paste0(folder, "/", files), ReadMorphNexus)

Ciampaglio <- function (m, p) {
  
  mtx <- m$matrix
  
  n_cells <- ceiling(length(m$matrix) * p)
  
  for (i in seq_len(n_cells)) {
    
    co <- sample(seq_len(ncol(mtx)), 1)
    
    ro <- sample(seq_len(nrow(mtx)), 1)
    
    mtx[ro, co] <- as.character(mean(as.numeric(mtx[-ro, co]), na.rm = TRUE))
    
  }
  
  m$matrix <- mtx
  
  return(m)
  
}

CSoV <- function (vectors) {
  
  return(sqrt(sum(apply(X = vectors, MARGIN = 2, FUN = var))))
  
}

CMPD <- function (vectors) {
  
  return(mean(dist(vectors, method = "euc")))
  
}

### .------

props <- seq(0.1, 0.8, 0.1)

reps <- 100

n_mtx <- length(files)

props_factor <- factor(x = as.character(c(rep(0, 8), rep(props, each = reps * 8))),
                       levels = as.character(seq(0, 0.8, 0.1)),
                       labels = as.character(seq(0, 0.8, 0.1)))

bigN <- length(props_factor) * n_mtx

method_factor <- factor(x = rep(c("Ciampaglio et al. (2001)", "Random replacement"), each = 4),
                        levels = c("Ciampaglio et al. (2001)", "Random replacement"),
                        labels = c("Ciampaglio et al. (2001)", "Random replacement"))

measure_factor <- factor(x = c("C_SoV", "C_MPD", "SoV", "MPD"),
                         levels = c("C_SoV", "C_MPD", "SoV", "MPD"),
                         labels = c("C_SoV", "C_MPD", "SoV", "MPD"))

df <- data.frame(matrix = rep(good_names, each = length(props_factor)),
                 propMD = rep(props_factor, n_mtx),
                 method = rep(method_factor, bigN / 8),
                 measure = rep(measure_factor, bigN / 4),
                 values = rep(0, length(props_factor) * n_mtx))

index <- 1

time_0 <- Sys.time()

for (i in seq_len(n_mtx)) {
  
  # Base case.
  
  m <- matrices[[i]]
  
  dm <- MorphDistMatrixFast(m, "GED")
  
  v <- L.PCoA(as.dist(dm), correction = "lingoes")$results$vectors
  
  df[index, "values"] <- CSoV(v)
  index <- index + 1
  
  df[index, "values"] <- CMPD(v)
  index <- index + 1
  
  df[index, "values"] <- sum_of_variance(v, axes = "all")
  index <- index + 1
  
  df[index, "values"] <- mpd(dm)
  index <- index + 1
  
  df[index, "values"] <- CSoV(v)
  index <- index + 1
  
  df[index, "values"] <- CMPD(v)
  index <- index + 1
  
  df[index, "values"] <- sum_of_variance(v, axes = "all")
  index <- index + 1
  
  df[index, "values"] <- mpd(dm)
  index <- index + 1
  
  for (j in seq_along(props)) {
    
    for (k in seq_len(reps)) {
      
      trep <- Sys.time()
      
      cm <- Ciampaglio(matrices[[i]], props[j])
      rrm <- missify(matrices[[i]], props[j])
      
      dcm <- MorphDistMatrixFast(cm, "GED")
      drrm <- MorphDistMatrixFast(rrm, "GED")
      
      cv <- L.PCoA(as.dist(dcm), correction = "lingoes")$results$vectors
      rrv <- L.PCoA(as.dist(drrm), correction = "lingoes")$results$vectors
      
      df[index, "values"] <- CSoV(cv)
      index <- index + 1
      
      df[index, "values"] <- CMPD(cv)
      index <- index + 1
      
      df[index, "values"] <- sum_of_variance(cv, axes = "all")
      index <- index + 1
      
      df[index, "values"] <- mpd(dcm)
      index <- index + 1
      
      df[index, "values"] <- CSoV(rrv)
      index <- index + 1
      
      df[index, "values"] <- CMPD(rrv)
      index <- index + 1
      
      df[index, "values"] <- sum_of_variance(rrv, axes = "all")
      index <- index + 1
      
      df[index, "values"] <- mpd(drrm)
      index <- index + 1
      
      time_out <- as.numeric(round(difftime(Sys.time(),
                                            time_0,
                                            units = "min"),
                                   digits = 2))
      
      trepfinal <- as.numeric(round(difftime(Sys.time(),
                                             trep,
                                             units = "sec"),
                                    digits = 2))
      
      message("Matrix ", i, "/", n_mtx, ". Prop ", props[j], ". ",
              "Replication ", k, "/", reps, ". ",
              time_out, " minutes elapsed. Last replication took ", trepfinal,
              " seconds.")
      
    }
    
  }
  
}

save(file = "Ciampaglio_data.R", df)

##### FIGURES. -----------------------------------------------------------------


# Clean the work space (commented for safety).

# rm(list = ls())

# Load required libraries.

library(dplyr)
library(ggplot2)

# Load the data.

load("Ciampaglio_data.R")

### DATA CLEANING. -------------------------------------------------------------

# Get the names of the matrices.

matricesNames <- levels(df$matrix)

## Change the labels of some of the factors.

# Change the labels of the proportion of missing data to percentages.

df$propMD <- factor(x = df$propMD,
                    labels = seq(0, 80, 10))

# Change the levels of the disparity measures (The MPDs were equal).

df$measure <- factor(x = df$measure,
                     levels = c("MPD", "C_MPD", "SoV", "C_SoV"),
                     labels = c("MPD",
                                "Cea2001 MPD",
                                "Sum of variances",
                                "Cea2001 sum of variances"))

## Create a data frame with the mean +- sd.

meanDataframe <- df %>%
  
  group_by(matrix, propMD, method, measure) %>% 
  
  summarize(meanValue = mean(values),
            lci = mean(values) - sd(values),
            uci = mean(values) + sd(values)) %>% 
  
  ungroup() %>% 
  
  mutate(lci = replace(x = lci, which(propMD == "0"), NA),
         uci = replace(x = uci, which(propMD == "0"), NA)) %>% 
  
  group_by(matrix, method, measure) %>% 
  
  mutate(scaledMeanValue = meanValue / max(c(meanValue, lci, uci), na.rm = TRUE),
         scaledLCI = lci / max(c(meanValue, lci, uci), na.rm = TRUE),
         scaledUCI = uci / max(c(meanValue, lci, uci), na.rm = TRUE))

### CREATE THE PLOTS. ----------------------------------------------------------

# Create a list for the plots.

ggplotList <- vector(mode = "list", length = length(matricesNames))

# Fill the list with the plots.

for (i in seq_along(matricesNames)) {
  
  plotDataframe <- meanDataframe %>% 
    filter(matrix == matricesNames[i])
  
  ggplotList[[i]] <- ggplot(plotDataframe) +
    
    geom_point(aes(x = propMD,
                   y = scaledMeanValue,
                   colour = measure)) +
    
    geom_errorbar(aes(x = propMD,
                      ymin = scaledLCI,
                      ymax = scaledUCI,
                      colour = measure),
                  width = 0.5) +
    
    facet_wrap(~method + measure,
               nrow = 2,
               ncol = 4) +
    
    ggtitle(label = matricesNames[i],
            subtitle = expression("Mean +- SD")) +
    
    xlab("Percentage of missing data added") +
    
    ylim(c(0, 1)) +
    
    ylab("Mean disparity") +
    
    theme_bw()
  
}

### CREATE THE FIGURE. ---------------------------------------------------------

dev.off()

pdf(file = "",
    width = 11,
    height = 6,
    paper = "a4r")

for (i in seq_along(ggplotList)) {
  
  plot(ggplotList[[i]])
  
  message(i, " / ", length(matricesNames), ".")
  
}

dev.off()