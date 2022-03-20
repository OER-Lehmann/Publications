### ANALYSIS 2 OF LEHMANN ET AL. (2018) ########################################
### See the associated README.txt file for the description of this analysis. ###

##### SET UP. ------------------------------------------------------------------

# Clean the workspace.

rm(list = ls())

# This script is based in version 0.2 of Claddis and the data structures and
# functions it provides. THIS SCRIPT WILL NOT WORK FOR VERSIONS > 0.2. In order
# to install version 0.2, use the following code (commented for safety).

# library("devtools")
# 
# install_github("graemetlloyd/Claddis",
#                ref = "ffa42d63d1f28eed387949b0aa247e9f2f49d6ee")

# Load required libraries.

library("Claddis")
library("ggplot2")
library("dplyr")

source("scripts/disp_functions.R")

##### DATA IMPORT. -------------------------------------------------------------

# Set the folder where the matrices are.

folder <- "Complete_matrices"

# Get the names of the .nex files and how many are there.

files <- dir(path = paste(getwd(), "/", folder, sep = ""), pattern = "\\.nex")
n_files <- length(files)

# Remove the .nex for the filenames.

good_names <- gsub(pattern = "\\.nex", replacement = "", x = files)
good_names <- gsub(pattern = "_", replacement = " ", x = good_names)

# Set up the value of the beta parameter for the Beta distribution.

beta_val <- (6 / 0.35) - 6

# Set the number of replications.

reps <- 200

# Create the vectors for the data.

mtx <- rep(good_names, each = (reps * 2))
repl <- rep(as.character(seq_len(reps)), each = (n_files * 2))
metric <- rep(c("GED", "MORD"), times = (n_files * reps))
p_value <- rep(0, (n_files * reps * 2))
corcoeff <- rep(0, (n_files * reps * 2))

index <- 1

## Start loop.

# Set random seed.

set.seed(25122015)

t0 <- Sys.time()

for (i in seq_len(n_files)) {
  
  ## LEVEL 1: READ MATRIX. ##
  
  for (j in seq_len(reps)) {
    
    m <- ReadMorphNexus(file = paste0(folder, "/", files[i]))
    
    ## LEVEL 2: REPLICATIONS PER MATRIX. ##
    
    # The 'check' object allows to... check... whether there are at least 3
    # characters scored for each taxon and at least 2 taxa scored for each
    # character.
    
    check <- TRUE
    
    while (check) {
      
      for (k in seq_len(nrow(m$matrix))) {
        
        ## LEVEL 3: MISSING ENTRIES PER TAXON. ##
        
        pMD <- rbeta(1, 6, beta_val)
        
        nc <- ncol(m$matrix)
        
        m$matrix[k, sample(seq_len(nc), size = floor(nc * pMD), replace = FALSE)] <- NA
        
      }
      
      notna_row <- apply(X = m$matrix, 
                         MARGIN = 1,
                         FUN = function (x) {sum(!is.na(x))})
      
      notna_col <- apply(X = m$matrix,
                         MARGIN = 2,
                         FUN = function (x) {sum(!is.na(x))})
      
      if (all(notna_row > 3) & all(notna_col > 2)) {
        
        check <- FALSE
        
      } else {
        
        check <- TRUE
        
      }
      
    }
    
    dM <- MorphDistMatrixFast(m, c("GED", "Max"))
    dMORD <- silent_TrimMorphDistMatrix(dM$max.dist.matrix)$dist.matrix
    dGED <- dM$GED.dist.matrix[rownames(dMORD), colnames(dMORD)]
    
    dMlist <- list(dGED, dMORD)
    
    for (k in seq_along(dMlist)) {
      
      ## LEVEL 3: ANALYSIS OF THE DISTANCE MATRICES. ##
      
      P <- pcoa(dMlist[[k]], "lingoes")
      
      if (length(P) == 7) {
        scores <- P$vectors.cor
      } else {scores <- P$vectors}
      
      centroid <- c(0, 0, 0)
      
      distances <- distanceToCentroid(vectors = scores, axes = "all")
      
      propMD <- rep(0, length(distances))
      
      for (l in seq_along(propMD)) {
        
        ## LEVEL 4: PROPORTION OF MISSING ENTRIES PER TAXON. ##
        
        propMD[l] <-  sum(is.na(m$matrix[l, ])) / ncol(m$matrix)
        
      }
      
      cormock <- cor.test(distances, propMD, method = "spearman")
      
      p_value[index] <- cormock$p.value
      corcoeff[index] <- cormock$estimate
      
      message("Completed ", index, " out of ", length(repl), ". ",
              round(difftime(Sys.time(), t0, units = "s"), 0), " seconds elapsed.")
      
      index <- index + 1
      
    }
    
  }
  
}


# Construct the data frame.

simulation_df <- data.frame(mtx, repl, metric, corcoeff, p_value)

# Save the data frame.

save(file = "Results/Analysis_2.R", simulation_df)


### FIGURE 4. ------------------------------------------------------------------

# Load data.

load("Results/Analysis_2.R")

# Load packages.

library("ggplot2")
library("reshape2")

# Create the discriminated data frame.

signifDataframe <- simulation_df %>% 
  
  group_by(mtx, metric) %>% 
  
  summarize(SN = sum(ifelse(test = p_value < 0.05,
                            yes = ifelse(test = corcoeff < 0,
                                         yes = TRUE,
                                         no = FALSE),
                            no = FALSE)),
            
            NS = sum(ifelse(test = p_value >= 0.05,
                            yes = TRUE,
                            no = FALSE)),
            
            SP = sum(ifelse(test = p_value < 0.05,
                            yes = ifelse(test = corcoeff > 0,
                                         yes = TRUE,
                                         no = FALSE),
                            no = FALSE)))

# Melt it.

signifDataframe <- melt(signifDataframe)

# Keep the desired matrices.

plotDataframe <- signifDataframe %>% 
  filter(mtx %in% c("Fiser et al. (2008)",
                    "Mayr (2011)",
                    "Perez (2017)"))

# Save Figure 4.

svg(filename = "Figure 4 raw.svg", width = 7, height = 3)

ggplot(plotDataframe) +
  
  geom_col(aes(x = metric, y = value / 2, fill = variable),
           colour = "black",
           position = position_dodge()) +
  
  scale_fill_manual(name = "Correlation",
                    labels = c("Significant and negative",
                               "Not significant",
                               "Significant and positive"),
                    values = c("black", "grey95", "darkgrey"),
                    guide = FALSE) +
  
  ylim(c(0, 100)) +
  
  ylab("Percentage of matrices") +
  
  xlab("Distance measure") +
  
  facet_wrap(~mtx) +
  
  theme_bw() +
  
  theme(panel.grid = element_blank())

dev.off()
