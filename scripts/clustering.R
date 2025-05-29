library(ggplot2)
library(doSNOW)
library(fpc)
library(dplyr)
library(QuantumClone)

base_dir <- getwd()

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args_full[grep("--file=", args_full)]))
script_dir <- dirname(script_path)

source(file.path(script_dir, "functions.R"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript clustering.R <output_dir> <sample_name>")
}

output_dir <- normalizePath(args[1])
sample_name <- args[2]
sample_dir <- file.path(output_dir, sample_name)

message("Sample directory:", sample_dir, "\n")

setwd(sample_dir)

snv_file <- list.files(sample_dir, pattern = "_SNVlist\\.txt$", full.names = TRUE)
freec_file <- list.files(sample_dir, pattern = "_freec\\.txt$", full.names = TRUE)

if (length(snv_file) == 0) {
  stop(paste("No SNVlist file found in", sample_dir))
}


vcf <- read.table(snv_file[1], header = TRUE, sep='\t')
vcf$Chr <- as.factor(vcf$Chr)
vcf$Start <- as.integer(vcf$Start)
vcf$Depth <- as.numeric(vcf$Depth)
vcf$Alt <- as.numeric(vcf$Alt)

SNV_list <- list(vcf)

if (length(freec_file) > 0) {
  freec <- read.table(freec_file, header = TRUE, sep = '\t')
  FREEC_list <- list(freec)
  Genotype_provided <- FALSE
} else {
  FREEC_list <- NULL
  Genotype_provided <- TRUE
}

#Sample_names <- lapply(SNV_list, function(z) z[1, 1])
Sample_names <- sample_name



contamination = 0 
nclone_range = 2:5
clone_priors = NULL
prior_weight = NULL
Initializations = 1
preclustering = "FLASH"
simulated = FALSE
epsilon = NULL
save_plot = FALSE
ncores = 1
restrict.to.AB = FALSE
output_directory = NULL
model.selection = "BIC"
optim = "default"
keep.all.models = FALSE
force.single.copy = FALSE



if (is.null(FREEC_list)) {
  message("FREEC_list is empty. Checking that there is a genotype column in all samples...")
  check <- TRUE
  for (i in 1:length(SNV_list)) {
    if (is.null(SNV_list[[i]][, "Genotype"])) {
      warning(paste("Sample", i, "does not have a Genotype column"))
      check <- FALSE
    }
  }
  if (!check) {
    stop("See warnings for more details")
  }
  else {
    message("Genotype is provided.")
    Genotype_provided <- TRUE
  }
} else {
  message("Checking that SNV_list and FREEC_list have the same number of samples...")
  if (length(SNV_list) != length(FREEC_list)) {
    stop("Incorrect input: number of samples different for SNV and FREEC")
  }
  message("Passed")
  Genotype_provided <- FALSE
}
if (length(model.selection) > 1) {
  stop("Model selection can only be one of BIC,AIC or a numeric value.")
} else {
  if (is.character(model.selection)) {
    if (model.selection != "BIC" && model.selection != 
        "AIC") {
      stop("Character value for model selection can only be BIC or AIC")
    }
  }
  else if (!is.numeric(model.selection)) {
    stop("Input model.selection is not numeric - and is not BIC or AIC. Please see documentation for model.selection.")
  }
}
if (is.null(epsilon)) {
  epsilon <- 1/mean(unlist(lapply(SNV_list, function(df) df$Depth)))
  message(paste("epsilon set to:", epsilon))
}
optims_accepted <- c("default", "optimx", "DEoptim", "compound")
if (length(optim) != 1) {
  stop("optim argument can only be of length 1")
} else if (!(optim %in% optims_accepted)) {
  stop(paste("optim can only be one of:", paste(optims_accepted, 
                                                collapse = ",")))
}
if (is.list(SNV_list) && length(contamination) != length(SNV_list)) {
  if (length(contamination) < length(SNV_list)) {
    warning("contamination and SNV_list have different lengths, will repeat contamination")
    contamination <- rep(contamination, times = length(SNV_list))
  }
  else if (length(contamination) < length(SNV_list)) {
    warning("contamination and SNV_list have different lengths, will use first arguments of contamination")
    contamination <- contamination[1:length(SNV_list)]
  }
} else if (is.data.frame(SNV_list) && length(contamination > 
                                           1)) {
  warning("length of contamination > 1 but only one sample, will use only first element")
  contamination <- contamination[1]
}
message(paste("Checking all possibilities for", sample_name))
Cell <- From_freq_to_cell(SNV_list = SNV_list, FREEC_list = FREEC_list, 
                          Sample_names = Sample_names, contamination = contamination, 
                          Genotype_provided = Genotype_provided, save_plot = save_plot, 
                          restrict.to.AB = restrict.to.AB, output_directory = output_directory, 
                          force.single.copy = force.single.copy)

message("Starting clustering for:", sample_name, "\n")

r <- Cluster_plot_from_cell(Cell = Cell, nclone_range = nclone_range, 
                            epsilon = epsilon, Sample_names = Sample_names, preclustering = preclustering, 
                            clone_priors = clone_priors, prior_weight = prior_weight, 
                            Initializations = Initializations, simulated = simulated, 
                            save_plot = save_plot, contamination = contamination, 
                            ncores = ncores, output_directory = output_directory, 
                            model.selection = model.selection, optim = optim, keep.all.models = keep.all.models)
if (length(SNV_list) == 1 & save_plot) {
  q <- One_D_plot(EM_out = r, contamination = contamination)
  if (is.null(output_directory)) {
    ggplot2::ggsave(plot = q, filename = paste(sample_name, 
                                               "/", "Density", sample_name, ".pdf", sep = ""), 
                    width = 6.04, height = 6.04)
  }
  else {
    ggplot2::ggsave(plot = q, filename = paste(output_directory, 
                                               "/", "Density", sample_name, ".pdf", sep = ""), 
                    width = 6.04, height = 6.04)
  }
}
if (length(unique(Cell[[1]]$id)) == length(Cell[[1]]$id)) {
  message("Post-processing output")
  if (keep.all.models) {
    print('all models')
    r <- lapply(r, FUN = function(z) {
      Tidy_output(r = z, Genotype_provided = Genotype_provided, 
                  SNV_list = SNV_list)
    })
  }
  else {
    r <- Tidy_output(r = r, Genotype_provided = Genotype_provided, 
                     SNV_list = SNV_list)
  }
} else {
  message("Keeping most likely state for each variant...")
  if (keep.all.models) {
    Crits <- as.numeric(as.character(unlist(lapply(r, 
                                                   function(z) {
                                                     z$Crit
                                                   }))))
    r <- r[[which.min(Crits)]]
  }
  r <- Tidy_output(r = r, Genotype_provided = Genotype_provided, 
                   SNV_list = SNV_list)
  message("Reclustering with most likely state of each variant...")
  r <- Cluster_plot_from_cell(Cell = r$filtered.data, nclone_range = nclone_range, 
                              epsilon = epsilon, Sample_names = Sample_names, preclustering = preclustering, 
                              clone_priors = clone_priors, prior_weight = prior_weight, 
                              Initializations = Initializations, simulated = simulated, 
                              save_plot = save_plot, contamination = contamination, 
                              ncores = ncores, output_directory = output_directory, 
                              model.selection = model.selection, optim = optim, 
                              keep.all.models = keep.all.models)
}
if (is.list(r$EM.output$centers)) {
  if (!keep.all.models) {
    normalized.centers <- list()
    for (i in 1:length(r$EM.output$centers)) {
      normalized.centers[[i]] <- r$EM.output$centers[[i]]/(1 - 
                                                             contamination[i])
    }
    r$EM.output$normalized.centers <- normalized.centers
  }
  else {
    for (mod in 1:length(r)) {
      normalized.centers <- list()
      for (i in 1:length(r$EM.output$centers)) {
        normalized.centers[[i]] <- r[[mod]]$EM.output$centers[[i]]/(1 - 
                                                                      contamination[i])
      }
      r[[mod]]$EM.output$normalized.centers <- normalized.centers
    }
  }
}

clustering <- r

cat("Clustering finished for:", sample_name, "\n")


df_clusters <- data.frame(
  Position = seq_along(clustering$cluster),
  Number = clustering$cluster)
write.csv(df_clusters, file.path(sample_dir, "clustering.csv"))

fik_df <- as.data.frame(clustering$EM.output$fik)
write.csv(fik_df, file.path(sample_dir, "FIK.csv"), row.names = FALSE)

weights_df <- data.frame(weights = clustering$EM.output$weights)
write.csv(weights_df, file.path(sample_dir, "weights.csv"))

val_df <- data.frame(val = clustering$EM.output$val)
write.csv(val_df, file.path(sample_dir, "VAL.csv"))

centers_df <- do.call(rbind, lapply(clustering$EM.output$centers, as.data.frame))
centers_df <- data.frame(centers_df)
write.csv(centers_df, file.path(sample_dir, "centers.csv"))

normalized_centers_df <- do.call(rbind, lapply(clustering$EM.output$normalized.centers, as.data.frame))
normalized_centers_df <- data.frame(normalized_centers_df)
write.csv(normalized_centers_df, file.path(sample_dir, "normalized_centers.csv"))

filtered_df <- do.call(rbind, lapply(clustering$filtered.data, as.data.frame))
filtered_df <- data.frame(filtered_df)
write.csv(filtered_df, file.path(sample_dir, "filtered.csv"))



#clustering <- One_step_clustering(vcf_list, FREEC_list = freec_list, contamination = 0, 
 #                             nclone_range = 2:5, clone_priors = NULL, prior_weight = NULL, 
  #                            Initializations = 1, preclustering = "FLASH", simulated = FALSE, 
   #                           epsilon = NULL, save_plot = TRUE, ncores = 1, restrict.to.AB = FALSE, 
    #                          output_directory = NULL, model.selection = "BIC", optim = "default", 
     #                         keep.all.models = FALSE, force.single.copy = FALSE)




