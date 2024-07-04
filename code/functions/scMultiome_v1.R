#!/usr/bin/env Rscript

# Helper functions to run scMultiome preprocessing pipeline
# @Author: Clara Savary (clara.savary@yahoo.fr, OmixAnalytics)


# Utility functions -----------------------------------------------------------

#' Import single-cell multiome data
#'
#' This function imports the gene and peak counts from Cell Ranger output files.
#'
#' @param path A string character indicating the path of the Cell Ranger output files.
#' @param sample A string character indicating the name of the sample. Important as it must match the Cell Ranger repository name.
#' @param min.cells A numeric value indicating that features must be detected in at least this many cells.
#' @param annotation A GRanges object containing genome annotation.
#' @param seqinfo A formal class Seqinfo containing genome annotation informations.
#'
#' @export
#' @return A Seurat object with three assays: RNA, ATAC and TF.
#' @examples
#' # Basic usage
#' import.scMultiome(path = "cellranger_out/data/", sample = "pbmc", min.cells = 1, annotation, seqinfo)
#' 
#' @author Clara Savary
import.scMultiome <- function(
    path,
    sample,
    min.cells = 1,
    annotation,
    seqinfo) {
  
  # Check arguments
  if (missing(path)) {
    stop("Argument 'path' is missing.")
  }
  
  if (!is.numeric(min.cells)) {
    stop("Argument 'min.cells' must be numeric.")
  }
  
  if (missing(annotation)) {
    stop("Argument 'annotation' is missing.")
  }
  
  if (missing(seqinfo)) {
    stop("Argument 'seqinfo' is missing.")
  }
  
  # Load cellranger-produced metadata with per-cell stats
  cr_metadata <- read.csv(
    file.path(path, "per_barcode_metrics.csv"),
    header    = TRUE,
    row.names = 1
  )
  
  # Load filtered matrices - contains both RNA and ATAC matrices
  counts <- Read10X_h5(file.path(path, "filtered_feature_bc_matrix.h5"))
  
  # list the data modalities for which counts are provided
  names(counts) %>% print
  
  # Initialize Seurat object with RNA
  seu <- CreateSeuratObject(
    counts    = counts$`Gene Expression`,
    assay     = "RNA",
    meta.data = cr_metadata,
    project   = sample,
    min.cells = min.cells
  )
    
  # Create the ATACseq assay
  seu[["ATAC"]] <- CreateChromatinAssay(
    counts    = counts$Peaks, 
    sep       = c(":", "-"),
    genome    = seqinfo,
    fragments = file.path(path, "atac_fragments.tsv.gz"),
    min.cells = min.cells
  )
  
  # Set gene annotations
  Annotation(seu[["ATAC"]]) <- annotation
  
  # Load filtered TF matrix
  tf_counts <- Read10X_h5(file.path(path, "analysis/tf_analysis/filtered_tf_bc_matrix.h5"))
  
  # Create the TF assay
  assay_TF <- CreateAssayObject(
    counts    = tf_counts,
    min.cells = min.cells
  )
  
  # Add TF assay to Seurat object
  seu[['TF']] <- assay_TF
  
  # Verbose
  message("Successful import of ", sample, " scMultiome data.")
  
  # Return Seurat object
  return(seu)

}


#' Add mitochondrial and ribosomal ratios
#'
#' This function add mitochondrial and ribosomal percentages to a Seurat object.
#'
#' @param seu A Seurat object.
#'
#' @export
#' @return A Seurat object with mitochondrial and ribosomal percentages.
#' @examples
#' # Basic usage
#' seu <- add_mito_ribo(seu)
#' 
#' @author Clara Savary
add_mito_ribo <- function(seu) {
  
  # Set default assay
  DefaultAssay(seu) <- "RNA"
  
  # Identify mitochondrial genes, which depends on the species, due to gene name differences
  if (species == "h_sapiens") {
    seu <- PercentageFeatureSet(seu, "^MT-", col.name = "percent_mito")
  } else if (species == "m_musculus") {
    seu <- PercentageFeatureSet(seu, "^mt-", col.name = "percent_mito")
  }
  
  # Identify ribosomal genes, which depends on the species, due to gene name differences
  if (species == "h_sapiens") {
    seu <- PercentageFeatureSet(seu, "^RPS|^RPL|^MRPS|^MRPL|^RP[SL]", col.name = "percent_ribo")
  } else if (species == "m_musculus") {
    seu <- PercentageFeatureSet(seu, "^Rps|^Rpl|^Mrps|^Mrpl|^Rp[sl]", col.name = "percent_ribo")
  }
  
  # Return
  seu
  
}


#' Add complexity metrics
#'
#' This function add complexity metrics to a Seurat object.
#'
#' @param seu A Seurat object.
#'
#' @export
#' @return A Seurat object with complexity metrics.
#' @examples
#' # Basic usage
#' seu <- add_complexity(seu)
#' 
#' @author Clara Savary
add_complexity <- function(seu) {
  
  # Set default assay
  DefaultAssay(seu) <- "RNA"
  
  # Add number of genes per UMI for each cell to metadata object
  seu[["log10nGene"]] <- log10(seu@meta.data$nFeature_RNA)
  
  # Add number of genes per UMI for each cell to metadata object
  seu[["log10nUMI"]] <- log10(seu@meta.data$nCount_RNA)
  
  # Add number of genes per UMI for each cell to metadata object
  seu[["log10GenesPerUMI"]] <- seu@meta.data$log10nGene/seu@meta.data$log10nUMI
  
  # Return
  seu
  
}


#' Add scATACseq metrics
#'
#' This function add scATACseq metrics to a Seurat object.
#'
#' @param seu A Seurat object.
#'
#' @export
#' @return A Seurat object with scATACseq metrics.
#' @examples
#' # Basic usage
#' seu <- add_atac_metrics(seu)
#' 
#' @author Clara Savary
add_atac_metrics <- function(seu) {
  
  # Set default assay
  DefaultAssay(seu) <- "ATAC"
  
  # Compute stats using Signac functions
  seu <- NucleosomeSignal(object = seu)
  seu <- TSSEnrichment(seu, fast = FALSE)
  
  # Return
  seu
  
}


#' Function to calculate thresholds of scMultiome metrics
#'
#' This function to calculate thresholds of scMultiome metrics of a Seurat object.
#'
#' @param seu A Seurat object.
#' @param sd_QC A numerical value to indicate number of standard deviation used for thresholds calculation.
#'
#' @export
#' @return A dataframe with thresholds of scMultiome metrics of a Seurat object.
#' @examples
#' # Basic usage
#' seu_thresholds(seu, sd_QC = 2)
#' 
#' @author Clara Savary
seu_thresholds <- function(seu, sd_QC) {
  
  # Sample name
  var_name <- seu@meta.data[["orig.ident"]] %>% levels()
  
  ### nFeature_RNA
  # the minimum number of features will be the greater of:
  # 400, or 2 standard deviations below the mean
  min_features <- max(
    config$min_nFeature_RNA,
    round(mean(seu@meta.data$nFeature_RNA) - sd_QC*sd(seu@meta.data$nFeature_RNA))
  )
  
  max_features <- round(
    mean(seu@meta.data$nFeature_RNA) + sd_QC*sd(seu@meta.data$nFeature_RNA)
  )
  
  ### percent_mito
  # by default
  min_mito <- 0
  
  # the max mitochondrial content will be the maximum of:
  # 5%, or 2 standard deviations above the mean
  # the parameter config$max_mito allows to set a hard upper threshold,
  # which takes precedence
  max_mito <- ifelse(
    !is.null(config$max_mito),
    config$max_mito,
    max(config$max_mito, round(mean(seu@meta.data$percent_mito) + 2*sd(seu@meta.data$percent_mito)))
  )

  ### nCount_RNA
  # set a max of 0 in case the value 2 standard deviations below the mean
  # is negative
  min_umi <- max(
    config$min_nCount_RNA,
    round(mean(seu@meta.data$nCount_RNA) - sd_QC*sd(seu@meta.data$nCount_RNA))
  )
  
  max_umi <- round(
    mean(seu@meta.data$nCount_RNA) + sd_QC*sd(seu@meta.data$nCount_RNA)
  )
  
  ### atac_peak_region_fragments
  min_peak <- max(
    config$min_peak,
    round(mean(seu@meta.data$atac_peak_region_fragments) - sd_QC*sd(seu@meta.data$atac_peak_region_fragments))
  )
  
  max_peak <- round(
    mean(seu@meta.data$atac_peak_region_fragments) + sd_QC*sd(seu@meta.data$atac_peak_region_fragments)
  )
  
  ### nucleosome_signal
  max_nucleosome_signal <- ifelse(
    !is.na(config$max_nucleosome_signal),
    config$max_nucleosome_signal,
    round(mean(seu@meta.data$nucleosome_signal) + 2*sd(seu@meta.data$nucleosome_signal))
    )
  
  ### TSS.enrichment
  max_TSS <- round(
    mean(seu@meta.data$TSS.enrichment) + sd_QC*sd(seu@meta.data$TSS.enrichment)
  )
  min_TSS <- round(
    mean(seu@meta.data$TSS.enrichment) - sd_QC*sd(seu@meta.data$TSS.enrichment)
  )  
  
  thresholds <- data.frame(
    "sample"                = var_name,
    "min_features"          = min_features,
    "max_features"          = max_features,
    "min_mito"              = min_mito,
    "max_mito"              = max_mito,
    "min_umi"               = min_umi,
    "max_umi"               = max_umi,
    "min_peak"              = min_peak,
    "max_peak"              = max_peak,
    "max_nucleosome_signal" = max_nucleosome_signal,
    "max_TSS"               = max_TSS,
    "min_TSS"               = min_TSS
  )
  
  # Return thresholds
  thresholds
  
}

#' Function to filter cells
#'
#' This function filters cells of a Seurat object.
#'
#' @param seu A Seurat object.
#' 
#' @export
#' @return A Seurat object with filtered cells.
#' @examples
#' # Basic usage
#' seu_filt <- filter_cells(seu)
#' 
#' @author Clara Savary
filter_cells <- function(seu) {
  
  sample <- seu@meta.data[["orig.ident"]] %>% levels()
  message(sample)
  
  # Filter cells
  seu <- subset(
    x      = seu,
    subset =
      nFeature_RNA > thresholds_all$min_features[thresholds_all$sample %in% sample] & 
      nFeature_RNA < thresholds_all$max_features[thresholds_all$sample %in% sample] &
      nCount_RNA > thresholds_all$min_umi[thresholds_all$sample %in% sample] &
      nCount_RNA < thresholds_all$max_umi[thresholds_all$sample %in% sample] &
      percent_mito > thresholds_all$min_mito[thresholds_all$sample %in% sample] &
      percent_mito < thresholds_all$max_mito[thresholds_all$sample %in% sample] &
      atac_peak_region_fragments > thresholds_all$min_peak[thresholds_all$sample %in% sample] &
      atac_peak_region_fragments < thresholds_all$max_peak[thresholds_all$sample %in% sample] &
      nucleosome_signal < thresholds_all$max_nucleosome_signal[thresholds_all$sample %in% sample] &
      TSS.enrichment > thresholds_all$min_TSS[thresholds_all$sample %in% sample] &
      TSS.enrichment < thresholds_all$max_TSS[thresholds_all$sample %in% sample]
  )
  
  # Return filtered Seurat object
  seu
  
}


#' Function to return filter cells criterion.
#'
#' This function return filter cells criterion of a Seurat object.
#'
#' @param criterion A character string indicating the feature or criterion name.
#' @param sample A character string indicating the name of the Seurat object.
#' 
#' @export
#' @return A vector object with criterion metrics.
#' @examples
#' # Basic usage
#' criterion(criterion = "nFeature_RNA", sample = "pbmc")
#' 
#' @author Selin Jessa
criterion <- function(criterion, sample) {
  
  seurat_prefilt <- seu_obj[[sample]]
  seurat         <- seu_obj_filt[[sample]]
  
  min_pre   <- round(min(seurat_prefilt@meta.data  %>% pull(criterion)), 2)
  mean_pre  <- mean(seurat_prefilt@meta.data %>% pull(criterion))
  max_pre   <- max(seurat_prefilt@meta.data  %>% pull(criterion))
  sd_pre    <- sd(seurat_prefilt@meta.data   %>% pull(criterion))
  
  min_post  <- min(seurat@meta.data  %>% pull(criterion))
  mean_post <- mean(seurat@meta.data %>% pull(criterion))
  max_post  <- max(seurat@meta.data  %>% pull(criterion))
  sd_post   <- sd(seurat@meta.data   %>% pull(criterion))
  
  return(c("min.preQC"   = min_pre,
           "mean.preQC"  = mean_pre,
           "max.preQC"   = max_pre,
           "sd.preQC"    = sd_pre,
           "min.postQC"  = min_post,
           "mean.postQC" = mean_post,
           "max.postQC"  = max_post,
           "sd.postQC"   = sd_post))
  
}


#' Function to return filter cells criterion metrics and warnings.
#'
#' This function return filter cells criterion metrics and warnings of a Seurat object.
#'
#' @param sample A character string indicating the name of the Seurat object.
#' 
#' @export
#' @return A table of filter cells criterion metrics and warnings.
#' @examples
#' # Basic usage
#' seu_criterion(sample = "pbmc")
#' 
#' @author Selin Jessa
seu_criterion <- function(sample) {
  
  seurat_prefilt <- seu_obj[[sample]]
  seurat         <- seu_obj_filt[[sample]]
  
  # Filtering criteria
  filtering_criteria <- c("nFeature_RNA", "nCount_RNA",  "percent_mito", "atac_peak_region_fragments", "nucleosome_signal", "TSS.enrichment")
  
  # Compute summary stats for each metric
  filtering_metrics <- sapply(filtering_criteria, criterion, sample = sample)
  
  # Round to 2 decimal places
  filtering_metrics <- apply(filtering_metrics, 2, round, 3)
  
  # Transform into a dataframe
  filtering_metrics <- filtering_metrics %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "criterion") %>%
    dplyr::select(criterion, min.preQC, min.postQC, max.preQC, max.postQC, mean.preQC, mean.postQC, sd.preQC, sd.postQC)
  
  # Add thresholds
  filtering_metrics$min.threshold <- c(
    thresholds_all$min_features[thresholds_all$sample %in% sample],
    thresholds_all$min_umi[thresholds_all$sample %in% sample],
    thresholds_all$min_mito[thresholds_all$sample %in% sample],
    thresholds_all$min_peak[thresholds_all$sample %in% sample],
    NA,
    thresholds_all$min_TSS[thresholds_all$sample %in% sample]
  )
  
  filtering_metrics$max.threshold <- c(
    thresholds_all$max_features[thresholds_all$sample %in% sample],
    thresholds_all$max_umi[thresholds_all$sample %in% sample],
    thresholds_all$max_mito[thresholds_all$sample %in% sample],
    thresholds_all$max_peak[thresholds_all$sample %in% sample],
    thresholds_all$max_nucleosome_signal[thresholds_all$sample %in% sample],
    thresholds_all$max_TSS[thresholds_all$sample %in% sample]
  )
  
  filtering_metrics %>% print()
  
  # Compute number of cells before and after filtering
  N_cells_metrics <- data.frame(
    "sample"         = sample,
    "N_cells_before" = dim(seurat_prefilt@meta.data)[1],
    "N_cells_after"  = dim(seurat@meta.data)[1]
  ) %>%
    mutate(Prop_kept = round(N_cells_after / N_cells_before, 2))
  
  N_cells_metrics %>% print()
  
  if (N_cells_metrics$N_cells_after < 1000) warnings$LOW_N_CELLS <- TRUE
  if (N_cells_metrics$Prop_kept < 0.6) warnings$HIGH_PROP_FILTERED <- TRUE
  if (filtering_metrics[filtering_metrics$criterion == "nCount_RNA", ]$mean.postQC < 2000) warnings$LOW_AVG_UMI <- TRUE
  if (filtering_metrics[filtering_metrics$criterion == "percent_mito", ]$max.postQC > 5) warnings$HIGH_MITO <- TRUE
  
  # Convert the list to a data frame
  warnings_df <- data.frame(
    WARNING = names(warnings),
    FLAG    = unname(unlist(warnings))
  )
  
  warnings_df$FLAG <- cell_spec(warnings_df$FLAG, background = ifelse(warnings_df$FLAG, "red", "green"))
  
  print(warnings_df %>% 
          kbl(escape = FALSE) %>% 
          kable_styling(position = "center")
  )
  cat("\n")
  
}


#' Function to perform peaks calling.
#'
#' This function perform peaks calling on a Seurat object.
#'
#' @param seu A Seurat object.
#' 
#' @export
#' @return A Seurat object with peaks calling.
#' @examples
#' # Basic usage
#' seu <- peak_call(seu)
#' 
#' @author Clara Savary
peak_call <- function(seu) {
  
  # Set up default assay
  DefaultAssay(seu) <- "ATAC"
  
  # Call peaks using MACS2
  # NOTE: this won't work in RStudio b/c the python version (3.7) needed for MACS2
  # conflicts with the python version used with jupyter (3.6)
  peaks <- CallPeaks(seu, macs2.path = config$macs2_path)
  
  # Remove peaks on non standard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  genome_id <- seqinfo@genome %>% unique()
  if(genome_id == "hg19"){peaks = subsetByOverlaps(x = peaks, ranges = Signac::blacklist_hg19, invert = TRUE)}
  if(genome_id == "hg38"){peaks = subsetByOverlaps(x = peaks, ranges = Signac::blacklist_hg38, invert = TRUE)}
  if(genome_id == "mm10"){peaks = subsetByOverlaps(x = peaks, ranges = Signac::blacklist_mm10, invert = TRUE)}
  
  # Quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(seu),
    features  = peaks,
    cells     = colnames(seu)
  )
  
  # Get fragment path
  sample <- seu@meta.data[["orig.ident"]] %>% levels()
  sample_path <- file.path(config$cellranger_dir, sample)
  
  # Create a new assay using the MACS2 peak set and add it to the Seurat object
  seu[["peaks"]] <- CreateChromatinAssay(
    counts    = macs2_counts,
    genome    = seqinfo,
    fragments = file.path(sample_path, "atac_fragments.tsv.gz"),
    min.cells = config$min_cells
  )
  
  # Set gene annotations
  Annotation(seu[["peaks"]]) <- annotation
  
  # Return Seurat object
  seu
  
}


#' Function to perform normalisation of a scRNA assay.
#'
#' This function return perform normalisation of a scRNA assay a Seurat object.
#'
#' @param seu A Seurat object.
#' 
#' @export
#' @return A Seurat object with normalised scRNA assay.
#' @examples
#' # Basic usage
#' seu <- rna_norm(sample = "pmbc")
#' 
#' @author Clara Savary
rna_norm <- function(seu) {
  
  # Set the default assay to RNA to run the log-normalization & scaling
  DefaultAssay(seu) <- "RNA"
  
  # Normalization 1: scale counts to 10000 UMIs per cell, and log2-transform the counts
  seu <- seu %>% 
    NormalizeData(normalization.method = "LogNormalize",
                  scale.factor         = 10000) %>% 
    # Identify variable genes
    FindVariableFeatures(selection.method = "vst",
                         nfeatures = config$n_features) %>%
    # Scale data
    ScaleData()
  
  # Normalization 2: using SCTransform
  # This command also identifies variable features and produces scaled values
  seu <- seu %>%
    SCTransform(method = "glmGamPoi", variable.features.n = config$n_features)
  
  # Run PCA at this step, based on scaled data
  seu <- RunPCA(seu,
                pc.genes        = VariableFeatures(.),
                npcs            = config$npcs,
                ndims.print     = 1:5,
                nfeatures.print = 5)
  
  # Return Seurat object with normalized RNA assay
  seu
  
}


#' Determine optimised number of pcs using ElbowPlot
#'
#' This function determine optimised number of pcs of a Seurat object.
#'
#' @param object A Seurat object.
#'
#' @export
#' @return Numeric value with number of optimised pcs.
#' @examples
#' # Basic usage
#' n_pcs <- Elbow_pcs(object = seu)
#' 
#' @author Clara Savary
Elbow_pcs <- function(object, ndims = 50, reduction = "pca"){
  
  # Perform Elbow plot for ranking principle components based on the percentage of variance explained by each one
  ElbowPlot(object, ndims = 50, reduction = reduction)
  
  ## First metric:
  # Determine percent of variation associated with each PC
  pct <- object[[reduction]]@stdev / sum(object[[reduction]]@stdev) * 100
  
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  ## Second metric:
  # Identify the PC where the percent change in variation between consecutive PCs is less than 0.1%:
  
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  
  ## Create a dataframe with values
  plot_df <- data.frame(
    pct  = pct,
    cumu = cumu,
    rank = 1:length(pct)
  )
  
  # Elbow plot to visualize 
  p <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw()
  
  # Plot results
  print(p)
  
  # Print results
  print(pcs)
  
  # Return results
  return(pcs)
  
}


#' Function to perform normalisation of a scATAC assay.
#'
#' This function return perform normalisation of a scATAC assay a Seurat object.
#'
#' @param seu A Seurat object.
#' 
#' @export
#' @return A Seurat object with normalised scATAC assay.
#' @examples
#' # Basic usage
#' seu <- atac_norm(seu)
#' 
#' @author Clara Savary
atac_norm <- function(seu) {
  
  # Set default assay
  DefaultAssay(seu) <- "ATAC"
  
  # Rename cell names using sample id
  sample <- seu@meta.data[["orig.ident"]] %>% levels()
  seu <- RenameCells(object = seu, add.cell.id = glue(sample, "_"))
  
  # Get ATACseq data with union peaks
  newcounts <- FeatureMatrix(fragments = Fragments(seu), features = combined_peaks, cells = colnames(seu))
  #seu[['peakunion']] <- CreateAssayObject(counts = newcounts)
  
  seu[["peakunion"]] <- CreateChromatinAssay(
    counts    = newcounts, 
    sep       = c(":", "-"),
    genome    = seqinfo,
    #fragments = file.path(path, "atac_fragments.tsv.gz"),
    min.cells = config$min_cells
  )
  
  DefaultAssay(seu) <- "peakunion"
  
  # Keep only standard chromosomes coordinates
  genome_id <- seqinfo@genome %>% unique()
  if(genome_id == "hg19"){standard_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)}
  if(genome_id == "hg38"){standard_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)}
  if(genome_id == "mm10"){standard_chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)}
  
  idx_standard_chroms <- which(as.character(seqnames(granges(seu[["peakunion"]]))) %in% standard_chroms)
  seu[["peakunion"]] <- subset(
    seu[["peakunion"]],
    features = rownames(seu[["peakunion"]])[idx_standard_chroms]
  )
  seqlevels(seu[["peakunion"]]@ranges) <- intersect(
    seqlevels(granges(seu[["peakunion"]])),
    unique(seqnames(granges(seu[["peakunion"]])))
  )
  
  # Normalisation and dimensionality reduction of chromatin assay
  seu <- RunTFIDF(seu, method = 1)
  seu <- FindTopFeatures(seu, min.cutoff = 'q0')
  seu <- RunSVD(object = seu, n = 50)
  
  n_pcs <- Elbow_pcs(seu, ndims = 50, reduction = "lsi")
  message(paste0("npcs selected is: ", n_pcs))
  
  seu <- RunUMAP(seu, reduction = "lsi", dims = 2:n_pcs)
  
  # Return Seurat object
  seu
  
}

#' Function to convert genomic coordinate string to a GRanges object.
#'
#' This function to convert genomic coordinate string to a GRanges object of a Seurat object.
#'
#' @param seu A Seurat object.
#' 
#' @export
#' @return A GRanges object with genomic coordinate of a Seurat object.
#' @examples
#' # Basic usage
#' peak_coord(seu)
#' 
#' @author Clara Savary
peak_coord <- function(seu) {
  
  DefaultAssay(seu) <- "ATAC"
  StringToGRanges(regions = rownames(seu))
  
}

# Plotting functions -----------------------------------------------------------

#' Plot density and boxplots of features
#'
#' This function plot density and boxplots of features of a Seurat object.
#'
#' @param seu A Seurat object.
#' @param var A character string indicating feature name to plot.
#' @param xint A numerical value indicating intercept on the x axis.
#'
#' @export
#' @return A plot displaying density and boxplots of features.
#' @examples
#' # Basic usage
#' geom_plot(seu, var = "nFeature_RNA", xint = 200)
#' 
#' @author Clara Savary
geom_plot <- function(seu, var, xint) {
  
  p1 <- seu@meta.data %>% 
    ggplot(aes_string(x = var, fill = "orig.ident")) + 
    geom_density(alpha = 0.1) + 
    scale_x_log10() + 
    geom_vline(xintercept = xint, linetype = "dashed", color = "red") +
    NoLegend()
  
  p2 <- seu@meta.data %>% 
    ggplot(aes_string(x = "orig.ident", y = var, fill = "orig.ident")) + 
    geom_violin() + 
    theme(axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.title      = element_text(hjust = 0.5, face = "bold"),
          axis.title.x    = element_blank(),
          legend.position = "none"
          ) + 
    geom_hline(yintercept = xint, linetype = "dashed", color = "red") 
  
  # Add a title to the combined plot
  var_name <- seu@meta.data[["orig.ident"]] %>% levels()
  
  title <- ggdraw() + 
    draw_label(
      glue(var, " | ", var_name),
      fontface = 'bold',
      x        = 0.5,
      hjust    = 0.5
    )
  
  
  # Plot genes detected per cell
  combined_plot <- cowplot::plot_grid(p1, p2, align = "hv", axis = "tb", rel_widths = c(0.7, 0.3))
  
  # Return combined plot with title
  plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
  
}


#' Plot histogramm of features
#'
#' This function plot histogramm of features of a Seurat object.
#'
#' @param seu A Seurat object.
#' @param var A character string to indicate which feature to plot.
#' @param xint A numerical value indicating intercept on the x axis.
#' 
#' @export
#' @return A plot displaying histogramm of RNA metrics.
#' @examples
#' # Basic usage
#' plot_histo(seu, var = "nFeature_RNA", xint = 200)
#' 
#' @author Clara Savary
plot_histo <- function(seu, var, xint) {
  ggplot(seu@meta.data, aes_string(var)) +
    geom_histogram(aes(fill = orig.ident), alpha = 0.8, bins = 50) +
    geom_vline(xintercept = xint, linetype = "dashed", size = 0.5)
}


#' Plot nGene and nUMI features
#'
#' This function plot nGene and nUMI features of a Seurat object.
#'
#' @param seu A Seurat object.
#' 
#' @export
#' @return A plot displaying nGene and nUMI features of a Seurat object.
#' @examples
#' # Basic usage
#' plot_nGene_nUMI(seu)
#' 
#' @author Clara Savary
plot_nGene_nUMI <- function(seu) {
  
  # Visualize the correlation between genes detected and
  # number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
  var_name <- seu@meta.data[["orig.ident"]] %>% levels()
  
  seu@meta.data %>% 
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent_mito)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    geom_vline(xintercept = config$min_nCount_RNA) +
    geom_hline(yintercept = config$min_nFeature_RNA) +
    ggtitle(glue("nUMI/nGene | ", var_name))
  
}


#' Plot histogramm of RNA metrics
#'
#' This function plot histogramm of RNA metrics of a Seurat object.
#'
#' @param seu A Seurat object.
#' 
#' @export
#' @return A plot displaying histogramm of RNA metrics of a Seurat object.
#' @examples
#' # Basic usage
#' plot_rna_complexity(seu)
#' 
#' @author Clara Savary
plot_rna_complexity <- function(seu) {
  
  var_name <- seu@meta.data[["orig.ident"]] %>% levels()
  
  p1 <- ggplot(seu@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_point() + geom_smooth(method = "lm") + ggtitle(glue("nUMI/nGene | ", var_name))
  p1 <- ggMarginal(p1, type = "histogram", fill = "lightgrey")
  
  p2 <- ggplot(seu@meta.data, aes(x = log10nUMI, y = log10nGene)) + geom_point() + geom_smooth(method = "lm") + ggtitle(glue("nUMI/nGene | ", var_name))
  p2 <- ggMarginal(p2, type = "histogram", fill="lightgrey")
  
  plot_grid(plotlist = list(p1, p2), ncol = 2, align = 'h', rel_widths = c(1, 1))
  
}


#' Plot nucleosome signal
#'
#' This function plot nucleosome signal of a Seurat object.
#'
#' @param seu A Seurat object.
#' 
#' @export
#' @return A plot displaying nucleosome signal of a Seurat object.
#' @examples
#' # Basic usage
#' plot_nucleosome(seu)
#' 
#' @author Clara Savary
plot_nucleosome <- function(seu) {
  
  # Violin plot of nucleosome signal
  p1 <- VlnPlot(seu, "nucleosome_signal", pt.size = 0.1, alpha = 0.1) +
    NoLegend() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
  
  # Deal with inf values
  sum(is.infinite(seu@meta.data$nucleosome_signal)) %>% print()
  seu@meta.data$nucleosome_signal[is.infinite(seu@meta.data$nucleosome_signal)] <- 0
  
  # Define nucleosome groups
  seu$nucleosome_group <- ifelse(seu$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  p2 <- FragmentHistogram(
    object   = seu,
    group.by = "nucleosome_group",
    region   = "chr1-1-10000000"
    ) +
    scale_fill_manual(values = c("forestgreen", "gray80"))
  
  # Add a title to the combined plot
  var_name <- seu@meta.data[["orig.ident"]] %>% levels()
  
  title <- ggdraw() + 
    draw_label(
      glue("Nucleosome signal | ", var_name),
      fontface = 'bold',
      x        = 0.5,
      hjust    = 0.5
    )
  
  # Plot genes detected per cell
  combined_plot <- cowplot::plot_grid(p1, p2, align = "hv", axis = "tb", rel_widths = c(0.3, 0.7))
  
  # Return combined plot with title
  plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
  
}


#' Plot TSS enrichment
#'
#' This function plot TSS enrichment of a Seurat object.
#'
#' @param seu A Seurat object.
#' 
#' @export
#' @return A plot displaying TSS enrichment of a Seurat object.
#' @examples
#' # Basic usage
#' plot_tss(seu)
#' 
#' @author Clara Savary
plot_tss <- function(seu) {
  
  # Plot TSS enrichment
  p1 <- VlnPlot(seu, "TSS.enrichment", pt.size = 0.1, alpha = 0.1) +
    NoLegend() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
  
  # Plot TSS groups
  seu$high.tss <- ifelse(seu$TSS.enrichment > 2, 'High', 'Low')
  p2 <- TSSPlot(seu, group.by = 'high.tss') +
    NoLegend() +
    scale_color_manual(values = c("forestgreen", "gray80"))
  
  # Add a title to the combined plot
  var_name <- seu@meta.data[["orig.ident"]] %>% levels()
  
  title <- ggdraw() + 
    draw_label(
      glue("TSS enrichment | ", var_name),
      fontface = 'bold',
      x        = 0.5,
      hjust    = 0.5
    )
  
  # Plot genes detected per cell
  combined_plot <- cowplot::plot_grid(p1, p2, rel_widths = c(0.3, 0.7))
  
  # Return combined plot with title
  plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
  
}


#' Custom violin plot
#'
#' This function plot features of a Seurat object using custom violin plot.
#'
#' @param object A Seurat object.
#' @param feature A character string to indicate which feature to plot.
#' @param y_int A numerical value to indicate y intercept.
#' @param pt_size A numerical value to indicate point size.
#' 
#' @export
#' @return A plot displaying custom violin plot of Seurat object.
#' @examples
#' # Basic usage
#' vln_custom(object = seu, feature = "nFeature_RNA")
#' 
#' @author Selin Jessa
vln_custom <- function(object, feature, y_int = NULL, pt_size = 0.1, alpha = 0.05) {
  
  p1 <- VlnPlot(object = object, features = feature, pt.size = pt_size, cols = "lightpink", alpha = alpha) +
    NoLegend() +
    ylab(NULL) +
    xlab(NULL) +
    theme(title = element_text(size = 8, face = "plain"),
          axis.text.x = element_blank()
    )
  
  if (!is.null(y_int)) p1 + geom_hline(yintercept = y_int, color = "red") 
  else p1
  
}


#' Plot features and associated thresholds of a Seurat object
#'
#' This function features and associated thresholds of a Seurat object.
#'
#' @param seu A Seurat object.
#' @param title A character string to indicate title of the plot.
#' 
#' @export
#' @return A plot displaying features and associated thresholds of a Seurat object.
#' @examples
#' # Basic usage
#' plot_thresholds(seu)
#' 
#' @author Clara Savary, adapted from Selin Jessa
plot_thresholds <- function(seu, title) {
  
  sample <- seu@meta.data[["orig.ident"]] %>% levels()
  
  title <- ggdraw() + 
    draw_label(
      glue(title, " | ", sample),
      fontface = 'bold',
      x        = 0.5,
      hjust    = 0.5
    )
  
  combined_plot <- plot_grid(
    vln_custom(seu, "nFeature_RNA", c(thresholds_all$min_features[thresholds_all$sample %in% sample], thresholds_all$max_features[thresholds_all$sample %in% sample]), pt_size = 0.1),
    vln_custom(seu, "nCount_RNA", c(thresholds_all$min_umi[thresholds_all$sample %in% sample], thresholds_all$max_umi[thresholds_all$sample %in% sample]), pt_size = 0.1),
    vln_custom(seu, "percent_mito", c(thresholds_all$min_mito[thresholds_all$sample %in% sample], thresholds_all$max_mito[thresholds_all$sample %in% sample]), pt_size = 0.1),
    vln_custom(seu, "nCount_ATAC", pt_size = 0.1),
    vln_custom(seu, "nFeature_ATAC", pt_size = 0.1),
    vln_custom(seu, "atac_peak_region_fragments", c(thresholds_all$min_peak[thresholds_all$sample %in% sample], thresholds_all$max_peak[thresholds_all$sample %in% sample]), pt_size = 0.1),
    vln_custom(seu, "TSS.enrichment", c(thresholds_all$min_TSS[thresholds_all$sample %in% sample], thresholds_all$max_TSS[thresholds_all$sample %in% sample]), pt_size = 0.1),
    vln_custom(seu, "nucleosome_signal", c(0, thresholds_all$max_nucleosome_signal[thresholds_all$sample %in% sample]), pt_size = 0.1),
    ncol  = 8,
    align = "h",
    axis  = "tb"
  )
  
  # Return combined plot with title
  plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
  
}


#' Plot sequencing depth correlation
#'
#' This function plot sequencing depth correlation of a Seurat object.
#'
#' @param seu A Seurat object.
#' @param seu A character string indicating assay of Seurat object.
#' 
#' @export
#' @return A plot displaying sequencing depth correlation of a Seurat object.
#' @examples
#' # Basic usage
#' DepthCor(seu)
#' 
#' @author Clara Savary
plot_depthcor <- function(seu, assay = "peaks") {
  
  DefaultAssay(seu) <- assay
  sample <- seu@meta.data[["orig.ident"]] %>% levels()
  DepthCor(seu) & ggtitle(sample)
  
}


#' Plot UMAP and gene expression levels
#'
#' This function plot UMAP and gene expression levels of a Seurat object.
#'
#' @param reduction A character string indicating which reduction to plot.
#' @param group.by A character string indicating which variable use to color samples.
#' 
#' 
#' @export
#' @return A plot displaying UMAP and gene expression levels of a Seurat object.
#' @examples
#' # Basic usage
#' umap_viz(reduction = "pca", group.by = "orig.ident")
#' 
#' @author Clara Savary
umap_viz <- function(reduction, group.by = "orig.ident"){
  
  p1 <- DimPlot(seu_combined, group.by = group.by, reduction = reduction) & NoAxes()
  p2 <- FeaturePlot(seu_combined, genes, reduction = reduction) & NoAxes() & NoLegend()
  
  p1 / p2
  
}
