# Script downloaded on the May 4th, 2023
# at https://zenodo.org/record/7221864#.ZFLtVuzP30o
# adapted by Clara Savary on the May 4th, 2023
# from the work of Selin Jessa | https://github.com/sjessa

### START
# Following the data structures vignette at
# https://github.com/timoast/signac/blob/master/vignettes/data_structures.Rmd

# Load libraries
library(here)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(Signac)
library(GenomicRanges)
library(AnnotationHub)
# BiocManager::install("biovizBase")

# Get seqinfo
seqinfo <- Seqinfo("mm10")

# List all EnsDb entries from AnnotationHub.
ah <- AnnotationHub()
query(ah, "EnsDb")

# Convert EnsDb to GRanges
gene.ranges <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

## Making a "short cut"
edb <- EnsDb.Mmusculus.v79

## Print some informations for this package
edb
# EnsDb for Ensembl:
# Backend: SQLite
# Db type: EnsDb
# Type of Gene ID: Ensembl Gene ID
# Supporting package: ensembldb
# Db created by: ensembldb package from Bioconductor
# script_version: 0.3.0
# Creation time: Thu May 18 13:38:26 2017
# ensembl_version: 79
# ensembl_host: localhost
# Organism: mus_musculus
# taxonomy_id: 10090
# genome_build: GRCm38
# DBSCHEMAVERSION: 2.0
# No. of genes: 43629.
# No. of transcripts: 104129.
# Protein data available.

# Convert to UCSC style
seqlevelsStyle(gene.ranges) <- "UCSC"
genome(gene.ranges) <- "mm10"

annotation <- gene.ranges

save(seqinfo, annotation, file = "data/", "genome_info/mm10_genome_info.Rda")

seqinfo <- Seqinfo(genome = "mm10")
save(seqinfo, annotation, file = "/Users/cs370/Documents/Bioinformatics/Github/scMultiome_pipelines/data/genome_info/mm10_seqinfo.Rda")

# save(seqinfo, annotation, file = "/scratch/rs55/cs4309/PPK_COMBO/scMultiome_pip/data/genome_info/mm10_genome_info.Rda")

# Promoter coordinates
gene.ranges <- genes(EnsDb.Mmusculus.v79)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

promoter.coords <- promoters(gene.ranges, upstream = 2500, downstream = 2500)

save(promoter.coords, file = here("data/", "genome_info/mm10_promoter_coordinates.Rda"))
# save(seqinfo, annotation, file = "/scratch/rs55/cs4309/PPK_COMBO/scMultiome_pip/data/genome_info/mm10_promoter_coordinates.Rda")

### END