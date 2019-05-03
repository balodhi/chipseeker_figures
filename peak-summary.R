#!/usr/bin/env Rscript

# Run the ChIPSeeker pipeline on .BED files


# ~~~~~ PACKAGES ~~~~~ # 
library("optparse")
library("tools")



# ~~~~~ FUNCTIONS ~~~~~ # 
msprintf <- function(fmt, ...) {
    message(sprintf(fmt, ...))
}


make_filename <- function (input_file, new_ext, out_dir = FALSE) {
    # Convert '/path/to/file.bed' to '/path/to/file_annotations.tsv'
    old_ext <- file_ext(input_file)
    filename_base <- gsub(pattern = sprintf('.%s$', old_ext), replacement = '', x = basename(input_file))
    filename_new <- sprintf('%s.%s', filename_base, new_ext)
    new_path <- file.path(dirname(input_file), filename_new)
    if(out_dir != FALSE){
        new_path <- file.path(out_dir, new_path)
        dir.create(path = dirname(new_path), recursive = TRUE, showWarnings = FALSE)
    }
    return(new_path)
}

check_numlines <- function(input_file, min_value = 0) {
    # make sure a file has >0 lines
    has_enough_lines <- FALSE
    if (length(readLines(input_file)) > min_value) has_enough_lines <- TRUE
    return(has_enough_lines)
}


validate_file <- function(input_file) {
    # make sure that all files are .bed, and that they have >0 lines
    # validation passes if all files are .bed
    all_exist <- all(file.exists(input_file))
    if ( ! isTRUE(all_exist)) {
        msprintf("WARNING: Input file do not exist:\n%s\nFile will not be processed\n\n", input_file)
        return(FALSE)
    }
    all_bed <- all(grepl(pattern = '*.bed$', x = basename(input_file)))
    if ( ! isTRUE(all_bed)) {
        msprintf("WARNING: Input file is not .bed:\n%s\nFile will not be processed\n\n", input_file)
        return(FALSE)
    }
    all_min_linenum <- all(sapply(input_file, check_numlines))
    if ( ! isTRUE(all_min_linenum)) {
        msprintf("WARNING: Input file does not have enough lines:\n%s\nFile will not be processed\n\n", input_file)
        return(FALSE)
    }
    return(TRUE)
}



find_all_beds <- function (input_dirs) {
    # find all .bed files in the supplied dirs
    return(dir(input_dirs, pattern = '.bed', full.names = TRUE, recursive = TRUE))
}


get_sampleID <- function(input_file, id_dirname = FALSE){
    # get the sample ID for a file
    # right now just use the basename but maybe some day do something fancier here
    sampleID <- basename(input_file)
    if(isTRUE(id_dirname)) sampleID <- basename(dirname(input_file))
    return(sampleID)
}

get_sample_outdir <- function(parent_outdir, sampleID, create = TRUE){
    # make a path for the sample's output directory
    output_path <- file.path(parent_outdir, sampleID)
    if(isTRUE(create)) dir.create(output_path, recursive = TRUE)
    return(output_path)
}

chipseeker_CompPipeline <- function(bed_file, sampleID, tss_dist, txdb, out_dir = FALSE, annoDb = "org.Ce.eg.db"){
	files = readLines(bed_file)
	NAMES  <- paste0("Peak",1:length(files))
	names(files) <- NAMES
	
	peaks_coverage_plot_file <- make_filename(input_file = bed_file, new_ext = 'Genomic_annotations.pdf', out_dir = out_dir)
  msprintf("Making Chrom Comp Coverages plot:\n%s\n\n", peaks_coverage_plot_file)
  sample_title <- paste0(sampleID, "Genomic Annotations")
  pdf(file = peaks_coverage_plot_file)
  peakAnnoList <- lapply(files,annotatePeak,tssRegion=c(-tss_dist, tss_dist),TxDb=txdb, annoDb = annoDb)
  print(plotAnnoBar(peakAnnoList)) # title = "Genomic Annotations"
  dev.off()
  
  peaks_coverage_plot_file <- make_filename(input_file = bed_file, new_ext = 'bingdings.pdf', out_dir = out_dir)
  msprintf("Making Chrom Comp Coverages plot:\n%s\n\n", peaks_coverage_plot_file)
  sample_title <- paste0(sampleID, "Distribution of Binding Sites among different ChIPseq data")
  pdf(file = peaks_coverage_plot_file)
  # peakAnnoList <- lapply(files,annotatePeak,tssRegion=c(-tss_dist, tss_dist),TxDb=txdb, annoDb = annoDb)
  print(plotDistToTSS(peakAnnoList)) # title = "test"
  dev.off()
  
  peaks_coverage_plot_file <- make_filename(input_file = bed_file, new_ext = 'pathway.pdf', out_dir = out_dir)
  msprintf("Making Chrom Comp Coverages plot:\n%s\n\n", peaks_coverage_plot_file)
  sample_title <- paste0(sampleID, "Pathway Enrichment Analysis")
  pdf(file = peaks_coverage_plot_file)
  # peakAnnoList <- lapply(files,annotatePeak,tssRegion=c(-tss_dist, tss_dist),TxDb=txdb, annoDb = annoDb)
  genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
  names(genes) = sub("_", "\n", names(genes))
  require(ReactomePA)
  comppathway <- compareCluster(geneCluster = genes, fun = "enrichPath", pvalueCutoff=1, organism = "celegans")
  print(dotplot(comppathway, showCategory = 20)) # title = "test"
  dev.off()
}

# org.Hs.eg.db
chipseeker_pipeline <- function(bed_file, sampleID, tss_dist, txdb, out_dir = FALSE, annoDb = "org.Ce.eg.db"){
    # the pipeline for ChIPSeeker peak annotations and plots
    # peak <- readPeakFile("/niv/ChIP-Seq/Advanced_Analysis/CE20 F0.vs.CE23 F0/DBA_20_F0.vs.23_F0.report_full.bed")
    # print(covplot(peak, weightCol = "V4", title = "adsfasdfasdfa"))
    msprintf("Reading peaks file...\n\n")
    peak <- readPeakFile(bed_file)
    
    peaks_coverage_plot_file <- make_filename(input_file = bed_file, new_ext = 'coverage.pdf', out_dir = out_dir)
    msprintf("Making Chrom Coverages plot:\n%s\n\n", peaks_coverage_plot_file)
    sample_title <- paste0(sampleID, " ChIP Peaks over Chromosomes")
    pdf(file = peaks_coverage_plot_file)
    print(covplot(peak, weightCol = "V5", title = sample_title)) # title = "ChIP Peaks over Chromosomes"
    dev.off()
    
    msprintf("Getting peak annotations...\n\n")
    peakAnno <- annotatePeak(peak, tssRegion = c(-tss_dist, tss_dist), 
                             TxDb = txdb, 
                             annoDb = annoDb)
    
    peak_anno_table_file <- make_filename(input_file = bed_file, new_ext = 'peak_anno.tsv', out_dir = out_dir)
    msprintf("Saving table:\n%s\n\n", peak_anno_table_file)
    write.table(peakAnno, quote=FALSE, sep="\t", row.names =FALSE, file=peak_anno_table_file)
    
    peak_anno_stats_file <- make_filename(input_file = bed_file, new_ext = 'peak_anno_stats.tsv', out_dir = out_dir)
    msprintf("Saving table:\n%s\n\n", peak_anno_stats_file)
    write.table(peakAnno@annoStat, quote=FALSE, sep="\t", row.names =FALSE, file=peak_anno_stats_file)
    
    tss_dist_file <- make_filename(input_file = bed_file, new_ext = 'tss_distance.txt', out_dir = out_dir)
    msprintf("Saving table:\n%s\n\n", tss_dist_file)
    cat(as.character(tss_dist), file = tss_dist_file)
    
    
    anno_piechart_plot_file <- make_filename(input_file = bed_file, new_ext = 'anno-piechart.pdf', out_dir = out_dir)
    msprintf("Making Peak Anno pie chart:\n%s\n\n", anno_piechart_plot_file)
    sample_title <- paste0("\n\n", sampleID, " Peak Types")
    pdf(file = anno_piechart_plot_file, height = 8, width = 8)
    print(plotAnnoPie(peakAnno, main = sample_title))
    dev.off()
    
    anno_bar_plot_file <- make_filename(input_file = bed_file, new_ext = 'anno-barplot.pdf', out_dir = out_dir)
    msprintf("Making Peak Anno pie chart:\n%s\n\n", anno_bar_plot_file)
    sample_title <- paste0("\n\n", sampleID, " Peak Types")
    pdf(file = anno_bar_plot_file, height = 4, width = 12)
    print(plotAnnoBar(peakAnno, main = sample_title))
    dev.off()

    msprintf("Making Upset plot...\n\n")
    # upset_plot_file <- file.path(output_directory, sprintf("upsetplot.pdf", sampleID))
    upset_plot_file <- make_filename(input_file = bed_file, new_ext = 'upsetplot.pdf', out_dir = out_dir)
    sample_title <- paste0(sampleID, " Peak Overlaps")
    pdf(file = upset_plot_file, width = 9, height = 4.5, onefile = F)
    print(upsetplot(peakAnno, vennpie=TRUE))
    text(x = 0, y = 1, sample_title) # add a title
    dev.off()
    
    
    extra_summarize(sampleID,peak, peakAnno, bed_file, out_dir, txdb)
}
extra_summarize <- function(sampleID,peak, peakAnno, bed_file, out_dir ,txdb) {
    suppressPackageStartupMessages(library(ReactomePA))
    pathway_file <- make_filename(input_file = bed_file, new_ext = 'pathway.pdf', out_dir = out_dir)
    gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance= 3000, TxDb=txdb)
    pathway <- enrichPathway(gene, organism = "celegans",pvalueCutoff=1)

    
    msprintf("Making pathway chart:\n%s\n\n", pathway_file)
    sample_title <- paste0("\n\n", sampleID, " Peak Types")
    pdf(file = pathway_file, height = 12, width = 12)
    print(dotplot(pathway), main = sample_title)
    dev.off()
    
 #   tss_file <- make_filename(input_file = bed_file, new_ext = 'tssregion.pdf', out_dir = out_dir)
 #   msprintf("Making pathway chart:\n%s\n\n", tss_file )
 #   sample_title <- paste0("\n\n", sampleID, " Peak Types")
 #   pdf(file = tss_file , height = 12, width = 12)
 #   print(peakHeatmap(bed_file, TxDb=txdb, upstream=3000, downstream=3000, color="red"), main = sample_title)
 #   dev.off()

    tfbindings_file <- make_filename(input_file = bed_file, new_ext = 'tfbindings.pdf', out_dir = out_dir)
    msprintf("Making tf bindings:\n%s\n\n", tfbindings_file )
    sample_title <- paste0("\n\n", sampleID, " Peak Types")
    pdf(file = tfbindings_file , height = 4, width = 10)
    print(plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS"), main = sample_title)
    dev.off()
    
 }

summarize_beds <- function(bed_files, tss_dist, id_dirname = FALSE, out_dir = FALSE) {
    # run the ChIPSeeker pipeline on all the .bed files
    
    
    # ~~~~~ VALIDATION ~~~~~ # 
    # check to make sure at least one files has >0 lines before we try to load data, because it takes a while to load
    any_min_linenum <- any(sapply(names(bed_files), check_numlines))
    if ( ! isTRUE(any_min_linenum)) {
        msprintf("ERROR: No input files have enough lines to be processed\nExiting...\n\n")
        quit()
    }
    
    # ~~~~~ LOAD DATA ~~~~~ # 
    message("\nLoading packages and data...\n")
    # source("http://bioconductor.org/biocLite.R")
    # biocLite("ChIPseeker")
    suppressPackageStartupMessages(library("ChIPseeker"))
    suppressPackageStartupMessages(library("clusterProfiler"))
    #suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
    suppressPackageStartupMessages(library("BSgenome.Celegans.UCSC.ce10"))
    suppressPackageStartupMessages(library("GenomicFeatures"))
    # txdb <- get("TxDb.Hsapiens.UCSC.hg19.knownGene")
    txdb <- makeTxDbFromUCSC("ce10", "refGene")
    # save(txdb, file = "TxDb.Hsapiens.UCSC.hg19.knownGene.Rdata")
    
    
    # ~~~~~ RUN ~~~~~ # 
    # iterate over bed files
    msprintf('\n------------------------------\n')
    msprintf('\n------------------------------\n')
    for(i in seq_along(bed_files)){
        bed_file <- names(bed_files[i])
        process_file <- bed_files[i] # TRUE or FALSE
        
        files_errors <- character()
        files_warnings <- character()
        #process_file = as.logical(TRUE)
        msprintf("Input File:\n%s\n\n\nFile will be processed:\n%s\n\n", bed_file, process_file)
        if(isTRUE(as.logical(process_file))){
            
            sampleID <- get_sampleID(input_file = bed_file, id_dirname = id_dirname)
            
            msprintf("Sample ID:\n%s\n\n\n", sampleID)
            
            result <- tryCatch(
                {
                    msprintf("Running ChIPSeeker pipeline for sample %s, file:\n%s\n\n", sampleID, bed_file)
                    chipseeker_pipeline(bed_file = bed_file, sampleID = sampleID, tss_dist = tss_dist, txdb = txdb, out_dir = out_dir,annoDb="org.Ce.eg.db")
                    #chipseeker_CompPipeline(bed_file = bed_file, sampleID = sampleID, tss_dist = tss_dist, txdb = txdb, out_dir = out_dir,annoDb="org.Ce.eg.db")
                },
                error = function(cond) {
                    msprintf("An error occured while running ChIPSeeker pipeline for sample %s, file:\n%s\n\n", sampleID, bed_file)
                    message("Original error message:")
                    message(cond)
                    return("error")
                },
                warning = function(cond) {
                    msprintf("An warning occured while running ChIPSeeker pipeline for sample %s, file:\n%s\n\n", sampleID, bed_file)
                    message("Original warning message:")
                    message(cond)
                    return("warning")
                },
                finally={
                    msprintf("Finished running ChIPSeeker pipeline for sample %s, file:\n%s\n\n", sampleID, bed_file)
                }
            )
            
            if(result == "error") files_errors <- c(files_errors, bed_file)
            if(result == "warning") files_warnings <- c(files_warnings, bed_file)
        }
        msprintf('\n------------------------------\n')
    }
    
    msprintf('The following files had errors:\n')
    msprintf('%s\n', files_errors)
    
    msprintf('The following files had warnings:\n')
    msprintf('%s\n', files_warnings)
    
    cat(files_errors, file = "file_errors.txt", append = TRUE)
    cat(files_warnings, file = "file_warnings.txt", append = TRUE)
}



# ~~~~~ SCRIPT ARGS ~~~~~ # 
option_list <- list(
    make_option(c("-d", "--dir"), action="store_true", default=FALSE,
                dest="dir_mode", help="Treat input items as directories to be searched for .bed files"),
    make_option(c("--id-dirname"), action="store_true", default=FALSE,
                dest="id_dirname", help="Take the sample ID from the .bed file's dirname, not its basename"),
    make_option(c("--out-dir"), type="character", default=FALSE,
                dest="out_dir", help="Path to the parent output directory. Will be created and appended to the input item's file path"),
    make_option(c("--tss-dist"), type="numeric", default=3000,
                dest = "tss_dist", help="TSS distance to use [default %default]",
                metavar="tss-dist")
)
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

dir_mode <- opt$options$dir_mode 
out_dir <- opt$options$out_dir
tss_dist <- opt$options$tss_dist
id_dirname <- opt$options$id_dirname

input_items <- opt$args

# get script dir
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
save.image(file.path(scriptPath, "loaded_args.Rdata"))
# quit()



# ~~~~~ RUN ~~~~~ #
# default output dir
if (isTRUE(dir_mode)) input_items <- find_all_beds(input_items)

msprintf('Input Items are:\n')
msprintf('%s\n', input_items)

validated_items <- sapply(input_items, validate_file)

summarize_beds(bed_files = validated_items, tss_dist = tss_dist, id_dirname = id_dirname, out_dir = out_dir)
