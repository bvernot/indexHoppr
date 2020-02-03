
# options(error=recover)
# options(error=NULL)



## I removed this code because it is nonstandard and also *really* slow
if (F) {
    ## install necessary packages
    cat('Installing new libraries if necessary (may take a few minutes)\n')
    list.of.packages <- c("devtools", "argparse")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("qvalue")
}

cat('Parsing arguments (may take a few minutes)\n')

library(argparse)


# create parser object
parser <- ArgumentParser()

parser$add_argument("-v", "--verbose", action="store_true", default=F,
                    help="Print extra output")
parser$add_argument("-local", "--local", action="store_true", default=F,
                    help="Load the library locally, mostly for testing new features.")
parser$add_argument("-nc", "--ncores", type='integer', default=1,
                    help="Number of cores to use.")
parser$add_argument("-sources", "--sources", type='integer', default=150,
                    help="Number of potential sources of contamination [~2x the number of expected libraries is a good value].")
parser$add_argument("-nlibs", "--num-libs", type='integer', default=NULL,
                    help="Number of libraries to search for contamination.")
parser$add_argument("-rcf", "--random-contam-factor", type='integer', default=5,
                    help="Number of contamination libraries to consider (as factor of -sources). Increase this if the smallest p-values are all 1/(number of contamination comparisons).")
parser$add_argument("-num-plot-libs", "--num-plot-libs", type='integer', default=10,
                    help="Number of potential contaminated libraries to plot.")
parser$add_argument("-splits", "--splits", required=T,
                    help="Splitsfile")
parser$add_argument("-libs", "--plot-libs", required=F, default=NULL,
                    help="One (or more, in the future) library to plot potential contaminations.")
parser$add_argument("-fk", "--fresh-kills", required=F,
                    help='Fresh kills json. You can download a new fresh_kills file with: curl "http://bioaps01:5984/default/_all_docs?include_docs=true" > freshkills.json')
parser$add_argument("-table", "--table", required=F, help="Name of output table. If not provided, no table is saved.")
parser$add_argument("-plots", "--plots", required=F, help="Name of output plots. If not provided, no plots are saved.")
parser$add_argument("-dump-data", "--dump-data", required=F, help="Name of file to save final data, in RDS. If not provided, no RDS is saved.")

if (interactive()) {
  args <- parser$parse_args(strsplit('--splits ~/GoogleDrive/debug_contamination_mt/splittingstats_190211_M06210_B24406_MTcapHuman -libs A17273 -nlibs 10 --plots test_plots_whatwhat.pdf --table test_table_whatwhat.pdf -nc 2 --sources 30', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('--splits ~/Documents/index_cross_contam/data/ludovic/splitstats_ludovic_orlando_001.myformat2.txt -nlibs 30 --plots test_plots_whatwhat.pdf --table test_table_whatwhat.pdf -nc 2 --sources 150 --local', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('--splits ~/GoogleDrive/debug_contamination_mt/splittingstats_190211_M06210_B24406_MTcapHuman.txt -nlibs 10 --plots test_plots_whatwhat.pdf --table test_table_whatwhat.pdf -nc 2 --sources 5 --local', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('--splits ~/GoogleDrive/debug_contamination_mt/splittingstats.190925_D00829_0282.lane2_B31560_MTcapHuman.txt -nlibs 50 --plots test_plots_whatwhat.pdf --table test_table_whatwhat.pdf -nc 10 --sources 200 --local', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('--splits ~/GoogleDrive/debug_contamination_mt/splittingstats.190925_D00829_0282.lane2_B31560_MTcapHuman.txt -nlibs 5 --plots test_plots_whatwhat.pdf --table test_table_whatwhat.pdf -nc 10 --sources 200 --local -fk ~/Downloads/freshkills_200123.json', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('--splits ~/Downloads/180817_B16399_MTcapAllmam_splittingstats.txt -nlibs 20 -lib A13792 --plots test_plots_whatwhat.pdf --table test_table_whatwhat.pdf -nc 10 --sources 200 --local -fk ~/Downloads/freshkills_200123.json', split = ' ')[[1]])
} else {
  # cat('hey\n')
  args <- parser$parse_args()
}

# Usage:
#
# Rscript debug_splitstats_by_indices.R ncores n_sources_of_contam splitsfile file_tag [fresh_kills_json]
#
# You can download a new fresh_kills file with:
# curl "http://bioaps01:5984/default/_all_docs?include_docs=true" > freshkills.json

cat('Loading libraries (may take a few minutes)\n')

library(devtools)

if (args$local) {
  # devtools::install()
  devtools::load_all()
} else {

    ## x = 'indexHoppr'
    ## x = 'index_cross_contam'
    ## if (!require(x,character.only = TRUE))
    ## {
    ##                                     # install.packages(x,dep=TRUE)
    ##     install_github('bvernot/index_cross_contam')
    ##     if(!require(x,character.only = TRUE)) stop("Package not found")
    ## }

    tryCatch({
        library(indexHoppr)
    },
    error = function(err) {
        cat('\n\n')
        cat('You must first install indexHoppr. Start R, and run:\n')
        cat('library(devtools)\ninstall_github("bvernot/index_cross_contam")\n')
        cat('\n\n')
        if (!interactive()) q()
        # stop(err)
    })
}
library(data.table)
setDTthreads(1)
library(doParallel)

## library(data.table)
## library(tidyverse)
## library(ggrepel)
## library(colorspace)
## library(foreach)
## library(doParallel)


# setwd('~/Google Drive/debug_contamination_mt/')
# basedir <- '~/Google Drive/debug_contamination_mt/'
# file_tag <- 'plate000'
# splitsfile <- '~/Google Drive/debug_contamination_mt/splittingstats_190211_M06210_B24406_MTcapHuman.txt'
# fresh_kills_json <- '~/Documents/soil_dna_capture/fresh_kills_2019_02_27.clean.json'





# ncores <- args$ncores
# n_sources_of_contam <- args$sources
# n.hits <- args$num_libs
# splitsfile <- args$splits
# file_tag <- args$prefix
# fresh_kills_json <- args$fresh_kills

if (is.null(args$fresh_kills)) {
  # use_fk <- F
  dt.fresh_kills <- NULL
} else {
                                        # use_fk <- T
    cat('Reading fresh kills database (may take a few minutes)\n')
  dt.fresh_kills <- load_fresh_kills(args$fresh_kills)
}


registerDoParallel(cores=args$ncores)
getDoParWorkers()


# print(args)
# print(length(args))



## allow the reading of multiple runs, and tag those runs (so we can see which wells come from
## which runs - but then have to deal w/ indices that repeat across runs..)
## q: why are some rows and/or columns seemingly swapped?
## q: the bands across the plate: are they the same indices across all plates?
## if so, then maybe it's some issue w/ contaminated indices from the beginning (like p5.431 or
## whatever is just present in all other p5 indices?)
## but then.. why is there this weird double band that breaks? wouldn't you expect that to also be
## present in all combinations?




# splits.mt.181009 <- load_splitstats('~/Google Drive/debug_contamination_mt/debug_contamination_mt/181009_B18927_MTcapHuman_splittingstats.txt', compute.swaps = T)
# splits.mt <- load_splitstats('~/Google Drive/debug_contamination_mt/splittingstats_190211_M06210_B24406_MTcapHuman.txt', compute.swaps = T)
# splits.mt.22401 <- load_splitstats('~/Google Drive/debug_contamination_mt/debug_contamination_mt/190208_M06210_B24401_MTcapHuman_splittingstats.txt', compute.swaps = T)
# splits.mt.22401 <- load_splitstats('~/Google Drive/debug_contamination_mt/debug_contamination_mt/190208_M06210_B24401_MTcapHuman_splittingstats.txt', compute.swaps = T)
# splits.shotgun <- load_splitstats('~/Google Drive/debug_contamination_mt/shotgun_190129_M02279_0409_000000000-D57YK_SN_A20417_split_stats.txt', compute.swaps = T)
#
#
# my.splits <- splits.mt.181009
# my.splits <- splits.mt
# my.splits <- splits.mt.22401
# my.splits <- splits.shotgun

# hey

cat('Reading split stats\n')

dt.idx <- fread(args$splits)
my.splits <- load_splitstats(dt.idx = dt.idx, n_contam_sources = args$sources, limit_search_libs = args$num_libs, compute.swaps = T, dt.fresh_kills = dt.fresh_kills)

cat('Computing corners (may take a long time)\n')
my.splits <- run_compute_swaps(my.splits, random_contam_factor = args$random_contam_factor)

# test.ids <- my.splits$test.ids
# target.ids <- my.splits$target.ids


file_tag <- 'whatwhat'



if (!is.null(args$table)) {
    ##table_file <- sprintf('test_table_%s.tsv', file_tag)
    full.filename = paste0(args$table, '.full.txt')
    cat('\n\nsaving full results table:', args$table, '\n')
    fwrite(my.splits$dt.swaps.test, args$table, sep = '\t')


    # summary.filename = paste0(args$table, '.summary.txt')
    summary.filename = args$table
    cat('\n\nsaving summary results table:', summary.filename, '\n')
    x.test = rbind(my.splits$dt.swaps.test[, .(test.stat, empirical.p, empirical.p.bonf, empirical.q, lib = id1.RG)],
                   my.splits$dt.swaps.test[, .(test.stat, empirical.p, empirical.p.bonf, empirical.q, lib = id2.RG)])
    x.test.summary = x.test[, .(test.stat = max(test.stat),
                                empirical.p = min(empirical.p),
                                empirical.p.bonf = min(empirical.p.bonf),
                                empirical.q = min(empirical.q)), lib]
    fwrite(x.test.summary, summary.filename, sep = '\t')
}

if (!is.null(args$plots)) {
    ##plot_file <- sprintf('test_plots_%s.pdf', file_tag)
    cat('\n\nsaving plot:', args$plots, '\n')
    pdf(args$plots, width=14, height=10)
    plot_debug_splits(my.splits = my.splits, n.hits = args$num_plot_libs, dt.fresh_kills = dt.fresh_kills, plot.libs = args$plot_libs)
    dev.off()
}

if (!is.null(args$dump_data)) {
    ##plot_file <- sprintf('test_plots_%s.pdf', file_tag)
    cat('\n\nsaving RDS:', args$dump_data, '\n')
    my.splits$args <- args
    saveRDS(my.splits, args$dump_data)
}
# pdf(sprintf('test_plots_%s_small.pdf', file_tag), width=7, height=5)
# plot_debug_splits(my.splits, args$num_libs, plot.libs = args$plot_libs)
# dev.off()




## make a plot w/ the flip of the matrix.. like if two samples are mixed together, you would expect
## to see p5a w/ p7b and p7b with p5a.  so then the x value is all the points on the upper triangle
## of the matrix, and the y axis is the corresponding point on the lower triangle.  If this is
## random, then there shouldn't be a correlation between the two.

## also, looks like there are some rows/cols that have low levels of pairing w/ all other p5 or p7
## indices. how to quantify that?
