### options
dynamic_require <- function(package) {
  if(eval(parse(text = paste("require(",package,")"))))
    return(T)
  
  install.packages(package)
  return(eval(parse(text=paste("require(",package,")"))))
}
suppressPackageStartupMessages(dynamic_require("argparse"))

# create parser object
parser <- ArgumentParser()

parser$add_argument("-v", "--verbose", action="store_true", default=F,
                    help="Print extra output")
parser$add_argument("-nc", "--ncores", type='integer', default=1,
                    help="Number of cores to use.")
parser$add_argument("-sources", "--sources", type='integer', default=150,
                    help="Number of potential sources of contamination [~2x the number of expected libraries is a good value].")
parser$add_argument("-nhits", "--num-hits", type='integer', default=30,
                    help="Number of potential contaminated libraries to plot")
parser$add_argument("-splits", "--splits", required=T,
                    help="Splitsfile")
parser$add_argument("-libs", "--plot-libs", required=F, default=NULL,
                    help="One (or more, in the future) library to plot potential contaminations.")
parser$add_argument("-fk", "--fresh-kills", required=F,
                    help='Fresh kills json. You can download a new fresh_kills file with: curl "http://bioaps01:5984/default/_all_docs?include_docs=true" > freshkills.json')
parser$add_argument("-prefix", "--prefix", required=T,
                    help="Prefix for output files.")

if (interactive()) {
  args <- parser$parse_args(strsplit('--splits ~/Documents/index_cross_contam/test_splits.txt -libs A17273 -nhits 10 --prefix what -nc 2 --sources 30', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('--splits ~/Documents/index_cross_contam/data/ludovic/splitstats_ludovic_orlando_001.myformat2.txt -nhits 30 --prefix what -nc 2 --sources 150', split = ' ')[[1]])
} else {
  args <- parser$parse_args()
}


# Usage:
# 
# Rscript debug_splitstats_by_indices.R ncores n_sources_of_contam splitsfile file_tag [fresh_kills_json]
# 
# You can download a new fresh_kills file with:
# curl "http://bioaps01:5984/default/_all_docs?include_docs=true" > freshkills.json



dynamic_require("data.table")
dynamic_require("tidyverse")
dynamic_require("ggrepel")
dynamic_require("colorspace")
dynamic_require("foreach")
dynamic_require("doParallel")

setDTthreads(1)

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



#basedir <- '/mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/debug_contamination_mt'

if (interactive()) {
  basedir <- '~/Documents/index_cross_contam/bin/'
} else {
    initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  basedir <- dirname(script.name)
}
cat('script', script.name, '\n')
cat('basedir', basedir, '\n')


ncores <- args$ncores
n_sources_of_contam <- args$sources
n.hits <- args$num_hits
splitsfile <- args$splits
file_tag <- args$prefix
fresh_kills_json <- args$fresh_kills
use_fk <- T
if (is.null(fresh_kills_json)) {
    use_fk <- F
    ## fresh_kills_json <- sprintf('%s/freshkills.json', basedir)
}


registerDoParallel(cores=ncores)
getDoParWorkers()


# print(args)
# print(length(args))


## loads functions and dt.fresh_kills - is slow
if (use_fk) source(sprintf('%s/fresh_kills_fns.R', basedir))
source(sprintf('%s/debug_splitstats_fns.R', basedir))

## allow the reading of multiple runs, and tag those runs (so we can see which wells come from 
## which runs - but then have to deal w/ indices that repeat across runs..)
## q: why are some rows and/or columns seemingly swapped?
## q: the bands across the plate: are they the same indices across all plates?  
## if so, then maybe it's some issue w/ contaminated indices from the beginning (like p5.431 or 
## whatever is just present in all other p5 indices?)
## but then.. why is there this weird double band that breaks? wouldn't you expect that to also be
## present in all combinations?




# rg.cats <- unique(splits.shotgun$dt.idx.top$RG.cat)
rg.cats <- c("expected", "PhiX", "unex_PhiX", "unknown", "in_fk")
rg.cat.colors <- qualitative_hcl(length(rg.cats), "Dark 3")
names(rg.cat.colors) <- rg.cats



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

my.splits <- load_splitstats(splitsfile, n_sources_of_contam, compute.swaps = T)


# test.ids <- my.splits$test.ids
# target.ids <- my.splits$target.ids




## start plotting
# setwd('~/Google Drive/debug_contamination_mt/')
# pdf(sprintf('test_plots_%s.pdf', file_tag), width=14, height=10)

plot_debug_splits <- function(my.splits, n.hits, plot.libs = NULL) {
  
  dt.swaps.bak <- data.table(my.splits$dt.swaps.bak)
  dt.swaps.test <- data.table(my.splits$dt.swaps.test)
  
  nseqs.hist = 10
  dt.idx.hist_plotting <- my.splits$dt.idx.top[RG.cat != 'expected' & nseqs > nseqs.hist]
  nseqs.label = dt.idx.hist_plotting[, min(tail(sort(nseqs),5))]
  if (use_fk) {
    dt.idx.hist_plotting[nseqs >= nseqs.label, desc := paste0(description_from_fresh_kills(roots_from_fresh_kills(p7ind, p5ind)), '_',
                                                              dates_in_fresh_kills(my.id = roots_from_fresh_kills(p7ind, p5ind)),
                                                              collapse = '\n'), .(p7ind, p5ind)]
    dt.idx.hist_plotting[nseqs >= nseqs.label & desc == '_', desc := paste0(p7seq, '_', p5seq)]
  } else {
    dt.idx.hist_plotting[nseqs >= nseqs.label, desc := paste0(p7seq, '_', p5seq)]
  }
  
  p1 <- ggplot(my.splits$dt.idx.top, aes(x=RG.fac, weight=nseqs/1e5, fill = RG.cat)) +
    geom_bar() + ylab('Read counts (100k)') + theme(axis.text.x = element_blank()) + xlab('Top Read Groups') +
    scale_fill_manual(values = rg.cat.colors) +
    # geom_vline(xintercept = 120) +
    NULL
  print(p1)
  
  p1 <- ggplot(my.splits$dt.idx.top[RG.full %in% tail(my.splits$top.RG, 200)], aes(x=RG.fac, weight=nseqs/1e5, fill = RG.cat)) +
    geom_bar() + ylab('Read counts (100k)') + theme(axis.text.x = element_blank()) + xlab('Top Read Groups') +
    scale_fill_manual(values = rg.cat.colors) +
    # geom_vline(xintercept = 120) +
    NULL
  print(p1)
  
  p1 <- ggplot(dt.idx.hist_plotting,
         aes(x=nseqs/1000, fill=RG.cat, weight=nseqs, group=interaction(p5ind, p7ind))) +
    geom_histogram(color=rgb(0,0,0,.4)) + scale_x_log10() + xlab('Number of sequences (x1000, log scale)') +
    ggtitle(sprintf('Unexpected index pairs w/ > %d reads', nseqs.hist)) + ylab('Number of sequences') +
    # geom_vline(data=dt.idx.hist_plotting[!is.na(desc)], aes(color = RG.cat, xintercept = nseqs/1000)) +
    # geom_point(data=dt.idx.hist_plotting[!is.na(desc)], aes(color = RG.cat, y=max(dt.idx.hist_plotting$nseqs))) +
    geom_text_repel(aes(label = desc, y=max(dt.idx.hist_plotting$nseqs)), box.padding = 2, force = 10, max.iter = 200000) +
    scale_fill_manual(values = rg.cat.colors) +
    scale_color_manual(values = rg.cat.colors)
  print(p1)
  
  
  p.tile <- ggplot(my.splits$dt.idx.top.p57, 
                   aes(x=p7.fac.maxseqs, y=p5.fac.maxseqs, fill=log(nseqs))) + 
    geom_tile() + scale_fill_distiller(type='div', palette = 2) +
    geom_abline(slope=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
    NULL
  print(p.tile)
  
  p1 <- p.tile + geom_point(data=dt.idx.hist_plotting[!is.na(desc)], color='black') +
    geom_point(data = my.splits$dt.idx.top[RG.cat == 'in_fk'], color='red') +
    geom_point(data = my.splits$dt.idx.top[RG.cat == 'expected'], color='blue') +
    geom_text_repel(data=dt.idx.hist_plotting, aes(label = desc), box.padding = 1) +
    NULL
  print(p1)  
  
  #### look at and plot swaps / squares
  
  dt.swaps.test[(x1.rg %like% 'A' | x2.rg %like% 'A')]
  dt.swaps.bak[(x1.rg %like% 'A' | x2.rg %like% 'A')]
  dt.swaps.bak[x > 8]
  
  dt.swaps.test[, min.nseqs := min(id1.nseqs, id2.nseqs), .(id1.nseqs, id2.nseqs)]
  dt.swaps.bak[, min.nseqs := min(id1.nseqs, id2.nseqs), .(id1.nseqs, id2.nseqs)]
  
  dt.swaps.hits <- dt.swaps.test[x-y > 5][order(x-y,decreasing = T)]
  
  
  p1 <- ggplot(dt.swaps.test, aes(x=x, y=y)) +
    geom_point(data=dt.swaps.bak, color='red', alpha=.2, size=3) +
    # geom_point(data=dt.swaps.bak[x1.rg %like% 'A' | x2.rg %like% 'A'], color='green', alpha=.2, size=3) +
    # geom_point(data=dt.swaps.test[x1.rg %like% 'A' | x2.rg %like% 'A'], color='green', alpha=.2, size=13, pch='x') +
    # geom_density_2d(data=dt.swaps4, color='black', h = c(5,5), n=80, bins=100) +
    geom_point(alpha=.2, color='black') +
    geom_point(data=dt.swaps.hits, color='black', alpha=1) +
    # geom_point(data=dt.swaps.hits[min.nseqs/(x1.nseqs+x2.nseqs) < 1], color='blue', alpha=.2, size=5) +
    # geom_point(data=dt.swaps2.noA[(x1.nseqs > 0 & x2.nseqs > 0) & prob.rg == T], color='blue') +
    # geom_point(data=dt.swaps2.noA[(my.id1 %like% 'A20371' | my.id2 %like% 'A20371') & (my.id1 %like% 'unknown.441.188' | my.id2 %like% 'unknown.441.188')],
    #            color='green', size=3) +
    # geom_smooth(data=dt.swaps4, color='black') +
    # geom_smooth(data=dt.swaps4, color='black', method='lm') +
    geom_abline(slope=1, intercept = c(-4,-5,-6,-7)) +
    xlab('log(nseqs in corner libraries)') +
    ylab('abs log(corner1/corner2)') +
    # scale_x_log10()
    NULL
  print(p1)

  p1 <- ggplot(dt.swaps.test, aes(x=x, y=y0)) +
      geom_point(data=dt.swaps.bak, color='red', alpha=.2, size=3) +
    geom_point(alpha=.2, color='black') +
    geom_point(data=dt.swaps.hits, color='black', alpha=1) +
    geom_abline(slope=1, intercept = c(-4,-5,-6,-7)) +
    geom_abline(slope=-1, intercept = c(4,5,6,7)) +
    xlab('log(nseqs in corner libraries)') +
    ylab('abs log(corner1/corner2)') +
    NULL
  print(p1)

  print(dt.swaps.test)
  # p1 <- ggplot(dt.swaps.test[, y0, .(score=round(x-y))][score > 4], aes(x=y0)) + geom_histogram() + facet_wrap(~score, scales = 'free_y') + xlab('log(corner1/corner2) [facets are rounded scores]')
  # print(p1)

  qqplot(dt.swaps.bak[, x-y],
         dt.swaps.test[, x-y],
         xlab = 'Randomized score distribution',
         ylab = 'Test score distribution')
  abline(b = 1, a=0)
  
  if (!is.null(plot.libs)) {
    dt.swaps.hits.plot <- dt.swaps.test[my.id1 %like% plot.libs | my.id2 %like% plot.libs][order(x-y,decreasing = T)]
    dt.swaps.hits.plot <- head(dt.swaps.hits.plot, n.hits)
  } else {
    dt.swaps.hits.plot <- head(dt.swaps.test[order(x-y,decreasing = T)], n.hits)
  }
  my.swap <- 1
  for (my.swap in dt.swaps.hits.plot[, .I]) {
    rg.square <- dt.swaps.hits.plot[my.swap, c(my.id1,my.id2,x1.rg,x2.rg)]
    my.id1 <- dt.swaps.hits.plot[my.swap, my.id1]
    my.id2 <- dt.swaps.hits.plot[my.swap, my.id2]
    my.id.subtitle <- 
      dt.swaps.hits.plot[my.swap, sprintf('y-x: %0.2f; nseqs: %d and %d; corners: %d and %d', 
                                     x-y, id1.nseqs, id2.nseqs, x1.nseqs, x2.nseqs)]
    dt.id.label <- my.splits$dt.idx[RG.full == my.id1 | RG.full == my.id2]
    
    if (use_fk) {
        dt.id.label[, desc := paste0(RG, '_', description_from_fresh_kills(roots_from_fresh_kills(p7ind, p5ind)), '_',
                                     dates_in_fresh_kills(my.id = roots_from_fresh_kills(p7ind, p5ind)),
                                     collapse = '\n'), .(p7ind, p5ind)]
        dt.id.label[desc == '__', desc := RG.full]
    } else {
        dt.id.label[, desc := RG.full]
    }

    a = cbind(my.splits$dt.idx.top.p57[RG.full == my.id1, .(xmin = p7.fac.maxseqs, ymin = p5.fac.maxseqs)],
              my.splits$dt.idx.top.p57[RG.full == my.id2, .(xmax = p7.fac.maxseqs, ymax = p5.fac.maxseqs, p7.fac.maxseqs, p5.fac.maxseqs)])
    p1 <- ggplot(my.splits$dt.idx.top.p57,
                 aes(x=p7.fac.maxseqs, y=p5.fac.maxseqs, fill=log(nseqs))) +
      geom_tile() + scale_fill_distiller(type='div', palette = 2) +
      # geom_point(data = my.splits$dt.idx.top[RG.cat == 'expected'], color='blue') +
      geom_abline(slope=1) +
      theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
      geom_rect(data=a, aes(xmin = xmin, ymin = ymin, 
                            xmax = xmax, ymax = ymax), 
                fill = NA, color='black') +
      geom_tile(data=my.splits$dt.idx.top.p57[RG.full %in% rg.square]) +
      geom_text_repel(data=dt.id.label, aes(label = desc), box.padding = 1, point.padding = 1, direction = 'y') +
      ggtitle(sprintf('Putative contaminated samples: %s and %s', my.id1, my.id2),
              subtitle = my.id.subtitle)
    print(p1)
  }
}


pdf(sprintf('test_plots_%s.pdf', file_tag), width=14, height=10)
plot_debug_splits(my.splits, n.hits, plot.libs = args$plot_libs)
dev.off()

pdf(sprintf('test_plots_%s_small.pdf', file_tag), width=7, height=5)
plot_debug_splits(my.splits, n.hits, plot.libs = args$plot_libs)
dev.off()




## make a plot w/ the flip of the matrix.. like if two samples are mixed together, you would expect
## to see p5a w/ p7b and p7b with p5a.  so then the x value is all the points on the upper triangle
## of the matrix, and the y axis is the corresponding point on the lower triangle.  If this is
## random, then there shouldn't be a correlation between the two.

## also, looks like there are some rows/cols that have low levels of pairing w/ all other p5 or p7
## indices. how to quantify that?
