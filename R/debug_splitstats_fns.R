



load_splitstats <- function(split_stats_file = NULL, dt.idx = NULL,
                            n_contam_sources, limit_search_libs = NULL,
                            p57.thresh = 500,
                            compute.swaps = F, filter_fully_unknown = T,
                            dt.fresh_kills = NULL) {
  if (is.null(dt.idx)) {
    dt.idx <- fread(split_stats_file)
  }
  # dt.idx <- fread(split_stats_file)
  setnames(dt.idx, '#seqs', 'nseqs', skip_absent=TRUE)

  if (is.null(dt.fresh_kills)) {
    use_fk <- F
  } else {
    use_fk <- T
  }

  ## add 'indices' for the p5 and p7 indexes
  ## these are strings
  if (!'p7ind' %in% names(dt.idx)) {
    dt.idx <- merge(dt.idx, dt.idx[, .N, p7seq][, .(p7seq, p7ind=as.character(.I))], by='p7seq')
  }
  if (!'p5ind' %in% names(dt.idx)) {
    dt.idx <- merge(dt.idx, dt.idx[, .N, p5seq][, .(p5seq, p5ind=as.character(.I))], by='p5seq')
  }

  # print(head(dt.idx))
  # print(is.data.table(dt.idx))
  # print(names(dt.idx))
  # print(head(dt.idx[, p7ind]))
  if (filter_fully_unknown) dt.idx <- dt.idx[p7ind != '-' & p5ind != '-']


  ## ISSUE : should I check to make sure that the indices are numeric/strings, and the seqs are characters?

  ######
  ## clean up the RG a bit - probably only works w/ MPI splitstats
  ## this also helps us to get the number of expected libraries
  dt.idx[, RG.cat := RG]
  dt.idx[RG == 'unexpected', RG := 'unknown']
  dt.idx[RG != 'PhiX' & (p7ind %like% 'Phi' | p5ind %like% 'Phi'), RG := 'unex_PhiX']
  dt.idx[RG != 'unknown' & RG != 'unexpected' & RG != 'PhiX' & RG != 'unex_PhiX', RG.cat := 'expected']
  dt.idx[p5ind == 'PhiX' & p7ind == 'PhiX', RG := 'PhiX']
  dt.idx[p5ind == 'PhiX' & p7ind == 'PhiX', RG.cat := 'PhiX']
  cat('Found N expected libraries:', dt.idx[RG.cat == 'expected', .N], '\n')

  ######
  ## if we don't specify the number of contaminants, use 2x the number of expected libraries
  if (is.null(n_contam_sources)) n_contam_sources <- 2 * dt.idx[RG.cat == 'expected', .N]

  ######
  ## filter down to just p5/p7 that are present in the top X rg
  ## so this is more than just the top X rg, this is all rg w/ p5/p7 that are found in the top X rg
  setorder(dt.idx, -nseqs)
  # dt.idx[, top.rg := .I <= n_target_libs]
  dt.idx[, top.rg := .I <= as.integer( 5 * max(n_contam_sources, limit_search_libs) )]
  cat('Plotting all libraries that share p5 or p7 with top N libs:', dt.idx[, sum(top.rg)], '\n')
  # top.rg <- dt.idx[top.rg == T, RG]
  ## this assumes that unique() preserves order
  top.p5 <- dt.idx[top.rg == T, unique(p5seq)]
  top.p7 <- dt.idx[top.rg == T, unique(p7seq)]
  dt.idx <- dt.idx[p5seq %in% top.p5 & p7seq %in% top.p7]

  ######
  ## make a dt with all possible pairs
  dt.idx.full <- dt.idx[, CJ(p7seq, p5seq, unique = T)]
  setnames(dt.idx.full, c('p7seq', 'p5seq'))
  dt.idx.full <- merge(dt.idx.full, unique(dt.idx[, .(p7seq, p7ind)]), by='p7seq')
  dt.idx.full <- merge(dt.idx.full, unique(dt.idx[, .(p5seq, p5ind)]), by='p5seq')
  dt.idx.full <- merge(dt.idx.full, dt.idx[, .(p7seq, p5seq, RG, RG.cat, nseqs, top.rg)], by=c('p7seq', 'p5seq'), all=T)
  dt.idx.full[is.na(top.rg), top.rg := F]
  dt.idx.full[is.na(nseqs), nseqs := 0]
  dt.idx.full[is.na(RG), RG := 'unknown']
  dt.idx <- dt.idx.full
  setorder(dt.idx, -nseqs)

  ## set the order of p5/p7 indices
  dt.idx[, p5seq.fac := factor(p5seq, levels = top.p5)]
  dt.idx[, p7seq.fac := factor(p7seq, levels = top.p7)]

  dt.idx[, p5.fac.maxseqs := factor(p5seq, levels = top.p5)]
  dt.idx[, p7.fac.maxseqs := factor(p7seq, levels = top.p7)]


  # ggplot(dt.idx, aes(x=p7seq.fac, y=p5seq.fac, fill=log(nseqs))) + geom_tile()

  ## make another RG column with more information, for reporting
  ## ISSUE : this is only unique for each possible RG if dt.idx is filtered to represent known indices
  ## i.e., if some reads have '-' and '50' as their indices, but the '-' represents two different
  ## indices (AAAAA and AAAAT), then they will both have the RG.full 'unknown.-.50'
  dt.idx[, RG.full := paste(RG, p7ind, p5ind, sep = '.')]
  dt.idx[, RG.full.fac := factor(RG.full, levels = RG.full)]


  ### ISSUE : NOT REALLY CORRECT - should look for contam in all "expected" libs?
  ## get the top n_contam_sources read groups / index pairs (even if they aren't supposed to be in the pool)
  ## these are the targetted libraries
  cat('n_contam_sources:', n_contam_sources, '\n')
  # target.RG.full <- dt.idx[, RG.full[.I < n_contam_sources]]
  # target.RG.full <- dt.idx[, head(RG.full, as.integer( 1.5 * (n_contam_sources + n_target_libs)) )]
  target.RG.full <- dt.idx[top.rg == T, RG.full]
  #print(dt.idx[top.rg == T, .(RG.full, nseqs)])
  # RG.levels <- dt.idx[, .(RG.nseqs = sum(nseqs)), .(RG.full)]
  # setkey(RG.levels, RG.nseqs)
  # RG.levels[, RG.fac := factor(RG.full, levels = rev(RG.full))]
  # top.RG = RG.levels[as.numeric(RG.fac) < n_target_libs, RG.full]
  #
  # ## order the p5 indices
  # p5.levels <- dt.idx[, .(p5.nseqs = sum(nseqs), p5.maxseqs = max(nseqs), RG.full = unique(RG.full)), .(p5seq)]
  # setkey(p5.levels, p5.nseqs)
  # p5.levels[, p5.fac := factor(p5seq, levels = rev(unique(p5seq)))]
  # setkey(p5.levels, p5.maxseqs)
  # p5.levels[, p5.fac.maxseqs := factor(p5seq, levels = rev(unique(p5seq)))]
  # top.p5 = p5.levels[as.numeric(p5.fac) < p57.thresh, unique(p5.fac)]
  # top.p5 = p5.levels[as.numeric(p5.fac.maxseqs) < p57.thresh, unique(p5.fac.maxseqs)]
  # top.p5.from_RG = p5.levels[RG.full %in% top.RG, unique(p5.fac)]
  # p5.levels[, p5seq := NULL]
  #
  # ## order the p7 indices
  # p7.levels <- dt.idx[, .(p7.nseqs = sum(nseqs), p7.maxseqs = max(nseqs), RG.full = unique(RG.full)), .(p7seq)]
  # setkey(p7.levels, p7.nseqs)
  # p7.levels[, p7.fac := factor(p7seq, levels = rev(unique(p7seq)))]
  # setkey(p7.levels, p7.maxseqs)
  # p7.levels[, p7.fac.maxseqs := factor(p7seq, levels = rev(unique(p7seq)))]
  # top.p7 = p7.levels[as.numeric(p7.fac) < p57.thresh, unique(p7.fac)]
  # top.p7 = p7.levels[as.numeric(p7.fac.maxseqs) < p57.thresh, unique(p7.fac.maxseqs)]
  # top.p7.from_RG = p7.levels[RG.full %in% top.RG, unique(p7.fac)]
  # p7.levels[, p7seq := NULL]

  ## decide if index pairs are in fresh kills (this is slow - could do just on RG.full %in% top.RG)
  cat('num in fresh kills\n')
  dt.idx[, in.fk := 0L]
  if (use_fk) {
    # print(head(dt.idx[RG.full %in% target.RG.full, .N, .(p7ind, p5ind)]))
    dt.idx[RG.full %in% target.RG.full, in.fk := num_in_fresh_kills(dt.fresh_kills, p7ind, p5ind), .(p7ind, p5ind)]
  } else {
      dt.idx[RG.full %in% target.RG.full, in.fk := 0L]
  }
  dt.idx[in.fk > 0 & (RG.cat == 'unknown' | RG.cat == 'unexpected'), RG.cat := 'in_fk']

  # print(head(dt.idx))

  ## cat('looking for den\n')
  ## dt.idx[RG.full %in% top.RG, fk_den := sum(description_from_fresh_kills(roots_from_fresh_kills(p7ind, p5ind)) %like% 'enisova'), .(p7ind, p5ind)]
  ## cat('doing nothing\n')
  ## print(dt.idx[fk_den > 0 & RG.cat == 'expected', .(RG, p7ind, p5ind)])
  ## dt.idx[fk_den > 0 & RG.cat == 'expected', description_from_fresh_kills(roots_from_fresh_kills(p7ind, p5ind)), .(p7ind, p5ind)]
  ## cat('done\n')

  ## ISSUES do we need this? it's setting factors for rg and indices, but.. indices already have levels?
  ## merge the data
  # dt.idx <- merge(dt.idx, RG.levels, by='RG.full')
  # dt.idx <- merge(dt.idx, p5.levels, by='RG.full')
  # dt.idx <- merge(dt.idx, p7.levels, by='RG.full')

  # dt.idx[, RG.rank := as.numeric(RG.fac)]

  top.RG.p5 <- dt.idx[RG.full %in% target.RG.full, p5seq.fac]
  top.RG.p7 <- dt.idx[RG.full %in% target.RG.full, p7seq.fac]

  ## assign colors for the different read categories - used only for histograms, I think?
  rg.cats <- c("expected", "PhiX", "unex_PhiX", "unknown", "in_fk", 'unexpected')
  rg.cat.colors <- qualitative_hcl(length(rg.cats), "Dark 3")
  names(rg.cat.colors) <- rg.cats

  ## look for cross-contamination only in the expected libraries
  test.ids <- dt.idx[RG.cat == 'expected', RG.full]
  if (!is.null(limit_search_libs)) {
    test.ids <- head(test.ids, limit_search_libs)
  }

  contam.ids <- head(target.RG.full, n_contam_sources)


  ## which of this data do we actually need to save?
  ret <- list(
    ## the "full" dataset [has already been trimmed]
    'dt.idx' = dt.idx,
    ## Top X RG, the real targets
    'top.RG' = target.RG.full,
    'test.ids' = test.ids,
    'contam.ids' = contam.ids,
    ##
    # 'top.p5.from_RG' = top.RG.p5,
    # 'top.p7.from_RG' = top.RG.p7,
    'top.p5.from_RG' = top.RG.p5,
    'top.p7.from_RG' = top.RG.p7,
    ## just the RG in the top X
    'dt.idx.top' = dt.idx[RG.full %in% target.RG.full],
    ## RG that share p5/p7 in the top X - much larger, the full matrix in which we might look for random targets
    'dt.idx.top.p57' = dt.idx[p5seq.fac %in% top.RG.p5 & p7seq.fac %in% top.RG.p7],
    'split_file' = split_stats_file,
    'rg.cat.colors' = rg.cat.colors)

  return(ret)
}

run_compute_swaps <- function(my.splits, random_contam_factor = 10) {

  ## look for cross-contamination only in the expected libraries
  test.ids <- my.splits$test.ids
  # contam.ids <- tail(my.splits$top.RG, n_contam_sources)
  contam.ids <- my.splits$contam.ids

  cat('looking for contamination IN   N libs:', length(test.ids), '\n')
  cat('looking for contamination FROM N libs:', length(contam.ids), '\n')

  ## this number - n_contam_sources - determines how many libraries are considered putative candidates as the sources of contamination
  ## this could be changed to get better stats e.g. for pools with a lower number of libraries in them, use a lower number
  ## could also be calculated on the fly?

    dt.swaps.null <- compute_swaps(my.splits, test.ids, contam.ids, rand_contam_ids = random_contam_factor * length(contam.ids))
    dt.swaps.test <- compute_swaps(my.splits, test.ids, contam.ids)
    
    setorder(dt.swaps.test, -test.stat)
    setorder(dt.swaps.null, -test.stat)

  x.p <- sapply(dt.swaps.test[, test.stat], function(test.x) dt.swaps.null[, (sum(test.x <= test.stat)+1) / (.N+1)])
  x.s <- sapply(dt.swaps.test[, test.stat], function(test.x) dt.swaps.null[, sum(test.x <= test.stat)])
  dt.swaps.test[, test.stat.rank := .I]
  dt.swaps.test[, empirical.p := x.p]
  dt.swaps.test[, empirical.n := x.s]
    dt.swaps.test[, empirical.p.bonf := pmin(x.p * .N, 1)]
    q.tmp <- tryCatch(dt.swaps.test[, qvalue(empirical.p)$qvalues],
                      error=function(cond) {
                          message(paste("q-value calculation failed. Could be due to low numbers of observations:", dt.swaps.test[, .N]))
                          ## dt.swaps.test[, rep(NA, .N)]
                          NA
                      })
    
    ## dt.swaps.test[, empirical.q := qvalue(empirical.p)$qvalues]
  dt.swaps.test[, empirical.q := q.tmp]
  dt.swaps.test[x.s == 0, empirical.p.flag := 'P_ISSUE_USE_RCF_FLAG']
  dt.swaps.test[x.s > 0, empirical.p.flag := '.']
    

  my.splits$dt.swaps.test <- dt.swaps.test
  my.splits$dt.swaps.null <- dt.swaps.null
  my.splits$test.ids <- test.ids
  my.splits$contam.ids <- contam.ids
  my.splits
}


compute_swaps <- function (my.splits, test_rg, contam.ids, rand_contam_ids = NULL, method = 'og') {
  if (method == 'og') compute_swaps_og(my.splits, test_rg, contam.ids, rand_contam_ids)
}


compute_swaps_og <- function (my.splits, test_rg, contam.ids, rand_contam_ids = NULL) {
  contam.p5 <- my.splits$dt.idx.top.p57[RG.full %in% contam.ids, p5seq]
  contam.p7 <- my.splits$dt.idx.top.p57[RG.full %in% contam.ids, p7seq]
  og.contam.ids <- contam.ids

  if (!is.null(rand_contam_ids))
      cat(sprintf('Randomly sampling %d from all [%d] RG to select sources of contamination\n', rand_contam_ids,
                  my.splits$dt.idx.top.p57[nseqs > 0 & p5seq %in% contam.p5 & p7seq %in% contam.p7 & !(RG.full %in% contam.ids), .N]))

  cat('|', paste0(rep('-', length(test_rg)), collapse = ''), '|\n  ')
  
  idx = 0
  ## first cycle through all of the potentially contaminated libraries (the libraries of interest, e.g. test libraries)
  dt.swaps <- foreach(my.id1 = test_rg, .combine = rbind) %dopar% {
    idx = idx + 1
    #cat('\n',my.id1, ':', idx, '/', length(test_rg))
    cat('.')

    # often the test libraries are also potential sources of contamination. we don't want to compute the full symmetric comparison
    # make sure that if there are rg in both sets, we don't do the test both ways
    rm_rg <- test_rg[my.id1 < test_rg]

    ## if no targets were provided, randomly sample from all RG (for getting background)
    ## should we restrict to only sampling from.. the indices present in the test targets?
    if (!is.null(rand_contam_ids)) {
      # ggplot(my.splits$dt.idx.top.p57, aes(x=p7.fac.maxseqs, y=p5.fac.maxseqs, fill=log(nseqs))) + geom_tile()
                                        # contam.ids <- my.splits$dt.idx.top.p57[nseqs > 0 & !(RG.full %in% contam.ids), sample(RG.full,rand_contam_ids)]
        ## sample random_contam_factor times more random contamination sources than real contamination sources, to help with p-values
      contam.ids <- my.splits$dt.idx.top.p57[nseqs > 0 & p5seq %in% contam.p5 & p7seq %in% contam.p7 & !(RG.full %in% contam.ids),
                                             sample(RG.full, min(.N, rand_contam_ids))]
    }

    foreach(my.id2 = contam.ids[!contam.ids %in% rm_rg], .combine = rbind) %do% {
      # cat('\n',my.id1, my.id2, ':', idx, '/', length(test_rg))
      if (my.id1 == my.id2) return(data.table())
      # cat('.')
      my.id1.p5 <- my.splits$dt.idx.top.p57[RG.full == my.id1, p5.fac.maxseqs]
      my.id1.p7 <- my.splits$dt.idx.top.p57[RG.full == my.id1, p7.fac.maxseqs]
      my.id2.p5 <- my.splits$dt.idx.top.p57[RG.full == my.id2, p5.fac.maxseqs]
      my.id2.p7 <- my.splits$dt.idx.top.p57[RG.full == my.id2, p7.fac.maxseqs]
      x1.rg = my.splits$dt.idx.top.p57[p7.fac.maxseqs == my.id1.p7 & p5.fac.maxseqs == my.id2.p5, RG.full]
      x2.rg = my.splits$dt.idx.top.p57[p7.fac.maxseqs == my.id2.p7 & p5.fac.maxseqs == my.id1.p5, RG.full]

      ## also test to see if either of the corner rg are 'expected' - then we shouldn't count them as evidence
      if ((x1.rg == my.id1 & x2.rg == my.id2) | (x1.rg == my.id2 & x2.rg == my.id1)) return(data.table())

      x1.nseqs <- my.splits$dt.idx.top.p57[p7.fac.maxseqs == my.id1.p7 & p5.fac.maxseqs == my.id2.p5, nseqs]
      x2.nseqs <- my.splits$dt.idx.top.p57[p7.fac.maxseqs == my.id2.p7 & p5.fac.maxseqs == my.id1.p5, nseqs]

      data.table(my.id1, my.id2, x1.nseqs, x2.nseqs,
                 id1.RG = my.splits$dt.idx.top.p57[RG.full == my.id1, RG],
                 id2.RG = my.splits$dt.idx.top.p57[RG.full == my.id2, RG],
                 id1.nseqs = my.splits$dt.idx.top.p57[RG.full == my.id1, nseqs],
                 id2.nseqs = my.splits$dt.idx.top.p57[RG.full == my.id2, nseqs],
                 x1.rg, x2.rg)
      ## also save both "swap" RG.full, and do the computation including the top X unexpected in_fk RG
    }
  }
  cat('\n')
  dt.swaps[, x := log((x1.nseqs+x2.nseqs+2)/2)]
  dt.swaps[, y := abs(log((x1.nseqs+1)/(x2.nseqs+1)))]
  dt.swaps[, y0 := log((x1.nseqs+1)/(x2.nseqs+1))]
  dt.swaps[, test.stat := x-y]
  dt.swaps
  # hey
}











## start plotting
# setwd('~/Google Drive/debug_contamination_mt/')
# pdf(sprintf('test_plots_%s.pdf', file_tag), width=14, height=10)

plot_all_libs_hist <- function(my.splits, nlibs = NULL) {
  ## plot a histogram of the read counts for expected and unexpected libraries
  if (!is.null(nlibs)) {
    dt.idx.top <- head(my.splits$dt.idx.top, nlibs)
  } else {
    dt.idx.top <- my.splits$dt.idx.top
  }
  p1 <- ggplot(dt.idx.top, aes(x=RG.full.fac, weight=nseqs/1e5, fill = RG.cat)) +
    geom_bar() +
    ylab('Read counts (100k)') +
    theme(axis.text.x = element_blank()) +
    xlab('Top Read Groups') +
    scale_fill_manual(values = my.splits$rg.cat.colors) +
    ggtitle(sprintf('Read counts for top %s index pairs', nrow(dt.idx.top))) +
    # geom_vline(xintercept = 120) +
    NULL
  print(p1)
}

add_description <- function(dt.idx.top_unexpected_index_combs, dt.fresh_kills = NULL) {
  # cat('adding desc:', nrow(dt.idx.top_unexpected_index_combs), '--', nrow(dt.fresh_kills), '\n')
  if (!is.null(dt.fresh_kills)) {
    dt.idx.top_unexpected_index_combs[, desc := paste0(libid_from_fresh_kills(dt.fresh_kills, roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind)), '__',
                                                       description_from_fresh_kills(dt.fresh_kills, my.id = roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind)), '__',
                                                   dates_in_fresh_kills(dt.fresh_kills, my.id = roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind)),
                                                   collapse = '\n'), .(p7ind, p5ind)]
    dt.idx.top_unexpected_index_combs[desc == '____', desc := paste0(RG, '_', p7ind, '_', p5ind)]
  } else {
    dt.idx.top_unexpected_index_combs[, desc := paste0(RG, '_', p7ind, '_', p5ind)]
  }
}


if (F) {

  dt.idx.tmp <- my.splits$dt.idx[RG == 'M5264']
  add_description(dt.idx.tmp, dt.fresh_kills)

  dt.idx.tmp[, libid_from_fresh_kills(dt.fresh_kills, roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind))]


  dt.idx.top_unexpected_index_combs[, roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind, this.debug = T), .(p7ind, p5ind)]

  dt.fresh_kills[doc.index_P7 == '320' & doc.index_P5 == '52']

  dt.idx.top_unexpected_index_combs[, description_from_fresh_kills(dt.fresh_kills,
                                                                   my.id = roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind)),
                                    .(p7ind, p5ind)]
  dt.idx.top_unexpected_index_combs[, dates_in_fresh_kills(dt.fresh_kills, my.id = roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind)), .(p7ind, p5ind)]

  dt.idx.top_unexpected_index_combs[, desc.x := paste0(description_from_fresh_kills(dt.fresh_kills, my.id = roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind)), '_',
                                                     dates_in_fresh_kills(dt.fresh_kills, my.id = roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind)),
                                                     collapse = '\n'), .(p7ind, p5ind)]

  dt.idx.top_unexpected_index_combs[, paste0(RG, '__', roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind), '__',
                                             description_from_fresh_kills(dt.fresh_kills, my.id = roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind)), '_',
                                                       dates_in_fresh_kills(dt.fresh_kills, my.id = roots_from_fresh_kills(dt.fresh_kills, p7ind, p5ind)),
                                                       collapse = '\n'), .(p7ind, p5ind)]


  if (!is.null(dt.fresh_kills)) {
    dt.id.label[, desc := paste0(RG, '_', description_from_fresh_kills(roots_from_fresh_kills(p7ind, p5ind)), '_',
                                 dates_in_fresh_kills(my.id = roots_from_fresh_kills(p7ind, p5ind)),
                                 collapse = '\n'), .(p7ind, p5ind)]
    dt.id.label[desc == '__', desc := RG.full]
  } else {
    dt.id.label[, desc := RG.full]
  }

}


plot_unexpected_libs_hist <- function(my.splits, dt.idx.top_unexpected_index_combs, nseqs.thresh = 10, dt.fresh_kills = NULL) {
  p1 <- ggplot(my.splits$dt.idx[nseqs > nseqs.thresh & RG.cat != 'expected'],
               aes(x=nseqs/100, fill=RG.cat, weight=nseqs, group=interaction(p5ind, p7ind))) +
    geom_histogram(color=rgb(0,0,0,.4)) +
    scale_x_log10() +
    ggtitle(sprintf('Unexpected index pairs w/ > %d reads', nseqs.thresh)) +
    xlab('Number of sequences (x100, log scale)') +
    ylab('Number of sequences') +
    geom_text_repel(data = dt.idx.top_unexpected_index_combs,
                    mapping = aes(label = desc,
                                  y=max(dt.idx.top_unexpected_index_combs$nseqs)),
                    box.padding = 2, force = 10, max.iter = 200000) +
    scale_fill_manual(values = my.splits$rg.cat.colors) +
    scale_color_manual(values = my.splits$rg.cat.colors)
  print(p1)
}
 # plot_unexpected_libs_hist(my.splits, dt.idx.top_unexpected_index_combs)

plot_simple_heatmap <- function(my.splits, dt.idx.top_unexpected_index_combs = NULL) {

  p.tile <- ggplot(my.splits$dt.idx.top.p57,
                   aes(x=p7.fac.maxseqs, y=p5.fac.maxseqs, fill=log(nseqs))) +
    geom_tile() + scale_fill_distiller(type='div', palette = 2) +
    geom_abline(slope=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
    NULL

  if (!is.null(dt.idx.top_unexpected_index_combs)) {
    ## restrict just to those libraries in the top p57
    ## if we don't do this, then it messes up the ordering of the axes in the heatmap, for some reason
    dt.idx.top_unexpected_index_combs.tmp <-
      dt.idx.top_unexpected_index_combs[!is.na(desc) & RG.full %in% my.splits$dt.idx.top.p57$RG.full]
    p.tile <- p.tile +
      geom_point(data=dt.idx.top_unexpected_index_combs.tmp, color='black') +
      geom_point(data = my.splits$dt.idx.top[RG.cat == 'in_fk'], color='red') +
      geom_point(data = my.splits$dt.idx.top[RG.cat == 'expected'], color='blue') +
      geom_text_repel(data=dt.idx.top_unexpected_index_combs.tmp, aes(label = desc), box.padding = 1) +
      NULL
  }

  # ggplot(my.splits$dt.idx.top.p57,
  #        aes(x=p7.fac.maxseqs, y=p5.fac.maxseqs, fill=log(nseqs))) +
  #   geom_tile() + scale_fill_distiller(type='div', palette = 2) +
  #   geom_point(data=, color='black') +
  #   # geom_point(data=head(my.splits$dt.idx.top.p57), color='black') +
  #   geom_abline(slope=1) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  #   NULL

  print(p.tile)
  p.tile
}
# plot_simple_heatmap(my.splits, dt.idx.top_unexpected_index_combs)

plot_simple_heatmap2 <- function(my.splits, dt.idx.labels = NULL, dt.idx.labels.colors = NULL, my.title = NULL, my.title.color = 'black') {

    p.tile <- ggplot(my.splits$dt.idx.top.p57,
                     aes(x=p7.fac.maxseqs, y=p5.fac.maxseqs, fill=log(nseqs))) +
        geom_tile() + scale_fill_distiller(type='div', palette = 2) +
        geom_abline(slope=1) +
        theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
        ggtitle(my.title) +
        theme(plot.title = element_text(color=my.title.color, face="bold")) +
        NULL

  if (!is.null(dt.idx.labels)) {
    if ('category' %in% colnames(dt.idx.labels)) {
      p.tile <- p.tile + geom_point(data=dt.idx.labels, aes(color = category))
    } else {
      p.tile <- p.tile + geom_point(data=dt.idx.labels, color = 'red')
    }
  }

  if (!is.null(dt.idx.labels) & 'desc' %in% colnames(dt.idx.labels)) {
    dt.idx.labels.text <-
      dt.idx.labels[!is.na(desc) & RG.full %in% my.splits$dt.idx.top.p57$RG.full]
    if ('category' %in% colnames(dt.idx.labels)) {
      p.tile <- p.tile + geom_text_repel(data=dt.idx.labels.text, aes(label = desc, color = category), box.padding = 1)
    } else {
      p.tile <- p.tile + geom_text_repel(data=dt.idx.labels.text, aes(label = desc), box.padding = 1)
    }
  }

  # ggplot(my.splits$dt.idx.top.p57,
  #        aes(x=p7.fac.maxseqs, y=p5.fac.maxseqs, fill=log(nseqs))) +
  #   geom_tile() + scale_fill_distiller(type='div', palette = 2) +
  #   geom_point(data=, color='black') +
  #   # geom_point(data=head(my.splits$dt.idx.top.p57), color='black') +
  #   geom_abline(slope=1) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  #   NULL

  print(p.tile)
  p.tile
}
if (F) {
  dt.idx.top_unexpected_index_combs <- head(my.splits$dt.idx[RG.cat != 'expected'],5)
  add_description(dt.idx.top_unexpected_index_combs, dt.fresh_kills)

  dt.idx.top_unexpected_index_combs[, category := sample(c('hey', 'what'), .N, replace = T)]
  plot_simple_heatmap2(my.splits, dt.idx.top_unexpected_index_combs)

  dt.idx.top_unexpected_index_combs <- my.splits$dt.idx[RG.cat == 'expected' | RG.cat == 'PhiX' | RG.full %in% my.splits$contam.ids]
  # add_description(dt.idx.top_unexpected_index_combs, dt.fresh_kills)

  dt.idx.top_unexpected_index_combs[, category := 'possible contaminants']
  plot_simple_heatmap2(my.splits, dt.idx.top_unexpected_index_combs, my.title = 'Possible Contaminants')
  plot_simple_heatmap2(my.splits, dt.idx.top_unexpected_index_combs)

}

plot_xy_null_and_test <- function(dt.swaps.test, dt.swaps.null, y.abs = T, test.stat.thresh = 5, dt.swaps.label = NULL) {

  dt.swaps.hits <- dt.swaps.test[test.stat > test.stat.thresh][order(test.stat,decreasing = T)]

  if (y.abs) {
    p1 <- ggplot(dt.swaps.test, aes(x=x, y=y))
  } else {
    p1 <- ggplot(dt.swaps.test, aes(x=x, y=y0))
  }

    if (nrow(dt.swaps.null) > 20000) {
        cat('Downsampling background points to avoid gigantic pdf:', nrow(dt.swaps.null), '\n')
        dt.swaps.null <- dt.swaps.null[sample(.N, 20000)]
    }

  p1 <- p1 +
    geom_point(data=dt.swaps.null, color='red', alpha=.2, size=3) +
    geom_point(alpha=.2, color='black') +
    geom_point(data=dt.swaps.hits, color='black', alpha=1) +
    geom_abline(slope=1, intercept = c(-4,-5,-6,-7)) +
    xlab('log(nseqs in corner libraries)') +
    NULL

  if (y.abs) {
    p1 <- p1 + ylab('abs log(corner1/corner2)')
  } else {
    p1 <- p1 + ylab('log(corner1/corner2)')
  }

  print(p1)
}


plot_qq_null_and_test <- function(dt.swaps.test, dt.swaps.null, plot.libs = NULL, test.stat.thresh = 5) {

  # qqplot(dt.swaps.null[, test.stat],
  #        dt.swaps.test[, test.stat],
  #        xlab = 'Randomized score distribution',
  #        ylab = 'Test score distribution')
  # abline(b = 1, a=0)

  dt.swaps.test.q <- data.table(dt.swaps.test)
  # dt.swaps.test.q[, test.stat := x-y]
  setorder(dt.swaps.test.q, test.stat)
  dt.swaps.test.q$test.stat.q <- quantile(dt.swaps.null[, test.stat], seq(0,1,1/(nrow(dt.swaps.test)-1)), names=T)
  # dt.swaps.test.q <- quantile(dt.swaps.test[, test.stat], seq(0,1,1/nrow(dt.swaps.test)))
  # dt.swaps.test.q <- quantile(dt.swaps.test[, test.stat], seq(0,1,.001))
  p1 <- ggplot(dt.swaps.test.q, aes(x=test.stat.q, y=test.stat)) +
    xlab('Background scores (random contaminates)') +
    ylab('Test scores (top X most abundant libraries as contaminants)')

  if (is.null(plot.libs)) {
    p1 <- p1 + geom_point(alpha=.5)
    p1 <- p1 + geom_text_repel(data=tail(dt.swaps.test.q[test.stat > test.stat.thresh], 10), aes(label=sprintf('%s-%s', id1.RG, id2.RG)), color='red')
  } else {
      dt.swaps.test.q.lab <- tail(dt.swaps.test.q[id1.RG %in% plot.libs | id2.RG %in% plot.libs][test.stat > test.stat.thresh], 10)
    p1 <- p1 + geom_point(aes(color = (id1.RG %in% plot.libs | id2.RG %in% plot.libs)), alpha=.5) + scale_color_manual('Lib of interest', values = c('black', 'red'))
    p1 <- p1 + geom_point(data=dt.swaps.test.q[id1.RG %in% plot.libs | id2.RG %in% plot.libs], alpha=.5, color='red')
    p1 <- p1 + geom_text_repel(data=dt.swaps.test.q.lab,
                               aes(label=sprintf('%s-%s', id1.RG, id2.RG)), color='red')
  }

  p1 <- p1 + geom_abline(slope = 1, lty=3)

  p1 <- p1 + ggtitle('q-q plot: real vs background scores')

  print(p1)
}

if (F) {
  plot_qq_null_and_test(dt.swaps.test, dt.swaps.null, plot.libs = 'A13602', test.stat.thresh = 5)
  plot_qq_null_and_test(dt.swaps.test, dt.swaps.null, test.stat.thresh = 5)
}



plot_putative_contam_squares_heatmap <- function(my.splits, dt.swaps.hits.plot, my.swap, dt.fresh_kills = NULL, title.prefix = '', title.color = 'black') {
  rg.square <- dt.swaps.hits.plot[my.swap, c(my.id1,my.id2,x1.rg,x2.rg)]
  my.id1 <- dt.swaps.hits.plot[my.swap, my.id1]
  my.id2 <- dt.swaps.hits.plot[my.swap, my.id2]
  #print(dt.swaps.hits.plot)
  my.id.subtitle <-
    dt.swaps.hits.plot[my.swap, sprintf('y-x: %0.2f; nseqs: %d and %d; corners: %d and %d    ',
                                        test.stat, id1.nseqs, id2.nseqs, x1.nseqs, x2.nseqs)]
  if ('empirical.p' %in% colnames(dt.swaps.hits.plot)) {
      my.id.subtitle <- paste0(my.id.subtitle,
                               dt.swaps.hits.plot[my.swap, sprintf('p: %0.4f;  p-bonf: %0.4f;  qval: %0.4f;  flag: %s, %d/%d;  rank: %d/%d',
                                                                   empirical.p, empirical.p.bonf, empirical.q,
                                                                   empirical.p.flag, empirical.n, my.splits$dt.swaps.null[, .N],
                                                                   test.stat.rank, my.splits$dt.swaps.test[, .N])])
  }
  # empirical.p empirical.p.bonf empirical.q
  dt.id.label <- my.splits$dt.idx[RG.full == my.id1 | RG.full == my.id2]

  add_description(dt.id.label, dt.fresh_kills)

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
    ggtitle(sprintf('%sPutative contaminated samples: %s and %s', title.prefix, my.id1, my.id2),
            subtitle = my.id.subtitle) +
      theme(plot.title = element_text(color=title.color, face="bold")) +
      NULL
  print(p1)
}


# dt <- rbind(dt.swaps.null[, .(test.stat, category = 'null')],
#             dt.swaps.test[, .(test.stat, category = 'test')])
# ggplot(dt, aes(x=test.stat, fill=category)) + geom_histogram()

############################################################
############################################################
#####
#####
##### MAIN PLOTTING FUNCTION - CALLS ALL THE OTHER FUNCTIONS
#####
#####
############################################################
############################################################

plot_debug_splits <- function(my.splits, n.hits, plot.libs = NULL, dt.fresh_kills = NULL) {

  # if (is.null(dt.fresh_kills)) {
  #   use_fk <- F
  # } else {
  #   use_fk <- T
  # }



  ## dt.swaps.test has the
  ## I modify these data tables, so make local copies of them
  dt.swaps.null <- data.table(my.splits$dt.swaps.null)
  dt.swaps.test <- data.table(my.splits$dt.swaps.test)

  ## big histogram of sequence counts for all libraries
  plot_all_libs_hist(my.splits)

  ## smaller histogram of sequence counts for top N libraries - easier to read
  plot_all_libs_hist(my.splits, nlibs = 200)


  ## get FK names for the top N unexpected index pairs
  ## these are added in a column
  dt.idx.top_unexpected_index_combs <- head(my.splits$dt.idx[RG.cat != 'expected' & RG.cat != 'PhiX'],5)
  add_description(dt.idx.top_unexpected_index_combs, dt.fresh_kills = dt.fresh_kills)
  plot_unexpected_libs_hist(my.splits, dt.idx.top_unexpected_index_combs)

  # dt.idx.top_unexpected_index_combs <- head(my.splits$dt.idx[RG.cat != 'expected'],5)
  # # dt.idx.top_unexpected_index_combs <- dt.idx.top_unexpected_index_combs[5]
  # add_description(dt.idx.top_unexpected_index_combs, dt.fresh_kills = dt.fresh_kills)
  # dt.idx.top_unexpected_index_combs

  plot_simple_heatmap2(my.splits, my.title = 'Number of sequenced reads per index pair')

  dt.idx.plot_index_combs <- my.splits$dt.idx[RG.cat == 'expected']
  plot_simple_heatmap2(my.splits, dt.idx.plot_index_combs, my.title = 'Expected Libraries')

  dt.idx.plot_index_combs <- my.splits$dt.idx[RG.full %in% my.splits$test.ids]
  plot_simple_heatmap2(my.splits, dt.idx.plot_index_combs, my.title = 'Searched for Contamination in these Libraries (usually all expected libs)')

  plot_simple_heatmap2(my.splits, dt.idx.top_unexpected_index_combs, my.title = 'Top 5 Unexpected Libraries')

  dt.idx.plot_index_combs <- my.splits$dt.idx[RG.full %in% my.splits$contam.ids]
  plot_simple_heatmap2(my.splits, dt.idx.plot_index_combs, my.title = 'Possible Contaminant Search Space')

  # plot_simple_heatmap(my.splits)
  # plot_simple_heatmap(my.splits, dt.idx.top_unexpected_index_combs)

  #### look at and plot swaps / squares

  # dt.swaps.test[(x1.rg %like% 'A' | x2.rg %like% 'A')]
  # dt.swaps.null[(x1.rg %like% 'A' | x2.rg %like% 'A')]
  # dt.swaps.null[x > 8]

  # dt.swaps.test[, min.nseqs := min(id1.nseqs, id2.nseqs), .(id1.nseqs, id2.nseqs)]
  # dt.swaps.null[, min.nseqs := min(id1.nseqs, id2.nseqs), .(id1.nseqs, id2.nseqs)]

  dt.swaps.hits <- dt.swaps.test[test.stat > 5][order(test.stat,decreasing = T)]
  # dt.swaps.label <- dt.swaps.test[my.id1 %like% plot.libs | my.id2 %like% plot.libs][order(test.stat,decreasing = T)]
  dt.swaps.label <- dt.swaps.test[id1.RG %in% plot.libs | id2.RG %in% plot.libs][order(test.stat,decreasing = T)]

  ## this is a plot that shows:
  # - x-axis: log(number of sequences on corner libraries)
  # - y-axis: log(corner1/corner2)
  ## x-y [test.stat] is the statistic used to identify putative contamination
  plot_xy_null_and_test(dt.swaps.test, dt.swaps.null, test.stat.thresh = 5)
  plot_xy_null_and_test(dt.swaps.test, dt.swaps.null, test.stat.thresh = 5, y.abs = F)

  plot_qq_null_and_test(dt.swaps.test, dt.swaps.null, test.stat.thresh = 5)
  plot_qq_null_and_test(dt.swaps.test, dt.swaps.null, plot.libs, test.stat.thresh = 5)

  hist(dt.swaps.test$empirical.p, col='black', breaks = seq(0,1,.01), main = 'p-value histogram')


  if (!is.null(plot.libs)) {
    dt.swaps.hits.plot <- dt.swaps.test[my.id1 %in% plot.libs | my.id2 %in% plot.libs][order(test.stat,decreasing = T)]
    dt.swaps.hits.plot <- head(dt.swaps.hits.plot, n.hits)
  } else {
    dt.swaps.hits.plot <- head(dt.swaps.test[order(test.stat,decreasing = T)], n.hits)
  }
  my.swap <- 1

  for (my.swap in dt.swaps.hits.plot[, .I]) {
    plot_putative_contam_squares_heatmap(my.splits, dt.swaps.hits.plot, my.swap, dt.fresh_kills)
  }

  # dt.swaps.hits.plot <- head(dt.swaps.null[order(test.stat,decreasing = T)], n.hits)
  dt.swaps.hits.plot <- head(dt.swaps.null[order(test.stat,decreasing = T)], 5)
  for (my.swap in dt.swaps.hits.plot[, .I]) {
    plot_putative_contam_squares_heatmap(my.splits, dt.swaps.hits.plot, my.swap, dt.fresh_kills, title.prefix = 'Background (top): ', title.color = 'red')
  }

  #dt.swaps.hits.plot <- tail(dt.swaps.null[order(test.stat,decreasing = T)], n.hits)
  dt.swaps.hits.plot <- tail(dt.swaps.null[order(test.stat,decreasing = T)], 5)
  for (my.swap in dt.swaps.hits.plot[, .I]) {
    plot_putative_contam_squares_heatmap(my.splits, dt.swaps.hits.plot, my.swap, dt.fresh_kills, title.prefix = 'Background (bottom): ', title.color = 'blue')
  }

}





## make a plot w/ the flip of the matrix.. like if two samples are mixed together, you would expect
## to see p5a w/ p7b and p7b with p5a.  so then the x value is all the points on the upper triangle
## of the matrix, and the y axis is the corresponding point on the lower triangle.  If this is
## random, then there shouldn't be a correlation between the two.

## also, looks like there are some rows/cols that have low levels of pairing w/ all other p5 or p7
## indices. how to quantify that?



