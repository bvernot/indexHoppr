 



load_splitstats <- function(split_stats_file, n_top_target,
                            rg.thresh = 500, p57.thresh = 500,
                            compute.swaps = F) {
  dt.idx <- fread(split_stats_file)
  
  if (!'p7ind' %in% names(dt.idx)) {
    dt.idx <- merge(dt.idx, dt.idx[, .N, p7seq][, .(p7seq, p7ind=.I)], by='p7seq')
  }
  if (!'p5ind' %in% names(dt.idx)) {
    dt.idx <- merge(dt.idx, dt.idx[, .N, p5seq][, .(p5seq, p5ind=.I)], by='p5seq')
  }
  
  dt.idx <- dt.idx[p7ind != '-' & p5ind != '-']
  dt.idx[p5ind == 'PhiX' & p7ind == 'PhiX', RG := 'PhiX']
  setnames(dt.idx, '#seqs', 'nseqs', skip_absent=TRUE)
  
  dt.idx.full <- dt.idx[, CJ(p7seq, p5seq, unique = T)]
  setnames(dt.idx.full, c('p7seq', 'p5seq'))
  dt.idx.full <- merge(dt.idx.full, unique(dt.idx[, .(p7seq, p7ind)]), by='p7seq')
  dt.idx.full <- merge(dt.idx.full, unique(dt.idx[, .(p5seq, p5ind)]), by='p5seq')
  dt.idx.full <- merge(dt.idx.full, dt.idx[, .(p7seq, p5seq, RG, nseqs)], by=c('p7seq', 'p5seq'), all=T)
  dt.idx.full[is.na(nseqs), nseqs := 0]
  dt.idx.full[is.na(RG), RG := 'unknown']
  dt.idx.full[RG == 'unexpected', RG := 'unknown']
  dt.idx.full[RG != 'PhiX' & (p7ind %like% 'Phi' | p5ind %like% 'Phi'), RG := 'unex_PhiX']
  dt.idx <- dt.idx.full
  
  dt.idx[, RG.full := paste(RG, p7ind, p5ind, sep = '.')]
  dt.idx[, RG.cat := RG]
  dt.idx[RG != 'unknown' & RG != 'unexpected' & RG != 'PhiX' & RG != 'unex_PhiX', RG.cat := 'expected']
  
  ## get the top 110 read groups / index pairs (even if they aren't supposed to be in the pool)
  RG.levels <- dt.idx[, .(RG.nseqs = sum(nseqs)), .(RG.full)]
  setkey(RG.levels, RG.nseqs)
  RG.levels[, RG.fac := factor(RG.full, levels = rev(RG.full))]
  top.RG = RG.levels[as.numeric(RG.fac) < rg.thresh, RG.full]
  
  ## order the p5 indices
  p5.levels <- dt.idx[, .(p5.nseqs = sum(nseqs), p5.maxseqs = max(nseqs), RG.full = unique(RG.full)), .(p5seq)]
  setkey(p5.levels, p5.nseqs)
  p5.levels[, p5.fac := factor(p5seq, levels = rev(unique(p5seq)))]
  setkey(p5.levels, p5.maxseqs)
  p5.levels[, p5.fac.maxseqs := factor(p5seq, levels = rev(unique(p5seq)))]
  top.p5 = p5.levels[as.numeric(p5.fac) < p57.thresh, unique(p5.fac)]
  top.p5 = p5.levels[as.numeric(p5.fac.maxseqs) < p57.thresh, unique(p5.fac.maxseqs)]
  top.p5.from_RG = p5.levels[RG.full %in% top.RG, unique(p5.fac)]
  p5.levels[, p5seq := NULL]
  
  ## order the p7 indices
  p7.levels <- dt.idx[, .(p7.nseqs = sum(nseqs), p7.maxseqs = max(nseqs), RG.full = unique(RG.full)), .(p7seq)]
  setkey(p7.levels, p7.nseqs)
  p7.levels[, p7.fac := factor(p7seq, levels = rev(unique(p7seq)))]
  setkey(p7.levels, p7.maxseqs)
  p7.levels[, p7.fac.maxseqs := factor(p7seq, levels = rev(unique(p7seq)))]
  top.p7 = p7.levels[as.numeric(p7.fac) < p57.thresh, unique(p7.fac)]
  top.p7 = p7.levels[as.numeric(p7.fac.maxseqs) < p57.thresh, unique(p7.fac.maxseqs)]
  top.p7.from_RG = p7.levels[RG.full %in% top.RG, unique(p7.fac)]
  p7.levels[, p7seq := NULL]
  
  ## decide if index pairs are in fresh kills (this is slow - could do just on RG.full %in% top.RG)
  cat('num in fresh kills\n')
  if (use_fk) {
      dt.idx[RG.full %in% top.RG, in.fk := num_in_fresh_kills(p7ind, p5ind), .(p7ind, p5ind)]
  } else {
      dt.idx[RG.full %in% top.RG, in.fk := 0]
  }
  dt.idx[in.fk > 0 & (RG.cat == 'unknown' | RG.cat == 'unexpected'), RG.cat := 'in_fk']
  
  ## cat('looking for den\n')
  ## dt.idx[RG.full %in% top.RG, fk_den := sum(description_from_fresh_kills(roots_from_fresh_kills(p7ind, p5ind)) %like% 'enisova'), .(p7ind, p5ind)]
  ## cat('doing nothing\n')
  ## print(dt.idx[fk_den > 0 & RG.cat == 'expected', .(RG, p7ind, p5ind)])
  ## dt.idx[fk_den > 0 & RG.cat == 'expected', description_from_fresh_kills(roots_from_fresh_kills(p7ind, p5ind)), .(p7ind, p5ind)]
  ## cat('done\n')
  
  ## merge the data
  dt.idx <- merge(dt.idx, RG.levels, by='RG.full')
  dt.idx <- merge(dt.idx, p5.levels, by='RG.full')
  dt.idx <- merge(dt.idx, p7.levels, by='RG.full')
  
  dt.idx[, RG.rank := as.numeric(RG.fac)]
  
  ret <- list('dt.idx' = dt.idx,
              'dt.idx.top' = dt.idx[RG.full %in% top.RG],
              'dt.idx.top.p57' = dt.idx[p5.fac %in% top.p5.from_RG & p7.fac %in% top.p7.from_RG],
              'top.RG' = top.RG,
              'top.p5.from_RG' = top.p5.from_RG,
              'top.p7.from_RG' = top.p7.from_RG,
              'split_file' = split_stats_file)
  
  if (compute.swaps) {
      ## look for cross-contamination only in the expected libraries
      test.ids <- ret$dt.idx[RG.cat == 'expected', RG.full]
      ## this number - 150 - determines how many libraries are considered putative candidates as the sources of contamination
      ## this could be changed to get better stats e.g. for pools with a lower number of libraries in them
      target.ids <- tail(ret$top.RG, n_top_target)
      dt.swaps.bak <- compute_swaps(ret, test.ids, target.ids, r.target = n_top_target)
      dt.swaps.test <- compute_swaps(ret, test.ids, target.ids)
      ret$dt.swaps.test <- dt.swaps.test
      ret$dt.swaps.bak <- dt.swaps.bak
      ret$test.ids <- test.ids
      ret$target.ids <- target.ids
  }
  return(ret)
}

compute_swaps <- function (my.splits, test_rg, target_rg, r.target = NULL) {
  target.p5 <- my.splits$dt.idx.top.p57[RG.full %in% target_rg, p5seq]
  target.p7 <- my.splits$dt.idx.top.p57[RG.full %in% target_rg, p7seq]
  og.target_rg <- target_rg
  
  if (!is.null(r.target))
    cat(sprintf('Randomly sampling %d from all RG to form background targets\n', r.target))
  
  idx = 0
  dt.swaps <- foreach(my.id1 = test_rg, .combine = rbind) %dopar% {
    idx = idx + 1
    cat('\n',my.id1, ':', idx, '/', length(test_rg))
    # make sure that if there are rg in both sets, we don't do the test both ways
    rm_rg <- test_rg[my.id1 < test_rg]
    
    ## if no targets were provided, randomly sample from all RG (for getting background)
    ## should we restrict to only sampling from.. the indices present in the test targets?
    if (!is.null(r.target)) {
      my.splits$dt.idx.top.p57[p5seq %in% target.p5 & p7seq %in% target.p7]
      target_rg <- my.splits$dt.idx.top.p57[nseqs > 0, sample(RG.full,r.target)]
    }
    
    foreach(my.id2 = target_rg[!target_rg %in% rm_rg], .combine = rbind) %do% {
      # cat('\n',my.id1, my.id2, ':', idx, '/', length(test_rg))
      if (my.id1 == my.id2) return(data.table())
      cat('.')
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
                 id1.nseqs = my.splits$dt.idx.top.p57[RG.full == my.id1, nseqs],
                 id2.nseqs = my.splits$dt.idx.top.p57[RG.full == my.id2, nseqs],
                 x1.rg, x2.rg)
      ## also save both "swap" RG.full, and do the computation including the top X unexpected in_fk RG
    }
  }
  dt.swaps[, x := log((x1.nseqs+x2.nseqs+2)/2)]
  dt.swaps[, y := abs(log((x1.nseqs+1)/(x2.nseqs+1)))]
  dt.swaps[, y0 := log((x1.nseqs+1)/(x2.nseqs+1))]
  dt.swaps
}
