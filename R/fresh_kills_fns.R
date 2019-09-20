library(data.table)
library(tidyverse)
# install.packages("rjson")
# install.packages('jsonlite')
library(jsonlite)

# ```curl 'http://bioaps01:5984/default/_all_docs?include_docs=true' > freshkills.json```
dt.fresh_kills <- fromJSON(fresh_kills_json, flatten = T)
# names(dt.fresh_kills)
dt.fresh_kills.rows <- data.table(dt.fresh_kills$rows)
# dt.fresh_kills.rows[doc.index_P7 == '417' & doc.index_P5 == '138', id]
# dt.fresh_kills.rows[id == 'G11097']
# dt.fresh_kills.rows[!is.na(doc.p7index.44)]

## clean it up a bit
dt.fresh_kills.rows <- dt.fresh_kills.rows[id != 'global_index_list']
dt.fresh_kills.rows <- dt.fresh_kills.rows[doc.type != 'run']
dt.fresh_kills.rows <- dt.fresh_kills.rows[!id %like% '_design']
## remove columns w/ all NA
dt.fresh_kills.rows <- dt.fresh_kills.rows[,which(unlist(lapply(dt.fresh_kills.rows, function(x)!all(is.na(x))))),with=F]

fetch_from_fresh_kills <- function(p7 = NULL, p5 = NULL, my.id = NULL) {
    if (!is.null(my.id)) {
        #cat('fkid:', my.id, '\n')
        return(dt.fresh_kills.rows[id %in% my.id])
    }
    else if (!is.null(p7) & !is.null(p5)) {
        #cat('fkp57:', p7, ':', p5, '\n')
        return(dt.fresh_kills.rows[doc.index_P7 == as.character(p7) & doc.index_P5 == as.character(p5)])
    }
    else {
        print('parent_from_fresh_kills: must give p7/p5 or id')
        return(NULL)
    }
}
fetch_from_fresh_kills(417, 138)
fetch_from_fresh_kills(my.id = 'A20322')

num_in_fresh_kills <- function(p7, p5) {
  dt.rows <- fetch_from_fresh_kills(p7, p5)
  dt.rows[, .N]
}
num_in_fresh_kills(417, 138)

dates_in_fresh_kills <- function(p7 = NULL, p5 = NULL, my.id = NULL) {
  if (!is.null(my.id))
    dt.rows <- fetch_from_fresh_kills(my.id = my.id)
  else if (!is.null(p7) & !is.null(p5)) 
    dt.rows <- fetch_from_fresh_kills(p7, p5)
  else {
    print('dates_in_fresh_kills: must give p7/p5 or id')
    return(NULL)
  }
  dt.rows[, doc.date]
}
dates_in_fresh_kills(417, 138)
dates_in_fresh_kills(my.id = 'A20322')

parent_from_fresh_kills <- function(p7 = NULL, p5 = NULL, my.id = NULL) {
  if (!is.null(my.id))
    dt.rows <- dt.fresh_kills.rows[id == my.id]
  else if (!is.null(p7) & !is.null(p5)) 
    dt.rows <- fetch_from_fresh_kills(p7, p5)
  else {
    print('parent_from_fresh_kills: must give p7/p5 or id')
    return(NULL)
  }
  return(dt.rows[, doc.parents] %>% unlist)
}
parent_from_fresh_kills(417, 138)
parent_from_fresh_kills(my.id = 'A20322')

root_from_fresh_kills <- function(my.id, lim.search = 10, this.debug = F) {
  parent <- my.id
  i <- 0
  #parent != ''
  while (T) {
    tmp <- parent_from_fresh_kills(my.id = parent)
    if (is.null(tmp)) break
    parent <- tmp
    if (i > lim.search) break
    if (this.debug) cat(my.id, i, parent, '\n')
    i = i+1
  }
  parent
}
root_from_fresh_kills(my.id = 'N2995')

roots_from_fresh_kills <- function(p7 = NULL, p5 = NULL, my.id = NULL, this.debug = F) {
  if (!is.null(my.id))
    dt.rows <- dt.fresh_kills.rows[id == my.id]
  else if (!is.null(p7) & !is.null(p5)) 
    dt.rows <- fetch_from_fresh_kills(p7, p5)
  else {
    print('roots_from_fresh_kills: must give p7/p5 or id')
    return(NULL)
  }
  sapply(dt.rows[, id], function(x) root_from_fresh_kills(x, this.debug = this.debug)) %>% unique
}
roots_from_fresh_kills(417, 138)
roots_from_fresh_kills(434, 119, this.debug = T)

description_from_fresh_kills <- function(my.id) {
  dt.fresh_kills.rows[id %in% my.id, doc.description]
}
description_from_fresh_kills(roots_from_fresh_kills(my.id='N2995'))
description_from_fresh_kills(roots_from_fresh_kills(434, 119))

description_from_fresh_kills(roots_from_fresh_kills(my.id='N2995'))
dates_in_fresh_kills(my.id = roots_from_fresh_kills(434, 119))
dates_in_fresh_kills(434, 119)

