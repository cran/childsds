##' ridits from ranks
##'
##' @title ridits from ranks
##' @param v vector of ranks
##' @param maxrank min possible rank
##' @param minrank max possible rank
##' @return vector of ridits
##' @author Mandy Vogel
##' @export
ridits.from.ranks <- function(v, maxrank = max(v), minrank = NULL){
  v <- stats::na.omit(v)
  if(is.null(minrank)) minrank <- min(c(min(v), 1))
  ranks <- minrank:maxrank
  v <- purrr::map_dbl(ranks, ~sum(. == v))
  rlang::set_names((cumsum(v) - 0.5 * v)/sum(v),
                   paste0("r",ranks))
}
##' propabilities
##'
##' @title propabilities
##' @param v vector of ranks
##' @param maxrank max possible rank
##' @param minrank min possible rank
##' @return 1-dim contingency table
##' @author Mandy Vogel
##' @export
props <- function(v, maxrank = max(v), minrank = NULL){
  v <- stats::na.omit(v)
  if(is.null(minrank)) minrank <- min(c(min(as.numeric(as.character(v))), 1))  
  prop.table(table(factor(v, levels = minrank:maxrank)))
}


##' mean ridits for vectors of paired observations
##'
##' @title mean ridits for vectors of paired observations
##' @param v vector 1
##' @param w vector 2
##' @param maxrank max possible rank
##' @param minrank min possible rank
##' @return 2-element named vector 
##' @author Mandy Vogel
##' @export
mean_ridits <- function(v, w, maxrank = max(v), minrank = NULL){
  v <- stats::na.omit(v)
  w <- stats::na.omit(w)
  if(is.null(minrank)) minrank <- min(c(min(v), 1))  
  props1 <- prop.table(table(factor(v, levels = minrank:maxrank)))
  ridits1 <- ridits.from.ranks(v, maxrank, minrank = minrank)
  props2 <- prop.table(table(factor(w, levels = minrank:maxrank)))
  ridits2 <- ridits.from.ranks(w, maxrank, minrank = minrank)  
  c(AY2 = sum(props2 * ridits1),
    AY1 = sum(props1 * ridits2))
}


##' calculate mid mean ranks
##'
##' @title calculate mid mean ranks
##' @param v1 vector of ranks
##' @param v2 vector of ranks
##' @param maxrank max possible rank
##' @param minrank min possible rank
##' @param conf.level confidence level alpha
##' @return list of marginal mean ranks and mid ranks
##' @author Mandy Vogel
##' @export
mid.mean.ranks <- function(v1, v2, maxrank = max(c(v2,v1)), minrank = NULL, conf.level = 0.05){
  if(length(v1) != length(v2)) stop("length of v1 and v2 differs")
  tmp <- data.frame(v1 = v1, v2 = v2) %>%
    stats::na.omit()
  if(nrow(tmp) < length(v1)) warning(paste(length(v1) - nrow(tmp),
					   " observations with missing values removed from data"))
  if(is.null(minrank)) minrank <- min(c(min(v), 1))
  tmp <- dplyr::mutate(tmp, dplyr::across(dplyr::where(is.factor), ~as.numeric(as.character(.))))
  v <- c(tmp$v1, tmp$v2)
  midranks <- purrr::map_dbl(rlang::set_names(minrank:maxrank), ~ sum(v < .) + 1 + sum(v == .)/2)
  props <- purrr::map(tmp, ~ prop.table(table(factor(., levels = minrank:maxrank))))
  marg.mean.rank = purrr::map(props, ~ sum(. * midranks))
  list(marginal.mean.ranks = unlist(marg.mean.rank),
       midranks = midranks)
}



