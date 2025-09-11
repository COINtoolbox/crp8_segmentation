filter_by_size <- function(L, min_size = 30L, max_size = Inf) {
  if (!any(L > 0)) return(L)
  tb   <- table(L[L > 0])
  keep <- as.integer(names(tb)[tb >= min_size & tb <= max_size])
  L[!(L %in% keep)] <- 0L
  compact_labels(L)
}