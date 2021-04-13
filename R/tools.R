has.col = function(x, col) {
  col %in% colnames(x)
}

is.true = function(x) {
  is.true = rep(FALSE,length(x))
  is.true[which(x==TRUE)] = TRUE
  is.true
}
