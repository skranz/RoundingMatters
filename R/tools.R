has.col = function(x, col) {
  col %in% colnames(x)
}

is.true = function(x) {
  is.true = rep(FALSE,length(x))
  is.true[which(x==TRUE)] = TRUE
  is.true
}

quick.df = function (...) {
    df = list(...)
    attr(df, "row.names") <- 1:length(df[[1]])
    attr(df, "class") <- "data.frame"
    df
}
