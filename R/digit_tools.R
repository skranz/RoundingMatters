# Some helper function to extract information
# related to the number of significant digits
# and rounding

examples.digit.tools = function(x) {
  x = c(1234567890,0.1234567890,10.10,0.123)
  significand(x)
  num.sig.digits(x)
  set.last.digit.zero(x)
  num.deci(x)
  rightmost.sig.digit(x, 1,1)
}


#' Get the significands of a numeric vector using
#' the character presentation of those numbers of R
#'
#' The significand is the integer of all significand digits, e.g.
#' the significand of 0.012 is 12.
#'
#' @param x a numeric vector.
#' @param num.deci If not NULL a vector that states the number reported decimal places for x. This can be used if we know that there were addtional trailing zeros.
#' @export
significand = function(x, num.deci=NULL) {
  options(scipen=999)
  if (!is.null(num.deci)) {
    return(round(x*10^num.deci))
  }

  #format(x,digits=max.digits,scientific = FALSE,justify = "none")
  str = as.character(x)
  str = gsub(".","",str,fixed=TRUE)
  #str = gsub("-","",str,fixed=TRUE)
  as.numeric(str)
}

#' Get the number of significand digits of a floating point number
#' using the character presentation of those numbers of R
#'
#' We assume that trailing zeros left of the decimal point are
#' significant digits while trailing zeros right of the decimal point
#' are not significant digits
#'
#' @param x a numeric vector
#' @export
num.sig.digits = function(x) {
  options(scipen=999)
  str = as.character(x)
  str = gsub(".","",str,fixed=TRUE)
  str = gsub("-","",str,fixed=TRUE)
  # Replace leading zeros
  str = gsub("^0+", "",str)
  nchar(str)
}

#' Get the number of significand digits of a floating point number
#' using the character presentation of those numbers of R
#'
#' @param x a numeric vector
#' @export
num.deci = function(x) {
  options(scipen=999)
  str = as.character(x)
  pos.dec = regexpr(pattern ='.',str, fixed=TRUE)
  res = nchar(str)-pos.dec
  res[pos.dec==-1] = 0
  as.integer(res)
}

str.right.of = function (str, pattern, ..., not.found = str) {
    pos = str.locate.first(str, pattern, ...)
    res = substring(str, pos[, 2] + 1, )
    rows = is.na(pos[, 2])
    res[rows] = not.found[rows]
    res
}

#' Get the last significant digit(s) of a floating point number
#'
#' @param x The vector of floating point numbers
#' @param r1 Starting position from right
#' @param r2 Ending position from right
#' @export
rightmost.sig.digit = function(x,r1=1,r2=1) {
  x = as.character(x)
  x = gsub(".","",x,fixed=TRUE)
  x = gsub("-","",x,fixed=TRUE)
  # replace leading zeros
  x = gsub("^0+", "",x)
  len = nchar(x)
  res = substring(x,len-r1+1,len-r2+1)
  rows = (len-r1+1 < 1)
  res[rows] = NA
  res
}



#' Convert numbers like 0.421 to 42.1\%
#'
#' @param x a vector of floating point numbers
#' @param digits to how many decimal digits shall the percentage be rounded?
#' @export
as.perc = function(x, digits=1) {
  paste0(round(x*100,digits),"%")
}

#' Sets the last digit of a number x to zero
#'
#' @param x a numeric vector
#'
#' @export
set.last.digit.zero = function(x) {
  str = as.character(x)
  nstr = paste0(substring(x,1,nchar(str)-1),"0")
  as.numeric(nstr)
}

