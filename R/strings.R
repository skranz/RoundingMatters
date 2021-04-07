str.left.of = function (str, pattern, ..., not.found = str) {
  pos = regexpr(pattern, str, fixed=TRUE)
  res = substring(str, 1, pos - 1)
  rows = (pos == -1)
  res[rows] = not.found[rows]
  res
}

str.right.of = function (str, pattern, ..., not.found = str) {
  pos = regexpr(pattern, str, fixed=TRUE)
  res = substring(str,pos+nchar(pattern))
  rows = (pos == -1)
  res[rows] = not.found[rows]
  res
}

str.between = function (str, start, end, ...)
{
    str.left.of(str.right.of(str, start, ...), end, ...)
}
