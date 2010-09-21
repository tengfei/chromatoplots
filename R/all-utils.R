mav <- function(x, na.rm = TRUE) median(abs(x), na.rm = na.rm) / 0.6745

greedy_match <- function(distances) {
  dist_col <- col(distances)
  dist_row <- row(distances)
  dist_order <- order(distances)  
  matching <- NULL
  tmp_dist <- distances
  sapply(dist_order, function(index) {
    if (!is.na(tmp_dist[index])) {
      a <- dist_col[index]
      b <- dist_row[index]
      matching <<- rbind(matching, cbind(a = a, b = b))
      tmp_dist[,a] <<- NA
      tmp_dist[b,] <<- NA
    }
  })
  return(matching)
}

which.median <- function(x) order(x)[ceiling(length(x)/2)]

# some stuff from XCMS (hidden by its namespace)

descendMin <- function(y, istart = which.max(y)) {
    if (!is.double(y)) y <- as.double(y)
    unlist(.C("DescendMin",
              y,
              length(y),
              as.integer(istart-1),
              ilower = integer(1),
              iupper = integer(1),
              DUP = FALSE, PACKAGE = "xcms")[4:5]) + 1
}
na.flatfill <- function(x) {

    realloc <- which(!is.na(x))
    if (realloc[1] > 1)
        x[1:(realloc[1]-1)] <- x[realloc[1]]
    if (realloc[length(realloc)] < length(x))
        x[(realloc[length(realloc)]+1):length(x)] <- x[realloc[length(realloc)]]

    x
}

## hack to get .local formals out of method definitions

getMethodLocal <- function(...) {
  eval(body(getMethod(...))[[2]])
  .local
}
