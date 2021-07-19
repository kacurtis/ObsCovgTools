# Hidden function to round up to specified significant digits
# Extended from code by JasonWang on stackoverflow at 
# https://stackoverflow.com/questions/37583715/round-up-values-to-a-specific-significant-figure-in-r
my_ceiling <- function(x, s){
  num_string <- format(x, scientific=TRUE)
  n <- strsplit(num_string, "e")
  n1 <- sapply(n, function(x) as.numeric(x[1]))
  n2 <- sapply(n, function(x) as.numeric(x[2]))
  ceiling(n1*10^(s-1))/(10^(s-1)) * 10^(n2)
}

# Hidden function to solve for upper one-tailed confidence limit (1-a) of bpue 
# when zero in observed effort n, given overdispersion d
solveucl <- function(a, d, n) {
  if (d==1) {
    bpue <- -1 * log(a)/n
  } else {
    bpue <- -1*(d-1)*log(a)/(n*log(d))
  }
  return(bpue)
}

