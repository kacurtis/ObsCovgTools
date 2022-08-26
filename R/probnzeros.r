#' Get probability of zero bycatch given effort
#' 
#' \code{probnzeros} returns probability of zero bycatch in a specified number 
#' of effort units, given bycatch per unit effort and dispersion index. 
#' 
#' @param n a vector of positive integers. Observed effort levels (in terms of 
#'   effort units, e.g., trips or sets) for which to calculate probability of 
#'   zero bycatch.
#' @param bpue a positive number. Bycatch per unit effort.
#' @param d a number greater than or equal to 1. Dispersion 
#'   index. The dispersion index corresponds to the variance-to-mean 
#'   ratio of effort-unit-level bycatch, so \code{d = 1} corresponds to Poisson-
#'   distributed bycatch, and \code{d > 1} corresponds to overdispersed bycatch.
#'   
#' @details
#' Calculated from the probability density at zero of the corresponding Poisson
#' (\code{d = 1}) or negative binomial (\code{d > 1}) distribution.
#' 
#' \strong{Caveat:} \code{probnzeros} assumes that (1) observer coverage is 
#' representative, (2) bycatch (\code{bpue}) is in terms of individuals (not 
#' weight) per unit effort, and (3) the specified dispersion index reflects 
#' the highest level of any hierarchical variance (e.g., using dispersion index 
#' at trip level if greater than that at set level). Violating these assumptions 
#' will likely result in negatively biased projections of the probability of 
#' observing zero bycatch at a given level of observer coverage. More conservative 
#' projections can be obtained by using a higher dispersion index \code{d}. Users 
#' may want to explore uncertainty in dispersion index and in bycatch per unit 
#' effort by varying those inputs.
#'   
#' @return Vector of same length as \code{n} with probabilities of zero bycatch. 
#' @return Returned invisibly
#' 
#' @export
probnzeros <- function(n, bpue, d) {
  
  # check input values
  if (any((ceiling(n) != floor(n)) | n<1)) stop("n must be a vector of positive integers")
  if (bpue<=0) stop("bpue must be > 0")
  if (d<1) stop("d must be >= 1")
  
  # calculate probability of observing zero bycatch in n units of effort
  ## one unit
  pz <- if(d==1) { exp(-1*bpue)
  } else { d^(-1*bpue/(d-1)) }
  ## n units
  pnz <- pz^n
  
  return(invisible(pnz))
}
