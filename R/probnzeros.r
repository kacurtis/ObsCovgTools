#' Get probability of zero bycatch given effort
#' 
#' \code{probnzeros} returns probability of zero bycatch in a specified number 
#' of trips or sets, given bycatch per unit effort and dispersion index. 
#' 
#' @param n a vector of positive integers. Observed effort levels (in terms of 
#'   trips or sets) for which to calculate probability of zero bycatch.
#' @param bpue a positive number. Bycatch per unit effort.
#' @param d a number greater than or equal to 1. Dispersion 
#'   index. The dispersion index corresponds to the variance-to-mean 
#'   ratio of effort-unit-level bycatch, so \eqn{d = 1} corresponds to Poisson-
#'   distributed bycatch, and \eqn{d > 1} corresponds to overdispersed bycatch.
#'   
#' @details
#' Calculated from the probability density at zero of the corresponding Poisson
#' (\eqn{d = 1}) or negative binomial (\eqn{d < 1}) distribution.
#' 
#' \strong{Caveat:} \code{probnzeros} assumes representative observer coverage 
#' and no hierarchical sources of variance (e.g., vessel- or trip-level variation). 
#' Violating these assumptions will likely result in negatively biased projections 
#' of the probability of observing zero bycatch at a given level of observer coverage. 
#' More conservative projections can be obtained by using higher-level units of effort 
#' (e.g., \code{bpue} as mean bycatch per trip instead of bycatch per set, and 
#' \code{n} as number of trips instead of number of sets). Landings (number or 
#' weight) do not represent independent sampling units unless an average independent 
#' fishing effort unit such as a trip lands less than or equal to one measurement 
#' unit of landings (e.g., <= one metric ton per trip).
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
