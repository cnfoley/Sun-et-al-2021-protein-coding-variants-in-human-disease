#' Inflation factor
#'
#' @param sig_z 
#' @param se_ratio 
#' @return Returned is denominator of IVW-uplift.
#' @export
inflate_factor = function(sig_z, se_ratio){
  return((1+sig_z)/sqrt(1+(se_ratio)^2))
}
