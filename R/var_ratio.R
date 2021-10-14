#' Ratio of variance of binomial RVs
#'
#' @param p1 
#' @param p2 
#' @return Returned is denominator of IVW-uplift.
#' @export
var_ratio = function(p1, p2){
  return(p1*(1-p1)/(p2*(1-p2)))
}
