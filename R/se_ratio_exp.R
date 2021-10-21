#' Expected value of standard error ratios
#'
#' @param p1 maf study 1
#' @param n1 sample size study 1
#' @param pi1_cond probability of disease conditional on expected value of genetic variant in study 1
#' @param p2 maf study 2
#' @param n2 sample size study 2 
#' @param pi2_cond probability of disease conditional on expected value of genetic variant in study 2
#' @return Return expected value of standard error ratios
#' @export
se_ratio_exp = function(p1, n1, pi1_cond, p2, n2, pi2_cond){
  return(sqrt(p2*(1-p2)*n2*pi1_cond*(1-pi1_cond)/(p1*(1-p1)*n1*pi2_cond*(1-pi2_cond))))
}
