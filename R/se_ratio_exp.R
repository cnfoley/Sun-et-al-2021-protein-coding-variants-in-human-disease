#' Expected value of standard error ratios
#'
#' @param p1 
#' @param n1 
#' @param pi1_cond 
#' @param p1 
#' @param n1 
#' @param pi1_cond 
#' @return Return expected value of standard error ratios
#' @export
se_ratio_exp = function(p1, n1, pi1_cond, p2, n2, pi2_cond){
  return(sqrt(p2*(1-p2)*n2*pi1_cond*(1-pi1_cond)/(p1*(1-p1)*n1*pi2_cond*(1-pi2_cond))))
}
