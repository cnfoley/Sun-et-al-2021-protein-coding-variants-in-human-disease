#' Exponent ratio
#'
#' @param p1_ast maf in cases only, study 1
#' @param n1_ast number of cases in study 1
#' @param pi1_cond probability of disease conditional on expected value of genetic variant in study 1
#' @param p2_ast maf in cases only, study 2
#' @param n2_ast number of cases in study 2
#' @param pi2_cond probability of disease conditional on expected value of genetic variant in study 2
#' @return Return expected value of standard error ratios
#' @export
exp_ratio = function(p1_ast, n1_ast, pi1_cond, p2_ast, n2_ast, pi2_cond){
  if(p1_ast==0 & p2_ast==0){
    return(n2_ast*pi1_cond*(1-pi1_cond)/(n1_ast*(pi2_cond*(1-pi2_cond))))
  }else{
    return(p2_ast*n2_ast*pi1_cond*(1-pi1_cond)/(p1_ast*n1_ast*(pi2_cond*(1-pi2_cond))))
  }
}
