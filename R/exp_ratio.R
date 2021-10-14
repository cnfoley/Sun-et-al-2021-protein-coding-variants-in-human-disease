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
exp_ratio = function(p1_ast, n1_ast, pi1_cond, p2_ast, n2_ast, pi2_cond){
  if(p1_ast==0 & p2_ast==0){
    return(n2_ast*pi1_cond*(1-pi1_cond)/(n1_ast*(pi2_cond*(1-pi2_cond))))
  }else{
    return(p2_ast*n2_ast*pi1_cond*(1-pi1_cond)/(p1_ast*n1_ast*(pi2_cond*(1-pi2_cond))))
  }
}
