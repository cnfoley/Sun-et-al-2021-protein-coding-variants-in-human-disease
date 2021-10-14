#' IVW z-score
#'
#' @param est1 
#' @param se1 
#' @param est2 
#' @param se2 
#' @return IVW z-score.
#' @export
ratio_sig_z = function(z1, se1, z2, se2){
  return(se1*z2/(se2*z1))
}