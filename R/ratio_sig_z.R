#' IVW z-score
#'
#' @param z1 z-score study 1
#' @param se1 standard error study 1
#' @param z2 a-score study 2
#' @param se2 standard error study 2 
#' @return IVW z-score.
#' @export
ratio_sig_z = function(z1, se1, z2, se2){
  return(se1*z2/(se2*z1))
}