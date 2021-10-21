#' IVW z-score
#'
#' @param est1 regression estimate study 1
#' @param se1 regression std error study 1
#' @param est2 regression estimate study 2
#' @param se2 regression std error study 2
#' @return IVW z-score.
#' @export
z_ivw = function(est1, se1, est2, se2){
  tmp1 = est1/se1^2 + est2/se2^2;
  tmp2 = 1/se1^2 + 1/se2^2;
  return(tmp1/sqrt(tmp2))
}
