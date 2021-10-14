#' Quantile computation
#'
#' @param data list taken from allelic heterogeneity simulation.
#' @param values numeric vector denoting the set of quantile summary values to record. Default is {0.1,0.5,0.9}, i.e., the 1st, 5th (median) and 9th deciles.
#' @return Returned are user defined quantile range of simulated and theoretically predicted values of IVW-uplift generated from "n_simulated_datasets" per combination of MAF and enrichment value.
#' @export
quantile_sums = function(data, values = c(0.1,0.5,0.9)){
  rm_obs = is.na(data[,2]) | data[,2]==Inf;
  data = data[!rm_obs,];
  maf = data$maf;
  maf_scale=data$maf_scale;
  obs = data[,1];
  exp = data[,2];
  obs_dta = data.frame(first = NA, median = NA, ninth = NA, maf = NA, maf_scale = NA)
  exp_dta = data.frame(first = NA, median = NA, ninth = NA, maf = NA, maf_scale = NA)
  count = 1;
  for(i in unique(maf)){
    for(j in unique(maf_scale)){
      obs_dta[count,] = c(quantile(obs[maf==i & maf_scale==j], values, na.rm = T), i, j)
      exp_dta[count,] = c(quantile(exp[maf==i & maf_scale==j], values, na.rm = T), i, j)
      count = count + 1;
    }
  }
  return(list(obsrved = obs_dta, expected = exp_dta))
}
