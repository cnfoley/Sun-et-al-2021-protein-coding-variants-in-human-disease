#' Allelic heterogeneity simulation and calculation of IVW uplift
#'
#' Simulation protocol to generate and estimate IVW-uplift owing to 
#' allelic-heterogeneity between 2 population cohorts.
#'
#' @param betas data.frame whose columns denote the regression coefficients from study 1 (UKB) and separately study 2 (FinGen) studies.
#' @param n1 numeric variable denoting the sample size of study 1. Default set to UKB sample size.
#' @param n2 numeric variable denoting the sample size of study 2. Default set to FG sample size.
#' @param maf_range numeric vector containing MAF values for study 1. Default considers a range of MAFs starting at the very rare MAF 0.01% and sequentially increasing to moderate rare 1%.
#' @param enrichment_value numeric vector denoting the allelic enrichment of study 2. Default considers a range of allelic enrichment values starting at no enrichment (i.e, a value of 1) to 50-fold enrichment.
#' @param disease_prevalence numeric vector denoting disease prevalence in study 1. Default is 0.5% which is the same order of mangnitude as the median disease prevalence across the 975 diseases in UKB.
#' @param disease_prevalence_scale_stdy2 numeric scalar denoting the scaled increase/decrease of disease prevalence in study 2. Default is 1.5, i.e,. 0.75%.  
#' @param n_simulated_datasets numeric scalar denoting the number of datasets to simulate for each value of MAF and enrichment value.  
#' @param effect_direction character denoting the direction of association (i.e., risk increasing or decreasing). Default is risk increasing 
#' @param meta_compare logical denoting whether we perform the primary meta-analysis (i.e., the meta-analysis which combines 2 studies from different populations) or to compare the relative effects of two meta-analyses: (i) 2 studies from different populations and (ii) 2 studies from the same population. Default is FALSE, i.e, to perform our primary meta-analysis.
#' @param quantile_range numeric vector denoting the set of quantile summary values to record. Default is {0.1,0.5,0.9}, i.e., the 1st, 5th (median) and 9th deciles.
#' @return Returned are results from the simulated and theoretically predicted values of IVW-uplift generated from "n_simulated_datasets" per combination of MAF and enrichment value.
#' @export
allelic_het_simulator <- function(betas, n1 = 392814, n2 = 260405, maf_range = c(1e-4,5e-4,1e-3,2.5e-3,5e-3, 1e-2),
                                   enrichment_value = c(1,5,10,20,30,50), disease_prevalence = 0.005, disease_prevalence_scale_stdy2 = 1.5,
                                  n_simulated_datasets = 1000, effect_direction = "increasing", meta_compare = FALSE, quantile_range = c(0.1,0.5,0.9)){

# trial 
set.seed(1234);

dta_inflate_factor = dta_inflate_factor_ukb_meta = data.frame(alpha_obs = NA, alpha_exp = NA, 
                               maf = NA, prev = NA, maf_scale = NA);
if(effect_direction=="increasing"){
        bta = betas[betas$ukb>0 & betas$fg>0,];
        #bta = bta[bta$ukb<=bta$fg,]
        #bta = bta[!bta$fg-bta$ukb>1,]
}else{
        bta = betas[betas$ukb<0 & betas$fg<0,];
        #bta = bta[bta$ukb>=bta$fg,]
        #bta = bta[!bta$fg-bta$ukb>1,]
        }
bta_length = dim(bta)[1];
count = 1;
iter = 1;
full_analysis = TRUE;
if(meta_compare){
  n1 = n1 - n2;
}
  for(i in maf_range){
    for(j in disease_prevalence){
      for(k in enrichment_value){
        for(t in 1:n_simulated_datasets){
          bta_loc = sample(x = 1:bta_length, size = 1, replace = T);
          bta_tmp = bta[bta_loc,1];
          bta_tmp2 = bta[bta_loc,2];#bta_tmp; #max(bta*i/maf2,0.1
          maf2 = i*(k);
          if(meta_compare){
            res3 = logit_run(n = n2, maf = i, pi_y = j, bta = bta_tmp);
          }
          res = logit_run(n = n1, maf = i, pi_y = j, bta = bta_tmp, full_analysis = full_analysis);
          res2 = logit_run(n = n2, maf = maf2, pi_y = j*disease_prevalence_scale_stdy2, bta = bta_tmp2, full_analysis = full_analysis);
          if(sign(res["Z"])==sign(res2["Z"])){
            if(!full_analysis){
                      tmp1 = ratio_sig_z(z1 = res["Z"], se1 = res["se"], z2 = res2["Z"], se2 = res2["se"]);
                      tmp2 = exp_ratio(p1_ast = res["pi_y_ast"], n1_ast = res["N_ast"], pi1_cond = res["pi_cond"],
                                p2_ast = res2["pi_y_ast"], n2_ast = res2["N_ast"], pi2_cond = res2["pi_cond"]);
                      tmp3 = z_ivw(est1 = res["theta"], se1 = res["se"], est2 = res2["theta"], se2 = res2["se"]);
                      #tmp4 = max(abs(res["Z"]), abs(res2["Z"]));
                      tmp0 = which.max(c(abs(res["Z"]), abs(res2["Z"])))
                      
                        tmp5 = (res["se"]/res2["se"]);
                        tmp6 = se_ratio_exp(p1 = i, n1 = n1, pi1_cond = res["pi_cond"], p2 = maf2, n2 = n2, pi2_cond = res2["pi_cond"]);
                        tmp7 = exp_ratio(p1_ast = res["pi_y_ast"], n1_ast = res["N_ast"], pi1_cond = 0.5,
                                         p2_ast = res2["pi_y_ast"], n2_ast = res2["N_ast"], pi2_cond = 0.5);
                        tmp_rhs = (tmp7*n1/n2)*res["pi_y_ast"]*(1-res["pi_y_ast"])/2;
                        tmp_lhs = res2["pi_y_ast"]*(1-res2["pi_y_ast"]);
                        inf_factor_exp = inflate_factor(sig_z = tmp2, se_ratio = tmp6);
                        inf_factor_obs = inflate_factor(sig_z = tmp1, se_ratio = tmp5);
                        dta_inflate_factor[count, ] = c(alpha_obs = inf_factor_obs, alpha_exp = inf_factor_exp, 
                                                        maf = i, prev = j, maf_scale = k)
            }else{
              tmp1 = ratio_sig_z(z1 = res["Z"], se1 = res["se"], z2 = res2["Z"], se2 = res2["se"]);
              tmp2 = (res["se"]/res2["se"]);
                inf_factor_obs = inflate_factor(sig_z = tmp1, se_ratio = tmp2);
              tmp3 = (res2["mn_y"]*(1-res2["mn_y"])/(res["mn_y"]*(1-res["mn_y"])))*ratio_sig_z(z1 = res["Zpred"], se1 = res["SEpred"], z2 = res2["Zpred"], se2 = res2["SEpred"]);
              tmp4 = (res2["mn_y"]*(1-res2["mn_y"])/(res["mn_y"]*(1-res["mn_y"])))*(res["SEpred"]/res2["SEpred"]);
              inf_factor_exp = inflate_factor(sig_z = tmp3, se_ratio = tmp4);
              
                if(meta_compare){
                  ### UKB meta with UKB
                  tmp1 = ratio_sig_z(z1 = res["Z"], se1 = res["se"], z2 = res3["Z"], se2 = res3["se"]);
                  tmp2 = (res["se"]/res3["se"]);
                  inf_factor_obs_ukb_meta = inflate_factor(sig_z = tmp1, se_ratio = tmp2);
                  tmp3 = (res3["mn_y"]*(1-res3["mn_y"])/(res["mn_y"]*(1-res["mn_y"])))*ratio_sig_z(z1 = res["Zpred"], se1 = res["SEpred"], z2 = res3["Zpred"], se2 = res3["SEpred"]);
                  tmp4 = (res3["mn_y"]*(1-res3["mn_y"])/(res["mn_y"]*(1-res["mn_y"])))*(res["SEpred"]/res3["SEpred"]);
                  inf_factor_exp_ukb_meta = inflate_factor(sig_z = tmp3, se_ratio = tmp4);
                  dta_inflate_factor_ukb_meta[count, ] = c(alpha_obs = inf_factor_obs_ukb_meta, alpha_exp = inf_factor_exp_ukb_meta, 
                                                  maf = i, prev = j, maf_scale = 1);
                }
              
              dta_inflate_factor[count, ] = c(alpha_obs = inf_factor_obs, alpha_exp = inf_factor_exp, 
                                              maf = i, prev = j, maf_scale = k)
              
            }
              
            count = count + 1;
          }
          if(t%%n_simulated_datasets==0){print(paste0("iteration ",iter," complete")); iter = iter+1}
        }
      }
    }
  }

rm_obs = is.na(dta_inflate_factor$alpha_exp) | dta_inflate_factor$alpha_exp==Inf 

if(meta_compare){
  #dta_inflate_factor_stored = dta_inflate_factor;
  dta_inflate_factor_new = dta_inflate_factor;
  dta_inflate_factor_new$alpha_obs = dta_inflate_factor$alpha_obs/dta_inflate_factor_ukb_meta$alpha_obs;
  dta_inflate_factor_new$alpha_exp = dta_inflate_factor$alpha_exp/dta_inflate_factor_ukb_meta$alpha_exp;
  dta_inflate_factor = dta_inflate_factor_new;
}
return(data = dta_inflate_factor)
}
    