rm(list=ls())
setwd("~/Desktop/sectors/pharma/biogen/projects/study_design/");
source("logit_power_fun.R");
source("prev_counts.R");
library(plotly)
library(speedglm)


# trial 
set.seed(1234);
n1 = 392814;
n2 = 260405;

#n2 = 1e6;
#n1 = 1.5*5e5;

dta_z = data.frame(out = NA, exp = NA, z_ivw = NA, z_1 = NA, z_2 = NA, 
                 max_obv = NA, maf = NA, prev = NA, maf_scale = NA);
dta_se = data.frame(se_ratio_out = NA, se_ratio_exp = NA, 
                   maf = NA, prev = NA, maf_scale = NA);
dta_maf_monitor = data.frame(p1_ast = NA, p2_ast = NA, 
                    maf = NA, prev = NA, maf_scale = NA);
dta_boost_monitor = data.frame(lhs = NA, rhs = NA, 
                    maf = NA, prev = NA, maf_scale = NA);
dta_inflate_factor = dta_inflate_factor_ukb_meta = data.frame(alpha_obs = NA, alpha_exp = NA, 
                               maf = NA, prev = NA, maf_scale = NA);

mafs = c(1e-4,5e-4,1e-3,2.5e-3,5e-3, 1e-2);#
scale_maf = c(1,5,10,20,30,50);#c(1e-2,1e-1),
prevs = seq(0.005,0.005,by = 2e-2);
fg_prev_scale = 1.5
bta = betas[betas$ukb>0 & betas$fg>0,];
bta = bta[bta$ukb<=bta$fg,]
bta = bta[!bta$fg-bta$ukb>1,]
bta_length = dim(bta)[1];
count = 1;
trials = 1e3;
iter = 1;
full_analysis = TRUE;
meta_compare = TRUE;
if(meta_compare){
  n1 = n1 - n2;
}
  for(i in mafs){
    for(j in prevs){
      for(k in scale_maf){
        for(t in 1:trials){
          bta_loc = sample(x = 1:bta_length, size = 1, replace = T);
          bta_tmp = bta[bta_loc,1];
          bta_tmp2 = bta[bta_loc,2];#bta_tmp; #max(bta*i/maf2,0.1
          maf2 = i*(k);
          if(meta_compare){
            res3 = logit_run(n = n2, maf = i, pi_y = j, bta = bta_tmp, full_analysis = full_analysis);
          }
          res = logit_run(n = n1, maf = i, pi_y = j, bta = bta_tmp, full_analysis = full_analysis);
          res2 = logit_run(n = n2, maf = maf2, pi_y = j*fg_prev_scale, bta = bta_tmp2, full_analysis = full_analysis);
          if(sign(res["Z"])==sign(res2["Z"])){
            if(!full_analysis){
                      tmp1 = ratio_sig_z(z1 = res["Z"], se1 = res["se"], z2 = res2["Z"], se2 = res2["se"]);
                      tmp2 = exp_ratio(p1_ast = res["pi_y_ast"], n1_ast = res["N_ast"], pi1_cond = res["pi_cond"],
                                p2_ast = res2["pi_y_ast"], n2_ast = res2["N_ast"], pi2_cond = res2["pi_cond"]);
                      tmp3 = z_ivw(est1 = res["theta"], se1 = res["se"], est2 = res2["theta"], se2 = res2["se"]);
                      #tmp4 = max(abs(res["Z"]), abs(res2["Z"]));
                      tmp0 = which.max(c(abs(res["Z"]), abs(res2["Z"])))
                      dta_z[count, ] = c(tmp1, tmp2, tmp3, res["Z"], res2["Z"], tmp0, maf = i, prev = j, maf_scale = k);
                      
                        tmp5 = (res["se"]/res2["se"]);
                        tmp6 = se_ratio_exp(p1 = i, n1 = n1, pi1_cond = res["pi_cond"], p2 = maf2, n2 = n2, pi2_cond = res2["pi_cond"]);
                        dta_se[count, ] = c(se_ratio_out = tmp5, se_ratio_exp = tmp6, maf = i, prev = j, maf_scale = k);
                        
                        dta_maf_monitor[count, ] = c(p1_ast = res["pi_y_ast"], p2_ast = res2["pi_y_ast"], 
                                                     maf = i, prev = j, maf_scale = k)
                        
                        tmp7 = exp_ratio(p1_ast = res["pi_y_ast"], n1_ast = res["N_ast"], pi1_cond = 0.5,
                                         p2_ast = res2["pi_y_ast"], n2_ast = res2["N_ast"], pi2_cond = 0.5);
                        tmp_rhs = (tmp7*n1/n2)*res["pi_y_ast"]*(1-res["pi_y_ast"])/2;
                        tmp_lhs = res2["pi_y_ast"]*(1-res2["pi_y_ast"]);
                        
                        dta_boost_monitor[count, ] = c(lhs = tmp_lhs, rhs = tmp_rhs, 
                                                     maf = i, prev = j, maf_scale = k)
                        
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
          if(t%%trials==0){print(paste0("iteration ",iter," complete")); iter = iter+1}
        }
      }
    }
  }

rm_obs = is.na(dta_inflate_factor$alpha_exp) | dta_inflate_factor$alpha_exp==Inf 
plot(1:dim(dta_inflate_factor[!rm_obs,])[1], dta_inflate_factor$alpha_obs[!rm_obs], ylim = c(0,50));
lines(1:dim(dta_inflate_factor[!rm_obs,])[1], dta_inflate_factor$alpha_exp[!rm_obs]);
cor(dta_inflate_factor$alpha_obs[!rm_obs], dta_inflate_factor$alpha_exp[!rm_obs] , method = "spearman");
round(dta_inflate_factor$alpha_obs- dta_z$z_ivw/dta_z$z_1,3);


if(meta_compare){
  #dta_inflate_factor_stored = dta_inflate_factor;
  dta_inflate_factor_new = dta_inflate_factor;
  dta_inflate_factor_new$alpha_obs = dta_inflate_factor$alpha_obs/dta_inflate_factor_ukb_meta$alpha_obs;
  dta_inflate_factor_new$alpha_exp = dta_inflate_factor$alpha_exp/dta_inflate_factor_ukb_meta$alpha_exp;
  res = quantile_sums(dta_inflate_factor_new)
  res = quantile_sums(dta_inflate_factor_ukb_meta)
}else{
  res = quantile_sums(dta_inflate_factor)
}



plot(1:dim(dta_z)[1], dta_z$out);
lines(1:dim(dta_z)[1], dta_z$exp);
dta_z
cor(dta_z$out, dta_z$exp , method = "spearman");
dta_z[dta_z$z_ivw<dta_z$z_strg,]

    plot(1:dim(dta_se)[1], dta_se$se_ratio_out);
    lines(1:dim(dta_se)[1], dta_se$se_ratio_exp);
    dta_se
    cor(dta_se$se_ratio_out, dta_se$se_ratio_exp , method = "spearman");

      plot(dta_maf_monitor$maf_scale, var_ratio(dta_maf_monitor$p1_ast, dta_maf_monitor$p2_ast));

    
    
    
    
    mat_obs = t(matrix(res$obsrved$median, nrow = 6, ncol = 6))
    #mat_obs = mat_obs[1:4,];
    fig_obs <- plot_ly(type = 'surface',
                       contours = list(
                         #x = list(show = TRUE, start = 1.5, end = 2, size = 0.04, color = 'white'),
                         y = list(show = TRUE, start = 1, end = 50, size = 5, color = '#2F4F4F'),
                         z = list(show = TRUE, start = 1, end = 10, size = 1, color = '#C0C0C0')),
                       x = unique(res$obsrved$maf), y = unique(res$obsrved$maf_scale), z = t(mat_obs))
    fig_obs = fig_obs %>% layout(autosize = T, 
                       yaxis = list(range = c(0, 50)), 
                       xaxis = list(range = c(1e-4, 1e-2)),
                       scene = list(
                         xaxis = list(title = "MAF",autotick = F,
                                      ticktext = list("1e-4", "5e-4","1e-3", "0.0025", "0.005", "0.01"), 
                                      tickvals = list(1e-4,5e-4,1e-3, 2.5e-3,5e-3, 1e-2),
                                      tickmode = "array",
                                      tickfont = list(size = 12),
                                      titlefont = list(size = 15)),
                         yaxis = list(title = "Fold-enrichment", autotick = F,
                                      ticktext = list("1", "5","10", "20", "30","50"), 
                                      tickvals = list(1,5,10,20, 30,50),
                                      tickmode = "array",
                                      tickfont = list(size = 12),
                                      titlefont = list(size = 15)),
                         zaxis = list(title = "IVW uplift",
                                      titlefont = list(size = 15)),
                         title = "Expected uplift",
                         camera = list(eye=list(x=-1.75,y=-1.75,z=.5))
                       )) 
    fig_obs
    
    mat_exp = t(matrix(res$expected$median, nrow = 6, ncol = 6))
    #mat_exp = mat_exp[1:4,]
    fig_exp <- plot_ly(type = 'surface',
                       contours = list(
                         #x = list(show = TRUE, start = 1e-4, end = 5e-4, size = 2.5e-4, color = 'white'),
                         y = list(show = TRUE, start = 1, end = 50, size = 5, color = '#2F4F4F'),
                         z = list(show = TRUE, start = 1, end = 10, size = 1, color = '#C0C0C0')),
                       x = unique(res$obsrved$maf), y = unique(res$obsrved$maf_scale), z = t(mat_exp))
    fig_exp = fig_exp %>% layout(autosize = T,                   yaxis = list(range = c(0, 50)), 
                       yaxis = list(range = c(0, 50)), 
                       xaxis = list(range = c(1e-4, 5e-3)),
                       scene = list(
                         xaxis = list(title = "MAF",autotick = F,
                                      ticktext = list("1e-4", "5e-4","1e-3", "0.0025", "0.005", "0.01"), 
                                      tickvals = list(1e-4,5e-4,1e-3, 2.5e-3,5e-3, 1e-2),
                                      tickmode = "array",
                                      tickfont = list(size = 12),
                                      titlefont = list(size = 15)),
                         yaxis = list(title = "Fold-enrichment", autotick = F,
                                      ticktext = list("1", "5","10", "20", "30","50"), 
                                      tickvals = list(1,5,10,20, 30,50),
                                      tickmode = "array",
                                      tickfont = list(size = 12),
                                      titlefont = list(size = 15)),
                         zaxis = list(title = "IVW uplift",
                                      titlefont = list(size = 15)),
                         title = "Expected uplift",
                         camera = list(eye=list(x=-1.75,y=-1.75,z=.5))
                       )) %>% hide_colorbar()  
    fig_exp;

    
    dta_inflate_factor_diff = dta_inflate_factor;
    dta_inflate_factor_diff$alpha_obs = abs(dta_inflate_factor_diff$alpha_obs - dta_inflate_factor_diff$alpha_exp)/abs(dta_inflate_factor_diff$alpha_obs);
    res_diff = quantile_sums(dta_inflate_factor_diff)
    mat_diff = t(matrix(res_diff$obsrved$median, nrow = 6, ncol = 6))*100
    #mat_exp = mat_exp[1:4,]
    fig_diff <- plot_ly(type = 'surface',
                       contours = list(
                         #x = list(show = TRUE, start = 1e-4, end = 5e-4, size = 2.5e-4, color = 'white'),
                         y = list(show = TRUE, start = 1, end = 50, size = 5, color = '#2F4F4F'),
                         z = list(show = TRUE, start = 1, end = 10, size = 1, color = '#C0C0C0')),
                       x = unique(res$obsrved$maf), y = unique(res$obsrved$maf_scale), z = t(mat_diff))
    fig_diff = fig_diff %>% layout(autosize = T,                   yaxis = list(range = c(0, 50)), 
                                 yaxis = list(range = c(0, 50)), 
                                 xaxis = list(range = c(1e-4, 5e-3)),
                                 scene = list(
                                   xaxis = list(title = "MAF",autotick = F,ticks = "outside",
                                                ticktext = list("1e-4", "5e-4","1e-3", "0.0025", "0.005", "0.01"), 
                                                tickvals = list(1e-4,5e-4,1e-3, 2.5e-3,5e-3, 1e-2),
                                                tickmode = "array",
                                                tickfont = list(size = 12),
                                                titlefont = list(size = 15)),
                                   yaxis = list(title = "Fold-enrichment", autotick = F,
                                                ticktext = list("1", "5","10", "20", "30","50"), 
                                                tickvals = list(1,5,10,20, 30,50), range=c(0,60),
                                                tickmode = "array",
                                                tickfont = list(size = 12),
                                                titlefont = list(size = 15)),
                                   zaxis = list(title = "Absolute relative error %", ticks = "outside", range=c(0,100),
                                                ticktext = list("", "20","40", "60", "80","100"), 
                                                tickvals = list(0,20,40,60,80,100),
                                                tickfont = list(size = 12),
                                                titlefont = list(size = 15)),
                                   title = "Absolute relative error (%)",
                                   camera = list(eye=list(x=0.25,y=-2.125,z=1))
                                 )) 
    fig_diff;
    
    library(gridExtra)
    library(grid)
    library(png)
    
    if(meta_compare){
      htmlwidgets::saveWidget(as_widget(fig_obs), "observed_ivw_uplift_ukb_vs_fg_meta.html")
      htmlwidgets::saveWidget(as_widget(fig_exp), "expected_ivw_uplift_ukb_vs_fg_meta.html")
      htmlwidgets::saveWidget(as_widget(fig_diff), "rel_err_ivw_uplift_ukb_vs_fg_meta.html")
      plot1 <- readPNG('submission_2/exp_rel_meta_uplift_maf_1e4_1e2_1000_trials.png')
      plot2 <- readPNG('submission_2/obs_rel_meta_uplift_maf_1e4_1e2_1000_trials.png')
      plot3 <- readPNG('submission_2/rel_err_meta_uplift_maf_1e4_1e2_1000_trials.png')
      plot4 <- readPNG("submission_2/log10_exp_obs_uplift_both_meta_maf_1e4_1e2_1000_trials")
      grid.arrange(rasterGrob(plot1),rasterGrob(plot2), rasterGrob(plot3), rasterGrob(plot4),
                   layout_matrix = rbind(c(1, 2), c(3, 4)),
                   ncol=2, widths = c(0.975,0.9), heights = c(5,6), 
                   top=textGrob("",
                                gp=gpar(fontsize=8,font=8)))
    }else{
      htmlwidgets::saveWidget(as_widget(fig_obs), "observed_ivw_uplift.html")
      htmlwidgets::saveWidget(as_widget(fig_exp), "expected_ivw_uplift.html")
      htmlwidgets::saveWidget(as_widget(fig_diff), "rel_err_ivw_uplift.html")
      plot1 <- readPNG('exp_uplift_maf_1e4_1e2_1000_trials.png')
      plot2 <- readPNG('obs_uplift_maf_1e4_1e2_1000_trials.png')
      plot3 <- readPNG('rel_err_maf_1e4_1e2_1000_trials.png')
      grid.arrange(rasterGrob(plot1),rasterGrob(plot2), rasterGrob(plot3),
                   layout_matrix = rbind(c(1, 2), c(3, 3)),
                   ncol=2, widths = c(0.975,1), heights = c(5,6), 
                   top=textGrob("MAF enrichment.\nTheoretical IVW uplift (left), simulated IVW uplift (right),\n% absolute relative error (bottom)",
                                gp=gpar(fontsize=8,font=8)))
    }
    
#    gridExtra::grid.arrange(fig_obs, fig_exp)
    
    df = quantile_sums(dta_inflate_factor_ukb_meta, rng = c(0.25,0.5,0.75));
    df_obs = df$obsrved
      df_obs$Result = "Reference UKB_1/UKB_2\n(simulated)"
    df_exp = df$expected
      df_exp$Result = "Reference UKB_1/UKB_2\n(theoretical)"  
    df_res = rbind(df_obs, df_exp);  
    df_res$shape = "no_enrichment";
    #df_res$first[1] = -4.99;
    
    df = quantile_sums(dta_inflate_factor_new, rng = c(0.25,0.5,0.75));
    df_obs = df$obsrved[df$obsrved$maf_scale==5,];
    df_obs$Result = "5-fold enrichment\n(simulated)"
    df_obs$shape = "Enriched";
    df_exp = df$expected[df$expected$maf_scale==5,]
    df_exp$Result = "5-fold enrichment\n(theoretical)"  
    df_exp$shape = "Enriched";
    df_res = rbind(df_res,df_obs, df_exp);  
    df_res = df_res[df_res$maf>=0.001,]
    df_res[,1:3] = log10(df_res[,1:3])
    
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    ggplot(data = df_res, aes(x=as.factor(maf), y=median, color=Result)) +
      geom_point(size = 2, position=position_dodge(width=0.6), aes(shape = Result)) +
      geom_errorbar(
        aes(ymin=first, ymax=ninth),
        width = 0.3,
        linetype = "solid",
        position=position_dodge(width=0.6)) +
      geom_hline(aes(yintercept = 0), linetype = "dashed") + 
      scale_y_continuous(limits = c(-0.7,0.7), 
                         breaks = c(-0.4, 0, 0.4)) + 
      xlab("MAF") + ylab("log10(FG/UKB meta Z-score ratio- \"IVW uplift\")") +
      scale_color_manual(values=cbPalette, name="Result", labels=c("5-fold enrichment\n(simulated)",
                                                                                      "5-fold enrichment\n(predicted)",
                                                                   "Reference analysis\nUKB_1 + UKB_2\n(simulated)",
                                                                   "Reference analysis\nUKB_1 + UKB_2\n(predicted)")) +
      scale_shape_manual(values=c(16,16,24,24), name  ="Result", labels=c("5-fold enrichment\n(simulated)",
                                                                          "5-fold enrichment\n(predicted)",
                                                                          "Reference analysis\nUKB_1 + UKB_2\n(simulated)",
                                                                          "Reference analysis\nUKB_1 + UKB_2\n(predicted)")) +
      theme_bw() +
      theme(legend.key.size = unit(1.5, 'cm'),
            text = element_text(size=10)) 
    
    ggsave("log10_exp_obs_uplift_both_meta_maf_1e4_1e2_1000_trials", device = "png", 
           dpi = 400, width = 5,height = 4.5, units = "in")
    
    #geom_segment(aes(x = 0.875, xend = 0.875, y = -4.85, yend = -5), colour = "#F8766D", arrow = arrow(length = unit(0.3, "cm")))