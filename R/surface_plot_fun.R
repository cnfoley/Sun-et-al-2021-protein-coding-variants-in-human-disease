#' Surface plots
#'
#' @param results data.frame computed from "allelic_heterogenetity_boost.R" 
#' @return Returned are surface plots of the observed (i.e., simulated) and expected (i.e., theoretically predicted) values for IVW-uplift.
#' @export
surface_plots = function(results){
  
  mat_obs = t(matrix(results$obsrved$median, nrow = 6, ncol = 6))
  #mat_obs = mat_obs[1:4,];
  fig_obs <- plot_ly(type = 'surface',
                     contours = list(
                       #x = list(show = TRUE, start = 1.5, end = 2, size = 0.04, color = 'white'),
                       y = list(show = TRUE, start = 1, end = 50, size = 5, color = '#2F4F4F'),
                       z = list(show = TRUE, start = 1, end = 10, size = 1, color = '#C0C0C0')),
                     x = unique(results$obsrved$maf), y = unique(results$obsrved$maf_scale), z = t(mat_obs))
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
  
  mat_exp = t(matrix(results$expected$median, nrow = 6, ncol = 6))
  #mat_exp = mat_exp[1:4,]
  fig_exp <- plot_ly(type = 'surface',
                     contours = list(
                       #x = list(show = TRUE, start = 1e-4, end = 5e-4, size = 2.5e-4, color = 'white'),
                       y = list(show = TRUE, start = 1, end = 50, size = 5, color = '#2F4F4F'),
                       z = list(show = TRUE, start = 1, end = 10, size = 1, color = '#C0C0C0')),
                     x = unique(results$obsrved$maf), y = unique(results$obsrved$maf_scale), z = t(mat_exp))
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
                      x = unique(results$obsrved$maf), y = unique(results$obsrved$maf_scale), z = t(mat_diff))
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
}