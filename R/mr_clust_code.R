#' MR-Clust mixture model fitting
#'
#' Assessment of clustered heterogeneity in Mendelian randomization analyses
#' using expectation-maximisation (EM) based model fitting of the MR-Clust
#' mixture model. Function output includes both data-tables and a visualisation
#' of the assingment of variants to clusters.
#'
#' @param  figure_4_panel Figure 4 panel to replicate - default is figure panel "a", i.e., Figure 4a.
#' @return Returned are: estimates of the putative number of clusters in the
#' sample, complete with allocation probabilities and summaries of the
#' association estimates for each variant; plots which visualise the allocation
#' of variants to clusters and; several summaries of the fitting process, i.e.
#' the BIC and likelihood estimates.
#' @export

figure_4_mr_clust = function(figure_4_panel = "a"){
  

      if(figure_4_panel == "a"){
        set.seed(123)
        theta = AllelicHeterogenitySimulator::ratio_est_a
        theta_se = AllelicHeterogenitySimulator::ratio_est_se_a
        bx = AllelicHeterogenitySimulator::bx_a
        by = AllelicHeterogenitySimulator::by_a
        bxse = AllelicHeterogenitySimulator::bxse_a
        byse =  AllelicHeterogenitySimulator::byse_a
        snp_names = AllelicHeterogenitySimulator::snp_names_a
      }else if(figure_4_panel == "b"){
        set.seed(1234)
        theta = AllelicHeterogenitySimulator::ratio_est_b
        theta_se = AllelicHeterogenitySimulator::ratio_est_se_b
        bx = AllelicHeterogenitySimulator::bx_b
        by = AllelicHeterogenitySimulator::by_b
        bxse = AllelicHeterogenitySimulator::bxse_b
        byse =  AllelicHeterogenitySimulator::byse_b
        snp_names = AllelicHeterogenitySimulator::snp_names_b
        
      }else if(figure_4_panel == "c"){
        set.seed(12345)
        theta = AllelicHeterogenitySimulator::ratio_est_c
        theta_se = AllelicHeterogenitySimulator::ratio_est_se_c
        bx = AllelicHeterogenitySimulator::bx_c
        by = AllelicHeterogenitySimulator::by_c
        bxse = AllelicHeterogenitySimulator::bxse_c
        byse =  AllelicHeterogenitySimulator::byse_c
        snp_names = AllelicHeterogenitySimulator::snp_names_c
        
      }else{
        set.seed(123456)
        theta = AllelicHeterogenitySimulator::ratio_est_d
        theta_se = AllelicHeterogenitySimulator::ratio_est_se_d
        bx = AllelicHeterogenitySimulator::bx_d
        by = AllelicHeterogenitySimulator::by_d
        bxse = AllelicHeterogenitySimulator::bxse_d
        byse =  AllelicHeterogenitySimulator::byse_d
        snp_names = AllelicHeterogenitySimulator::snp_names_d
        
      }
      

      res_em = mr_clust_em_jnk_optional(theta = theta, theta_se = theta_se, junk_mixture = FALSE,
                                    bx = bx, by = by, bxse = bxse, byse = byse, obs_names = snp_names,
                                    max_iter = 100000, min_clust_search = 10)
  
  if(figure_4_panel == "a"){
    # relabel clusters to match colours
    relab = c(3,1,4,2);
    clsts = res_em$plots$two_stage$data$cluster_class;
    clst_class = rep("Null", length(clsts));
    clst_class_num = rep(5, length(clsts));
    for(i in unique(clsts)){
      if(i!="Null"){
        clst_class[clsts==as.character(i)] = as.character(relab[as.numeric(i)]);
        clst_class_num[clsts==as.character(i)]= relab[as.numeric(i)];
      }
    }
    res_em$results$best$cluster_class = clst_class;
    res_em$results$best$cluster = clst_class_num;
    res_em$plots$two_stage = two_stage_plot(res_em$results$best, bx, by, bxse, byse, snp_names)
    res_em$plots$two_stage$data = res_em$plots$two_stage$data %>% 
    filter(probability >=0.7) %>%
      # select(-cluster_class) %>%
      # left_join(res_em1_junk$plots$two_stage$data %>% select(observation, cluster_class)) %>% 
      dplyr::mutate(text_labels = ifelse(observation == "HCN4", "HCN4", 
                                  ifelse(observation == "SCN5A", "SCN5A",
                                         ifelse(observation == "SCN10A", "SCN10A",""))))
  }else if(figure_4_panel == "b"){
        res_em$plots$two_stage = two_stage_plot(res_em$results$best, bx, by, bxse, byse, snp_names)
        res_em$plots$two_stage$data = res_em$plots$two_stage$data %>%
          filter(probability >=0.7) %>%
          dplyr::mutate(text_labels = ifelse(observation == "HCN4", paste0(""),  
                                       ifelse(observation == "HCN4_gwa", paste0("HCN4"), 
                                              ifelse(observation == "SCN5A", paste0(""),
                                                     ifelse(observation == "SCN10A", paste0("SCN10A,SCN5A"),"")))))
        nudges_x = nudges_y = rep(0,19);
        nudges_x[c(3,4,8,18)] = c(-0.125,0.025,0.1,-0.025);
        nudges_y[c(3,4,8,18)] = c(0.04,0.095,0.025,-0.05);
  }else if (figure_4_panel == "c"){
    
    clsts = res_em$plots$two_stage$data$cluster_class;
    ord = order(-res_em$plots$two_stage$data$cluster_mean)
    dta = res_em$plots$two_stage$data[ord,]
    clsts = unique(clsts[ord])
    clst_class = rep("", dim(dta)[1]);
    clst_class_num = rep(0, dim(dta)[1]);
    count = 1;
      for(i in unique(clsts)){
        clst_class[dta$cluster_class==i] = as.character(count);
        clst_class_num[dta$cluster==as.numeric(i)]= count;
        count = count+1;
      }
    
    dta$cluster_class = clst_class;
    dta$cluster = clst_class_num;
    
    res_em$results$best = dta;
    res_em$plots$two_stage = two_stage_plot(res_em$results$best, bx[ord], by[ord], bxse[ord], byse[ord], snp_names[ord])
    res_em$plots$two_stage$data = res_em$plots$two_stage$data %>% filter(probability >=0.7)
    text_labels = NULL;
  }else{
    res_em$plots$two_stage$data = res_em$plots$two_stage$data %>% filter(probability >=0.7)
    text_labels = NULL;
  }


res_em$plots$two_stage$layers[[1]]$aes_params$shape = 20
res_em$plots$two_stage$layers[[1]]$mapping$size = NULL
res_em$plots$two_stage$layers[[1]]$aes_params$alpha = 0.7
res_em$plots$two_stage$layers[[1]]$aes_params$size = 0.05

plot_out = res_em$plots$two_stage + 
  # ggplot(aes(color = ~cluster_class,  size = ~probability)) + 
  ggplot2::xlim(c(0, 1)) + 
  ggplot2::ylim(c(-0.6, 0.15)) + 
  ggplot2::xlab("Genetic association with AF") + 
  ggplot2::ylab("Genetic association with pulse") + ggplot2::ggtitle("") + 
  ggplot2::theme(axis.text=element_text(size=12),axis.title=element_text(size=16),
                 legend.text=element_text(size=12), legend.title=element_text(size=16),legend.key.size = unit(1.5, 'lines'))
  if(figure_4_panel == "a"){plot_out = plot_out +ggplot2::geom_point(aes(text = observation, color = cluster_class))+
                            ggrepel::geom_text_repel(mapping = aes(label = text_labels), size = 5,nudge_y = 0.025, nudge_x = -0.125, min.segment.length = 1)}
  if(figure_4_panel == "b"){plot_out = plot_out + ggplot2::geom_point(shape=c(rep(10,9),rep(20,10)), size = 2.25, aes(text = observation, color = cluster_class)) +
                            ggrepel::geom_text_repel(mapping = aes(label = text_labels), size = 5,nudge_y = 0.025, nudge_x = -0.125, min.segment.length = 1)}
  if(figure_4_panel == "c"){plot_out = plot_out +ggplot2::geom_point(aes(text = observation, color = cluster_class))+
                            ggplot2::ylim(c(-0.1, 0.15)) + ggplot2::xlim(c(0, 0.5))}
  if(figure_4_panel == "d"){plot_out = plot_out +ggplot2::geom_point(aes(text = observation, color = cluster_class))}

return(plot_out)
}