
  logit = function(pi){
    log(pi/(1-pi));
  }
  expit = function(a,b,g){
    tmp = a + b*g;
    1/(1+exp(-tmp));
  }

  logit_run = function(n, maf, pi_y, bta, full_analysis = T){
    
    genos = rbinom(n = n, size = 2, prob = maf);
    tmp_pi_y = pi_y;
    pi_cond = expit(logit(pi_y),bta, 2*maf); 
    pi_y = expit(logit(pi_y),bta, genos); 
    y = rbinom(n = n, size = 1, prob = pi_y);
    if(full_analysis){
      dta = data.frame(Y = y, X = (genos-mean(genos)));
      #model = glm(formula = Y~X, family = binomial(link = "logit"), data = dta);
      model = speedglm(formula = Y~X, data = dta, family = binomial(link = "logit"));
      #smy = summary.glm(model);
      smy = summary(model)
      coefs = smy$coefficients
      est = coefs["X","Estimate"];
      se = coefs["X","Std. Error"];
      z = coefs["X","z value"];
      
      n_ast = sum(y==1);    
      p_ast = sum((genos-mean(genos))[y==1])/2/n_ast;
      Znum = sum((genos-mean(genos))[y==1]);
      Zden = sqrt(sum(pi_y*((genos-mean(genos))^2)));
      Z_pred = Znum/Zden;
      SE_pred = (Zden)/(sum((genos-mean(genos))^2));
      return(c(theta = est, se = se, Z = z, Zpred = Z_pred, SEpred = SE_pred, mn_y = mean(y)))
      
    }else{
      dta = data.frame(Y = y, X = (genos-mean(genos)));
      #model = glm(formula = Y~X, family = binomial(link = "logit"), data = dta);
      model = speedglm(formula = Y~X, data = dta, family = binomial(link = "logit"));
      #smy = summary.glm(model);
      smy = summary(model)
      coefs = smy$coefficients
      est = coefs["X","Estimate"];
      se = coefs["X","Std. Error"];
      z = coefs["X","z value"];
      
      n_ast = sum(y==1);    
      p_ast = sum((genos-mean(genos))[y==1])/2/n_ast;
      return(c(theta = est, se = se, Z = z, N_ast = n_ast, pi_y_ast = p_ast, pi_cond = pi_cond))
    }
  }


  ratio_sig_z = function(z1, se1, z2, se2){
    return(se1*z2/(se2*z1))
  }
  
  exp_ratio = function(p1_ast, n1_ast, pi1_cond, p2_ast, n2_ast, pi2_cond){
    if(p1_ast==0 & p2_ast==0){
      return(n2_ast*pi1_cond*(1-pi1_cond)/(n1_ast*(pi2_cond*(1-pi2_cond))))
    }else{
    return(p2_ast*n2_ast*pi1_cond*(1-pi1_cond)/(p1_ast*n1_ast*(pi2_cond*(1-pi2_cond))))
    }
  }
  
  se_ratio_exp = function(p1, n1, pi1_cond, p2, n2, pi2_cond){
    return(sqrt(p2*(1-p2)*n2*pi1_cond*(1-pi1_cond)/(p1*(1-p1)*n1*pi2_cond*(1-pi2_cond))))
  }
  
  
  z_ivw = function(est1, se1, est2, se2){
              tmp1 = est1/se1^2 + est2/se2^2;
              tmp2 = 1/se1^2 + 1/se2^2;
          return(tmp1/sqrt(tmp2))
  }
  
  var_ratio = function(p1, p2){
    return(p1*(1-p1)/(p2*(1-p2)))
  }
  
  inflate_factor = function(sig_z, se_ratio){
    return((1+sig_z)/sqrt(1+(se_ratio)^2))
  }
  
  quantile_sums = function(dta, rng = c(0.1,0.5,0.9)){
    rm_obs = is.na(dta[,2]) | dta[,2]==Inf;
    dta = dta[!rm_obs,];
    maf = dta$maf;
    maf_scale=dta$maf_scale;
    obs = dta[,1];
    exp = dta[,2];
    obs_dta = data.frame(first = NA, median = NA, ninth = NA, maf = NA, maf_scale = NA)
    exp_dta = data.frame(first = NA, median = NA, ninth = NA, maf = NA, maf_scale = NA)
    count = 1;
    for(i in unique(maf)){
      for(j in unique(maf_scale)){
        obs_dta[count,] = c(quantile(obs[maf==i & maf_scale==j], rng, na.rm = T), i, j)
        exp_dta[count,] = c(quantile(exp[maf==i & maf_scale==j], rng, na.rm = T), i, j)
        count = count + 1;
      }
    }
    return(list(obsrved = obs_dta, expected = exp_dta))
  }
  
  
  