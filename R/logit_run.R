#' Simulated disease endpoints
#'
#' @param n sample size
#' @param maf minor allele frequency
#' @param pi_y disease prevalence
#' @param bta genetic effect size
#' @return Return simulated disease endpoint.
#' @export
  logit_run = function(n, maf, pi_y, bta){
    
    logit = function(pi){
      log(pi/(1-pi));
    }
    expit = function(a,b,g){
      tmp = a + b*g;
      1/(1+exp(-tmp));
    }
    full_analysis = TRUE;
    genos = rbinom(n = n, size = 2, prob = maf);
    tmp_pi_y = pi_y;
    pi_cond = expit(logit(pi_y),bta, 2*maf); 
    pi_y = expit(logit(pi_y),bta, genos); 
    y = rbinom(n = n, size = 1, prob = pi_y);
    if(full_analysis){
      dta = data.frame(Y = y, X = (genos-mean(genos)));
      #model = glm(formula = Y~X, family = binomial(link = "logit"), data = dta);
      model = speedglm::speedglm(formula = Y~X, data = dta, family = binomial(link = "logit"));
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
      model = speedglm::speedglm(formula = Y~X, data = dta, family = binomial(link = "logit"));
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

  
  
  


  