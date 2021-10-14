library(readxl)
n1 = 392814;
n2 = 260405;

# Disease prevalance
ukb_fg_dta = read_xlsx("~/Desktop/sectors/pharma/biogen/projects/FinnGen CWAS paper/Supplementary_Tables.xlsx", col_names = T, sheet = 1)
ukb_fg_dta = ukb_fg_dta[c(-1),];
colnames(ukb_fg_dta) = ukb_fg_dta[1,];
ukb_fg_dta = ukb_fg_dta[c(-1),];
head(ukb_fg_dta)


ukb_cases = as.numeric(ukb_fg_dta$`UKB cases`)
summary(ukb_cases)
prev_ukb = quantile(ukb_cases, c(0.1,0.5,0.9))/n1;
prev_ukb
  fg_cases = as.numeric(ukb_fg_dta$`FG cases (R6)`)
  summary(fg_cases)
  prev_fg = quantile(fg_cases, c(0.1,0.5,0.9), na.rm = T)/n2;
  prev_fg
    fg_cases_2 = as.numeric(ukb_fg_dta$`FG cases (R5)`)
    summary(fg_cases_2)
    prev_fg_2 = quantile(fg_cases_2, c(0.1,0.5,0.9), na.rm = T)/n2;
    prev_fg_2
    
prev_ratio =  fg_cases/ukb_cases;
summary(prev_ratio)

alpha_ukb = logit(prev_ukb);
alpha_fg = logit(prev_fg);
alpha_ukb;
alpha_fg;

# MAFs
ukb_fg_maf = read_xlsx("~/Desktop/sectors/pharma/biogen/projects/FinnGen CWAS paper/Supplementary_Tables.xlsx", col_names = T, sheet = 3)
ukb_fg_maf = ukb_fg_maf[c(-1,-2),];
colnames(ukb_fg_maf) = ukb_fg_maf[1,];
ukb_fg_maf = ukb_fg_maf[c(-1),];
head(ukb_fg_maf)

ukb_maf = as.numeric(ukb_fg_maf$`A1 freq (UKB)`)
summary(ukb_maf);
ukb_maf[ukb_maf>0.5] = 1- ukb_maf[ukb_maf>0.5]
maf_ukb = quantile(ukb_maf[ukb_maf<0.05], c(0.1,0.5,0.9));
maf_ukb
  fg_maf = as.numeric(ukb_fg_maf$`A1 freq (FG)`)
  summary(fg_maf);
  fg_maf[fg_maf>0.5] = 1- fg_maf[fg_maf>0.5]
  maf_fg = quantile(fg_maf[fg_maf<0.05], c(0.1,0.5,0.9));
  maf_fg
  
  
  
  
  # Beta estimates
  ukb_fg_betas = read_xlsx("~/Desktop/sectors/pharma/biogen/projects/FinnGen CWAS paper/Supplementary_Tables.xlsx", col_names = T, sheet = 2)
  colnames(ukb_fg_betas) = ukb_fg_betas[3,];
  ukb_fg_betas = ukb_fg_betas[-c(1:4),];
  head(ukb_fg_betas)
  tmp_ukb = as.numeric(as.data.frame(ukb_fg_betas[,"MAF (UKB)"])$`MAF (UKB)`);
  tmp_fg = as.numeric(as.data.frame(ukb_fg_betas[,"MAF_FG"])$`MAF_FG`);
  maf_uplift = tmp_ukb<= tmp_fg;
  bta_ukb = as.data.frame(ukb_fg_betas[maf_uplift,"BETA_UKB"]);
  bta_fg = as.data.frame(ukb_fg_betas[maf_uplift,"BETA_FG"]);
  betas = data.frame(ukb = as.numeric(tolower(bta_ukb$BETA_UKB)), fg = as.numeric(tolower(bta_fg$BETA_FG)));
  #mafs_ukb_fg = 
  