# Empirical Analysis III - Pset 2 
# Script created by: Yixin Sun, Nadia Lucas, Camilla Schneier
# With minor edits from Franco Calle.

packages = c("tidyverse", "haven", "knitr", 
             "MatchIt", "sandwich", "stargazer", 
             "gridExtra")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

root <- "C:/Users/franc/Dropbox/Franco First Year/Empirical Analysis III/problem_sets/pset2/question2"

# Reading in data and setting up control variables
lalondeDataFrame <- 
  read_dta(file.path(root, "lalonde2.dta")) %>%
  mutate(treated = if_else(treated == 1, 1, 0), 
         sample = case_when(sample == 1 ~ "NSW", 
                            sample == 2 ~ "CPS", 
                            sample == 3 ~ "PSID")) 

nsw_df <- filter(lalondeDataFrame, !is.na(treated))

controls_indiv <- 
  c("age", "educ", "black", "married", "nodegree", "re74", #"re75", 
    "hisp", "kids18", "kidmiss") 

controls <- paste(controls_indiv, collapse = " + ")

f_full <- 
  paste("treated", controls, sep = " + ") %>%
  paste("re78", ., sep = " ~ ") %>%
  formula


# ============================================================
# Functions Used Throughout
# ============================================================ 
# function for performing t-test over all control variables
bal_tab <- function(cov, treat, data){
  tt_formula <- paste0(cov, " ~ ", treat)
  
  tt <- t.test(formula = formula(tt_formula), data = data)
  ests <- tt$estimate
  names(ests) <- sub("mean in group ()", "\\1",names(ests))
  counts <- xtabs(tt_formula[c(1,3)],data)
  names(counts) <- paste0("n",names(counts))
  cbind(          	
    coef = cov,
    as.list(ests),
    data.frame(pvalue = tt$p.value, stringsAsFactors = FALSE)) %>%
    format(digits = 3)
}

# Function for creating desired latex tables
print_output <- function(m, outpath){
  if(class(m) != "list"){m <- list(m)}
  
  cov <- map(m, vcovHC, type = "HC1")
  robust_se <- map(cov, function(x) sqrt(diag(x)))
  stargazer(m, keep.stat = "n", keep = "treated", se = robust_se,
            type = "latex", dep.var.labels = c("No Controls", "Controls"), 
            star.cutoffs = NA, omit.table.layout = "n", style = "qje", 
            out = outpath, digits = 2)
}

# function for plotting the supports of continuous variables
plot_support <- function(var, name, d){
  d %>%
    mutate(treated = as.factor(treated)) %>%
    ggplot(aes(x = get(var), y = ..density.., fill = treated, color = treated)) + 
    geom_histogram(position = "identity", alpha = .5, bins = 20) + 
    theme_minimal() + scale_fill_manual(values = c("#69b3a2", "#404080")) +
    scale_color_manual(values = c("#69b3a2", "#404080"))  + 
    theme(legend.justification=c(1, 1), legend.position=c(1, 1), 
          legend.key.width=unit(0.3,"cm"),legend.key.height=unit(0.3,"cm")) + 
    xlab(var) +
    ggtitle(name)
}

# Function for producing mirror density plots
mirror_density <- function(data){
  data_treated <- filter(data, treated == 1)
  data_control <- filter(data, treated == 0)
  
  ggplot(data = data_treated) +
    geom_histogram(aes(x = propensity_score, y = ..density.., fill = "Control"), 
                   fill = "#404080", color = "#404080", alpha = .5) + 
    theme_minimal() + 
    geom_histogram(data = data_control, 
                   aes(x = propensity_score, y = -..density.., fill = "Treated"), 
                   fill = "#69b3a2", color = "#69b3a2", alpha = .5) + 
    xlab("Propensity Score") 
}

# ============================================================
# Part A: Is the data is consistent with randomization of the treatment?
# ============================================================
# calculate F-statistic
lm(formula(paste("treated", controls, sep = "~")), data = nsw_df) %>%
  summary

# perform t.test for each control variable
# bal_tab() is a user written function - see beginning of code
balance <- 
  map_df(controls_indiv, bal_tab, 'treated', nsw_df) %>%
  kable(format = "latex", 
        col.names = c("Variable", "Control", "Treatment", "P-Value"), 
        escape = F, 
        digits = 2, 
        booktabs = TRUE)

sum(nsw_df$treated)
sum(1-nsw_df$treated)

# output the balance table to latex
writeLines(balance, file.path(root, "3a_balance.tex"))

# ============================================================
# Part B - Estimate effect using NSW sample
# ============================================================
nsw_lm <- lm(formula(re78 ~ treated), data = nsw_df)
nsw_full_lm <- lm(f_full, data = nsw_df)

# print_output() is a user written function for printing pretty regression tables 
  # see beginning of code file
print_output(list(nsw_lm, nsw_full_lm), file.path(root, "3b_nsw_reg.tex"))

# ============================================================
# Part C - Estimate effect using NSW sample + CPS
# ============================================================
# create new dataset of NSW and CPS variables
# created new treated variable that is 1 if it is NSW treatment, and
  # 0 for CPS

nsw_cps_df <- 
  lalondeDataFrame %>%
  filter((sample == "NSW" & treated ==1 | sample == "CPS")) %>%
  mutate(treated = if_else(is.na(treated), 0, 1))


nsw_cps_lm <- lm(formula(re78 ~ treated), data = nsw_cps_df)
nsw_cps_full_lm <- lm(f_full, data = nsw_cps_df)

list(nsw_cps_lm, nsw_cps_full_lm) %>%
  print_output(file.path(root, "3c_nsw_cps_reg.tex"))

# ============================================================
# Part D -  Investigate covariate balancing and support 
# between the treated and the CPS sample.
# ============================================================ 
# Run t test on variables in CPS and NSW dataset
balance_cps <- 
  map_df(controls_indiv, bal_tab, 'treated', nsw_cps_df) %>%
  kable(format = "latex", 
        col.names = c("Variable", "Control", "Treatment", "P-Value"), 
        escape = F, 
        digits = 2, 
        booktabs = TRUE)
writeLines(balance_cps, file.path(root, "3d_balance_cps.tex"))

sum(nsw_cps_df$treated)
sum(1-nsw_cps_df$treated)

# applying the plotting function over the continuous variables and arranging 
# them in a grid for better display
# plot_support() is a user written function - see beginning of code
balance_plots <- 
  list(c("age", "Age (Years)"), 
       c("educ", "Years of Schooling"), 
       c("kids18" ,"Kids Younger Than 18"), 
       c("re74", "Real Earnings, 1974")) %>%
  map(function(x) plot_support(x[1], x[2], nsw_cps_df))

balance_plots_all <- grid.arrange(balance_plots[[1]], balance_plots[[2]], 
                                  balance_plots[[3]], balance_plots[[4]], ncol =2)
ggsave(balance_plots_all, file = file.path(root, "3d_balance_plots.png"), width = 13, height = 10)

# ============================================================
# Part E -  Estimate the effect using 1 nearest neighbor propensity score matching.
# ============================================================ 
# first calculate propensity score by running glm
match_f <- formula(paste("treated", controls, sep = "~"))

nsw_cps_df <- 
  nsw_cps_df %>%
  dplyr::select(treated, all_of(controls_indiv), re78) %>%
  na.omit() %>%
  mutate(treated = as.factor(treated))

# check support of propensity score via graph -------------------------------
sample_ps <- glm(match_f, family = binomial, data = nsw_cps_df)
nsw_cps_df <- cbind(nsw_cps_df, propensity_score = sample_ps$fitted.values)

# mirror_density() is a user written function for plotting the support
# see beginning of code
mirror_density(nsw_cps_df)
ggsave(file = file.path(root, "3e_propensity_score.png"), 
       width = 10, height = 6)

# trimming the support to not include <.02 ----------------------------------
nsw_cps_trimmed_df <- filter(nsw_cps_df, propensity_score > .02)
mirror_density(nsw_cps_trimmed_df)
ggsave(file = file.path(root, "3e_propensity_score_trimmed.png"), 
       width = 10, height = 6)

# find nearest neighbor match ---------------------------------------------- 
matches <- matchit(match_f, method = "nearest", data = nsw_cps_trimmed_df)
matches_df <- match.data(matches)

nn_lm <- lm(formula(re78 ~ treated), data = matches_df)
nn_full_lm <- lm(f_full, data = matches_df)

list(nn_lm, nn_full_lm) %>%
  print_output(file.path(root, "3e_nn_reg.tex"))

# matching with replacement
matches_replace <- matchit(match_f, method = "nearest", data = nsw_cps_trimmed_df, replace = TRUE)
matches_replace_df <- match.data(matches_replace)

nn_replace_lm <- lm(formula(re78 ~ treated), data = matches_replace_df)
nn_replace_full_lm <- lm(f_full, data = matches_replace_df)

list(nn_replace_lm, nn_replace_full_lm) %>%
  print_output(file.path(root, "3e_nn_replace_reg.tex"))


# assessing match quality -------------------------------------------------- 
# estimate bias reduction from matching
bias <- 
  summary(matches_replace)$reduction %>%
  dplyr::select(bias = "Mean Diff.") %>%
  mutate(coef = rownames(.))

# compute t-test for covariates and print balance + bias table 
balance_matches <- 
  map_df(controls_indiv, bal_tab, 'treated', matches_replace_df) %>%
  left_join(bias) %>%
  kable(format = "latex", 
        col.names = c("Variable", "Control", "Treatment", "P-Value", "% Reduction in Bias"), 
        escape = F, 
        digits = 2, 
        booktabs = TRUE)
writeLines(balance_matches, file.path(root, "3e_balance_matches.tex"))


# ============================================================
# Part F -  Estimate the effct using the propensity score and local linear regression.
# ============================================================ 
# function that runs the local linear regression on treated and untreated sample
loclin <- function(data, t){
  data <- filter(data, treated == t)
  lp <- loess(re78 ~ propensity_score, data = data, degree = 1)
  data$ll_estimate <- predict(lp)
  return(data)
}

# apply local linear function to the unmatched dataset
unmatched_loclin <- 
  bind_rows(loclin(nsw_cps_df, 1), loclin(nsw_cps_df, 0)) %>%
  t.test(ll_estimate ~ treated, .)

# apply local linear function to the matched dataset
matched_loclin <- 
  bind_rows(loclin(matches_replace_df, 1), loclin(matches_replace_df, 0)) %>%
  t.test(ll_estimate ~ treated, .)

# print results of ATT for matched dataset
loclin_output <-
  matched_loclin$estimate %>%
  t() %>%
  as_tibble %>%
  dplyr::select(mean_control = 1, mean_treat = 2) %>%
  mutate(ATT = mean_treat - mean_control, 
         SE = matched_loclin$stderr) %>%
  kable(format = "latex", 
        col.names = c("Mean Control", "Mean Treatment", "ATT", "Std. Error"), 
        escape = F, 
        digits = 2, 
        booktabs = TRUE)

writeLines(loclin_output, file.path(root, "3f_loclin.tex"))
