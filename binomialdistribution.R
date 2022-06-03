# Binomial Distribution for Summary Allele data
# gnomAD v2.1

# Libraries
library(tidyverse)
library(plotly)
library(nloptr)
library(Summix)

# Load merged genomic data
gdat = read.csv("~/gnomad_genome_ANmerge.txt.gz", sep="")

# Filter test frame, sample 100k snps
samplesnps = 100000
test_frame = gdat %>% 
  select(ref_AF_eur_1000G, ref_AF_afr_1000G, gnomad_AF_afr, gnomad_AFR_AN) %>%
  sample_n(samplesnps) %>% 
  mutate(obs_afr = floor(gnomad_AF_afr * gnomad_AFR_AN))

# initialize Log-likelihood test frame
ll_sguess = data.frame(AFR_proportion = seq(0.01, .99, 0.01)) %>% 
  mutate(EUR_proportion = 1 - AFR_proportion,
         LL_value = numeric(length(AFR_proportion)))

# Evaluate log-likelihood over grid search
for (i in 1:length(ll_sguess$AFR_proportion)){
  ll_val = dbinom(test_frame$obs_afr, test_frame$gnomad_AFR_AN, ((ll_sguess$AFR_proportion[i]*test_frame$ref_AF_afr_1000G) + (ll_sguess$EUR_proportion[i]*test_frame$ref_AF_eur_1000G)))
  ll_val[ll_val == 0] = min(ll_val[ll_val != 0])
  ll_sguess$LL_value[i] = sum(log(ll_val))
}

# Maximum log likelihood value
max_LL = ll_sguess[which.max(ll_sguess$LL_value),]

# Create AFR potting line
AFR_line = data.frame(x = rep(max_LL[[1]], 2),
                      y = c(0, 1),
                      z = c(0,0))
# EUR plotting line
EUR_line = data.frame(x = c(0, 1),
                      y = rep(max_LL[[2]], 2),
                      z = c(0,0))
# Vertical/Solution plotting line
vert_line = data.frame(x = rep(max_LL[[1]], 2),
                       y = rep(max_LL[[2]], 2),
                       z = c(-min(ll_sguess$LL_value), max(ll_sguess$LL_value)))

# 3D plot of log-likelihood and the maximum value
plot_ly(ll_sguess, x = ~AFR_proportion, y = ~EUR_proportion, z = ~-LL_value, type = "scatter3d", mode = "lines") %>% 
  add_trace(AFR_line, x = AFR_line$x, y = AFR_line$y, z = AFR_line$z) %>% 
  add_trace(EUR_line, x = EUR_line$x, y = EUR_line$y, z = EUR_line$z) %>% 
  add_trace(vert_line, x = vert_line$x, y = vert_line$y, z = vert_line$z)

# Initiliaze Summix (Least squares) ancestry estimation method
sum_sguess = data.frame(AFR_proportion = seq(0.001, .999, 0.001)) %>% 
  mutate(EUR_proportion = 1 - AFR_proportion,
         sum_value = numeric(length(AFR_proportion)))
# Evaluate summix over grid search
for (i in 1:length(sum_sguess$AFR_proportion)){
  sum_sguess$sum_value[i] = sum((test_frame$gnomad_AF_afr - (sum_sguess$AFR_proportion[i]*test_frame$ref_AF_afr_1000G) -  (sum_sguess$EUR_proportion[i]*test_frame$ref_AF_eur_1000G))^2)
}
# Minimum Summix value
min_sum = sum_sguess[which.min(sum_sguess$sum_value),]
# Summix plottng lines
AFR_sumline = data.frame(x = rep(min_sum[[1]], 2),
                         y = c(0, 1),
                         z = rep(min(sum_sguess, 2)))
EUR_sumline = data.frame(x = c(0, 1),
                         y = rep(min_sum[[2]], 2),
                         z = rep(min(sum_sguess, 2)))
vert_sumline = data.frame(x = rep(min_sum[[1]], 2),
                          y = rep(min_sum[[2]], 2),
                          z = c(min(sum_sguess), max(sum_sguess)))
# Plot Summix least squares minimization algorithm
plot_ly(sum_sguess, x = ~AFR_proportion, y = ~EUR_proportion, z = ~sum_value, type = "scatter3d", mode = "lines") %>% 
  add_trace(AFR_sumline, x = AFR_sumline$x, y = AFR_sumline$y, z = AFR_sumline$z) %>% 
  add_trace(EUR_sumline, x = EUR_sumline$x, y = EUR_sumline$y, z = EUR_sumline$z) %>% 
  add_trace(vert_sumline, x = vert_sumline$x, y = vert_sumline$y, z = vert_sumline$z)

# Verify estimates using grid descent (SQP) method on inverted log likelihood eq.
# Starting guesses
starting = c(0.5, 0.5)
# Minimization EQ
fn.ancmix.ll = function(x){
  mf = dbinom(test_frame$obs_afr, test_frame$gnomad_AFR_AN, ((x[1]*test_frame$ref_AF_afr_1000G) + (x[2]*test_frame$ref_AF_eur_1000G)))
  mf[mf == 0] = min(mf[mf != 0])
  mf = -1 * sum(log(mf))
}
# Equalit EQ
heq.ancmix.ll = function(x){
  equality = x[1] + x[2] - 1
}
# SQP Minimization method
S = slsqp(starting,
          fn.ancmix.ll,
          heq = heq.ancmix.ll)
