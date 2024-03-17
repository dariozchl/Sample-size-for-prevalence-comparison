
##########################################################################################
# Simulation to determine the number of persons examined by the gold standard to achieve a certain power to reject the hypothesis of two prevalences being the same
# Author: Dario Zocholl
# latest changes: March 15, 2024
# accompanying software code for Manteltext "Bayesian approaches to medical studies under limited information"

library(tidyverse)
library(rstan)
library(foreach)
library(doParallel)

prev_model_hier <- stan_model("compare_prevalence_hierarchical.stan")

nsim = 1e3

prop_gold = 0.01
true_pos = c(0.05, 0.02)
true_Se = c(0.9,0.9)
true_Sp = c(0.95,0.95)
counter <- 0
sign <- FALSE
power <- c()

while(sum(sign)/nsim < 0.99){
  
  counter <- counter+1
  prop_gold <- seq(0,0.3,by=0.01)[counter]
  print(prop_gold)
  
  
  ncores <- detectCores()-1; 
  cl <- makeCluster(ncores); registerDoParallel(cl);
  start.time <- Sys.time()
  
  sign <- foreach(i=1:nsim, .combine=rbind, .packages=c("rstan")) %dopar% {
    
    n <- c(10000, 10000)
    n_Se <- n_Sp <- y_Se <- y_Sp <- n_gold <- c()
    
    true_data_1 <- rbinom(n=n[1],size=1,prob=true_pos[1])
    test_data_1 <- rbinom(n=n[1],size=1,prob=c(1-true_Sp[1], true_Se[1])[true_data_1+1]) # a positive test has probability 1-Sp if true status = 0 and probability Se if true status = 1
    person_gold_1 <- sample(1:n[1], size=round(n[1]*prop_gold)) # sample participants for gold standard 
    true_data_gold_1 <- true_data_1[person_gold_1]
    test_data_gold_1 <- test_data_1[person_gold_1]
    n_Se[1] <- sum(true_data_1[person_gold_1]==1)
    y_Se[1] <- sum(test_data_gold_1[which(true_data_1[person_gold_1]==1)] == 1)
    n_Sp[1] <- sum(true_data_1[person_gold_1]==0)
    y_Sp[1] <- sum(test_data_gold_1[which(true_data_1[person_gold_1]==0)] == 0)
    
    true_data_2 <- rbinom(n=n[2],size=1,prob=true_pos[2])
    test_data_2 <- rbinom(n=n[2],size=1,prob=c(1-true_Sp[2], true_Se[2])[true_data_2+1]) # a positive test has probability 1-Sp if true status = 0 and probability Se if true status = 1
    person_gold_2 <- sample(1:n[2], size=round(n[2]*prop_gold)) # sample participants for gold standard 
    true_data_gold_2 <- true_data_2[person_gold_2]
    test_data_gold_2 <- test_data_2[person_gold_2]
    n_Se[2] <- sum(true_data_2[person_gold_2]==1)
    y_Se[2] <- sum(test_data_gold_2[which(true_data_2[person_gold_2]==1)] == 1)
    n_Sp[2] <- sum(true_data_2[person_gold_2]==0)
    y_Sp[2] <- sum(test_data_gold_2[which(true_data_2[person_gold_2]==0)] == 0)
    
    fit_raw <- sampling(prev_model_hier, 
                        data = list(mu=c(1.81,1.84), Sigma=matrix(c(2, 0, 0, 2), nrow=2), 
                                    prev_beta=c(1,1), n=n, y=c(sum(test_data_1), sum(test_data_2)), n_Se=n_Se, y_Se=y_Se, n_Sp=n_Sp, y_Sp=y_Sp),
                        warmup = 1000, iter = 6000, chains = 4, cores = 1, refresh=1000, control=list(adapt_delta=0.95, max_treedepth=15))
    
    sign <- (sum(extract(fit_raw)$diff > 0) / 2e4) > 0.95
    return(sign)
  }
  stopCluster(cl) # end parallel computing
  print(Sys.time() - start.time)
  print(sum(sign)/nsim)
  power[counter] <- sum(sign)/nsim
}

data.frame(power, "prop.gold"=seq(0,0.3,by=0.01)[1:length(power)])
data.frame(power, "prop.gold"=seq(0,0.3,by=0.01)[1:length(power)]) %>% 
  mutate(lower=qbeta(0.025, power*nsim, nsim-power*nsim+1),
         upper=qbeta(0.975, power*nsim+1, nsim-power*nsim)) %>% 
  ggplot(.) + geom_line(aes(x=prop.gold, y=power), linewidth=0.2) + geom_ribbon(aes(ymin=lower, ymax=upper, x=prop.gold), alpha=0.2) + 
  xlab("Proportion of clinical interviews") + ylab("Power") +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) + scale_x_continuous(limits=c(0,0.27), breaks=seq(0,0.25,by=0.05), expand=c(0,0)) + 
  theme_bw()

ggsave("power_prevalence_test.eps", device=cairo_ps, width=5, height=2)




