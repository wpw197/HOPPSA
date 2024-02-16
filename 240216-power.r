## ------------------------------------------------------------------------
## 
## Script name: FSH in HOPPSA Power Calculation
##
## Purpose of script: Simulate sample size and power calculations for the 
## proposed extension of the HOPPSA study
##
## Authors: Emily Lane, Adam Brentnall
##
## Last Update: 15/02/2024
##
## ------------------------------------------------------------------------
##
## Notes: The following script performs the sample size / power calculation 
## to detect a difference in FSH levels between participants undergoing 
## hysterectomy and opportunistic salpingectomy and participants undergoing
## hysterectomy.
##
## Assumptions:
## 1. The age of menopause is normally distributed with mean age 47 and sd 4.4
## amongst those who have undergone hysterectomy with opportunistic 
## salpingectomy and normally distributed with mean age 49 and sd 4.4 in those
## who have undergone hysterectomy. 
## 2. log(FSH) given age and menopausal status is normally distributed 
## 3. Models for the distribution of FSH given age and menopausal status
## were built using data collected in the Study of Women's Health Across 
## the Nation (SWAN).
## 4. An increase in the mean log(FSH) of 0.1 has been applied to post-menopausal
## women to adjust for differences in BMI between the SWAN and HOPPSA cohorts.
## This was chosen because ... XXXXXXX
## ------------------------------------------------------------------------
## load up the packages we will need: 
packages <- c("tidyverse", "parallel")
##if needed install
##install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE, quietly = TRUE)

# Models ------------------------------------------------------------------
## Pre-menopausal mean and sd for the distribution of FSH given age FSH is 
## measured:
## Model1: log(FSH) ~ alpha1 + beta1*Age
## Model1 residuals: r1
## Model2: r1^2 ~ alpha2 + beta2*Age
## mean1 = mean(FSH) (by age)
## Model3: log(mean1) ~ alpha3 + beta3*Age
## sd(log(fsh)) = sqrt(alpha2 +beta2*Age)
## Post-menopausal mean and sd for the distribution of FSH given age FSH is 
## measured:
## Model4: log(FSH) ~ alpha4 + beta4*Age
## Model4 residuals: r2
## Model5: r2^2 ~ alpha5 + beta5*Age
## mean2 = mean(FSH) (by age)
## Model6: log(mean2) ~ alpha6 + beta6*Age
## sd(log(fsh)) = sqrt(alpha5 +beta5*Age)

# Functions  --------------------------------------------------------------

## One simulation replicate

fn.sim<-function(idx, params){
    ## Input parameters
    ## idx simulation replicate (not used)
    ## params = list of parameters (see below)
    
    ## output df
    sims <- data.frame(
        test = 1:params$n_tests,
        mean_fsh_surgery = rep(NA,params$n_tests),
        mean_fsh_no_surgery = rep(NA,params$n_tests),
        p_val= rep(NA,params$n_tests), 
        sd_fsh_surgery = rep(NA,params$n_tests),
        sd_fsh_no_surgery = rep(NA, params$n_tests), 
        mean_age_surgery = rep(NA,params$n_tests),
        mean_age_no_surgery = rep(NA,params$n_tests), 
        surgery_p_postmeno = rep(NA, params$n_tests),
        no_surgery_p_postmeno = rep(NA,params$n_tests)
    )

    ## this will hold the age and follow up, probability of post-menopausal status, post-menopausal status, and FSH in those who have 
    ## had hysterectomy and opportunistic salpingectomy 
    surgery_p <- cbind(sample_size = rep(params$n_samp, params$n_samp),age = sample(params$ages, size = params$n_samp, replace = TRUE), 
                       follow_up = sample(params$follow_up, size = params$n_samp, replace = TRUE), p1 = NA, outcome1 = NA, fsh1 = NA)
 
    ## same for those undergoing just hysterectomy 
    no_surgery_p <-  surgery_p

    ## first test
    j=1

    ##loop over number people
    for(l in 1:params$n_samp){
        
    # calculates the probability that the patient is post menopausal at their first fsh test given they were pre-menopausal at study entry
        surgery_p[l,((3*j)+1)] <- (pnorm((surgery_p[l,2]+surgery_p[l,3]),mean=params$mean_age_surgery, sd=params$sd_age) - pnorm(surgery_p[l,2],mean=params$mean_age_surgery, sd=params$sd_age))/ 
      (1-pnorm(surgery_p[l,2], mean=params$mean_age_surgery, sd=params$sd_age))
        
        no_surgery_p[l,((3*j)+1)] <- (pnorm((no_surgery_p[l,2]+no_surgery_p[l,3]),mean=params$mean_age_no_surgery, sd=params$sd_age) - pnorm(no_surgery_p[l,2],mean=params$mean_age_no_surgery, sd=params$sd_age))/ 
            (1-pnorm(no_surgery_p[l,2],mean=params$mean_age_no_surgery, sd=params$sd_age))
        
        ## uses the probability of being post-menopausal to simulate whether or not the participant is post-menopausal at their first test
        surgery_p[l, ((3*j)+2)]<-rbinom(1,1,surgery_p[l,((3*j)+1)])
        no_surgery_p[l,((3*j)+2)]<-rbinom(1,1,no_surgery_p[l,((3*j)+1)])
    
        ## given the menopausal status and age of the participants a value for FSH is simulated 
        surgery_p[l,((3*j)+3)] <- ifelse(surgery_p[l, ((3*j)+2)]==1, rnorm(1, (0.1+params$alpha6 + params$beta6*(surgery_p[l,2]+surgery_p[l,3])),
                                                                           sqrt(params$alpha5+params$beta5*(surgery_p[l,2]+surgery_p[l,3]))),
                                         rnorm(1, params$alpha3 + params$beta3*(surgery_p[l,2]+surgery_p[l,3]),
                                               sqrt(params$alpha2+params$beta2*(surgery_p[l,2]+surgery_p[l,3]))))
        
        no_surgery_p[l,((3*j)+3)] <- ifelse(no_surgery_p[l, ((3*j)+2)]==1, rnorm(1, (0.1+params$alpha6 + params$beta6*(no_surgery_p[l,2]+no_surgery_p[l,3])),
                                                                                 sqrt(params$alpha5+params$beta5*(no_surgery_p[l,2]+no_surgery_p[l,3]))),
                                            rnorm(1, params$alpha3 + params$beta3*(no_surgery_p[l,2]+no_surgery_p[l,3]),
                                                  sqrt(params$alpha2+params$beta2*(no_surgery_p[l,2]+no_surgery_p[l,3]))))
    
    }
    
    ## adding results to the sims data frame 
        sims[sims$test == j, 2] <- mean(surgery_p[,((3*j)+3)])
        sims[sims$test == j, 3] <- mean(no_surgery_p[,((3*j)+3)])
        test <- wilcox.test(surgery_p[,((3*j)+3)], no_surgery_p[,((3*j)+3)])
        sims[sims$test == j, 4]  <- test$p.value
        sims[sims$test == j, 5]  <- sd(surgery_p[,((3*j)+3)])
        sims[sims$test == j, 6]  <- sd(no_surgery_p[,((3*j)+3)])
        sims[sims$test == j, 7] <- mean(surgery_p[,2])
        sims[sims$test == j, 8] <-  mean(no_surgery_p[,2])
        sims[sims$test == j, 9] <- (sum(surgery_p[, ((3*j)+2)])/length(surgery_p[,1]))*100
        sims[sims$test == j, 10] <- (sum(no_surgery_p[, ((3*j)+2)])/length(no_surgery_p[,1]))*100
        ## saving the first set of results
        no_surgery_p1 <- no_surgery_p
        surgery_p1 <- surgery_p

    
    ## further tests (assume 2+)
    for(j in 2:params$n_tests){
    
        ## adding columns for the next probabilities and tests
        surgery_p <- data.frame(surgery_p, p = NA, outcome = NA, fsh = NA)
        no_surgery_p <- data.frame(no_surgery_p, p = NA, outcome = NA, fsh = NA)
        colnames(surgery_p)[((3*j)+1):((3*j)+3)] <- c(paste0('p', j), paste0('outcome',j), paste0('fsh',j))
        colnames(no_surgery_p)[((3*j)+1):((3*j)+3)] <- c(paste0('p', j), paste0('outcome',j), paste0('fsh',j))
        
        for(l in 1:dim(surgery_p)[1]){
            
            ## probability of post-menopausal given status at last test
            ## if they were post-menopausal at the last test their probability is 1
            surgery_p[l,((3*j)+1)] <- ifelse(surgery_p[l,((3*(j-1))+2)]==1, 1,
            (pnorm((surgery_p[l,2]+surgery_p[l,3]+(j-1)),mean=params$mean_age_surgery, sd=params$sd_age) - pnorm(surgery_p[l,2]+surgery_p[l,3]+(j-2),mean=params$mean_age_surgery, sd=params$sd_age))/ 
            (1-pnorm(surgery_p[l,2]+surgery_p[l,3]+(j-2),mean=params$mean_age_surgery, sd=params$sd_age)))
            
            no_surgery_p[l,((3*j)+1)] <- ifelse(no_surgery_p[l,((3*(j-1))+2)]==1, 1,
            (pnorm((no_surgery_p[l,2]+no_surgery_p[l,3]+(j-1)),mean=params$mean_age_no_surgery, sd=params$sd_age)-
             pnorm(no_surgery_p[l,2]+no_surgery_p[l,3]+(j-2),mean=params$mean_age_no_surgery, sd=params$sd_age))/ 
            (1-pnorm(no_surgery_p[l,2]+no_surgery_p[l,3]+(j-2),mean=params$mean_age_no_surgery, sd=params$sd_age)))
            
            ## generate menopausal status
            surgery_p[l, ((3*j)+2)]<-rbinom(1,1,surgery_p[l,((3*j)+1)])
            no_surgery_p[l,((3*j)+2)]<-rbinom(1,1,no_surgery_p[l,((3*j)+1)])
      
            ## generating FSH values
            surgery_p[l,((3*j)+3)] <- ifelse(surgery_p[l, ((3*j)+2)]==1, rnorm(1, (0.1+params$alpha6 + params$beta6*(surgery_p[l,2]+surgery_p[l,3]+(j-1))),
                                                                               sqrt(params$alpha5+params$beta5*(surgery_p[l,2]+surgery_p[l,3]+(j-1)))),
                                             rnorm(1, params$alpha3 + params$beta3*(surgery_p[l,2]+surgery_p[l,3]+(j-1)),
                                                   sqrt(params$alpha2+params$beta2*(surgery_p[l,2]+surgery_p[l,3]+(j-1)))))
            
            
            no_surgery_p[l,((3*j)+3)] <- ifelse(no_surgery_p[l, ((3*j)+2)]==1, rnorm(1, (0.1+params$alpha6 + params$beta6*(no_surgery_p[l,2]+no_surgery_p[l,3]+(j-1))),
                                                                                     sqrt(params$alpha5+params$beta5*(no_surgery_p[l,2]+no_surgery_p[l,3]+(j-1)))),
                                                rnorm(1, params$alpha3 + params$beta3*(no_surgery_p[l,2]+no_surgery_p[l,3]+(j-1)),
                                                      sqrt(params$alpha2+params$beta2*(no_surgery_p[l,2]+no_surgery_p[l,3]+(j-1)))))
            
        }
    
        ## saving results
        surg_get_names <- paste0('surgery_p', 1:(j-1))
        no_surg_get_names <- paste0('no_surgery_p', 1:(j-1))
        surgery_fsh <- surgery_p[,((3*j)+3)]
        no_surgery_fsh <- no_surgery_p[,((3*j)+3)]

        for(i in 1:(j-1)){

            surgery_fsh_m <- cbind(surgery_fsh, get(surg_get_names[i])[,((3*(i-1))+6)])
            
            no_surgery_fsh_m <- cbind(no_surgery_fsh,get(no_surg_get_names[i])[,((3*(i-1))+6)])
            
        }
        
        ## getting the mean FSH for each patient
        surgery_fsh <- rowMeans(surgery_fsh_m)
        no_surgery_fsh <- rowMeans(no_surgery_fsh_m)

        ## saving mean FSH in each group
        sims[sims$test == j, 2] <- mean(surgery_fsh)
        sims[sims$test == j, 3] <- mean(no_surgery_fsh)

        ## testing for a difference in FSH between groups 
        test <- wilcox.test(surgery_fsh, no_surgery_fsh)

        ## saving the p value 
        sims[sims$test == j, 4]  <- test$p.value

        ## saving standard deviation of FSH
        sims[sims$test == j, 5]  <- sd(surgery_fsh)
        sims[sims$test == j, 6]  <- sd(no_surgery_fsh)

        ## mean age at surgery
        sims[sims$test == j, 7] <- mean(surgery_p[,2])
        sims[sims$test == j, 8] <-  mean(no_surgery_p[,2])

        ## proportion post-menopausal in each group
        sims[sims$test == j, 9] <- (sum(surgery_p[, ((3*j)+2)])/length(surgery_p[,1]))*100
        sims[sims$test == j, 10] <- (sum(no_surgery_p[, ((3*j)+2)])/length(no_surgery_p[,1]))*100

        ## saving results for this test
        assign(paste0('surgery_p',j), surgery_p)
        assign(paste0('no_surgery_p',j), no_surgery_p)

             
    }

    sims
}

fn.summary<-function(inres){

    sims<-bind_rows(inres)
 
    ## performing power calculation 
    ## power = proportion of simulations where p value in the wilcoxon test is < 0.05
    out <- sims %>% 
        group_by(test) %>% 
        summarise(age = mean(mean_age_surgery),fsh_surgery = mean(mean_fsh_surgery), fsh_nosurgery = mean(mean_fsh_no_surgery),
                  sdfsh_surgery =mean(sd_fsh_surgery), sdfsh_nosurgery = mean(sd_fsh_no_surgery),
                  prop_no_surg = mean(no_surgery_p_postmeno), prop_sur = mean(surgery_p_postmeno), power = sum(p_val < 0.05)/n())

    out

}

# Parameters --------------------------------------------------------------

# Input parameters for the simulation


## all input parameters
params <- list(
  
    ## fsh model parameters - see above for description
    alpha2 = -1.15123373,
    beta2 = 0.03330157, 
    alpha3 = 0.60678151,
    beta3 = 0.05137627,
    alpha5 = 1.80283572, 
    beta5 = -0.02719853,
    alpha6 = 2.69220927,
    beta6 = 0.03432768, 
    
    ## simulation parameters
    ages =  seq(40,50,1), # ages of study entry
#    ages = agedist,
    
    follow_up = seq(1,6,1), # years of follow-up between study entry and first FSH test
    
    n_sims = 1000, # number of simulations to use for power calculation   
    n_tests = 3, # number of yearly FSH tests
    
    n_samp = 700, # number of participants per treatment arm
    
    mean_age_surgery = 47,
    
    mean_age_no_surgery = 49,
    
    sd_age = 4.4
  
)

# Run simulation  --------------------------------------------------------------
## random seed
myseed<-240216

##number cores
ncore<-detectCores()-1 # ncores -1

## random number generator for parallel computing
RNGkind("L'Ecuyer-CMRG")

##base scenario
## set seed reproducibility
set.seed(myseed)
## run sims
myres0 <- mclapply(1:params$n_sims, function(i){fn.sim(i, params)}, mc.cores = ncore)  
res0<-fn.summary(myres0)
res0$power

### alternative age distributiob
## truncated normal
agedist<- floor(rnorm(100000,45,6))
agedist<-agedist[agedist>38]
agedist<-agedist[agedist<55]
summary(agedist)
sd(agedist)
hist(agedist, breaks=seq(min(agedist-0.5), max(agedist)+0.5))
params$ages<-agedist
set.seed(myseed)
myres1 <- mclapply(1:params$n_sims, function(i){fn.sim(i, params)}, mc.cores = ncore)  
res1<-fn.summary(myres1)
res1$power


############################################################
## if had uniform broader age range U(37-53)
params$ages<-seq(37,53)
## set seed reproducibility
set.seed(myseed)
myres2 <- mclapply(1:params$n_sims, function(i){fn.sim(i, params)}, mc.cores = ncore)  
res2<-fn.summary(myres2)
res2$power

############################################################
## if had more narrow age range U(42-47)
params$ages<-seq(42,47)
## set seed reproducibility
set.seed(myseed)
myres3 <- mclapply(1:params$n_sims, function(i){fn.sim(i, params)}, mc.cores = ncore)  
res3<-fn.summary(myres3)
res3$power

