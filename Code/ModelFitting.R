################################################
# Fit and compare nonlinear mixed-effects models
################################################
rm(list=ls())

library(nlme)
library(lme4)
library(bbmle) # for *ICtab()

# lme4 has superceded nlme.
# https://stats.stackexchange.com/questions/5344/how-to-choose-nlme-or-lme4-r-library-for-mixed-effects-models
# We use both and compare.
# Conclusion: same pt estimates, tighter errors using lme4

######################
# Import prepared data
######################
load(file='../Results/ProcessedData.Rdata')

#############################
# Mean initial specimen value
#############################
IestC <- mean(dat$d13C[dat$EStor==0], na.rm = T)
IestN <- mean(dat$d15N[dat$EStor==0], na.rm = T)

#########################################################################################
#########################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (Random) intercept-only model (i.e. no effect of fixation and storage)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~
#~~~~ Nitrogen ~~~~~
#~~~~~~~~~~~~~~~~~~~
fit0N<-nlme(d15N ~ I,
            fixed=I~1,
            random=I~1|SpecimenID,
            start=c(I=IestN),
            data=dat, verbose=F,
            na.action = na.omit)
summary(fit0N)

#~~~~~~~~~~~~~~~~~
#~~~~ Carbon ~~~~~
#~~~~~~~~~~~~~~~~~
fit0C<-nlme(d13C ~ I,
           fixed=I~1,
           random=I~1|SpecimenID,
           start=c(I=IestC),
           data=dat, verbose=F)
summary(fit0C)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~  Repeat null model using lme4:nlmer ~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod0 <- function(I){ I }
mod0gr <- deriv(body(mod0)[[2]],
                namevec = c('I'),
                function.arg=mod0)

#~~~~ Nitrogen ~~~~~
starts <- apply(coef(fit0N),2,mean)

fit0Nr <- nlmer(
  # Response
  d15N ~ 
    # Fixed effects
    mod0gr(I) ~ 
    # Random effects
    (I | SpecimenID), 
  # Data
  data = dat, 
  start = starts,
  na.action = na.omit)

summary(fit0Nr)

#~~~~ Carbon ~~~~~
starts <- apply(coef(fit0C),2,mean)

fit0Cr <- nlmer(
  # Response
  d13C ~ 
    # Fixed effects
    mod0gr(I) ~ 
    # Random effects
    (I | SpecimenID), 
  # Data
  data = dat, 
  start = starts)

summary(fit0Cr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Monophasic model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~
#~~~~ Nitrogen ~~~~~
#~~~~~~~~~~~~~~~~~~~
starts <- c(I=IestN, A=-0.03,l1=0.001)

fit1N<-nlme(d15N ~ I + A*(1-exp(l1*Pres)),
            fixed=I+A+l1~1,
            random=I~1|SpecimenID,
            start=starts,
            data=dat, verbose=F,
            na.action = na.omit)
summary(fit1N)

#~~~~~~~~~~~~~~~~~
#~~~~ Carbon ~~~~~
#~~~~~~~~~~~~~~~~~
starts <- c(I=IestC,A=-0.9,l1=-.001)

fit1C<-nlme(d13C ~ I + A*(1-exp(l1*Pres)),
            fixed=I+A+l1~1,
            random=I~1|SpecimenID,
            start=starts,
            data=dat, verbose=F)
summary(fit1C)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~  Repeat monophasic model using lme4:nlmer ~~~~~~~~~
# ~~~~~~      in order to get better estimates of  ~~~~~~~~~
# ~~~~~~                 uncertanties              ~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod1 <- function(Pres, I, A, l1){
  I + A*(1-exp(l1*Pres))
}
mod1gr <- deriv(body(mod1)[[2]],
                namevec = c('I','A','l1'),
                function.arg=mod1)

#~~~~ Carbon ~~~~~
starts <- apply(coef(fit1C),2,mean)
starts <- list(nlpars= apply(coef(fit1C),2,mean), theta=vcov(fit1C))

fit1Cr <- nlmer(
  # Response
  d13C ~ 
    # Fixed effects
    mod1gr(Pres=Pres, I, A, l1) ~ 
    # Random effects
    (I | SpecimenID), 
  # Data
  data = dat, 
  start = starts,
  verbose=F,
  control=nlmerControl(tolPwrss=1E-4))

summary(fit1Cr)

##############################################################################
##############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Biphasic decay model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~
#~~~~ Nitrogen ~~~~~
#~~~~~~~~~~~~~~~~~~~
starts <- c(apply(coef(fit1N),2,mean),l2=0.001)

fit2N<-nlme(d15N ~ I + A*(1-exp(l1*Fix)*exp(l2*EStor)),   
            fixed=I+A+l1+l2~1,
            random=I~1|SpecimenID,
            start=starts,
            data=dat, 
            verbose=F,
            control = nlmeControl(pnlsTol = 0.07, msVerbose = TRUE),
            na.action = na.omit)
summary(fit2N)

#~~~~~~~~~~~~~~~~~
#~~~~ Carbon ~~~~~
#~~~~~~~~~~~~~~~~~
starts <- c(apply(coef(fit1C),2,mean),l2=-0.001)

fit2C<-nlme(d13C ~ I + A*(1-exp(l1*Fix)*exp(l2*EStor)),
           fixed=I+A+l1+l2~1,
           random=I~1|SpecimenID,
           start=starts,
           data=dat, verbose=F)
summary(fit2C)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~  Repeat biophasic model using lme4:nlmer ~~~~~~~~~
# ~~~~~~     in order to get better estimates of  ~~~~~~~~~
# ~~~~~~                uncertanties              ~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod2 <- function(Fix, EStor, I, A, l1, l2){
    I + A*(1-exp(l1*Fix)*exp(l2*EStor)) # equivalent to exp( l1*Fix + l2*EStor )
}
mod2gr <- deriv(body(mod2)[[2]],
                namevec = c('I','A','l1','l2'),
                function.arg=mod2)

#~~~~ Carbon ~~~~~
starts <- apply(coef(fit2C),2,mean)

fit2Cr <- nlmer(
  # Response
  d13C ~ 
  # Fixed effects
  mod2gr(Fix=Fix, EStor=EStor, I, A, l1, l2) ~ 
  # Random effects
  (I | SpecimenID), 
  # Data
  data = dat, 
  start = starts)

summary(fit2Cr)

##############################################################################
##############################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Information-theoretic model comparison
#        of all three models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~
#~~~~ Nitrogen ~~~~~
#~~~~~~~~~~~~~~~~~~~
bbmle::BICtab(fit0N,fit1N,fit2N, logLik=TRUE, base=TRUE,sort=FALSE, weights=TRUE)
bbmle::AICctab(fit0N,fit1N,fit2N, logLik=TRUE, base=TRUE,sort=FALSE, weights=TRUE)

#~~~~~~~~~~~~~~~~~
#~~~~ Carbon ~~~~~
#~~~~~~~~~~~~~~~~~
bbmle::BICtab(fit0C,fit1C,fit2C, logLik=TRUE, base=TRUE,sort=FALSE, weights=TRUE)
bbmle::AICctab(fit0C,fit1C,fit2C, logLik=TRUE, base=TRUE,sort=FALSE, weights=TRUE)

bbmle::BICtab(fit0Cr,fit1Cr,fit2Cr, logLik=TRUE, base=TRUE,sort=FALSE, weights=TRUE)
bbmle::AICctab(fit0Cr,fit1Cr,fit2Cr, logLik=TRUE, base=TRUE,sort=FALSE, weights=TRUE)

##############################################################################
##############################################################################
# Save models and fits
save(file='../Results/ModelFits.Rdata',
     fit0N,fit1N,fit2N,
     fit0C,fit1C,fit2C,
     fit1Cr,fit2Cr,
     mod1,mod1gr,
     mod2,mod2gr)

#############################################################################
#############################################################################
# Comment on model-fitting methods:
# lme4:nlmer producted tighter uncertainty bounds around point estimates
#  Why?  Haven't bothered to find out, but it's the replacement for nmle 
#  so perhaps there are indeed some improvements.
#  It could also be because we're providing the gradient to nlmer.
#  Note also the nlmer doesn't count degrees of freedom
#############################################################################
#############################################################################
#############################################################################
#############################################################################





