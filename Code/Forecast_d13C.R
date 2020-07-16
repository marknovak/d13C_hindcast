##############################################################################
# Forecast into the future (ignoring individual variation) from time = 0 forward
# (i.e. predict population mean, ignoring random effect)
##############################################################################
rm(list=ls())

library(MASS) # for mvnorm

# Load models and fits
load(file='../Results/ModelFits.Rdata')
# Load processed data
load(file='../Results/ProcessedData.Rdata')

#~~~~~~~~~~~~~~~~~~~~~
# Specify some options
#~~~~~~~~~~~~~~~~~~~~~
# Sample information
t_Fix <- 14 #  days
t_EStor <- 365*4 #  years
# How many samples to draw
samps<-9999
# Number of sample EStorage time-series to make predictions for (for plotting)
n<-200
#~~~~~~~~~~~~~~~~~~~~

Fix <- c(seq(0,t_Fix,length=n/2),rep(t_Fix,n/2))
EStor <-   c(rep(0,n/2),seq(1,t_EStor,length=n/2))
Pres <- Fix+EStor
ndat <- data.frame(Fix,EStor,Pres)

##############################################################################
##################
# Monophasic model
##################
# Extract parameter estimates and (co)variances
fixef <- summary(fit1Cr)$coefficients[,1]
vcov <- summary(fit1Cr)$vcov

#Resample parameter estimates from fixed effects and cov matrix
pars.resamp <- mvrnorm(samps, mu=fixef , Sigma=vcov)

# Predict y value for each parameter sample given new 'data'
preds <- apply(pars.resamp, 1, function(x){
  mod1(ndat$Pres, 
       x['I'], 
       x['A'], 
       x['l1'])})

#determine 2.5 and 97.5% quantiles
dfCI <- data.frame(t(apply(preds, 1, quantile, c(0.025,0.5,0.975))) )
dfM <- apply(preds,1,mean)

#Combine output in data frame
fit1Cr.Predictions <- data.frame( ndat,
                                  Pred=mod1(ndat$Pres, 
                                            fixef['I'], 
                                            fixef['A'], 
                                            fixef['l1']),
                                  mu=dfM,
                                  med50=dfCI$"X50.",
                                  l95=dfCI$"X2.5.",
                                  u95=dfCI$"X97.5.")

##############################################################################
################
# Biphasic model
################
# Extract parameter estimates and (co)variances for biphasic model
fixef <- summary(fit2Cr)$coefficients[,1]
vcov <- summary(fit2Cr)$vcov

# Resample parameter estimates from fixed effects and cov matrix
pars.resamp <- mvrnorm(samps, mu=fixef , Sigma=vcov)

# Predict y value for each parameter sample given new 'data'
preds <- apply(pars.resamp, 1, function(x){
  mod2(ndat$Fix, 
       ndat$EStor, 
       x['I'], 
       x['A'], 
       x['l1'], 
       x['l2'])})

#determine 2.5 and 97.5% quantiles
dfCI <- data.frame(t(apply(preds, 1, quantile, c(0.025,0.5,0.975))) )
dfM <- apply(preds,1,mean)

#Combine output in data frame
fit2Cr.Predictions <- data.frame( ndat,
                                  Pred=mod2(ndat$Fix, 
                                            ndat$EStor, 
                                            fixef['I'], 
                                            fixef['A'], 
                                            fixef['l1'], 
                                            fixef['l2']),
                                  mu=dfM,
                                  med50=dfCI$"X50.",
                                  l95=dfCI$"X2.5.",
                                  u95=dfCI$"X97.5.")


###############################################################################
#######
# Plots
#######
pdf('../Results/Figures/Fig2_alt-Forecasts_d13C_monophasic.pdf',height=4,width=3)
par(mfcol=c(2,1),cex=0.8,cex.axis=0.9,cex.lab=1.1,
    tcl=-0.2,las=1,mgp=c(2.5,0.3,0))

  # ylims for forecasts
  ylims_f <-range(c(fit1Cr.Predictions$u95,fit1Cr.Predictions$l95))
  
  # ylims for uncertainties
  ylims_u <- c(0,max(c(fit1Cr.Predictions$u95-fit1Cr.Predictions$l95)))
  
  par(mar=c(1,4,2.5,0.75))
  out <- fit1Cr.Predictions
  ylims<-range(c(out$l95,
                 out$u95))
  plot(out$Pres,out$Pred, 
       type='n',ylim=ylims_f,
       ylab=expression(paste(delta^13,'C')),
       xlab='')
    polygon(c(out$Pres,rev(out$Pres)), 
            c(out$l95,rev(out$u95)), border=NA,col='grey80')
    lines(out$Pres,
          out$Pred,
          col='black',lwd=2)

  
  # Intercept parameter estimate
  points(0,fixef['I'],pch=19)
  
  # Intercept + Asymptotic effect parameter estimates
  # Remember not to think of 'A' as the asymptote itself.
  # The asymptotic value is I+A.
  # points(max(Pres),(fixef['I']+fixef['A']),pch=22,bg='black')
  abline(h=fixef['I']+fixef['A'],lty=3)
  box(lwd=1)
  legend('topright',legend='A',bty='n',cex=1.5,inset=0)

  par(mar=c(2.5,4,1,0.75))
  out$uncert <- out$u95-out$l95
  ylims<-c(0,max(out$uncert))
  plot(out$Pres,
       out$uncert, 
       ylab=expression(paste(delta^13,'C uncertainty')),
       xlab='',
       ylim=ylims_u,
       type='l',lwd=2,col='black')
  title(xlab='Total preservation time (days)',line=1.5)
  box(lwd=1)
  legend('topright',legend='B',bty='n',cex=1.5,inset=0)

dev.off()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf('../Results/Figures/Fig2-Forecasts_d13C_mono_and_biphasic.pdf',height=4,width=5.5)
par(mfcol=c(2,2),cex=0.8,cex.axis=0.9,cex.lab=1.1,
    tcl=-0.2,las=1,mgp=c(2.5,0.3,0))

  # ylims for forecasts
  ylims_f <-range(c(fit1Cr.Predictions$u95,fit1Cr.Predictions$l95,
                    fit2Cr.Predictions$u95,fit2Cr.Predictions$l95))
  
  # ylims for uncertainties
  ylims_u <- c(0,max(c(fit1Cr.Predictions$u95-fit1Cr.Predictions$l95,
                     fit2Cr.Predictions$u95-fit2Cr.Predictions$l95)))
  
    # ~~~~~~~~~~~~~~
    # Monophasic model
    # ~~~~~~~~~~~~~~
  par(mar=c(1,4,2.5,0))
    out <- fit1Cr.Predictions
    ylims<-range(c(out$l95,
                   out$u95))
    plot(out$Pres,out$Pred, 
         type='n',ylim=ylims_f,
         ylab=expression(paste(delta^13,'C')),
         xlab='')
      polygon(c(out$Pres,rev(out$Pres)), 
              c(out$l95,rev(out$u95)), border=NA,col='grey80')
      lines(out$Pres,
            out$Pred,
            col='black',lwd=2)

      # Intercept parameter estimate
      points(0,fixef['I'],pch=19)
      
      # Intercept + Asymptotic effect parameter estimates
      # Remember not to think of 'A' as the asymptote itself.
      # The asymptotic value is I+A.
      # points(max(Pres),(fixef['I']+fixef['A']),pch=22,bg='black')
      abline(h=fixef['I']+fixef['A'],lty=3)
      box(lwd=1)
      legend('topright',legend='A',bty='n',cex=1.5,inset=0)
  
  par(mar=c(2.5,4,1,0))
    out$uncert <- out$u95-out$l95
    ylims<-c(0,max(out$uncert))
    plot(out$Pres,
         out$uncert, 
         ylab=expression(paste(delta^13,'C uncertainty')),
         xlab='',
         ylim=ylims_u,
         type='l',lwd=2,col='black')
    title(xlab='Total preservation time (days)',line=1.5)
    box(lwd=1)
    legend('topright',legend='B',bty='n',cex=1.5,inset=0)
    
    # ~~~~~~~~~~~~~~~~
    # Biphasic model
    # ~~~~~~~~~~~~~~~~
    par(mar=c(1,3,2.5,1))
    out <- fit2Cr.Predictions
    plot(out$Pres,out$Pred, 
         type='n',ylim=ylims_f,
         ylab='',
         xlab='')
    polygon(c(out$Pres,rev(out$Pres)), 
            c(out$l95,rev(out$u95)), border=NA,col='grey80')
    lines(out$Pres,
          out$Pred,
          col='black',lwd=2)

    # Intercept parameter estimate
    points(0,fixef['I'],pch=19)
    
    # Intercept + Asymptotic effect parameter estimates
    # Remember not to think of 'A' as the asymptote itself.
    # The asymptotic value is I+A.
    # points(max(Pres),(fixef['I']+fixef['A']),pch=22,bg='black')
    abline(h=fixef['I']+fixef['A'],lty=3)
    box(lwd=1)
    legend('topright',legend='C',bty='n',cex=1.5,inset=0)
    
    
    par(mar=c(2.5,3,1,1)) 
    out$uncert <- out$u95-out$l95
    ylims<-c(0,max(out$uncert))
    plot(out$Pres,
         out$uncert, 
         ylab='',
         xlab='',
         ylim=ylims_u,
         type='l',lwd=2,col='black')
    title(xlab='Total preservation time (days)',line=1.5)
    box(lwd=1)
    legend('topright',legend='D',bty='n',cex=1.5,inset=0)  

dev.off()

# Uncertainty increases rapidly because it's the fixation effect (l1) that's
# poorly contrained.  It subsequently decreases both because the EStorage 
# effect (l2) is (relatively) well-constrained (i.e. it's clearly negative),
# but also because the asymptote (I+A, which is also well-constrained) 
# is being approached (i.e. if you EStore the sample long enough then its
# final value is more easily predicted).
################################################################################
################################################################################
################################################################################
