######################################################################
# Apply Hindcast model to data
# Below, this is done twice:
# 1) to create a time-series for a single hypothetical specimen
# 2) to infer the initial values from the observed (final) values of a set of specimens
######################################################################
rm(list=ls())

# load functions for hindcasting
source(file='Hindcast_functions.R')

# load processed data
load(file='../Results/ProcessedData.Rdata')

# load fit models
load(file='../Results/ModelFits.Rdata')
# The parameter estimates of the fits are hard-coded into the above functions
# when a fit is not provided, so the fits don't necessarily need to be loaded

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get means and SDs of last specimens from dataset
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract last sample for each specimen
lasts <- as.vector(by(dat, dat$SpecimenID,
                      function(x){rownames(x)[which.max(x$Pres)]})) 
ldat <- dat[lasts,] 

# Remove experimental specimens
ldat <- subset(ldat,Source=='Field')

ldat.mean <- mean(ldat$d13C)
ldat.sd <- as.numeric(names(sort(table(ldat$d13C_sd),decreasing=TRUE)[1]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Apply time-series creation method to a single hypothetical sample
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ts.preds <- Hindcast_d13C(obs_d13C = ldat.mean,
                           precision.dp=2,
                           obs_d13C_SD = ldat.sd,
                           t_Fix = 10,
                           t_EStor = 365*10,
                           model = 'Monophasic',
                           fit = fit1Cr,
                           createTimeSeries = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Infer initial value from final (observed) value for many specimens
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make hindcasts from final observed value
preds <- Hindcast_d13C(specimenID=ldat$SpecimenID,
                       obs_d13C = ldat$d13C,
                       obs_d13C_SD = ldat$d13C_sd,
                       t_Fix = ldat$Fix,
                       t_EStor = ldat$EStor,
                       fit = fit1Cr,
                       model = 'Monophasic')

# Compare to first sample for each specimen
firsts <- as.vector(by(dat, dat$SpecimenID,
                       function(x){rownames(x)[which.min(x$Pres)]})) 
fdat <- dat[firsts,] 
fdat <- subset(fdat,Source=='Field')

# confidence interval for observed values
int95 <- 1.96 * fdat$d13C_sd/sqrt(fdat$d13C_sd_n)
fdat$CI.97.5 <- fdat$d13C + int95
fdat$CI.2.5 <- fdat$d13C - int95

# Combine predictions with first observed value
fdat.preds <- merge(
  fdat[,c('SpecimenID','d15N','d13C','Fix','EStor','Pres','CI.2.5','CI.97.5')],
  preds,
  by.x='SpecimenID',
  by.y='specimenID')



##################################################################
# Combine time-series and multi-specimen hindcasts into one figure
##################################################################

fill.cols.95 <- c('grey80','grey70')
fill.cols.68 <- c('grey60','grey45')

pdf('../Results/Figures/Fig3-Hindcasts_d13C_monophasic.pdf',height=4,width=3)
par(mfcol=c(2,1),
  cex=0.8,cex.axis=0.9,cex.lab=1.1,
    tcl=-0.2,las=1,mgp=c(2.5,0.3,0))
par(mar=c(3,4,0.5,0.5))
  ylims<-range(c(ldat.mean,ts.preds$pred.CI.2.5,ts.preds$pred.CI.97.5))
  plot(ts.preds$t_Pres, ts.preds$pred.mean,
       ylab=expression(paste(delta^13,'C')),
       xlab='',
       type='n', ylim=ylims)
    polygon(c(ts.preds$t_Pres, rev(ts.preds$t_Pres)) , 
            c(ts.preds$pred.CI.2.5, rev(ts.preds$pred.CI.97.5)),
            border=NA,col=fill.cols.95[2])
    polygon(c(ts.preds$t_Pres, rev(ts.preds$t_Pres)) , 
            c(ts.preds$pred.CI.16, rev(ts.preds$pred.CI.84)),
            border=NA,col=fill.cols.68[2])
    trans<-ts.preds$t_Pres<=ts.preds$t_Pres_mx[1] # to indicate asymptote precision transition
    # trans<-ts.preds$t_Fix<max(ts.preds$t_Fix) # to indicate fixation transition
    polygon(c(ts.preds$t_Pres[trans], rev(ts.preds$t_Pres[trans])) ,
            c(ts.preds$pred.CI.2.5[trans], rev(ts.preds$pred.CI.97.5[trans])),
            border=NA,col=fill.cols.95[1])
    polygon(c(ts.preds$t_Pres[trans], rev(ts.preds$t_Pres[trans])) ,
            c(ts.preds$pred.CI.16[trans], rev(ts.preds$pred.CI.84[trans])),
            border=NA,col=fill.cols.68[1])
    lines(ts.preds$t_Pres, ts.preds$pred.point,
          lwd=2)
    points(ts.preds$t_Pres[1], ts.preds$pred.point[1],
           pch=19,lwd=2)
    title(xlab='Total preservation time (days)',line=1.5)
    legend('topright',legend='A',bty='n',cex=1.5,inset=0)

par(mar=c(3,4,0.5,0.5),mgp=c(2,0.3,0))
par(pty='s',las=1)
  xlims <- range(c(fdat.preds$d13C, fdat.preds$pred.point))
  plot(fdat.preds$d13C, fdat.preds$pred.point, 
       xlim=xlims, ylim=xlims, type='n',
       xlab='', 
       ylab=expression(paste('Predicted initial ',delta^13,'C')))
    abline(0,1,lty=2,col='grey60')
    arrows(fdat.preds$d13C, fdat.preds$pred.CI.2.5, 
           fdat.preds$d13C, fdat.preds$pred.CI.97.5, 
           length=0, code=3, angle=90)
    arrows(fdat.preds$CI.2.5, fdat.preds$pred.point, 
           fdat.preds$CI.97.5, fdat.preds$pred.point, 
           length=0, code=3, angle=90)
    points(fdat.preds$d13C, fdat.preds$pred.point,
           pch=21, bg='grey',cex=0.3)
    title(xlab=expression(paste('Actual initial ',delta^13,'C')),line=1.5)
    legend('bottomright',legend='B',bty='n',cex=1.5,inset=0)
    box(lwd=1)

dev.off()

##################################################################
##################################################################
##################################################################