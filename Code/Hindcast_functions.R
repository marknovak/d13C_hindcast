##############################################################################
##############################################################################
# Hindcast time series to time = 0 given an specimen's reading at t=x
##############################################################################
##############################################################################
# NOTE: 
# The naive approach works fine as long as the length of 
# fixation and EStorage aren't very long.  Given the empirical parameter 
# estimates and the assumed uncertainty/precison of the sample, 
# once total duration gets past about a year the hindcast
# predictions become nonsensical (not accounting for asymptotic effect).

# The naive approach simply reverses the forward model.  It therefore
# allows values to increase to infinity, failing to consider the asymptotic
# EStorage effect.

# The naive approach also fails to consider the fact that fixation
# occurred before EStorage.  In running backwards, it's therefore possible
# for the EStorage effect to "use up" the whole asymptotic offset when the 
# EStorage time is long enough such that no fixation even "occurs".

# The "non-naive" approach fixes these issues.

##############################################################################
# Define functions
# See "Hindcast_d13C.R" for application to data.
##############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to determine the precision of a number
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate the maximum possible periods of time before either
# fixation or the combination of fixation and EStorage could have proceeded
# before their asymptotic effect was achieved within 
# one estimated standard deviation of the sample's measured value
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Max_times <- function(t_Fix=NULL, t_Pres=NULL, obs_d13C, precision.dp=NULL, A, l1, l2=NULL){
  
  if(is.null(t_Fix) & is.null(t_Pres) | 
     !is.null(t_Fix) & !is.null(t_Pres) | 
     is.null(obs_d13C) & is.null(precision.dp) | 
     is.null(A) | is.null(l1)){
    stop('Too much or too little information provided to "Max_times()" function.')
  }
  
  if(!is.null(precision.dp)){
    # if precision is provided (in numbers of decimal places)
    # we cap at 1E-16 which is typical computer "machine" precision
    pr <- 10^(-precision.dp) - 1E-16
  }else{
    # if not, use precision of provided observation
    # (Note that this fails if data is provided as, for example, "obs_d13C=-24.00"
    # because R rounds these number to -24 automatically)
    if(decimalplaces(obs_d13C)==0){
      warning('Precision assumed to be 0 d.p. due to provided value of "obs_d13C" (e.g., -24.00 is read as -24 rather than being at 2 d.p.).  Provide precision using "precision.dp" instead.',
              call. = FALSE,
              immediate. = TRUE)
      }
    pr <- (10^(-decimalplaces(obs_d13C))) - 1E-16
  }
  
  names(A) <- names(l1) <- names(l2) <- NULL
  
  # For monophasic model
  if(!is.null(t_Pres)){
    # The maximum time before preservation (fixation + storage) effect reaches asymptotic effect
    # to within precision
    t_Pres_mx <- log(abs(pr/A))/l1
    out<-c('t_Pres_mx'=t_Pres_mx)
    return(out)
  }
  
  # For biphasic model
  if(!is.null(t_Fix)){
    # The maximum time before fixation effect alone reaches asymptotic effect
    # to within precision
    t_Fix_mx <- log(abs(pr/A))/l1
    # The predicted amount of the asymptotic effect remaining at the end of actual fixation
    Fix_pred <- A*exp(l1*t_Fix)
    # The maximum time after fixation before ethanol storage effect reaches asymptotic effect
    # to within precision
    t_EStor_mx <- log(abs(pr/Fix_pred))/l2
    if(is.infinite(t_EStor_mx) | t_EStor_mx<0){  t_EStor_mx <- 0  }
    out<-c('t_Fix_mx'=t_Fix_mx, 
           't_EStor_mx'=t_EStor_mx)
    return(out)
  }
}

#~~~~~~~~~~~~~~~~~~
# Test the function
#~~~~~~~~~~~~~~~~~~
# Max_times(t_Pres=100, obs_d13C=-25.5, A=-1, l1=-0.01)

# Max_times(t_Pres=100, obs_d13C=-25.5, precision.dp=1, A=-1, l1=-0.01)

# Max_times(t_Fix=10, obs_d13C=-25.5, A=-1, l1=-0.01, l2=-0.01)

# Max_times(t_Fix=10, obs_d13C=-25.5, precision.dp=1, A=-1, l1=-0.01, l2=-0.01)

# Max_times(t_Fix=10, obs_d13C=-25.5, precision.dp=2, A=-1, l1=-0.01, l2=-0.01)

# Max_times(t_Fix=200, obs_d13C=-25.5, A=-1, l1=-0.01, l2=-0.01)

# Max_times(t_Fix=250, obs_d13C=-25.5, A=-1, l1=-0.01, l2=-0.01)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to hindcast sample(s) given (hardcoded) fits from models
# (If needed, this could be sped up by using an 'apply' function instead of 
# looping through each of the samples.)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Hindcast_d13C<-function(obs_d13C,
                        specimenID=NULL,
                        precision.dp=NULL,
                        obs_d13C_SD=NULL, 
                        obs_d13C_SE=NULL, 
                        obs_d13C_n=NULL,
                        t_Pres=NULL,
                        t_Fix=NULL, 
                        t_EStor=NULL,
                        model='Monophasic',
                        fit=NULL,
                        naive=FALSE,
                        reSamps=999*2,
                        CI.quantiles=c(0.025, 0.16, 0.84, 0.975),
                        createTimeSeries=FALSE,
                        n_TimeSeries=100,
                        progressBar=TRUE){
  
  require(MASS)
  
  if(is.null(specimenID)){specimenID <- 1:length(obs_d13C)}
  
  if(!missing(model) & length(model)>1) stop("Only one 'model' allowed at a time.")
  
  model<-match.arg(model,choices=c('Biphasic','Monophasic'))
  
  if(createTimeSeries & length(obs_d13C)>1){
    stop(paste('A time series can only be created for single sample but you specified ',length(obs_d13C)))
  }
  
  if(all(c(is.null(obs_d13C_SD),is.null(obs_d13C_SE)))){
    stop("Provide either Standard Deviation or Standard Error (and sample size) of lab standard or control samples in order apply it to focal sample.")
  }
  
  if(is.null(obs_d13C_SD)){
    obs_d13C_SD <- obs_d13C_SE * sqrt(obs_d13C_n)
  }
  
  if(length(obs_d13C_SD)==1 & length(obs_d13C)>1){
    warning('Using same Standard Deviation for all samples.  If not desired, specify SD or Standard Error (and sample size) of control samples for each sample.',
            call. = FALSE,
            immediate. = TRUE)
    obs_d13C_SD <- rep(obs_d13C_SD, length(obs_d13C))
  }
  
  if(is.null(t_Pres)){
    t_Pres <- t_Fix + t_EStor
  }
  
  if(any(t_Pres != t_Fix + t_EStor)){
    warning("The durations of fixation and ethanol storage don't all add up to the provided total preservation time.",
            call. = FALSE,
            immediate. = TRUE)
  }
  
  if(length(obs_d13C)!=length(t_Pres) |
     length(obs_d13C)!=length(t_Fix) |
     length(obs_d13C)!=length(t_EStor)){
    warning("Using same time periods of total preservation, fixation and storage for all samples. If not desired, provide them for each sample.",
            call. = FALSE,
            immediate. = TRUE)
    t_Pres <- rep(t_Pres, length(obs_d13C))
    t_Fix <- rep(t_Fix, length(obs_d13C))
    t_EStor <- rep(t_EStor, length(obs_d13C))
  }
  
  if(is.null(precision.dp)){
    warning("No precison for obs_d13C provided, so using decimal places of obs_d13C.",
            call. = FALSE,
            immediate. = TRUE)
  }
  
  # Define hind-casting model
  if(model=='Monophasic'){
    hmod <- function(t_Pres, obs_d13C, precision.dp, A, l1, naive){
      mx_times <- Max_times(t_Pres=t_Pres, 
                            obs_d13C=obs_d13C, 
                            precision.dp=precision.dp, 
                            A=A, l1=l1)
      if(!naive & t_Pres >= mx_times['t_Pres_mx']){
        n_t_Pres <- mx_times['t_Pres_mx']
      }else{
        n_t_Pres <- t_Pres
      }
      pred <-  obs_d13C - A*(exp(l1*n_t_Pres)) 
      return(pred)
    }
  }
  if(model=='Biphasic'){
    hmod <- function(t_Fix, t_EStor, obs_d13C, precision.dp, A, l1, l2, naive){
      mx_times <- Max_times(t_Fix=t_Fix, 
                            obs_d13C=obs_d13C, 
                            precision.dp=precision.dp,
                            A=A, l1=l1, l2=l2)
      if(!naive & t_Fix >= mx_times['t_Fix_mx']){
        n_t_Fix <- mx_times['t_Fix_mx']
        n_t_EStor <- 0
      }
      if(!naive & t_Fix < mx_times['t_Fix_mx'] & t_EStor >= mx_times['t_EStor_mx']){
        n_t_Fix <- t_Fix
        n_t_EStor <- mx_times['t_EStor_mx']
      }  
      if(naive | t_Fix < mx_times['t_Fix_mx'] & t_EStor < mx_times['t_EStor_mx']){
        n_t_Fix <- t_Fix
        n_t_EStor <- t_EStor
      }
      pred <- obs_d13C - A*(exp(l1*(n_t_Fix))*exp(l2*(n_t_EStor)))
      return(pred)
    }
  }
  
  # Define point estimates and Variance-Covariance matrix
  # Extract directly from the fit if it is provided, otherwise use hard-coded estimates
  # Note that estimates for 'I' are removed as these are what is being predicted.
  if(!is.null(fit)){
      fixef <- summary(fit)$coefficients[-1,1]
      vcov <- summary(fit)$vcov[-1,-1]
      if(model=='Monophasic' & length(fixef)!=2){
        stop("The number of parameters estimated in the provided model fit
             is not 3 (as assumed the monophasic model)")}
      if(model=='Biphasic' & length(fixef)!=3){
        stop("The number of parameters estimated in the provided model fit
             is not 4 (as assumed under the biphasic model)")}
  }
  if(is.null(fit)){ # hard-coded estimates from Taylor et al. data
    if(model=='Monophasic'){
      fixef <- c( A = -0.923724613, 
                  l1 = -0.003135802)
      
      vcov <- matrix(c( 4.964753e-03, -2.770482e-05,
                        -2.770482e-05,  1.173272e-06),
                     ncol=length(fixef),byrow=T,
                     dimnames=list( names(fixef), names(fixef) ))
    }
    if(model=='Biphasic'){
      fixef <- c( A = -0.922786952, 
                  l1 = -0.000449796, 
                  l2 = -0.003208472)
      
      vcov <- matrix(c(  8.343273e-05,  7.522671e-08, -4.282089e-07,
                         7.522671e-08,  5.162026e-05, -6.940043e-07,
                        -4.282089e-07, -6.940043e-07,  1.496840e-06),
                     ncol=length(fixef),byrow=T,
                     dimnames=list( names(fixef), names(fixef) ))
    }
  }
  
  # Resample parameter estimates from fixed effects and cov matrix
  pars.resamp.e <- mvrnorm(reSamps, mu=fixef, Sigma=vcov)
  
  # The originals are going to be over-written, so keep them
  t_Fix_obs <- t_Fix
  t_EStor_obs <- t_EStor
  t_Pres_obs <- t_Pres
  obs_d13C_obs <- obs_d13C
  obs_d13C_SD_obs <- obs_d13C_SD
  
  if(createTimeSeries){
    n <- n_TimeSeries
    t_Fix <- rev(c(seq(0,t_Fix,length=n/2),rep(t_Fix,n/2)))
    t_EStor <- rev(c(rep(0,n/2),seq(1,t_EStor,length=n/2)))
    t_Pres <- t_Fix + t_EStor
    obs_d13C <- rep(obs_d13C,n)
    obs_d13C_SD <- rep(obs_d13C_SD,n)
    
    t_Fix_obs <- rep(t_Fix_obs,n)
    t_EStor_obs <- rep(t_EStor_obs,n)
    t_Pres_obs <- rep(t_Pres_obs,n)
    specimenID <- rep(specimenID,n)
  }
  
  if(!createTimeSeries){
    t_Fix <- rep(0,length(obs_d13C))
    t_EStor <- rep(0,length(obs_d13C))
    t_Pres <- t_Fix + t_EStor
  }
  
  pred.point <- pred.mean <- pred.median <- t_Pres_mx <- rep(NA,length(obs_d13C))
  t_FixEStor_mx <- array(NA, dim=c(length(obs_d13C),2))
  pred.CI <- array(NA, dim=c(length(obs_d13C),length(CI.quantiles)))
  
  if(progressBar){
    # create progress bar
    pb <- txtProgressBar(min = 0, max = length(obs_d13C), style = 3)
  }
  for(i in 1:length(obs_d13C)){ # this loop could be converted to apply-style 
    pars.resamp <- cbind(obs_d13C=rnorm(reSamps, obs_d13C[i], obs_d13C_SD[i]), 
                         pars.resamp.e)
    
    # Predict value for each parameter sample
    if(model=='Monophasic'){
      preds <- apply(pars.resamp, 1, function(x){ hmod(t_Pres=t_Pres[i], 
                                                       obs_d13C=x['obs_d13C'], 
                                                       precision.dp=precision.dp,
                                                       A=x['A'],
                                                       l1=x['l1'],
                                                       naive=naive)})
      preds.adj <- apply(pars.resamp, 1, function(x){ x['obs_d13C'] -
                                                      hmod(t_Pres=t_Pres_obs[i],
                                                       obs_d13C=x['obs_d13C'],
                                                       precision.dp=precision.dp,
                                                       A=x['A'],
                                                       l1=x['l1'],
                                                       naive=naive)})
      preds <- preds + preds.adj

      pred.point[i] <- hmod( t_Pres=t_Pres[i],
                              obs_d13C=obs_d13C[i],
                              precision.dp=precision.dp,
                              A=fixef['A'],
                              l1=fixef['l1'],
                              naive=naive)
      pred.point.adj <- hmod( t_Pres=t_Pres_obs[i],
                              obs_d13C=obs_d13C[i],
                              precision.dp=precision.dp,
                              A=fixef['A'],
                              l1=fixef['l1'],
                              naive=naive)
      pred.point[i] <- pred.point[i] + obs_d13C[i] - pred.point.adj
        
      t_Pres_mx[i] <- Max_times( t_Pres = t_Pres[i], 
                                 obs_d13C = obs_d13C[i], 
                                 precision.dp=precision.dp,
                                 A = fixef['A'], 
                                 l1 = fixef['l1'])
    }
    
    if(model=='Biphasic'){
      preds <- apply(pars.resamp, 1, function(x){ hmod(t_Fix=t_Fix[i], 
                                                       t_EStor=t_EStor[i], 
                                                       obs_d13C=x['obs_d13C'], 
                                                       precision.dp=precision.dp,
                                                       A=x['A'], 
                                                       l1= x['l1'], 
                                                       l2=x['l2'],
                                                       naive=naive)}) 
      preds.adj <- apply(pars.resamp, 1, function(x){ x['obs_d13C'] -
                                                      hmod(t_Fix=t_Fix_obs[i], 
                                                           t_EStor=t_EStor_obs[i], 
                                                           obs_d13C=x['obs_d13C'], 
                                                           precision.dp=precision.dp,
                                                           A=x['A'], 
                                                           l1= x['l1'], 
                                                           l2=x['l2'],
                                                           naive=naive)}) 
      preds <- preds + preds.adj

      pred.point[i] <- hmod( t_Fix=t_Fix[i],
                              t_EStor=t_EStor[i],
                              obs_d13C=obs_d13C[i], 
                              precision.dp=precision.dp,
                              A=fixef['A'], 
                              l1=fixef['l1'],
                              l2=fixef['l2'],
                              naive=naive)
      pred.point.adj <- hmod( t_Fix=t_Fix_obs[i],
                              t_EStor=t_EStor_obs[i],
                              obs_d13C=obs_d13C[i], 
                              precision.dp=precision.dp,
                              A=fixef['A'], 
                              l1=fixef['l1'],
                              l2=fixef['l2'],
                              naive=naive)
      pred.point[i] <- pred.point[i] + obs_d13C[i] - pred.point.adj
      
      t_FixEStor_mx[i,] <- Max_times(t_Fix = t_Fix[i],
                                     obs_d13C = obs_d13C[i], 
                                     precision.dp=precision.dp,
                                     A = fixef['A'], 
                                     l1 = fixef['l1'],
                                     l2 = fixef['l2'])
    }
    
    # Determine mean, quantiles, and median
    pred.mean[i] <- mean(preds)
    pred.CI[i,] <- quantile(preds, CI.quantiles)
    pred.median[i] <- quantile(preds,0.5)
    
    if(progressBar){
      setTxtProgressBar(pb, i)
    }
    
  }
  
  if(progressBar){ close(pb) }
  
  # Package results
  colnames(pred.CI) <- 100*CI.quantiles

  Out <- data.frame(specimenID = specimenID,
                    obs_d13C = obs_d13C,
                    obs_d13C_SD = obs_d13C_SD,
                    t_Fix = t_Fix,
                    t_EStor = t_EStor,
                    t_Pres = t_Pres,
                    t_Fix_mx = t_FixEStor_mx[,1],
                    t_EStor_mx = t_FixEStor_mx[,2],
                    t_Pres_mx = t_Pres_mx,
                    model = rep(model, length(obs_d13C)),
                    naive = rep(naive, length(obs_d13C)),
                    pred.point = pred.point,
                    pred.mean = pred.mean,
                    pred.median = pred.median,
                    pred.CI = pred.CI,
                    change.point = pred.point - obs_d13C,
                    change.mean = pred.mean - obs_d13C,
                    change.median = pred.median - obs_d13C
  )
  return(Out)
  
}
################################################################################
################################################################################
################################################################################
