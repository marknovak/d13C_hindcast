#######################################################
# Script to pre-process the data prior to model-fitting
#######################################################
rm(list = ls())

odat <-
  read.csv(
    '../Data/Taylor_FormalinFixation.csv',
    header = T,
    stringsAsFactors = FALSE
  )

# Convert storage years to days
odat$Storage_days <- 365 * odat$Storage_years

# Combine to total preservation
odat$Preservation_days <- odat$Fixation_days + odat$Storage_days


# Range of variation for in-house standards
range(c(odat$d13C_sd_Run1, odat$d13C_sd_Run2), na.rm=TRUE)
range(c(odat$d15N_sd_Run1, odat$d15N_sd_Run2), na.rm=TRUE)


# Due to an isotope analyzer failure, many of the 'experimental' (i.e. 'lab-processed')
# specimen sample had to be run a second time.   Only the d15N readings were affected,
# and for some specimens there was not enough material remaining to permit a second run.
# Therefore, for the 'experimental' specimens, we:
# (i) take the average d13C of both Run 1 and 2, using only Run 1 when there is no Run 2
# (ii) use Run 2 for d15N when both Run 1 and 2 were available
# (iii) use Run 1 for d15N when only Run 1 was available

odat$d13C <- odat$d13C_sd <- odat$d13C_sd_n <- NA
odat$d15N <- odat$d15N_sd <- odat$d15N_sd_n <- NA

# Which have both Run 1 and Run 2 d13C values?  Take their average.
# Those that have only a Run 1 value will use that value (because na.rm=TRUE)
odat$d13C <-
  apply(cbind(odat$d13C_Run1, odat$d13C_Run2), 1, mean, na.rm = TRUE)
odat$d13C_sd <-
  apply(cbind(odat$d13C_sd_Run1, odat$d13C_sd_Run2), 1, mean, na.rm = TRUE)
odat$d13C_sd_n <-
  apply(cbind(odat$d13C_sd_n_Run1, odat$d13C_sd_nRun2),
        1,
        mean,
        na.rm = TRUE)

# Identify 'experimental' specimens
expmt <- odat$Source == 'Experiment'

# For 'field' specimens, use Run 1 d15N
odat$d15N[!expmt] <- odat$d15N_Run1[!expmt]
odat$d15N_sd[!expmt] <- odat$d15N_sd_Run1[!expmt]
odat$d15N_sd_n[!expmt] <- odat$d15N_sd_n_Run1[!expmt]

# For 'experiment' specimens....
# Use Run 2 for d15N when both Run 1 and 2 were available or when only Run 2 was available
sel <- expmt & !is.na(odat$d15N_Run2)
odat$d15N[sel] <- odat$d15N_Run2[sel]
odat$d15N_sd[sel] <- odat$d15N_sd_Run2[sel]
odat$d15N_sd_n[sel] <- odat$d15N_sd_n_Run2[sel]

# Use Run 1 for d15N when only Run 1 was available
sel <- expmt & is.na(odat$d15N_Run2) & !is.na(odat$d15N_Run1)
odat$d15N[sel] <- odat$d15N_Run1[sel]
odat$d15N_sd[sel] <- odat$d15N_sd_Run1[sel]
odat$d15N_sd_n[sel] <- odat$d15N_sd_n_Run1[sel]

# This many specimens did not have a d15N value in either Run 1 or Run 2
sum(is.na(odat$d15N))


# Select columns to keep
dat <- odat[, c(
  'Collector_Number',
  'Source',
  'Fixation_days',
  'Storage_days',
  'Preservation_days',
  'd13C',
  'd13C_sd',
  'd13C_sd_n',
  'd15N',
  'd15N_sd',
  'd15N_sd_n'
)]

# Abbreviate predictors
colnames(dat) <-
  sub('Collector_Number', 'SpecimenID', colnames(dat))
colnames(dat) <- sub('Fixation_days', 'Fix', colnames(dat))
colnames(dat) <- sub('Storage_days', 'EStor', colnames(dat))
colnames(dat) <- sub('Preservation_days', 'Pres', colnames(dat))

# Export processed data for subsequent importing
save(file = '../Results/ProcessedData.Rdata', dat)


# Reviewer request:  What was mean difference between runs?
mean(odat$d13C_Run1 - odat$d13C_Run2, na.rm=TRUE)
sd(odat$d13C_Run1 - odat$d13C_Run2, na.rm=TRUE)
mean(odat$d15N_Run1 - odat$d15N_Run2, na.rm=TRUE)
sd(odat$d15N_Run1 - odat$d15N_Run2, na.rm=TRUE)
#######################################################
#######################################################
#######################################################