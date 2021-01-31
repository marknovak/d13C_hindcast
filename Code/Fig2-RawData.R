####################################
# Visualize the raw, un-modeled data
####################################
rm(list = ls())
######################
# Import prepared data
######################
load(file = '../Results/ProcessedData.Rdata')

#######
# Setup
#######
dat <-
  dat[order(dat$EStor), ] # order by EStor to avoid need to do so later

pdf('../Results/Figures/Fig2-RawData.pdf', height = 5, width = 6)

p<-layout(matrix(c(1, 2, 3, 4, 6, 5),
            ncol = 2, byrow = FALSE),
            heights = c(1, 0.2, 1))
# layout.show(p)
dmar <- c(4, 4, 2, 1)
par(
  cex = 0.8,
  las = 1,
  mgp = c(1.8, 0.3, 0),
  tcl = -0.2,
  cex.lab = 1.1,
  cex.axis = 0.9,
  mar = dmar
)

tdat <- subset(dat, Source == 'Field')
specimens <- unique(tdat$SpecimenID)

tdat$jit.EStor <- jitter(tdat$EStor, 0.3)
tdat$cols <- gray.colors(max(tdat$Fix) + 1)[1 + tdat$Fix]

leg.txt <- sort(unique(tdat$Fix))
leg.col <- gray.colors(max(tdat$Fix) + 1)[leg.txt + 1]


plot(
  tdat$jit.EStor,
  # setup empty plot
  tdat$d13C,
  type = 'n',
  xlab = 'Time in ethanol (days)',
  ylab = expression(paste(delta ^ 13, 'C')),
  main = 'Field-processed'
)
for (i in 1:length(specimens)) {
  # plot connecting lines
  sdat <- subset(tdat, SpecimenID == specimens[i])
  points(sdat$jit.EStor[!is.na(sdat$d13C)],
         sdat$d13C[!is.na(sdat$d13C)],
         type = 'l',
         col = 'grey')
}
points(tdat$jit.EStor, # add points over lines
       tdat$d13C,
       pch = 21,
       bg = tdat$cols)
mtext(
  'A',
  side = 3,
  at = 2150,
  line = -1.5,
  cex = 1.2
)

# ~~~ Legend
par(mar = c(0, 2.5, 0, 0))
plot(1,
     1,
     type = 'n',
     ann = FALSE,
     axes = FALSE)
legend(
  'bottom',
  legend = leg.txt,
  pch = 21,
  pt.bg = leg.col,
  horiz = TRUE,
  title = 'Fixation time (days)',
  inset = 0,
  xjust = 0,
  x.intersp = 0.5,
  text.width = 0.03,
  cex=0.9
)

# ~~~ 15N
par(mar = dmar)
plot(
  tdat$jit.EStor,
  # setup empty plot
  tdat$d15N,
  type = 'n',
  xlab = 'Time in ethanol (days)',
  ylab = expression(paste(delta ^ 15, 'N'))
)
for (i in 1:length(specimens)) {
  # plot connecting lines
  sdat <- subset(tdat, SpecimenID == specimens[i])
  points(sdat$jit.EStor[!is.na(sdat$d15N)],
         sdat$d15N[!is.na(sdat$d15N)],
         type = 'l',
         col = 'grey')
}
points(tdat$jit.EStor, # add points over lines
       tdat$d15N,
       pch = 21,
       bg = tdat$cols)
mtext(
  'B',
  side = 3,
  at = 2150,
  line = -1.5,
  cex = 1.2
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tdat <- subset(dat, Source == 'Experiment')
specimens <- unique(tdat$SpecimenID)

tdat$jit.Fix <- jitter(tdat$Fix, 0.3)
tdat$cols <- gray.colors(max(tdat$Fix) + 1)[1 + tdat$Fix]

leg.txt <- sort(unique(tdat$Fix))
leg.col <- gray.colors(max(tdat$Fix) + 1)[leg.txt + 1]


plot(
  tdat$jit.Fix,
  # setup empty plot
  tdat$d13C,
  type = 'n',
  xlab = 'Fixation time (days)',
  ylab = expression(paste(delta ^ 13, 'C')),
  main = 'Lab-processed'
)

for (i in 1:length(specimens)) {
  # plot connecting lines
  sdat <- subset(tdat, SpecimenID == specimens[i])
  points(sdat$jit.Fix[!is.na(sdat$d13C)],
         sdat$d13C[!is.na(sdat$d13C)],
         type = 'l',
         col = 'grey')
}
points(tdat$jit.Fix, # add points over lines
       tdat$d13C,
       pch = 21,
       bg = tdat$cols)
mtext(
  'C',
  side = 3,
  at = 13.25,
  line = -1.5,
  cex = 1.2
)


# ~~~ 15N
par(mar = dmar)
plot(
  tdat$jit.Fix,
  # setup empty plot
  tdat$d15N,
  type = 'n',
  xlab = 'Fixation time (days)',
  ylab = expression(paste(delta ^ 15, 'N'))
)
for (i in 1:length(specimens)) {
  # plot connecting lines
  sdat <- subset(tdat, SpecimenID == specimens[i])
  points(sdat$jit.Fix[!is.na(sdat$d15N)],
         sdat$d15N[!is.na(sdat$d15N)],
         type = 'l',
         col = 'grey')
}
points(tdat$jit.Fix, # add points over lines
       tdat$d15N,
       pch = 21,
       bg = tdat$cols)
mtext(
  'D',
  side = 3,
  at = 13.25,
  line = -1.5,
  cex = 1.2
)


dev.off()
