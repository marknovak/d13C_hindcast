rm(list = ls())

mod3 <- function(x) {
  Fix_dur = pmin(x, tau)
  Stor = ifelse(x <= tau, 0, pmax(0, x - tau))
  I + A * (1 - exp(l1 * Fix_dur) * exp(l2 * Stor)) # equivalent to exp( l1*Fix_dur + l2*Stor )
}

mod3bar <- function(x) {
  Fix_dur = pmin(x, tau)
  Stor = ifelse(x <= tau, 0, pmax(0, x - tau))
  out = I + A * (1 - exp(l1 * Fix_dur) * exp(l2 * Stor)) # equivalent to exp( l1*Fix_dur + l2*Stor )
  out = pmax(out, Abar)
}


pdf('../Results/Figures/Fig1-Conceptual.pdf',
    height = 5,
    width = 3)
par(
  mfrow = c(3, 1),
  mar = c(2, 4.5, 1, 1),
  cex = 0.7,
  tcl = -0.3,
  mgp = c(1, 0.5, 0)
)

I  = -20
A  = -1
l1 = -0.008
l2 = -0.002
tau = 150
tmax = 700

# ~~~~~~~~~~~~

tend = 700
ybuf = 0.2
yadj = 0
Abar = I + A + yadj

ylims = c(I + A - ybuf + yadj, I + ybuf)

plot(
  1,
  1,
  type = 'n',
  xlim = c(0, tmax),
  ylim = ylims,
  axes = F,
  ann = F
)
title(ylab = expression(paste(delta ^ 13, C)), mgp = c(0.5, 0, 0))

axis(1, at = c(0, tau, tmax), labels = c(NA, NA, NA))
mtext(
  'Fixation',
  side = 1,
  line = 0.3,
  at = tau / 2,
  cex = 0.7
)
mtext(
  'Storage time',
  side = 1,
  line = 0.3,
  at = tau + (tmax - tau) / 2,
  cex = 0.7
)

axis(
  2,
  at = c(I, I + A + yadj),
  labels = c('Initial', 'Asymptote'),
  las = 2
)
abline(h = I, lty = 3)
abline(h = I + A + yadj, lty = 3)

curve(mod3,
      0,
      tend,
      add = T,
      lwd = 3,
      lty = 3)
curve(mod3bar,
      0,
      tend,
      add = T,
      lwd = 3,
      lty = 1)

points(
  c(0, tend),
  mod3bar(c(0, tend)),
  pch = 21,
  bg = c('white', 'black'),
  cex = 2,
  lwd = 2
)

mtext(
  'A',
  side = 3,
  at = tmax,
  line = -1,
  cex = 1.2
)

box(lwd = 2, bty = 'l')
# ~~~~~~~~~~~~

tend = 700
ybuf = 0.2
yadj = 0.2
Abar = I + A + yadj

ylims = c(I + A - ybuf + yadj, I + ybuf)

plot(
  1,
  1,
  type = 'n',
  xlim = c(0, tmax),
  ylim = ylims,
  axes = F,
  ann = F
)
title(ylab = expression(paste(delta ^ 13, C)), mgp = c(0.5, 0, 0))

axis(1, at = c(0, tau, tmax), labels = c(NA, NA, NA))
mtext(
  'Fixation',
  side = 1,
  line = 0.3,
  at = tau / 2,
  cex = 0.7
)
mtext(
  'Storage time',
  side = 1,
  line = 0.3,
  at = tau + (tmax - tau) / 2,
  cex = 0.7
)

axis(
  2,
  at = c(I, I + A + yadj),
  labels = c('Initial', 'Asymptote'),
  las = 2
)
abline(h = I, lty = 3)
abline(h = I + A + yadj, lty = 3)

curve(
  mod3,
  0,
  tend,
  add = T,
  lwd = 3,
  lty = 1,
  col = 'grey'
)
curve(mod3,
      0,
      tend,
      add = T,
      lwd = 2,
      lty = 2)
curve(mod3bar,
      0,
      tend,
      add = T,
      lwd = 3,
      lty = 1)

points(
  c(0, tend),
  mod3bar(c(0, tend)),
  pch = 21,
  bg = c('white', 'black'),
  cex = 2,
  lwd = 2
)

mtext(
  'B',
  side = 3,
  at = tmax,
  line = -1,
  cex = 1.2
)

box(lwd = 2, bty = 'l')
# ~~~~~~~~~~~~

tend = 700
ybuf = 0.2
yadj = 0.45
Abar = I + A + yadj

ylims = c(I + A - ybuf + yadj, I + ybuf)

plot(
  1,
  1,
  type = 'n',
  xlim = c(0, tmax),
  ylim = ylims,
  axes = F,
  ann = F
)
title(ylab = expression(paste(delta ^ 13, C)), mgp = c(0.5, 0, 0))

axis(1, at = c(0, tau, tmax), labels = c(NA, NA, NA))
mtext(
  'Fixation',
  side = 1,
  line = 0.3,
  at = tau / 2,
  cex = 0.7
)
mtext(
  'Storage time',
  side = 1,
  line = 0.3,
  at = tau + (tmax - tau) / 2,
  cex = 0.7
)

axis(
  2,
  at = c(I, I + A + yadj),
  labels = c('Initial', 'Asymptote'),
  las = 2
)
abline(h = I, lty = 3)
abline(h = I + A + yadj, lty = 3)

curve(
  mod3,
  0,
  tend,
  add = T,
  lwd = 3,
  lty = 1,
  col = 'grey'
)
curve(mod3,
      0,
      tend,
      add = T,
      lwd = 2,
      lty = 2)
curve(mod3bar,
      0,
      tend,
      add = T,
      lwd = 3,
      lty = 1)

points(
  c(0, tend),
  mod3bar(c(0, tend)),
  pch = 21,
  bg = c('white', 'black'),
  cex = 2,
  lwd = 2
)

mtext(
  'C',
  side = 3,
  at = tmax,
  line = -1,
  cex = 1.2
)

box(lwd = 2, bty = 'l')

dev.off()