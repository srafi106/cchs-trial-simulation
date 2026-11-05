# Functions for Part B: Bayesian analysis of a two-arm trial with a binary
# endpoint, using a robust mixture prior for the control event probability.

# Control prior: w0 * Beta(a1, b1) + (1 - w0) * Beta(1, 1)
# Treatment prior: Beta(1, 1)
# Returns P(p_treatment < p_control | data) for each simulated trial, the
# posterior weight on the informative component, and the posterior mean
# risk difference (treatment minus control).
post_success = function(yc, nc, yt, nt, w0, a1, b1) {
  R = length(yc)
  nc = rep_len(nc, R)
  nt = rep_len(nt, R)
  l1 = log(w0) + lbeta(a1 + yc, b1 + nc - yc) - lbeta(a1, b1)
  l2 = log(1 - w0) + lbeta(1 + yc, 1 + nc - yc) - lbeta(1, 1)
  wp = 1 / (1 + exp(l2 - l1))
  g = seq(0.0002, 0.5998, length.out = 1500)
  dg = g[2] - g[1]
  ps = numeric(R)
  for (s in seq(1, R, by = 1000)) {
    ii = s:min(R, s + 999)
    Rk = length(ii)
    G = rep(g, each = Rk)
    m = length(g)
    Ft = matrix(pbeta(G, rep(1 + yt[ii], m), rep(1 + nt[ii] - yt[ii], m)), nrow = Rk)
    f1 = matrix(dbeta(G, rep(a1 + yc[ii], m), rep(b1 + nc[ii] - yc[ii], m)), nrow = Rk)
    f2 = matrix(dbeta(G, rep(1 + yc[ii], m), rep(1 + nc[ii] - yc[ii], m)), nrow = Rk)
    ps[ii] = rowSums((wp[ii] * f1 + (1 - wp[ii]) * f2) * Ft) * dg
  }
  mean_c = wp * (a1 + yc) / (a1 + b1 + nc) + (1 - wp) * (1 + yc) / (2 + nc)
  mean_t = (1 + yt) / (2 + nt)
  list(p_success = ps, wpost = wp, rd = mean_t - mean_c)
}

# One-sided pooled z test for two proportions
freq_success = function(yc, nc, yt, nt, alpha_one) {
  pc = yc / nc
  pt = yt / nt
  pb = (yc + yt) / (nc + nt)
  z = (pc - pt) / sqrt(pb * (1 - pb) * (1 / nc + 1 / nt))
  z > qnorm(1 - alpha_one)
}
