# Functions for Part A: virtual participants, randomization schemes, analysis.

prepare_population = function(pop) {
  pop$prov_id = as.integer(factor(pop$prov))
  pop$age65 = as.integer(pop$age_grp == "A65plus")
  pop$female = as.integer(pop$sex == "Female")
  pop$inc = as.integer(sub("Q", "", pop$inc_q))
  list(pop = pop, np = max(pop$prov_id), cell_p = pop$N / sum(pop$N))
}

# Draw n participants with replacement from the cells, proportional to cell size
draw_participants = function(P, n) {
  pop = P$pop
  i = sample.int(nrow(pop), n, replace = TRUE, prob = P$cell_p)
  list(prov = pop$prov_id[i], age65 = pop$age65[i], female = pop$female[i],
       inc = pop$inc[i], p0 = pop$p_prev[i])
}

dummies = function(x, lv) outer(x, lv, "==") * 1

# Binary baseline indicators used to measure balance
make_B = function(p, np) {
  cbind(p$female, p$age65, dummies(p$inc, 1:5), dummies(p$prov, 1:np))
}

# Covariates for the adjusted analysis
make_Xad = function(p, np) {
  cbind(p$age65, p$female, dummies(p$inc, 2:5), dummies(p$prov, 2:np))
}

# Level indices for minimization, unique across factors
make_F = function(p, np) {
  cbind(p$prov, np + 1 + p$age65, np + 3 + p$female, np + 4 + p$inc)
}

rand_simple = function(n) rbinom(n, 1, 0.5)

# Permuted blocks within strata
rand_pbr = function(strata, block = 4) {
  arm = integer(length(strata))
  for (s in unique(strata)) {
    idx = which(strata == s)
    ns = length(idx)
    nb = ceiling(ns / block)
    sq = as.vector(replicate(nb, sample(rep(0:1, block / 2))))
    arm[idx] = sq[seq_len(ns)]
  }
  arm
}

# Pocock-Simon minimization on marginal imbalance, biased coin with prob p_bc
rand_min = function(Fm, p_bc = 0.8) {
  n = nrow(Fm)
  arm = integer(n)
  cnt = matrix(0L, nrow = max(Fm), ncol = 2)
  for (i in seq_len(n)) {
    idx = Fm[i, ]
    d0 = sum(abs((cnt[idx, 1] + 1L) - cnt[idx, 2]))
    d1 = sum(abs(cnt[idx, 1] - (cnt[idx, 2] + 1L)))
    if (d0 == d1) {
      a = rbinom(1, 1, 0.5)
    } else {
      best = if (d0 < d1) 0L else 1L
      a = if (runif(1) < p_bc) best else 1L - best
    }
    arm[i] = a
    cnt[cbind(idx, a + 1L)] = cnt[cbind(idx, a + 1L)] + 1L
  }
  arm
}

schemes = c("Simple", "Province blocks", "Province x age x sex blocks",
            "Minimization")

assign_arm = function(s, p, Fm, block_size = 4, p_bc = 0.8) {
  switch(s,
    "Simple" = rand_simple(length(p$prov)),
    "Province blocks" = rand_pbr(p$prov, block_size),
    "Province x age x sex blocks" =
      rand_pbr(p$prov * 4 + p$age65 * 2 + p$female, block_size),
    "Minimization" = rand_min(Fm, p_bc))
}

# Largest absolute standardized mean difference over the columns of B
max_smd = function(B, arm) {
  B1 = B[arm == 1, , drop = FALSE]
  B0 = B[arm == 0, , drop = FALSE]
  s = sqrt((apply(B1, 2, var) + apply(B0, 2, var)) / 2)
  ok = s > 0
  max((abs(colMeans(B1) - colMeans(B0)))[ok] / s[ok])
}

# Logistic regression for the treatment effect; returns estimate and standard error
fit_arm = function(y, arm, Xad = NULL) {
  X = cbind(1, arm)
  if (!is.null(Xad)) {
    keep = colSums(abs(Xad)) > 0
    X = cbind(X, Xad[, keep, drop = FALSE])
  }
  f = suppressWarnings(glm.fit(X, y, family = binomial()))
  if (!isTRUE(f$converged) || f$rank < ncol(X)) return(c(NA_real_, NA_real_))
  V = chol2inv(qr.R(f$qr))
  c(unname(f$coefficients[2]), sqrt(V[2, 2]))
}
