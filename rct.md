Fusion
================
2026-06-15

``` r
set.seed(20260607)

if (!requireNamespace("ranger", quietly = TRUE)) {
  install.packages("ranger", repos = "https://cloud.r-project.org")
}
library(ranger)

logistic <- function(u) 1 / (1 + exp(-u))
```

## 1. DGP

``` r
simulate_data_v3 <- function(n_rct = 1000, n_ec = 3000) {
  X_rct <- matrix(rnorm(n_rct * 5), n_rct, 5)
  X_ec  <- matrix(rnorm(n_ec * 5),  n_ec,  5)
  # Selection shift on X1..X3 ONLY (X4, X5 identically distributed in both) (it means in grids we have same bias not idiosyncratic)
  X_ec[, 1] <- X_ec[, 1] + 0.5
  X_ec[, 2] <- X_ec[, 2] + 0.3
  X_ec[, 3] <- X_ec[, 3] - 0.3
  colnames(X_rct) <- colnames(X_ec) <- paste0("X", 1:5)

  m0_fun    <- function(X) 1 + 0.5 * X[, 1] - 0.3 * X[, 2] + 0.4 * X[, 3]
  tau_fun   <- function(X) 3 + 1.0 * X[, 2]                # heterogeneous CATE
  idx_fun   <- function(X) X[, 4] + 0.6 * X[, 5]           # bias index (off-selection)
  delta_fun <- function(X) -2 * logistic(4 * (idx_fun(X) - 1))

  A_rct <- rbinom(n_rct, 1, 0.5)
  Y_rct <- m0_fun(X_rct) + tau_fun(X_rct) * A_rct + rnorm(n_rct, sd=0.3)
  Y_ec  <- m0_fun(X_ec) + delta_fun(X_ec) + rnorm(n_ec, sd=0.3)

  # NCO: unaffected by treatment, but hit by the SAME bias mechanism (scaled).
  NCO_rct <- 0.2 * X_rct[, 2] + rnorm(n_rct, sd = 0.5)
  NCO_ec  <- 0.2 * X_ec[, 2] + 0.5 * logistic(4 * (idx_fun(X_ec) - 1)) +
             rnorm(n_ec, sd = 0.5)

  dat <- data.frame(
    Y = c(Y_rct, Y_ec), NCO = c(NCO_rct, NCO_ec),
    A = c(A_rct, rep(0, n_ec)), S = c(rep(1, n_rct), rep(0, n_ec)),
    rbind(X_rct, X_ec)
  )
  attr(dat, "true_ate_S1") <- 3
  attr(dat, "delta_fun")   <- delta_fun
  attr(dat, "idx_fun")     <- idx_fun
  dat
}

dat0 <- simulate_data_v3()
cat("=== Data ===\n")
```

    ## === Data ===

``` r
cat("n_RCT =", sum(dat0$S == 1), " n_EC =", sum(dat0$S == 0),
    " | true ATE in S=1 pop =", attr(dat0, "true_ate_S1"), "\n\n")
```

    ## n_RCT = 1000  n_EC = 3000  | true ATE in S=1 pop = 3

``` r
K <- 5
Xm0    <- as.matrix(dat0[, paste0("X", 1:5)])
folds0 <- sample(rep(1:K, length.out = nrow(dat0)))
```

## 2. Cross-fitted regression helper

``` r
cf_predict <- function(y, X, train_filter, folds,
                       num.trees = 300, min.node.size = 10, weights = NULL) {
  n <- length(y)
  out <- rep(NA_real_, n)
  for (k in sort(unique(folds))) {
    idx <- which(folds != k & train_filter)
    if (length(idx) < 30) next
    if (is.null(weights)) {
      f <- ranger(y = y[idx], x = X[idx, , drop = FALSE],
                  num.trees = num.trees, min.node.size = min.node.size)
    } else {
      f <- ranger(y = y[idx], x = X[idx, , drop = FALSE],
                  case.weights = weights[idx],
                  num.trees = num.trees, min.node.size = min.node.size)
    }
    te <- which(folds == k)
    out[te] <- predict(f, data = X[te, , drop = FALSE])$predictions
  }
  out
}
```

``` r
cf_enrollment <- function(S, X, folds) {
  n <- length(S)
  e <- numeric(n)
  for (k in sort(unique(folds))) {
    tr <- folds != k
    g <- glm(S ~ ., data = data.frame(S = S[tr], X[tr, , drop = FALSE]),
             family = binomial())
    te <- which(folds == k)
    e[te] <- predict(g, newdata = data.frame(X[te, , drop = FALSE]),
                     type = "response")
  }
  e
}

e_pre <- cf_enrollment(dat0$S, Xm0, folds0) #propensity score e(x) = p(s=1|x)
keep  <- (dat0$S == 1) | (dat0$S == 0 & e_pre >= 0.1 & e_pre <= 0.9) #buna gerek yok bence
cat("EC kept after trimming:", sum(keep & dat0$S == 0), "of",
    sum(dat0$S == 0), "\n\n")
```

    ## EC kept after trimming: 2806 of 3000

``` r
dat   <- dat0[keep, ]
Xm    <- Xm0[keep, , drop = FALSE]
n     <- nrow(dat)
folds <- sample(rep(1:K, length.out = n))
piA   <- 0.5
```

## 3. Doubly robust pseudo-outcome

Doubly robust pseudo-outcome for the bias function delta(.) psi_i has
E\[psi \| Z\] = delta(Z) = E\[Y\|Z,EC\] - E\[Y\|Z,RCT ctrl\] (full-Z
nuisances; smoothing axis chosen later, per Lee et al.)

``` r
make_pseudo <- function(yvec, S, A, e_hat, mu_ec, mu_rct_ctrl,
                        piA = 0.5, clip = 0.05) {
  e <- pmin(pmax(e_hat, clip), 1 - clip)
  ec_ctrl  <- (S == 0)                  # all EC are controls
  rct_ctrl <- (S == 1) & (A == 0)
  (mu_ec - mu_rct_ctrl) +
    ifelse(ec_ctrl,  (yvec - mu_ec)       / (1 - e),         0) -
    ifelse(rct_ctrl, (yvec - mu_rct_ctrl) / (e * (1 - piA)), 0)
}

e_hat <- cf_enrollment(dat$S, Xm, folds) #the enrollment (selection) propensity, e(X) = P(S = 1 | X).
mu01  <- cf_predict(dat$Y,   Xm, dat$S == 1 & dat$A == 0, folds)  # E[Y | X]   trained on RCT controls (S=1, A=0), predicted for all
mu00  <- cf_predict(dat$Y,   Xm, dat$S == 0,              folds)  # E[Y | X]   trained on external controls (S=0), predicted for all
mu11  <- cf_predict(dat$Y,   Xm, dat$S == 1 & dat$A == 1, folds)  # [Y | X]   trained on RCT treated   (S=1, A=1), predicted for all
mn1   <- cf_predict(dat$NCO, Xm, dat$S == 1 & dat$A == 0, folds)  # E[NCO | X] trained on RCT controls
mn0   <- cf_predict(dat$NCO, Xm, dat$S == 0,              folds)  # E[NCO | X] trained on external controls

psi_d <- make_pseudo(dat$Y,   dat$S, dat$A, e_hat, mu00, mu01, piA)
psi_n <- make_pseudo(dat$NCO, dat$S, dat$A, e_hat, mn0,  mn1,  piA)
```

## 4. Lee-Okui-Whang style band on a chosen 1-D axis V

- local linear ghat(x) of psi on V
- Lemma-3 linearization: ghat(x)-g(x) ~ (1/n) sum_i l_i(x), l_i(x) =
  K_h(V_i - x) \* U_i / (h \* fhat(x)), U_i = psi_i - ghat(V_i)
- multiplier bootstrap sup for the critical value (their Remark 6)

``` r
band_lee <- function(V, psi, n_grid = 80, h = NULL, qlim = c(0.05, 0.95)) {
  # V    = axis values v(X_i) to smooth along (one covariate or an index)
  # psi  = the DR pseudo-outcomes psi_delta,i from eq.(2), computed upstream
  # n_grid = how many v-points to evaluate the band at
  # h    = bandwidth; NULL means "build it below"
  # qlim = trim the grid to this quantile range of V (interior only)
  
  ok <- is.finite(V) & is.finite(psi)        # TRUE for rows where BOTH V and psi are finite (drop NA/Inf)
  V <- V[ok]; psi <- psi[ok]                 # keep only clean rows in both vectors (they stay aligned)
  n <- length(V)                             # n = number of usable observations
  glim <- as.numeric(quantile(V, qlim))      # grid endpoints = 5th & 95th percentiles of V (trim tails)
  grid <- seq(glim[1], glim[2], length.out = n_grid)  # 80 equally spaced v-values where the band lives
  
  if (is.null(h)) {                          # no bandwidth given -> construct one
    h_pilot <- 1.06 * sd(V) * n^(-1/5)       # Silverman pilot, rate n^(-1/5)  (LOW prefer the RSW plug-in b_hat)
    h <- h_pilot * n^(1/5 - 2/7)             # undersmooth: pilot x n^(1/5-2/7)  =>  final h ~ n^(-2/7)
  }                                          # (so smoothing bias is negligible vs band width)
  
  ghat <- numeric(n_grid)    # storage for ghat(v): local-linear estimate of E[delta|V=v]
  for (g in seq_len(n_grid)) {   # loop over grid points to fill ghat
    d <- V - grid[g]                         # centered distances d_i = V_i - v  (v = this grid point)
    k <- dnorm(d / h)  # Gaussian kernel weights K((V_i - v)/h)
    Sw <- sum(k); Sx <- sum(k * d); Sxx <- sum(k * d^2)   # weighted moments: Sum k, Sum k*d, Sum k*d^2
    Sy <- sum(k * psi); Sxy <- sum(k * d * psi)           # weighted cross-moments: Sum k*psi, Sum k*d*psi
    den <- Sw * Sxx - Sx^2                   # determinant of the 2x2 weighted Gram matrix (>= 0 always)
    ghat[g] <- if (den > 1e-10) (Sxx * Sy - Sx * Sxy) / den   # local-linear INTERCEPT at d=0 (Cramer's rule) = ghat(v) burada cramer's rule var
    else Sy / max(Sw, 1e-10)      # fallback to local-constant mean Sum(k*psi)/Sum(k) if window degenerate
  }                                          # end ghat loop
  
  g_at_obs <- approx(grid, ghat,  # interpolate ghat from the grid back onto every observation's V_i
                     xout = pmin(pmax(V, glim[1]), glim[2]), rule = 2)$y  # clamp V_i into [glim]; rule=2 = flat extrapolation
  U <- psi - g_at_obs  # residuals U_i = psi_delta,i - ghat(V_i)
  
  #the following is creating uncertainty around the band
  Lmat <- matrix(0, n, n_grid)               # storage: n x n_grid matrix of influence terms l_i(v)
  sdv  <- numeric(n_grid)                    # storage: pointwise SD s_n(v) at each grid point
  for (g in seq_len(n_grid)) {               # second loop: build influence terms and their SD
    k <- dnorm((V - grid[g]) / h)            # kernel weights at this grid point (same h as ghat)
    fhat <- max(mean(k) / h, 1e-10)          # kernel density estimate f_V(v) = (1/nh) Sum k, floored > 0
    li <- k * U / (h * fhat)                 # influence term l_i(v) = K((V_i-v)/h) * U_i / (h * f_V(v))
    li <- li - mean(li)                      # center l_i to mean 0 (stabilizes bootstrap; extra step vs paper)
    Lmat[, g] <- li                          # store this column of influence terms
    sdv[g] <- sqrt(mean(li^2) / n)           # s_n(v) = sqrt( (1/n) * mean(l^2) ) = estimated SD of ghat(v)
  }                                          # end variance loop 
  
  sdv <- pmax(sdv, 1e-10)                     # floor the SDs (avoid divide-by-zero when studentizing later)
  list(grid = grid, ghat = ghat, sd = sdv,    # return ingredients: grid, ghat(v), s_n(v),
       L = Lmat, h = h, n = n, keep = ok)      # influence matrix L, bandwidth, n, and the finite-row mask
}  

V_e   <- e_hat
V_idx <- Xm[, 4] + 0.6 * Xm[, 5]

bd_e   <- band_lee(V_e,   psi_d)
bd_X   <- lapply(1:5, function(j) band_lee(Xm[, j], psi_d))
bd_idx <- band_lee(V_idx, psi_d)
```

NOTE: no bootstrap / critical value / band formed here -\> done
downstream One bootstrap sup shared across a LIST of bands -\>
simultaneous validity across axes. All bands must be built on the same
observations (same order).

``` r
mb_critical <- function(band_list, B = 1000, alpha = 0.05) {
  # band_list = a LIST of band_lee() outputs, one per axis; one shared bootstrap covers them all
  # B         = number of bootstrap draws
  # alpha     = error level -> returns the (1-alpha) simultaneous critical value c_{1-alpha}
  
  n <- band_list[[1]]$n                       # sample size, read from the first axis's band object
  stopifnot(all(vapply(band_list, function(b) b$n, 1L) == n))   # require EVERY axis to share the same n
  Lall  <- do.call(cbind, lapply(band_list, function(b) b$L))   # cbind all influence matrices -> n x (ALL grid points)
  sdall <- do.call(c,     lapply(band_list, function(b) b$sd))  # concatenate all SD vectors, aligned to Lall's columns
  M <- numeric(B)                             # storage for the B bootstrap sup-statistics
  for (b in seq_len(B)) {                     # bootstrap loop
    xi <- rnorm(n)                            # draw ONE N(0,1) per observation, shared across every column
    M[b] <- max(abs(colSums(xi * Lall)) / (n * sdall))   # sup over all grid pts of |(1/n) Sum_i xi_i*l_i(v)| / s_n(v)
  }                                           #   (xi recycles down each column; colSums = Sum_i xi_i*l_i per grid point)
  as.numeric(quantile(M, 1 - alpha))          # c_{1-alpha} = (1-alpha) empirical quantile of those maxima
}

c_e   <- mb_critical(list(bd_e),  B = 1000)
c_X   <- mb_critical(bd_X,        B = 1000)   # simultaneous across X1..X5
c_idx <- mb_critical(list(bd_idx), B = 1000)
```

Equivalence-based classification: the POSITIVE transportable region is
where the whole band sits inside \[-eps, +eps\]. “Band covers 0” is only
absence of evidence and is reported separately.

``` r
band_regions <- function(bd, ccrit, eps) {
  # bd    = one band_lee() output (a single axis)
  # ccrit = the simultaneous critical value c_{1-alpha} from mb_critical()
  # eps   = tolerance on the outcome scale, defining the equivalence interval [-eps, eps]
  
  Lb <- bd$ghat - ccrit * bd$sd               # lower band edge:  ghat(v) - c_{1-alpha} * s_n(v)
  Ub <- bd$ghat + ccrit * bd$sd               # upper band edge:  ghat(v) + c_{1-alpha} * s_n(v)
  list(lower = Lb, upper = Ub,                # return both edges, plus two logical maps over bd$grid:
       excludes0      = (Lb > 0) | (Ub < 0),          # FALSIFY: band wholly above 0 OR wholly below 0 = detectable bias
       transportable  = (Lb >= -eps) & (Ub <= eps))   # CERTIFY (eq.3 set T_eps): the whole band lies within [-eps, eps]
}

eps_tol <- 0.3   # tolerance defining "negligible bias"; set substantively

r_e   <- band_regions(bd_e,   c_e,   eps_tol)
r_X4  <- band_regions(bd_X[[4]], c_X, eps_tol)
r_idx <- band_regions(bd_idx, c_idx, eps_tol)


bd_nco <- band_lee(Xm[, 4], psi_n)              # NCO band joins the family

family  <- c(list(bd_e), bd_X, list(bd_idx), list(bd_nco)) #burada neden join ettik

c_joint <- mb_critical(family, B = 1000)
cat(sprintf("Joint 95%% critical value over %d axes (incl. NCO): c = %.2f\n",
            length(family), c_joint))
```

    ## Joint 95% critical value over 8 axes (incl. NCO): c = 3.70

``` r
r_gate <- function(eps) band_regions(bd_idx, c_joint, eps)   # regions on the GATE axis, because selection happens on bd_idx, this is actually providing conservative region
r_nco  <- band_regions(bd_nco, c_joint, eps_tol)

excl_any <- mean(do.call(c, lapply(family, function(b)
  band_regions(b, c_joint, eps_tol)$excludes0)))
cat(sprintf("Falsification: bands exclude 0 on %.0f%% of the pooled grid; NCO on %.0f%%.\n\n",
            100 * excl_any, 100 * mean(r_nco$excludes0)))
```

    ## Falsification: bands exclude 0 on 65% of the pooled grid; NCO on 9%.

## 6. From a certified REGION (on the grid) to a certified UNIT

Conservative lookup: a unit is certified only if it lies INSIDE the grid
range and BOTH grid points flanking it are certified. Outside the grid
(the trimmed tails) there is no certificate -\> never borrowed.

``` r
certified_units <- function(V, bd, region) {
  ng   <- length(bd$grid) #mapping
  out  <- rep(FALSE, length(V)) # number of grid points
  inb  <- is.finite(V) & V >= bd$grid[1] & V <= bd$grid[ng] # in bounds; true only for units whose axis value sits inside the grid range
  j    <- pmin(pmax(findInterval(V, bd$grid), 1L), ng - 1L) # eger grid j ve j+1 arasina duserse j ye assign ediliyor
  out[inb] <- region$transportable[j[inb]] & region$transportable[j[inb] + 1L]
  out
}
```

## 7. The fused estimator

``` r
#(same as before; included so this file self-contains)

fused_ate_v3 <- function(dat, e_hat, pi_borrow, mu11, mu0_fused,
                         piA = 0.5, clip = 0.02) {
  S <- dat$S; A <- dat$A; Y <- dat$Y
  n  <- length(Y); P1 <- mean(S)
  e  <- pmin(pmax(e_hat, clip), 1 - clip)
  w   <- e / (e * (1 - piA) + pi_borrow * (1 - e))
  rho <- ifelse(S == 1 & A == 0, 1, ifelse(S == 0, pi_borrow, 0))
  plug <- mean((S / P1) * (mu11 - mu0_fused))
  ctrt <- mean(ifelse(S == 1 & A == 1, (Y - mu11) / (P1 * piA), 0))
  cctl <- -mean((w * rho / P1) * ifelse(A == 0, Y - mu0_fused, 0))
  psi  <- plug + ctrt + cctl
  phi <- (S / P1) * (mu11 - mu0_fused - psi) +
    ifelse(S == 1 & A == 1, (Y - mu11) / (P1 * piA), 0) -
    (w * rho / P1) * ifelse(A == 0, Y - mu0_fused, 0)
  list(psi = psi, se = sqrt(var(phi) / n),
       m = w * pi_borrow * (1 - e))      # m = borrowed mass per row (for budget)
}
```

## 8. Guarantee machinery

1)  bias-aware critical value: the 1-alpha quantile of \|N(t,1)\|, t =
    B/se
2)  lambda-layer: orthogonal pseudo-outcome whose mean/P1 estimates the
    bias ACTUALLY injected (eq. 20), assumption-free

``` r
cv_bias_aware <- function(t, alpha = 0.05) {
  t <- abs(t)
  f <- function(c) pnorm(c - t) - pnorm(-c - t) - (1 - alpha)
  uniroot(f, c(1e-6, t + qnorm(1 - alpha) + 8))$root} #this finds critical value

m_fun <- function(e, pi_x, piA) {
  D <- e * (1 - piA) + pi_x * (1 - e)
  ifelse(pi_x > 0, pi_x * e * (1 - e) / D, 0)}

mprime_fun <- function(e, pi_x, piA) {
  D <- e * (1 - piA) + pi_x * (1 - e); N <- pi_x * e * (1 - e)
  ifelse(pi_x > 0, (pi_x * (1 - 2 * e) * D - N * ((1 - piA) - pi_x)) / D^2, 0)}
borrowed_bias_ci <- function(psi_d, delta_hat, S, e_hat, pi_x,
                             piA = 0.5, clip = 0.02, alpha = 0.05) {
  e  <- pmin(pmax(e_hat, clip), 1 - clip)
  psi_lam <- m_fun(e, pi_x, piA) * psi_d +
    delta_hat * mprime_fun(e, pi_x, piA) * (S - e)
  P1 <- mean(S); B <- mean(psi_lam) / P1
  se <- sd((psi_lam - B * S) / P1) / sqrt(length(S))
  z  <- qnorm(1 - alpha / 2)
  list(B = B, se = se, lo = B - z * se, hi = B + z * se) #burada yeni bir bound olustu idea is from https://rdrr.io/github/kolesarm/RDHonest/man/CVb.html
}
EPS_GATE <- eps_tol                      # the tolerance the decision is made at
EPS_GRID <- sort(unique(c(eps_tol, 0.5, 0.75, 1.0))) #decision bunlarla nasil degisirdi
ec <- dat$S == 0

delta_Z <- cf_predict(psi_d, Xm, rep(TRUE, n), folds,
                      num.trees = 400, min.node.size = 20)
res_rct <- fused_ate_v3(dat, e_hat, rep(0, n), mu11, mu01, piA)
```

## 9. Band-gated borrowing across a tolerance grid, with guarantees

``` r
cat(sprintf("RCT-only AIPW benchmark: %.3f (SE %.3f)\n\n", res_rct$psi, res_rct$se))
```

    ## RCT-only AIPW benchmark: 2.999 (SE 0.038)

``` r
cat("=== Band-gated borrowing: pi(x) = 1{v(x) in certified region} ===\n")
```

    ## === Band-gated borrowing: pi(x) = 1{v(x) in certified region} ===

``` r
cat(sprintf("%-5s %-10s %-14s %-16s %-8s %-20s\n",
            "eps", "%grid cert", "EC borrowed", "psi (SE)", "budget", "bias-aware 95% CI"))
```

    ## eps   %grid cert EC borrowed    psi (SE)         budget   bias-aware 95% CI

``` r
store <- list()
for (eps in EPS_GRID) {
  reg  <- r_gate(eps)
  pig  <- as.numeric(certified_units(V_idx, bd_idx, reg))     # the DECISION
  if (sum(pig[ec]) > 0) {
    rho <- ifelse(dat$S == 1 & dat$A == 0, 1, ifelse(dat$S == 0, pig, 0))
    mu0 <- cf_predict(dat$Y, Xm, dat$A == 0 & rho > 0, folds, weights = rho)
  } else mu0 <- mu01                                          # nothing borrowed
  res <- fused_ate_v3(dat, e_hat, pig, mu11, mu0, piA)
  budget <- eps * mean(res$m) / mean(dat$S)                   # (G2) bound
  hw  <- cv_bias_aware(budget / res$se) * res$se
  cat(sprintf("%-5.2f %-10s %-14s %-16s %-8.3f [%.3f, %.3f]\n",
              eps, sprintf("%.0f%%", 100 * mean(reg$transportable)),
              sprintf("%d (%.0f%%)", sum(pig[ec]), 100 * mean(pig[ec])),
              sprintf("%.3f (%.3f)", res$psi, res$se),
              budget, res$psi - hw, res$psi + hw))
  store[[as.character(eps)]] <- list(pi = pig, res = res, budget = budget, reg = reg)
}
```

    ## 0.30  45%        1258 (45%)     3.001 (0.037)    0.110    [2.831, 3.171]
    ## 0.50  62%        1665 (59%)     3.002 (0.036)    0.245    [2.698, 3.306]
    ## 0.75  69%        1883 (67%)     3.026 (0.036)    0.417    [2.550, 3.501]
    ## 1.00  72%        1982 (71%)     3.044 (0.035)    0.585    [2.400, 3.687]

``` r
pick <- store[[as.character(EPS_GATE)]]
bb <- borrowed_bias_ci(psi_d, delta_Z, dat$S, e_hat, pick$pi, piA)
cat(sprintf("(G3) assumption-free check -- bias actually injected: B-hat = %.3f, 95%% CI [%.3f, %.3f]\n",
            bb$B, bb$lo, bb$hi))
```

    ## (G3) assumption-free check -- bias actually injected: B-hat = -0.013, 95% CI [-0.032, 0.006]

``` r
if (!is.null(attr(dat0, "delta_fun"))) {
  cat(sprintf("[oracle, sim only] realized injected bias = %.3f (should sit in the CI)\n",
              mean(pick$res$m * attr(dat0, "delta_fun")(Xm)) / mean(dat$S)))
}
```

    ## [oracle, sim only] realized injected bias = -0.007 (should sit in the CI)

``` r
medeps <- median(abs(bd_idx$ghat) + (r_gate(1)$upper - bd_idx$ghat))
cat(sprintf("Power note: median minimum-certifiable eps on this axis = %.2f;\n", medeps))
```

    ## Power note: median minimum-certifiable eps on this axis = 0.33;

## 10. Picture: the decision drawn ON the band

``` r
png("band_gated_borrowing.png", width = 1100, height = 650)
reg <- pick$reg
plot(bd_idx$grid, bd_idx$ghat, type = "l", lwd = 2,
     ylim = range(c(reg$lower, reg$upper, EPS_GATE, -EPS_GATE)),
     xlab = "v(X) = X4 + 0.6 X5", ylab = expression(bar(delta)(v)),
     main = sprintf("Borrowing decided by the band (eps = %.2f, joint c = %.2f)",
                    EPS_GATE, c_joint))
polygon(c(bd_idx$grid, rev(bd_idx$grid)), c(reg$lower, rev(reg$upper)),
        col = rgb(.6, .6, 1, .35), border = NA)
lines(bd_idx$grid, bd_idx$ghat, lwd = 2)
abline(h = c(-EPS_GATE, 0, EPS_GATE), lty = c(3, 2, 3),
       col = c("darkgreen", "black", "darkgreen"))
if (any(reg$transportable))
  rug(bd_idx$grid[reg$transportable], col = "darkgreen", lwd = 2)
if (sum(pick$pi[ec]) > 0)
  rug(V_idx[ec][pick$pi[ec] == 1], side = 3, col = "blue", lwd = 1)
legend("bottomleft",
       c("band (joint c)", "tolerance +/- eps", "certified region (gate)",
         "borrowed EC units (top rug)"),
       lty = c(1, 3, 1, 1), lwd = c(2, 1, 2, 1),
       col = c("black", "darkgreen", "darkgreen", "blue"), bty = "n", cex = .85)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
cat("\nSaved plot to band_gated_borrowing.png\n")
```

    ## 
    ## Saved plot to band_gated_borrowing.png
