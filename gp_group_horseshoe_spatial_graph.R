## ================================================================
## GP-Group Horseshoe Spatial Graph Estimation
## Spatially varying cell-cell interaction networks
##
## Key design:
##   - HSGP with Matern(nu=3/2) kernel for spatial smoothness
##   - Spectral prior on basis weights w_k ~ N(0, S_matern(omega_k))
##   - Group horseshoe delta_{s'} for edge selection (one per neighbor)
##   - Global tau for overall graph sparsity
##
## Model for node s:
##   y_s(l) = sum_{s'} delta_{s'} * phi(l) %*% w_{s'} + eps(l)
##   w_{s',k} ~ N(0, S_matern(omega_k, nu, rho))   [GP smoothness]
##   delta_{s'} ~ Horseshoe(tau)                    [edge sparsity]
##   tau ~ Half-Cauchy(0,1)                         [global sparsity]
##   eps ~ N(0, sigma_sq)
##
## Dependencies: MASS
## ================================================================


rm(list = ls())
library(MASS)

## 2. Build HSGP Basis + Spectral Weights
##    coords  : n x 2 spatial coordinates
##    m       : basis functions per dimension (m^2 total in 2D)
##    nu      : Matern smoothness parameter
##    rho     : length-scale (default: median pairwise distance)
##    alpha   : GP marginal std dev
##    L_factor: domain boundary expansion
##
## Returns:
##   phi    : n x m^2 basis matrix (fixed, deterministic)
##   S_diag : m^2 spectral densities (prior variance for w_k)
##   omega_sq: m^2 eigenfrequencies
## ----------------------------------------------------------------
build_hsgp_basis <- function(coords, m = 5, nu = 1.5,
                             rho = NULL, alpha = 1, L_factor = 1.5) {
  n  <- nrow(coords)
  m2 <- m^2
  
  # Normalize coordinates to [-1, 1]
  coords_norm <- apply(coords, 2, function(x) {
    r <- range(x)
    2 * (x - r[1]) / (r[2] - r[1]) - 1
  })
  
  L <- L_factor  # boundary for normalized domain [-L, L]
  
  # Default rho: median pairwise distance on normalized coords
  if (is.null(rho)) {
    D   <- as.matrix(dist(coords_norm))
    rho <- median(D[upper.tri(D)])
    cat("  Auto length-scale rho =", round(rho, 4), "\n")
  }
  
  # 1D Laplacian eigenfunctions and eigenvalues
  # phi_j(x) = (1/sqrt(L)) * sin(j*pi*(x + L) / (2L))
  # lambda_j = (j*pi / (2L))^2
  phi_1d <- function(x, m) {
    sapply(1:m, function(j) sin(j * pi * (x + L) / (2 * L)) / sqrt(L))
  }
  lambda_1d <- function(m) sapply(1:m, function(j) (j * pi / (2 * L))^2)
  
  phi_x <- phi_1d(coords_norm[, 1], m)   # n x m
  phi_y <- phi_1d(coords_norm[, 2], m)   # n x m
  lam_x <- lambda_1d(m)                  # m eigenvalues (x)
  
  # 2D tensor product: phi_{ij}(l) = phi_i(x) * phi_j(y)
  # omega^2_{ij} = lambda_i + lambda_j  (2D Laplacian eigenvalue)
  phi    <- matrix(0, n, m2)
  S_diag <- rep(0, m2)
  omega_sq_vec <- rep(0, m2)
  
  idx <- 1
  for (i in 1:m) {
    for (j in 1:m) {
      phi[, idx]        <- phi_x[, i] * phi_y[, j]
      omega_sq_vec[idx] <- lam_x[i] + lam_y[j]
      S_diag[idx]       <- matern_spectral_density(omega_sq_vec[idx],
                                                   nu = nu, rho = rho,
                                                   alpha = alpha)
      idx <- idx + 1
    }
  }
  
  list(phi = phi, S_diag = S_diag, omega_sq = omega_sq_vec,
       m = m, m2 = m2, rho = rho, nu = nu)
}


## ----------------------------------------------------------------
## 3. Group Horseshoe + GP Spectral Prior Gibbs Sampler
##
## Model (for one node s):
##   y = Z %*% theta + eps,  eps ~ N(0, sigma_sq * I)
##
##   Parameterization:
##     theta_{j,k} = delta_j * w_{j,k}
##     Z[, group_j] = x_j * phi   (n x m2 block)
##
##   Priors:
##     theta_{j,k} ~ N(0, lambda_j^2 * tau^2 * S_diag[k])
##     lambda_j    ~ Half-Cauchy(0,1)   [group horseshoe]
##     tau         ~ Half-Cauchy(0,1)   [global shrinkage]
##     sigma_sq    ~ sparseIG (shape = 0.1, rate = 0.1)
##
##   Full conditionals (CORRECTED):
##     theta       | rest ~ MVN(A^{-1} Z'y / sigma_sq, A^{-1})
##                          A = Z'Z/sigma_sq + D^{-1}
##                          D = blkdiag(lambda_j^2 tau^2 diag(S_diag))
##
##     lambda_j^{-2} | rest ~ Gamma((m2+1)/2,  1/nu_j + ||theta_j||^2_{S^{-1}} / (2*tau^2))
##       [sigma_sq does NOT appear — it is not in the theta prior]
##
##     tau^{-2}    | rest ~ Gamma((q*m2+1)/2,  1/xi + sum_j ||theta_j||^2_{S^{-1}} / (2*lambda_j^2))
##       [sigma_sq does NOT appear — same reason]
##
##     sigma^{-2}  | rest ~ Gamma((n+0.2)/2,   (||y - Z theta||^2 + 0.2) / 2)
##       [only likelihood and hyperprior contribute — theta prior has no sigma_sq]
##
## Arguments:
##   y       : n-vector response (scaled)
##   X_s     : n x q matrix of neighbor expressions (scaled)
##   phi     : n x m2 GP basis matrix
##   S_diag  : m2 spectral densities (GP prior variances)
##   nmc, burn, thin: MCMC settings
## ----------------------------------------------------------------
gp_group_horseshoe_sampler <- function(y, X_s, phi, S_diag,
                                       nmc = 3000, burn = 1000, thin = 5) {
  n  <- length(y)
  q  <- ncol(X_s)       # number of neighbors
  m2 <- ncol(phi)       # basis dimension
  P  <- q * m2          # total coefficients
  
  effective_ss <- floor((nmc - burn) / thin)
  
  ## --- Build blocked design matrix Z (n x P) ---
  ## Z[, ((j-1)*m2+1):(j*m2)] = x_j * phi
  Z <- matrix(0, n, P)
  for (j in 1:q) {
    idx <- ((j - 1) * m2 + 1):(j * m2)
    Z[, idx] <- X_s[, j] * phi
  }
  
  ZtZ <- crossprod(Z)   # P x P  (precompute)
  Zty <- crossprod(Z, y)
  
  ## --- Initialize ---
  theta     <- rep(0, P)
  sigma_sq  <- var(y)
  tau_sq    <- 1
  lambda_sq <- rep(1, q)
  nu_g      <- rep(1, q)
  xi        <- 1
  
  ## --- Storage ---
  theta_store  <- matrix(0, P, effective_ss)
  lambda_store <- matrix(0, q, effective_ss)
  tau_store    <- numeric(effective_ss)
  sigma_store  <- numeric(effective_ss)
  
  iter_saved <- 0
  
  for (iter in 1:nmc) {
    
    ## ----------------------------------------------------------
    ## 1. Update theta | y, lambda, tau, sigma_sq
    ##    Prior precision for theta_{j,k}: 1 / (lambda_sq[j] * tau_sq * S_diag[k])
    ##    D = blkdiag(lambda_j^2 * tau^2 * diag(S_diag))  for j = 1,...,q
    ## ----------------------------------------------------------
    prior_prec <- rep(0, P)
    for (j in 1:q) {
      idx <- ((j - 1) * m2 + 1):(j * m2)
      prior_prec[idx] <- 1 / (lambda_sq[j] * tau_sq * S_diag)
    }
    
    A     <- ZtZ / sigma_sq + diag(prior_prec)
    A_inv <- tryCatch(
      solve(A, tol = 1e-60),
      error = function(e) solve(A + diag(1e-8, P), tol = 1e-60)
    )
    theta_mu <- A_inv %*% Zty / sigma_sq
    theta    <- as.vector(mvrnorm(10, mu = theta_mu, Sigma = A_inv))
    
    ## ----------------------------------------------------------
    ## 2. Update lambda_sq[j] | theta, tau
     ##    FC: lambda_j^{-2} ~ Gamma((m2+1)/2,  1/nu_j + ||theta_j||^2_{S^{-1}} / (2*tau^2))
    ##    where ||theta_j||^2_{S^{-1}} = sum_k theta_{j,k}^2 / S_diag[k]
    ## ----------------------------------------------------------
    for (j in 1:q) {
      idx  <- ((j - 1) * m2 + 1):(j * m2)
      th_j <- theta[idx]
      
      # Auxiliary variable update for Half-Cauchy prior on lambda_j
      nu_g[j]   <- 1 / rgamma(1, 1, rate = 1 + 1 / lambda_sq[j])
      
      # GP-adjusted weighted norm (spectral structure)
      gp_weighted_norm <- sum(th_j^2 / S_diag)
      
      # CORRECTED rate: no sigma_sq
      s_lambda     <- 1 / nu_g[j] + gp_weighted_norm / (2 * tau_sq)
      lambda_sq[j] <- 1 / rgamma(1, (m2 + 1) / 2, rate = s_lambda)
    }
    
    ## ----------------------------------------------------------
    ## 3. Update tau_sq | theta, lambda
    ##    CORRECTED: sigma_sq removed from rate.
    ##    FC: tau^{-2} ~ Gamma((q*m2+1)/2,  1/xi + sum_j ||theta_j||^2_{S^{-1}} / (2*lambda_j^2))
    ## ----------------------------------------------------------
    xi <- 1 / rgamma(1, 1, rate = 1 + 1 / tau_sq)
    
    
    s_tau  <- 1 / xi + global_norm / 2
    tau_sq <- 1 / rgamma(1, (P + 1) / 2, rate = s_tau)
    
    ## ----------------------------------------------------------
    ## 4. Update sigma_sq | y, theta
    ##    FC: sigma^{-2} ~ Gamma((n+0.2)/2,  (||y - Z*theta||^2 + 0.2) / 2)
    ## ----------------------------------------------------------
    resid    <- y - Z %*% theta
    E1       <- sum(resid^2)
    sigma_sq <- 1 / rgamma(1, (n + 0.2) / 2, rate = (E1 + 0.2) / 2)
    
    if (iter %% 500 == 0) cat("    iter:", iter, "/ sigma:", round(sqrt(sigma_sq), 3),
                              "/ tau:", round(sqrt(tau_sq), 4), "\n")
    
    ## --- Store ---
    if (iter > burn && (iter - burn) %% thin == 0) {
      iter_saved <- iter_saved + 1
      theta_store[, iter_saved]  <- theta
      lambda_store[, iter_saved] <- lambda_sq
      tau_store[iter_saved]      <- tau_sq
      sigma_store[iter_saved]    <- sigma_sq
    }
  }
  
  list(
    theta_store  = theta_store[, 1:iter_saved],
    lambda_store = lambda_store[, 1:iter_saved],
    tau_store    = tau_store[1:iter_saved],
    sigma_store  = sigma_store[1:iter_saved]
  )
}

gp_node <- function(y, X_s, coords, rho, alpha = 1,
                           sigma_sq_init = NULL, ci_majority = 0.5) {
  n  <- nrow(coords)
  q  <- ncol(X_s)

  ## Precompute kernel matrix once — shared across neighbors
  D   <- as.matrix(dist(coords))
  K   <- matern32_kernel(D, rho = rho, alpha = alpha)
  K   <- K + diag(1e-6, n)   # numerical jitter

  ## Estimate sigma_sq from residual variance (simple init)
  ## Use OLS residuals as a fast proxy
  if (is.null(sigma_sq_init)) {
    ols_fit  <- tryCatch(
      lm.fit(cbind(1, X_s), y),
      error = function(e) list(residuals = y)
    )
    sigma_sq <- max(var(ols_fit$residuals), 1e-4)
  } else {
    sigma_sq <- sigma_sq_init
  }

  beta_mean   <- matrix(0, n, q)
  beta_sd     <- matrix(0, n, q)
  edge_active <- integer(q)

  for (j in 1:q) {
    xj <- X_s[, j]   # n-vector

    ## Effective noise for this neighbor:
    ## y_j_eff = y / x_j treated as noisy obs of beta_j
    ## This is an approximation — marginalizing over other betas
    ## is intractable; we use the "backfitting" residual instead.
    ## Residual after removing other neighbors (OLS approximation):
    r_j <- tryCatch({
      if (q > 1) {
        Xother <- X_s[, -j, drop = FALSE]
        b_other <- lm.fit(Xother, y)$coefficients
        b_other[is.na(b_other)] <- 0
        y - Xother %*% b_other
      } else {
        y
      }
    }, error = function(e) y)

    ## Locations where xj is non-negligible
    ## GP posterior: beta_j | r_j ~ GP
    ## Observation model: r_j(l) = x_j(l) * beta_j(l) + noise
    ## => diag noise variance: sigma_sq / x_j(l)^2  (heteroscedastic)
    ## To avoid /0, add small floor
    xj_safe <- pmax(abs(xj), 1e-3) * sign(xj + 1e-10)

    ## Noise variance per location (heteroscedastic)
    noise_var <- sigma_sq / xj_safe^2   # n-vector

    ## GP posterior with heteroscedastic noise:
    ## Sigma_noise = diag(noise_var)
    ## K_post_inv  = (K + Sigma_noise)^{-1}
    Sigma_n    <- diag(noise_var)
    K_noise    <- K + Sigma_n

    ## Cholesky solve — O(n^3) — the main cost
    L <- tryCatch(
      chol(K_noise),
      error = function(e) chol(K_noise + diag(1e-4, n))
    )
    ## Posterior mean: K (K + Sigma_n)^{-1} r_j / x_j
    ## Effective observations: z_j = r_j / x_j
    z_j      <- r_j / xj_safe
    alpha_vec <- backsolve(L, forwardsolve(t(L), z_j))
    post_mean <- as.vector(K %*% alpha_vec)

    ## Posterior variance (diagonal): diag(K) - diag(K (K+Sn)^{-1} K)
    ## Efficient: v_k = L^{-1} K[,k], then var_k = K[k,k] - ||v_k||^2
    V <- backsolve(L, t(K))   # n x n
    post_var <- pmax(diag(K) - colSums(V^2), 1e-8)
    post_sd  <- sqrt(post_var)

    beta_mean[, j] <- post_mean
    beta_sd[, j]   <- post_sd

    ## Edge selection: 95% CI excludes 0 at majority of locations
    ci_lo <- post_mean - 1.96 * post_sd
    ci_hi <- post_mean + 1.96 * post_sd
    excludes_zero  <- (ci_lo > 0) | (ci_hi < 0)
    edge_active[j] <- as.integer(mean(excludes_zero) > ci_majority)
  }

  list(
    beta_mean   = beta_mean,    # n x q
    beta_sd     = beta_sd,      # n x q
    edge_active = edge_active,  # q-vector
    sigma_sq    = sigma_sq
  )
}


## ----------------------------------------------------------------
## 4. Recover Spatial Edge Maps from Posterior
##    fit       : output of gp_group_horseshoe_sampler
##    phi       : n x m2 basis matrix
##    S_diag    : m2 spectral densities
##    q         : number of neighbors
##    m2        : basis dimension
##    nbr_names : character vector of neighbor names
##    edge_threshold: kappa cutoff for edge detection
##
## Returns per-node:
##   edge_maps    : n x q matrix, posterior mean spatial edge strength
##   edge_maps_sd : n x q matrix, posterior sd of edge strength
##   kappa_mean   : q-vector, group-level shrinkage
##   edge_summary : data.frame with edge-level diagnostics
## ----------------------------------------------------------------
summarize_node <- function(fit, phi, S_diag, q, m2, nbr_names,
                           edge_threshold = 0.1) {
  n  <- nrow(phi)
  SS <- ncol(fit$theta_store)
  
  ## Group kappa: kappa_j = 1/(1 + lambda_j)
  ## kappa near 1 -> fully shrunk (no edge), near 0 -> active edge
  kappa_samples <- 1 / (1 + fit$lambda_store)   # q x SS
  kappa_mean    <- rowMeans(kappa_samples)
  kappa_sd      <- apply(kappa_samples, 1, sd)
  
  ## Recover spatial fields: beta_j(l) = phi %*% theta_j (per posterior sample)
  edge_maps    <- matrix(0, n, q)
  edge_maps_sd <- matrix(0, n, q)
  
  for (j in 1:q) {
    idx           <- ((j - 1) * m2 + 1):(j * m2)
    th_j_samples  <- fit$theta_store[idx, , drop = FALSE]   # m2 x SS
    field_samples <- phi %*% th_j_samples                    # n x SS
    
    edge_maps[, j]    <- rowMeans(field_samples)
    edge_maps_sd[, j] <- apply(field_samples, 1, sd)
  }
  
  ## Edge summary
  edge_active <- kappa_mean < 0.5 #threshold Kappa's
  
  edge_summary <- data.frame(
    neighbor    = nbr_names,
    kappa_mean  = round(kappa_mean, 4),
    kappa_sd    = round(kappa_sd, 4),
    edge_active = edge_active,
    map_max     = round(apply(edge_maps, 2, max), 4),
    map_min     = round(apply(edge_maps, 2, min), 4),
    map_range   = round(apply(edge_maps, 2, function(x) diff(range(x))), 4),
    map_sd      = round(apply(edge_maps, 2, sd), 4)
  )
  
  list(
    edge_maps    = edge_maps,
    edge_maps_sd = edge_maps_sd,
    kappa_mean   = kappa_mean,
    kappa_sd     = kappa_sd,
    edge_summary = edge_summary
  )
}


## ----------------------------------------------------------------
## 5. Main Wrapper
##    Y        : n x p normalized expression matrix (spots x cell types)
##    coords   : n x 2 spatial coordinates
##    m        : HSGP basis functions per dimension
##    nu       : Matern smoothness (0.5, 1.5, 2.5)
##    rho      : GP length-scale (NULL = auto from median pairwise dist)
##    symmetry : "AND" (conservative) or "OR" (liberal)
##    n_cores  : number of parallel workers (default: detectCores() - 1)
##               set to 1 to run sequentially (useful for debugging)
## ----------------------------------------------------------------
gp_group_horseshoe_graph <- function(Y, coords,
                                     m = 5, nu = 1.5, rho = NULL,
                                     nmc = 3000, burn = 1000, thin = 5,
                                     edge_threshold = 0.1,
                                     symmetry = "AND",
                                     n_cores = NULL) {
  library(parallel)
  
  n  <- nrow(Y)
  p  <- ncol(Y)
  cell_types <- colnames(Y)
  if (is.null(cell_types)) cell_types <- paste0("CT", 1:p)
  
  ## Resolve number of cores
  max_cores <- detectCores(logical = FALSE)
  if (is.null(n_cores)) n_cores <- max(1L, max_cores - 1L)
  n_cores <- min(n_cores, p)
  cat("Using", n_cores, "core(s) for", p, "nodewise regressions.\n\n")
  
  ## Build HSGP basis once — shared across all workers
  cat("Building HSGP basis (m =", m, ", m^2 =", m^2, ", nu =", nu, ")...\n")
  hsgp   <- build_hsgp_basis(coords, m = m, nu = nu, rho = rho)
  phi    <- hsgp$phi
  S_diag <- hsgp$S_diag
  m2     <- hsgp$m2
  cat("  Spectral density range: [", round(min(S_diag), 6),
      ",", round(max(S_diag), 4), "]\n\n")
  
  ## Scale Y once
  Y_scaled <- scale(Y)
  
  ## Per-node worker function
  run_node <- function(s) {
    y_s       <- as.vector(Y_scaled[, s])
    X_s       <- Y_scaled[, -s, drop = FALSE]
    nbr_names <- cell_types[-s]
    q         <- ncol(X_s)
    
    fit  <- gp_group_horseshoe_sampler(
      y = y_s, X_s = X_s, phi = phi, S_diag = S_diag,
      nmc = nmc, burn = burn, thin = thin
    )
    summ <- summarize_node(fit, phi, S_diag, q, m2, nbr_names, edge_threshold)
    
    list(
      s            = s,
      nbr_idx      = (1:p)[-s],
      edge_active  = as.integer(summ$edge_summary$edge_active),
      kappa_mean   = summ$kappa_mean,
      edge_maps    = summ$edge_maps,
      edge_summary = summ$edge_summary
    )
  }
  
  ## Run nodewise regressions
  if (n_cores == 1L) {
    node_results <- lapply(1:p, function(s) {
      cat("=== Node", s, ":", cell_types[s], "===\n")
      res <- run_node(s)
      cat("  Active edges:", sum(res$edge_active), "/", p - 1, "\n\n")
      res
    })
    
  } else {
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
    
    clusterExport(cl, varlist = c(
      "Y_scaled", "cell_types", "phi", "S_diag", "m2", "p",
      "nmc", "burn", "thin", "edge_threshold",
      "gp_group_horseshoe_sampler", "summarize_node",
      "matern_spectral_density", "build_hsgp_basis"
    ), envir = environment())
    
    clusterEvalQ(cl, library(MASS))
    
    cat("Dispatching", p, "nodes to", n_cores, "workers...\n")
    node_results <- parLapply(cl, 1:p, run_node)
    cat("All nodes complete.\n\n")
  }
  
  ## Assemble results
  adj_matrix        <- matrix(0,  p, p, dimnames = list(cell_types, cell_types))
  kappa_matrix      <- matrix(NA, p, p, dimnames = list(cell_types, cell_types))
  edge_maps_list    <- vector("list", p)
  edge_summary_list <- vector("list", p)
  names(edge_maps_list) <- names(edge_summary_list) <- cell_types
  
  for (res in node_results) {
    s <- res$s
    adj_matrix[s, res$nbr_idx]   <- res$edge_active
    kappa_matrix[s, res$nbr_idx] <- res$kappa_mean
    edge_maps_list[[s]]          <- res$edge_maps
    edge_summary_list[[s]]       <- res$edge_summary
  }
  
  ## Symmetrize
  if (symmetry == "AND") {
    adj_sym <- (adj_matrix * t(adj_matrix))
  } else {
    adj_sym <- ((adj_matrix + t(adj_matrix)) > 0) * 1
  }
  diag(adj_sym) <- 0
  
  list(
    adj            = adj_sym,
    adj_raw        = adj_matrix,
    kappa          = kappa_matrix,
    edge_maps      = edge_maps_list,
    edge_summaries = edge_summary_list,
    phi            = phi,
    S_diag         = S_diag,
    coords         = coords,
    hsgp_params    = list(m = m, m2 = m2, nu = nu, rho = hsgp$rho)
  )
}


## ----------------------------------------------------------------
## 6. Visualization
## ----------------------------------------------------------------

## 6a. Spatially varying edge strength map for edge (s -> t)
plot_edge_map <- function(result, node_s, node_t, plot_sd = FALSE) {
  cell_types <- rownames(result$adj)
  coords     <- result$coords
  
  s_idx   <- which(cell_types == node_s)
  t_idx   <- which(cell_types == node_t)
  nbr_pos <- which((1:length(cell_types))[-s_idx] == t_idx)
  
  edge_strength <- result$edge_maps[[node_s]][, nbr_pos]
  is_active     <- result$adj[node_s, node_t] == 1
  
  df <- data.frame(x = coords[, 1], y = coords[, 2],
                   strength = edge_strength)
  
  title_str <- paste0("Edge: ", node_s, " -> ", node_t,
                      ifelse(is_active, "  [ACTIVE]", "  [shrunk]"))
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    ggplot(df, aes(x = x, y = y, color = strength)) +
      geom_point(size = 1.5) +
      scale_color_gradient2(low = "steelblue", mid = "white",
                            high = "firebrick", midpoint = 0,
                            name = "Edge\nStrength") +
      labs(title = title_str, x = "X", y = "Y") +
      theme_bw(base_size = 12)
  } else {
    pal <- colorRampPalette(c("steelblue", "white", "firebrick"))(100)
    br  <- cut(edge_strength, breaks = 100, labels = FALSE)
    plot(coords, col = pal[br], pch = 16, cex = 0.8,
         main = title_str, xlab = "X", ylab = "Y")
  }
}

## 6b. All edge maps for one node as a panel
plot_all_edges <- function(result, node_s) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    stop("ggplot2 and tidyr required for panel plot.")
  }
  library(ggplot2); library(tidyr)
  
  cell_types <- rownames(result$adj)
  s_idx      <- which(cell_types == node_s)
  nbr_names  <- cell_types[-s_idx]
  coords     <- result$coords
  maps       <- result$edge_maps[[node_s]]   # n x q
  active     <- result$adj[node_s, -s_idx]
  
  df <- as.data.frame(maps)
  colnames(df) <- paste0(nbr_names, ifelse(active == 1, "*", ""))
  df$x <- coords[, 1]; df$y <- coords[, 2]
  
  df_long <- pivot_longer(df, -c(x, y), names_to = "neighbor",
                          values_to = "strength")
  
  ggplot(df_long, aes(x = x, y = y, color = strength)) +
    geom_point(size = 0.8) +
    facet_wrap(~neighbor) +
    scale_color_gradient2(low = "steelblue", mid = "white",
                          high = "firebrick", midpoint = 0) +
    labs(title = paste("All edge maps for:", node_s),
         subtitle = "* = active edge (AND rule)") +
    theme_bw(base_size = 10)
}

## 6c. Global adjacency heatmap with kappa
plot_adj_heatmap <- function(result) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("reshape2", quietly = TRUE)) {
    image(result$adj, main = "Adjacency"); return(invisible())
  }
  library(ggplot2); library(reshape2)
  
  signal <- 1 - result$kappa
  diag(signal) <- NA
  
  df     <- melt(signal, na.rm = TRUE)
  adj_df <- melt(result$adj, na.rm = FALSE)
  colnames(df)     <- c("From", "To", "signal")
  colnames(adj_df) <- c("From", "To", "active_val")
  
  adj_df <- adj_df[adj_df$From != adj_df$To, ]
  df <- merge(df, adj_df, by = c("From", "To"))
  df$active <- df$active_val == 1
  
  ggplot(df, aes(x = To, y = From, fill = signal)) +
    geom_tile(color = "gray85") +
    geom_tile(data = subset(df, active),
              fill = NA, color = "black", linewidth = 1.2) +
    scale_fill_gradient(low = "white", high = "darkred",
                        name = "1 - kappa\n(signal)") +
    labs(title = "Cell-Cell Interaction Graph",
         subtitle = "Black border = active edge | Color = signal strength") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

## 6d. Spectral density diagnostic
plot_spectral_density <- function(result) {
  df <- data.frame(
    basis_idx = 1:length(result$S_diag),
    S         = result$S_diag
  )
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    ggplot(df, aes(x = basis_idx, y = S)) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      scale_y_log10() +
      labs(title = paste0("Matern(nu=", result$hsgp_params$nu,
                          ") Spectral Densities"),
           subtitle = "High S = smooth (low-freq) basis allowed | Low S = rough basis suppressed",
           x = "Basis function index (ordered by frequency)",
           y = "S(omega_k) [log scale]") +
      theme_bw(base_size = 12)
  }
}


## ================================================================
## EXAMPLE USAGE
## ================================================================
if (FALSE) {
  set.seed(42)
  n <- 600; p <- 12
  coords <- matrix(runif(n * 2, 0, 10), n, 2)
  Y <- matrix(rnorm(n * p), n, p)
  colnames(Y) <- paste0("CT", 1:p)
  
  result <- gp_group_horseshoe_graph(
    Y      = Y,
    coords = coords,
    m      = 4,
    nu     = 1.5,
    rho    = NULL,
    nmc    = 3000,
    burn   = 1000,
    thin   = 5,
    edge_threshold = 0.1,
    symmetry = "AND"
  )
  
  print(result$adj)
  print(result$edge_summaries[["CT1"]])
  plot_edge_map(result, "CT1", "CT3")
  plot_all_edges(result, "CT1")
  plot_adj_heatmap(result)
  plot_spectral_density(result)
}
