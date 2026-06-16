#!/usr/bin/env Rscript

options(warn = 1)

suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)
  library(RGCCA)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
out_dir <- if (length(script_file) > 0) dirname(normalizePath(script_file)) else getwd()

n <- 200
n_locs <- 101
n_comp <- 3
n_groups <- 4
sigma_noise <- 1

norm_l2 <- function(x) sqrt(as.numeric(t(x) %*% x))

a_gen <- function(x, id) {
  if (id == 0) return(0 * x)
  if (id == 1) return(exp(-80 * (x - 1 / 4)^2))
  if (id == 2) return(exp(-80 * (x - 3 / 4)^2))
  if (id == 3) return(exp(-80 * (x - 1 / 2)^2))
  stop("unknown loading id")
}

generate_data <- function(seed = 0) {
  C <- matrix(c(
    0, 1, 1, 1,
    1, 0, 0, 1,
    1, 0, 0, 1,
    1, 1, 1, 0
  ), ncol = 4, byrow = TRUE)

  locs_D <- seq(0, 1, length.out = n_locs)

  A_locs <- H <- vector("list", n_groups)

  rho <- 0.9
  set.seed(seed)
  Sca <- diag(1, 12, 12)
  Cor <- diag(1, 12, 12)
  Sca[1:4, 1:4] <- 1.5 * Sca[1:4, 1:4]
  Sca[5:8, 5:8] <- 1.2 * Sca[5:8, 5:8]
  Cor[1, 3] <- Cor[3, 1] <- rho
  Cor[2, 4] <- Cor[4, 2] <- -rho
  Cor[4 + 1, 4 + 2] <- Cor[4 + 2, 4 + 1] <- -rho * 0.9
  Cor[4 + 3, 4 + 4] <- Cor[4 + 4, 4 + 3] <- rho * 0.9
  Cor[8 + 1, 8 + 4] <- Cor[8 + 4, 8 + 1] <- rho * 0.8

  HH <- mvrnorm(n = n, mu = rep(0, 12), Sigma = Sca %*% Cor %*% Sca, empirical = TRUE)
  HH[, 8 + 2] <- 0 * HH[, 8 + 2]
  HH[, 8 + 3] <- 0 * HH[, 8 + 3]

  H1 <- HH[, 1:4]
  H2 <- HH[, 5:8]
  H3 <- HH[, 9:12]
  H4 <- 0 * HH[, 9:12]

  A_locs[[1]] <- cbind(a_gen(locs_D, 1), a_gen(locs_D, 2), a_gen(locs_D, 3), a_gen(locs_D, 0))[, 1:n_comp]
  H[[1]] <- cbind(H1[, 1], H2[, 1], H3[, 1], H4[, 1])[, 1:n_comp]

  A_locs[[2]] <- cbind(a_gen(locs_D, 2), a_gen(locs_D, 1), a_gen(locs_D, 0), a_gen(locs_D, 0))[, 1:n_comp]
  H[[2]] <- cbind(H1[, 2], H2[, 2], H3[, 2], H4[, 2])[, 1:n_comp]

  A_locs[[3]] <- cbind(a_gen(locs_D, 1), a_gen(locs_D, 3), a_gen(locs_D, 0), a_gen(locs_D, 0))[, 1:n_comp]
  H[[3]] <- cbind(H1[, 3], H2[, 3], H3[, 3], H4[, 3])[, 1:n_comp]

  A_locs[[4]] <- cbind(a_gen(locs_D, 1), a_gen(locs_D, 3), a_gen(locs_D, 2), a_gen(locs_D, 0))[, 1:n_comp]
  H[[4]] <- cbind(H1[, 4], H2[, 4], H3[, 4], H4[, 4])[, 1:n_comp]

  X_locs <- vector("list", n_groups)
  for (g in 1:n_groups) {
    X_locs[[g]] <- H[[g]] %*% t(A_locs[[g]])
    for (h in 1:n_comp) {
      norm_H <- var(H[[g]][, h])
      norm_H <- if (is.na(norm_H) || norm_H == 0) 1 else sqrt(norm_H)
      norm_A <- norm_l2(A_locs[[g]][, h])
      norm_A <- if (is.na(norm_A) || norm_A == 0) 1 else norm_A
      H[[g]][, h] <- H[[g]][, h] / norm_H
      A_locs[[g]][, h] <- A_locs[[g]][, h] / norm_A
    }
  }

  noise <- vector("list", n_groups)
  set.seed(seed + 1000)
  for (g in 1:n_groups) {
    EE <- rnorm(n_locs * n, sd = sigma_noise)
    EE <- matrix(EE, nrow = n)
    noise[[g]] <- scale(EE, scale = FALSE)
  }

  X <- vector("list", n_groups)
  for (g in 1:n_groups) {
    X[[g]] <- X_locs[[g]] + noise[[g]]
    rownames(X[[g]]) <- paste0("n", 1:n)
    colnames(X[[g]]) <- paste0("g", g, "p", seq_len(n_locs))
  }
  names(X) <- paste0("X", seq_len(n_groups))

  list(X = X, C = C)
}

write_inputs <- function(X) {
  for (g in seq_along(X)) {
    write.csv(format(X[[g]], digits = 16), file = file.path(out_dir, paste0("X", g, ".csv")))
  }
}

write_vector_mtx <- function(x, path) {
  x <- as.numeric(x)
  nz <- which(x != 0)
  lines <- c(
    "%%MatrixMarket matrix coordinate real general",
    paste(length(x), 1, length(nz))
  )
  if (length(nz) > 0) {
    values <- formatC(x[nz], digits = 17, format = "fg")
    lines <- c(lines, paste(nz, 1, values))
  }
  writeLines(lines, path)
}

tau_as_component_major <- function(tau) {
  if (is.null(dim(tau))) {
    tau <- as.numeric(tau)
    if (length(tau) == 1) {
      tau_matrix <- matrix(tau, nrow = n_comp, ncol = n_groups)
    } else if (length(tau) == n_groups) {
      tau_matrix <- matrix(rep(tau, each = n_comp), nrow = n_comp, ncol = n_groups)
    } else if (length(tau) == n_comp * n_groups) {
      tau_matrix <- matrix(tau, nrow = n_comp, ncol = n_groups)
    } else {
      stop("unexpected tau length")
    }
  } else {
    tau_matrix <- as.matrix(tau)
  }
  as.numeric(t(tau_matrix))
}

write_case <- function(data, case_dir, tau) {
  path <- file.path(out_dir, case_dir)
  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  fit <- rgcca(
    data$X,
    method      = "rgcca",
    scale_block = FALSE,
    scale       = FALSE,
    bias        = TRUE,
    connection  = data$C,
    init        = "svd",
    superblock  = FALSE,
    tau         = tau,
    ncomp       = n_comp,
    scheme      = "factorial",
    comp_orth   = TRUE,
    sparsity    = 1,
    verbose     = FALSE,
    quiet       = TRUE,
    tol         = 1e-8
  )

  for (g in seq_len(n_groups)) {
    block_name <- paste0("X", g)
    for (h in seq_len(n_comp)) {
      write_vector_mtx(fit$a[[g]][, h], file.path(path, paste0("ref_weights_", block_name, "_comp", h, ".mtx")))
      write_vector_mtx(fit$astar[[g]][, h], file.path(path, paste0("ref_weights_star_", block_name, "_comp", h, ".mtx")))
      write_vector_mtx(fit$Y[[g]][, h], file.path(path, paste0("ref_components_", block_name, "_comp", h, ".mtx")))
    }
  }

  tau_used <- tau_as_component_major(fit$call$tau)
  write_vector_mtx(tau_used, file.path(path, "ref_tau.mtx"))
}

data <- generate_data(seed = 0)
write_inputs(data$X)
invisible(write_case(data, "cov", 1))
invisible(write_case(data, "cor", 0))
invisible(write_case(data, "rgcca", "optimal"))

cat("RGCCA CRAN references written to", out_dir, "\n")
