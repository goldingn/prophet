# fitting functions for greta models - uses same inputs as stan

#' @noRd
#' @import greta
greta_fit <- function (type = c("linear", "logistic"),
                       data,
                       init,
                       iter,
                       ...) {
  
  type <- match.arg(type)
  
  if (type != "linear")
    stop ("not implemented")
  
  T         <- data$T
  K         <- data$K
  t         <- data$t
  y         <- data$y
  S         <- data$S
  A         <- data$A
  t_change  <- data$t_change
  X         <- data$X
  sigmas    <- data$sigmas
  tau       <- data$tau
  
  unlist_tf <- greta:::.internals$utils$samplers$unlist_tf
  
  # parameters & their priors
  k <- normal(0, 5)
  m <- normal(0, 5)
  delta <- laplace(0, tau, dim = length(t_change))
  sigma_obs <- normal(0, 0.5, truncation = c(0, Inf))
  beta <- normal(0, sigmas)
  
  gamma <- -t_change * delta
  
  eta <- (k + A %*% delta) * t + (m + A %*% gamma) + X %*% beta
  distribution(y) <- normal(eta, sigma_obs)
  
  mod <- model(k, m, delta, sigma_obs, beta, gamma)
  
  # set up initial values; flatten init in the correct order
  
  # prophet inits
  inits <- lapply(init(), as.array)
  
  # rearrange inits to match greta's order
  greta_names <- names(mod$dag$parameters_example)
  tf_names <- vapply(names(inits),
                     function(name) mod$dag$tf_name(get(name)$node),
                     "")
  order <- match(greta_names, tf_names)
  inits <- unlist_tf(inits[order])
  
  # optimise with BFGS
  fn <- function(par) {
    mod$dag$send_parameters(par)
    -mod$dag$log_density()
  }
  
  gr <- function(par) {
    mod$dag$send_parameters(par)
    -mod$dag$gradients()
  }
  
  o <- optim(inits, fn, gr,
             method = "BFGS",
             control = list(maxit = iter))
  
  mod$dag$send_parameters(o$par)
  par <- mod$dag$trace_values()
  
  # convert par into a named list
  nm <- names(mod$target_greta_arrays)
  values <- lapply(nm,
                   function (name) {
                     pattern <- paste0("^", name)
                     idx <- grep(pattern, names(par))
                     val <- as.vector(par[idx])
                     names(val) <- NULL
                     val
                   })
  names(values) <- nm
  list(par = values)
  
}

