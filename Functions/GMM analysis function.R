#' Box–Cox logit transformation
#'
#' Applies a Box–Cox–type transformation to the logit of a proportion.
#' This transformation maps values in (0, 1) to the real line and allows
#' flexible control of tail behavior via the parameter \code{lambda}.
#'
#' When \code{lambda = 0}, the transformation reduces to the standard logit:
#' \deqn{\log(x / (1 - x))}.
#'
#' @param x Numeric vector of proportions strictly between 0 and 1.
#'
#' @param lambda Numeric scalar controlling the Box–Cox transformation.
#'   If \code{lambda = 0}, the logit transformation is used.
#'
#' @return Numeric vector of transformed values on the real line.
#'
#' @details
#' The transformation is defined as:
#' \deqn{
#' f(x) =
#' \begin{cases}
#' \log\left(\frac{x}{1-x}\right), & \lambda = 0 \\
#' \frac{(x/(1-x))^\lambda - 1}{\lambda}, & \lambda \neq 0
#' \end{cases}
#' }
#'
#' This transformation is commonly used to stabilize variance and improve
#' Gaussianity when modeling proportions with flexible skewness.
#'
#' @seealso \code{\link{InvBoxCoxLogit}}
#'
#' @export
BoxCoxLogit <- function(x, lambda = -0.07) {
  if (lambda == 0) {
    log(x / (1 - x))
  } else {
    ((x / (1 - x))^lambda - 1) / lambda
  }
}

#' Inverse Box–Cox logit transformation
#'
#' Computes the inverse of the Box–Cox logit transformation, mapping values
#' from the real line back to the (0, 1) interval.
#'
#' When \code{lambda = 0}, this function reduces to the inverse logit
#' (logistic) transformation.
#'
#' @param y Numeric vector of real-valued transformed data.
#'
#' @param lambda Numeric scalar controlling the inverse Box–Cox transformation.
#'   Must be identical to the value used in \code{\link{BoxCoxLogit}}.
#'
#' @return Numeric vector of proportions in the interval (0, 1).
#'
#' @details
#' The inverse transformation is defined as:
#' \deqn{
#' f^{-1}(y) =
#' \begin{cases}
#' \frac{e^y}{1 + e^y}, & \lambda = 0 \\
#' \frac{(\lambda y + 1)^{1/\lambda}}
#' {1 + (\lambda y + 1)^{1/\lambda}}, & \lambda \neq 0
#' \end{cases}
#' }
#'
#' This function is used to recover proportions after modeling in the
#' transformed space.
#'
#' @seealso \code{\link{BoxCoxLogit}}
#'
#' @export
InvBoxCoxLogit <- function(y, lambda = -0.07) {
  if (lambda == 0) {
    exp(y) / (1 + exp(y))
  } else {
    z <- (lambda * y + 1)^(1 / lambda)
    z / (1 + z)
  }
}

#' Bayesian GMM analysis of WL ratio series using Stan
#'
#' Fits a Bayesian Gaussian Mixture Model (GMM) to a univariate series of
#' woody-to-leafy (WL) ratios using a Stan model. The WL ratios are transformed
#' into odd-ratios prior to inference. Posterior distributions of mixture
#' parameters are estimated via MCMC sampling using CmdStan.
#'
#' Optionally produces diagnostic and distribution plots and returns either
#' posterior summaries only or posterior samples for selected parameters.
#'
#' @param WLseries Numeric vector of WL ratios. Values must be strictly between
#'   0 and 1, with no missing values.
#'
#' @param stan_file Character string giving the path to the Stan model file
#'   implementing the Gaussian mixture model.
#'
#' @param graph.plot Logical. If TRUE, produces MCMC trace plots and density
#'   plots of the fitted distributions.
#'
#' @param N_samples Integer. Number of posterior samples to randomly extract
#'   for each parameter (mu, sigma, lambda). If set to 0 (default), no posterior
#'   samples are returned and only posterior summaries are provided.
#'
#' @param study_site Character string used as a label in plots (e.g. site name).
#'
#' @param chains Integer. Number of MCMC chains to run.
#'
#' @param parallel_chains Integer. Number of chains to run in parallel.
#'
#' @param iter_warmup Integer. Number of warmup (burn-in) iterations per chain.
#'
#' @param iter_sampling Integer. Number of post-warmup sampling iterations per
#'   chain.
#'
#' @param refresh Integer. Frequency of progress reporting from CmdStan. Use
#'   0 for silent sampling.
#'
#' @param show_messages Logical. If TRUE, prints Stan sampler messages.
#'
#' @param mc.cores Integer. Number of CPU cores to use for parallel execution.
#'
#' @param check_toolchain Logical. If TRUE, checks that the CmdStan toolchain
#'   is correctly installed before sampling.
#'
#' @return
#' If \code{N_samples == 0}, returns a named numeric vector containing posterior
#' quantiles (25%, 50%, 75%) of all GMM parameters.
#'
#' If \code{N_samples > 0}, returns a list with:
#' \describe{
#'   \item{fit_parameters}{Named numeric vector of posterior quantiles.}
#'   \item{muH, muW}{Posterior samples of the means of herbaceous and woody components.}
#'   \item{sigmaH, sigmaW}{Posterior samples of the standard deviations of each component.}
#'   \item{lambda}{Posterior samples of the Box–Cox logit transformation parameter.}
#' }
#'
#' @details
#' The WL ratios are transformed using a Box–Cox logit transformation inside
#' the Stan model. The mixture consists of two Gaussian components interpreted
#' as herbaceous and woody charcoal sources.
#'
#' Posterior predictive densities are visualized both in the original WL space
#' and in the transformed space when \code{graph.plot = TRUE}.
#'
#' @seealso
#' \code{\link[cmdstanr]{cmdstan_model}},
#' \code{\link[bayesplot]{mcmc_trace}}
#'
#' @author
#' Fiona Cornet, Charly Favier
#'
#' @export
analyze_WLseries_GMM <- function(
    WLseries,
    stan_file       = "Stan/Herb_Woody_Charcoal_GMM.stan",
    graph.plot      = TRUE,
    N_samples       = 0,
    study_site      = "",
    chains          = 5,
    parallel_chains = 5,
    iter_warmup     = 1000,
    iter_sampling   = 2000,
    refresh         = 0,
    show_messages   = FALSE,
    mc.cores        = 5,  
    check_toolchain = TRUE 
) {

  # ------------------------------------------------------------------------
  # Load required libraries
  # ------------------------------------------------------------------------
  required_pkgs <- c(
    "ggplot2", "cmdstanr", "bayesplot", "dplyr",
    "patchwork", "ggpubr"
  )
  
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Required package '%s' is not installed. Please install it first.", pkg))
    }
    # Chargement silencieux
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  # ------------------------------------------------------------------------
  # Optional: check CmdStan toolchain
  # ------------------------------------------------------------------------
  if (check_toolchain) {
    if (is.null(cmdstan_path())) {
      message("CmdStan n'est pas installé. Exécutez `cmdstanr::install_cmdstan()` pour l'installer.")
    } else {
      check_cmdstan_toolchain()
    }
  }
  
  # ------------------------------------------------------------------------
  # Global cores option for parallel execution
  # ------------------------------------------------------------------------
  if (!is.numeric(mc.cores) || length(mc.cores) != 1 || mc.cores <= 0) {
    stop("'mc.cores' must be a single positive numeric value.")
  }
  options(mc.cores = mc.cores)
  # This option is commonly used to control how many cores are used for parallel chains in Stan interfaces. [web:22][web:23][web:26]
  
  # ------------------------------------------------------------------------
  # Input checks
  # ------------------------------------------------------------------------
  if (!is.character(stan_file) || length(stan_file) != 1)
    stop("'stan_file' must be a single character string giving the path to a Stan file.")
  if (!file.exists(stan_file))
    stop(sprintf("Stan file not found: %s", stan_file))
  
  if (!is.numeric(WLseries))
    stop("Argument 'WLseries' must be a numeric vector.")
  if (any(is.na(WLseries)))
    stop("Argument 'WLseries' contains missing values (NA).")
  if (any(WLseries <= 0 | WLseries >= 1))
    stop("All values in 'WLseries' must be strictly between 0 and 1.")
  if (length(WLseries) < 10)
    warning("Few observations (<10) for this site: parameter estimation may be unreliable.")
  
  if (!is.logical(graph.plot) || length(graph.plot) != 1)
    stop("'graph.plot' must be a single logical value (TRUE/FALSE).")
  
  # Basic checks for MCMC control arguments
  for (arg_name in c("chains", "parallel_chains", "iter_warmup", "iter_sampling", "refresh")) {
    val <- get(arg_name)
    if (!is.numeric(val) || length(val) != 1 || val < 0)
      stop(sprintf("Argument '%s' must be a single non-negative numeric value.", arg_name))
  }
  if (!is.logical(show_messages) || length(show_messages) != 1)
    stop("'show_messages' must be a single logical value (TRUE/FALSE).")

  if (!is.numeric(N_samples) || length(N_samples) != 1 || N_samples < 0 || N_samples != floor(N_samples)) {
    stop("'N_samples' must be a single non-negative integer.")
  }
  
  if (N_samples > chains * iter_sampling) {
    stop(sprintf(
      "'N_samples' (%d) must be <= total post-warmup draws (%d = chains × iter_sampling).",
      N_samples, chains * iter_sampling
    ))
  }
  
  message(sprintf("Using Stan model file: %s", stan_file))
  
  # ------------------------------------------------------------------------
  # Compile/load CmdStan model
  # ------------------------------------------------------------------------
  BayesGMM_charcoal <- cmdstan_model(stan_file)
  
  # ------------------------------------------------------------------------
  # Prepare data for Stan
  # ------------------------------------------------------------------------
  data_charcoal_GMM <- list(
    N = length(WLseries),
    WL_odd_ratios = WLseries / (1 - WLseries)
  )
  
  # ------------------------------------------------------------------------
  # Fit the Stan model (all control parameters now come from function arguments)
  # ------------------------------------------------------------------------
  fit <- BayesGMM_charcoal$sample(
    data           = data_charcoal_GMM,
    chains         = chains,
    parallel_chains = parallel_chains,
    iter_warmup    = iter_warmup,
    iter_sampling  = iter_sampling,
    refresh        = refresh,
    show_messages  = show_messages
  )
  
  posterior_draws <- fit$draws(variables = c("theta", "mu", "sigma", "lambda"))
  
  if (graph.plot) plot(mcmc_trace((posterior_draws)))
 
  # ------------------------------------------------------------------------
  # Compute parameter quantiles
  # ------------------------------------------------------------------------
  matrix_param <- apply(posterior_draws, 3, function(param_array) {
    stats::quantile(param_array, probs = c(0.25, 0.5, 0.75))
  })
  
  # ------------------------------------------------------------------------
  # Distributions and visualization
  # ------------------------------------------------------------------------
  if (graph.plot) {
    xx <- seq(1e-4, to = 1 - 1e-4, length.out = 512)
    dx <- diff(xx)[1]
    
    distributions <- data.frame(
      x = xx[-1] - dx/2,
      herbaceous = matrix_param[2,1] * diff(stats::pnorm(BoxCoxLogit(xx, lambda = matrix_param[2,7]),
                                                  matrix_param[2,3], matrix_param[2,5])) / dx,
      woody = matrix_param[2,2] * diff(stats::pnorm(BoxCoxLogit(xx, lambda = matrix_param[2,7]),
                                             matrix_param[2,4], matrix_param[2,6])) / dx
    )
    distributions$total <- distributions$herbaceous + distributions$woody
    
    p1 <- ggplot() +
      geom_density(aes(x = WLseries), color = "#1e88e5", fill = "#1e88e5", alpha = 0.3, adjust=.7) +
      geom_line(data = distributions, aes(x = x, y = herbaceous),
                         color = "#e01e37", linewidth = 1) +
      geom_line(data = distributions, aes(x = x, y = woody),
                         color = "#a7c957", linewidth = 1) +
      geom_line(data = distributions, aes(x = x, y = total),
                         color = "black", linetype = "dashed", size = 1) +
      labs(title = study_site, x = "WL Ratio", y = "Density") +
      theme_pubr()
    
    p2 <- ggplot() +
      geom_density(aes(x = BoxCoxLogit(WLseries, lambda = matrix_param[2,7])), color = "#1e88e5", fill = "#1e88e5", alpha = 0.3, adjust=.7) +
      stat_function(
        fun = function(x) matrix_param[2,1] * dnorm(x, matrix_param[2,3], matrix_param[2,5]),
        color = "#e01e37", size = 1) +
      stat_function(
        fun = function(x) matrix_param[2,2] * dnorm(x, matrix_param[2,4], matrix_param[2,6]),
        color = "#a7c957", size = 1) +
      stat_function(
        fun = function(x) matrix_param[2,1] * dnorm(x, matrix_param[2,3], matrix_param[2,5]) +
          matrix_param[2,2] * dnorm(x, matrix_param[2,4], matrix_param[2,6]),
        color = "black", linetype = "dashed", size = 1) +
      labs(title=study_site,x = "BC logit-transformed WL Ratio", y = "Density")+
      scale_x_continuous(limits = c(-4.5, 4.5), expand = c(0, 0)) +
      theme_pubr() 
    
    plot_combine <- p1 / p2 + plot_annotation(tag_levels = "A")
    plot(plot_combine)
  }
  
  # ------------------------------------------------------------------------
  # Extract and format output parameters
  # ------------------------------------------------------------------------
  matrix_param <- cbind(matrix_param[,1:4], 
                        InvBoxCoxLogit(matrix_param[,3], lambda = matrix_param[2,7]),
                        InvBoxCoxLogit(matrix_param[,4], lambda = matrix_param[2,7]),
                        matrix_param[,5:7])
  colnames(matrix_param)[5:6] <- c("WL[1]", "WL[2]")
  
  param_names <- colnames(matrix_param)
  quantile_names <- rownames(matrix_param) %>% sub("^\\.", "", .)
  
  formatted_colnames <- as.vector(sapply(param_names, function(p) {
    paste0(p, quantile_names)
  }))
  
  fit_parameters <- setNames(as.vector(matrix_param), formatted_colnames)
  
  if(N_samples > 0){
    posterior_draws_matrix = data.frame(apply(posterior_draws, 3, c))
    muH = sample(posterior_draws_matrix$mu.1., N_samples,replace = F)
    muW = sample(posterior_draws_matrix$mu.2., N_samples,replace = F)
    sigmaH = sample(posterior_draws_matrix$sigma.1., N_samples,replace = F)
    sigmaW = sample(posterior_draws_matrix$sigma.2., N_samples,replace = F)
    lambda = sample(posterior_draws_matrix$lambda, N_samples,replace = F)
    return(list(fit_parameters = fit_parameters,
                muH = muH, muW = muW,
                sigmaH = sigmaH, sigmaW = sigmaW,
                lambda = lambda))
   } else{
  return(fit_parameters)
  }
}

#' Bayesian Gaussian Mixture Model analysis of WL charcoal ratios along a sediment core
#'
#' Performs a depth-resolved Bayesian Gaussian Mixture Model (GMM) analysis of
#' woody-to-leafy (WL) charcoal ratios along a sediment core. The core is
#' discretized into depth bins, within which a two-component Gaussian mixture
#' model is fitted using Stan. Posterior samples of mixture parameters are
#' extracted for each depth bin and summarized using posterior quantiles.
#'
#' The function provides:
#' \itemize{
#'   \item depth binning of the sediment core (regular or adaptive),
#'   \item Bayesian inference of a two-component GMM per depth bin,
#'   \item posterior sampling of mixture parameters,
#'   \item depth-resolved uncertainty quantification via posterior quantiles.
#' }
#'
#' This function is intended as a first step in a two-stage workflow. Its output
#' can be passed to \code{\link{analyze_WLcore_proportions}} to reconstruct
#' continuous woody and herbaceous fuel proportions along an age model.
#'
#' @param WLRatio_series Numeric vector of woody-to-leafy (WL) ratios measured
#'   along the sediment core. Values must be strictly between 0 and 1.
#'
#' @param depth_series Numeric vector of depths corresponding to
#'   \code{WLRatio_series}.
#'
#' @param age_series Numeric vector of ages corresponding to
#'   \code{WLRatio_series}.
#'
#' @param n_bins Integer. Number of depth bins used to discretize the core.
#'   Must be smaller than the number of unique depth values.
#'
#' @param sampling_type Character string specifying the depth binning strategy.
#'   Either \code{"regular"} for equal-width depth bins or \code{"adaptive"}
#'   for quantile-based bins with equal numbers of observations.
#'
#' @param stan_file Character string giving the path to the Stan model file
#'   implementing the Gaussian mixture model.
#'
#' @param graph.plot Logical. If TRUE, enables diagnostic or summary plots
#'   produced by \code{\link{analyze_WLseries_GMM}}.
#'
#' @param N_samples Integer. Number of posterior samples to extract per depth bin.
#'
#' @param chains Integer. Number of MCMC chains used for each GMM fit.
#'
#' @param parallel_chains Integer. Number of chains run in parallel.
#'
#' @param iter_warmup Integer. Number of warmup (burn-in) iterations per chain.
#'
#' @param iter_sampling Integer. Number of post-warmup sampling iterations per
#'   chain.
#'
#' @param refresh Integer. Frequency of progress reporting from CmdStan.
#'
#' @param show_messages Logical. If TRUE, prints Stan sampler messages.
#'
#' @param mc.cores Integer. Number of CPU cores used for parallel computation.
#'
#' @param progress Logical. If TRUE, displays a progress bar during depth-bin
#'   GMM fitting.
#'
#' @param check_toolchain Logical. If TRUE, checks that the CmdStan toolchain
#'   is properly installed before model fitting.
#'
#' @return
#' A list containing depth-resolved posterior draws and summaries:
#' \describe{
#'   \item{charcoal_data}{Input data with depth bin assignments.}
#'   \item{depth_bins}{Vector of depth-bin centers.}
#'   \item{muH_series, muW_series}{Posterior samples of component means per depth bin.}
#'   \item{sigmaH_series, sigmaW_series}{Posterior samples of component standard deviations per depth bin.}
#'   \item{lambda_series}{Posterior samples of mixture weights per depth bin.}
#'   \item{n_bins}{Number of depth bins used.}
#'   \item{parameter_quantiles}{Data frame of posterior quantiles (0.05, 0.5, 0.95)
#'         for all GMM parameters along depth.}
#' }
#'
#' @details
#' Gaussian mixture models are fitted independently within each depth bin using
#' \code{\link{analyze_WLseries_GMM}}. No smoothing or temporal reconstruction
#' is performed at this stage; posterior samples are returned in raw
#' depth-resolved form.
#'
#' @seealso
#' \code{\link{analyze_WLseries_GMM}},
#' \code{\link{analyze_WLcore_proportions}}
#'
#' @author
#' Fiona Cornet, Charly Favier
#'
#' @export
analyze_WLcore_GMM <- function(
    WLRatio_series,
    depth_series,
    age_series,
    n_bins,
    sampling_type = c("regular", "adaptive"),
    stan_file       = "Stan/Herb_Woody_Charcoal_GMM.stan",
    graph.plot      = TRUE,
    N_samples     = 1000,
    chains          = 5,
    parallel_chains = 5,
    iter_warmup     = 1000,
    iter_sampling   = 2000,
    refresh         = 0,
    show_messages   = FALSE,
    mc.cores        = 5,
    progress        = TRUE,
    check_toolchain = TRUE
){
  
  # ------------------------------------------------------------------------
  # Libraries
  # ------------------------------------------------------------------------
  # ------------------------------------------------------------------------
  # Load required libraries
  # ------------------------------------------------------------------------
  required_pkgs <- c(
     "dplyr", "progress"
  )
  
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Required package '%s' is not installed. Please install it first.", pkg))
    }
    # Chargement silencieux
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  # ------------------------------------------------------------------------
  # Input checks
  # ------------------------------------------------------------------------
  
  ## ---- WLRatio series ----
  if (!is.numeric(WLRatio_series))
    stop("'WLRatio_series' must be a numeric vector.")
  
  if (anyNA(WLRatio_series))
    stop("'WLRatio_series' contains NA values.")
  
  if (any(WLRatio_series < 0 | WLRatio_series > 1))
    stop("'WLRatio_series' must contain values strictly between 0 and 1.")
  
  ## ---- Depth & Age series ----
  if (!is.numeric(depth_series) || !is.numeric(age_series))
    stop("'depth_series' and 'age_series' must be numeric vectors.")
  
  if (length(depth_series) != length(WLRatio_series) ||
      length(age_series)   != length(WLRatio_series))
    stop("'WLRatio_series', 'depth_series' and 'age_series' must have the same length.")
  
  if (anyNA(depth_series) || anyNA(age_series))
    stop("'depth_series' or 'age_series' contains NA values.")

  ## ---- n_bins ----
  
  # Valeur par défaut dépendante de WLRatio_series
  if (is.null(n_bins)) {
    n_bins <- floor(length(WLRatio_series) / 200)
  }
  
  if (!is.numeric(n_bins) || length(n_bins) != 1 || n_bins < 2 || n_bins %% 1 != 0)
    stop("'n_bins' must be a single integer >= 2.")
  
  if (n_bins > length(unique(depth_series)))
    stop("'n_bins' must be smaller than the number of unique depth values.")
  
  ## ---- sampling_type ----
  sampling_type <- match.arg(sampling_type)
  
  ## ---- Stan file ----
  if (!is.character(stan_file) || length(stan_file) != 1)
    stop("'stan_file' must be a single character string.")
  
  if (!file.exists(stan_file))
    stop(sprintf("Stan file not found: %s", stan_file))
  
  ## ---- Graphical / logical flags ----
  for (arg in c("graph.plot", "progress", "show_messages", "check_toolchain")) {
    if (!is.logical(get(arg)) || length(get(arg)) != 1)
      stop(sprintf("'%s' must be a single logical value.", arg))
  }
  
  ## ---- MCMC settings ----
  for (arg in c("chains", "parallel_chains", "iter_warmup", "iter_sampling", "refresh")) {
    val <- get(arg)
    if (!is.numeric(val) || length(val) != 1 || val < 0 || val %% 1 != 0)
      stop(sprintf("'%s' must be a single non-negative integer.", arg))
  }
  
  ## ---- N_samples ----
  if (!is.numeric(N_samples) || length(N_samples) != 1 ||
      N_samples < 0 || N_samples %% 1 != 0)
    stop("'N_samples' must be a non-negative integer.")
  
  if (N_samples > chains * iter_sampling)
    stop(sprintf(
      "'N_samples' (%d) must be <= total post-warmup draws (%d = chains × iter_sampling).",
      N_samples, chains * iter_sampling
    ))
  
  ## ---- mc.cores ----
  if (!is.numeric(mc.cores) || length(mc.cores) != 1 || mc.cores < 1)
    stop("'mc.cores' must be a positive integer.")

  # ------------------------------------------------------------------------
  # 1. Preprocessing
  # ------------------------------------------------------------------------
  charcoal_data <- data.frame(
    Depth   = depth_series,
    Age     = age_series,
    WLRatio = WLRatio_series
  )
  
  charcoal_data <- charcoal_data %>%
    dplyr::mutate(WL.gt.0.5 = as.integer(WLRatio > 0.5))
  
  if (sampling_type == "regular") {
    breaks  <- seq(min(depth_series), max(depth_series), length.out = n_bins + 1)
  } else {
    breaks  <- stats::quantile(depth_series, probs = seq(0, 1, length.out = n_bins + 1))
  }
  centers <- (breaks[-1] + breaks[-length(breaks)]) / 2
  
  charcoal_data <- charcoal_data %>%
    dplyr::mutate(
      Depthbin        = cut(Depth, breaks = breaks, include.lowest = TRUE),
      Depthbin_center = centers[as.integer(Depthbin)]
    )
  
  # ------------------------------------------------------------------------
  # 2. GMM per depth bin
  # ------------------------------------------------------------------------
  depth_bins <- sort(unique(charcoal_data$Depthbin_center))
  n_depthbins <- length(depth_bins)
  
  muH_series    <- matrix(NA, n_depthbins, N_samples)
  muW_series    <- matrix(NA, n_depthbins, N_samples)
  sigmaH_series <- matrix(NA, n_depthbins, N_samples)
  sigmaW_series <- matrix(NA, n_depthbins, N_samples)
  lambda_series <- matrix(NA, n_depthbins, N_samples)
  
  if (progress) {
    pb <- progress::progress_bar$new(
      total = n_depthbins,
      format = "GMM bins [:bar] :percent eta: :eta"
    )
  }
  
  for (k in seq_along(depth_bins)) {
    
    charcoal_extr <- charcoal_data %>%
      dplyr::filter(
        Depthbin_center == depth_bins[k],
        WLRatio > 0.04,
        WLRatio < 0.96
      )
    
    output <- analyze_WLseries_GMM(
      WLseries        = charcoal_extr$WLRatio,
      stan_file       = stan_file,
      N_samples       = N_samples,
      study_site      = depth_bins[k],
      chains          = chains,
      parallel_chains = parallel_chains,
      iter_warmup     = iter_warmup,
      iter_sampling   = iter_sampling,
      refresh         = refresh,
      show_messages   = show_messages,
      mc.cores        = mc.cores,
      check_toolchain = check_toolchain
    )
    
    muH_series[k, ]    <- output$muH
    muW_series[k, ]    <- output$muW
    sigmaH_series[k, ] <- output$sigmaH
    sigmaW_series[k, ] <- output$sigmaW
    lambda_series[k, ] <- output$lambda
    
    if (progress) pb$tick()
  }
  
  # ------------------------------------------------------------------------
  # 3. Quantiles
  # ------------------------------------------------------------------------
  get_quantiles <- function(mat, name, depth) {
    q <- t(apply(mat, 1, quantile, probs = c(0.05, 0.5, 0.95)))
    data.frame(
      param = name,
      depth = depth,
      q0.05 = q[, 1],
      q0.5  = q[, 2],
      q0.95 = q[, 3]
    )
  }
  
  parameter_quantiles <- dplyr::bind_rows(
    get_quantiles(muH_series,    "muH",    depth_bins),
    get_quantiles(muW_series,    "muW",    depth_bins),
    get_quantiles(sigmaH_series, "sigmaH", depth_bins),
    get_quantiles(sigmaW_series, "sigmaW", depth_bins),
    get_quantiles(lambda_series, "lambda", depth_bins)
  )
  
  # ------------------------------------------------------------------------
  # Return
  # ------------------------------------------------------------------------
  return(list(
    charcoal_data       = charcoal_data,
    depth_bins          = depth_bins,
    muH_series          = muH_series,
    muW_series          = muW_series,
    sigmaH_series       = sigmaH_series,
    sigmaW_series       = sigmaW_series,
    lambda_series       = lambda_series,
    n_bins              = n_bins,
    parameter_quantiles = parameter_quantiles
  ))
}

#' Reconstruction of woody vs herbaceous fuel proportions from GMM posterior draws
#'
#' Reconstructs continuous woody and herbaceous fuel proportions along an age
#' model by propagating posterior draws of Gaussian Mixture Model (GMM)
#' parameters obtained from \code{\link{analyze_WLcore_GMM_quantiles}}.
#'
#' Posterior samples of GMM parameters (means and standard deviations of woody
#' and herbaceous components) are first smoothed along depth using B-spline
#' regression. These smoothed parameters are then combined with an empirical
#' estimate of woody fuel occurrence to reconstruct continuous fuel proportions
#' with full uncertainty propagation.
#'
#' The function performs:
#' \itemize{
#'   \item spline-based smoothing of GMM posterior parameters along depth,
#'   \item projection of posterior parameters onto a target depth/age grid,
#'   \item empirical estimation of woody fuel occurrence using local likelihood,
#'   \item reconstruction of woody and herbaceous fuel proportions per posterior draw,
#'   \item computation of uncertainty envelopes from posterior quantiles.
#' }
#'
#' @param gmm_draws List returned by \code{\link{analyze_WLcore_GMM_quantiles}},
#'   containing posterior samples of GMM parameters per depth bin as well as
#'   the original charcoal data.
#'
#' @param depth_estimation Numeric vector of depths at which GMM parameters
#'   should be evaluated and fuel proportions reconstructed.
#'
#' @param age_estimation Numeric vector of ages corresponding to
#'   \code{depth_estimation}. Must have the same length and define the
#'   reconstruction axis.
#'
#' @param nn_locfit Numeric in (0, 1]. Nearest-neighbour fraction used by
#'   \code{\link[locfit]{locfit}} to estimate the empirical probability of
#'   woody fuel occurrence (WLRatio > 0.5) as a function of age.
#'
#' @param h_locfit Positive numeric value giving the bandwidth (in age units)
#'   used by \code{\link[locfit]{locfit}}.
#'
#' @param progress Logical. If TRUE, displays a progress bar during posterior
#'   reconstruction across samples.
#'
#' @return
#' A data frame containing reconstructed woody and herbaceous fuel proportions
#' along the \code{age_estimation} grid, including uncertainty envelopes:
#' \describe{
#'   \item{Age}{Reconstruction age grid.}
#'   \item{P_W}{Median reconstructed woody fuel proportion.}
#'   \item{P_W_lower, P_W_upper}{Lower and upper credible bounds (5\% and 95\%) for woody fuel proportion.}
#'   \item{P_H}{Median reconstructed herbaceous fuel proportion.}
#'   \item{P_H_lower, P_H_upper}{Lower and upper credible bounds for herbaceous fuel proportion.}
#' }
#'
#' @details
#' For each posterior draw, GMM parameters are smoothed along depth using
#' B-spline regression and evaluated at the target depth grid. These parameters
#' define a probabilistic mapping between observed woody occurrence and
#' underlying fuel proportions.
#'
#' An empirical probability of woody fuel occurrence is estimated from the
#' original charcoal data using local likelihood regression on the binary
#' indicator \code{WLRatio > 0.5}. This empirical curve is combined with the
#' posterior GMM parameters to reconstruct continuous fuel proportions with
#' uncertainty propagation.
#'
#' This function represents the second step of a two-stage workflow:
#' \enumerate{
#'   \item \code{\link{analyze_WLcore_GMM_quantiles}} — depth-resolved Bayesian GMM inference
#'   \item \code{analyze_WLcore_proportions} — continuous fuel reconstruction along an age model
#' }
#'
#' @seealso
#' \code{\link{analyze_WLcore_GMM_quantiles}},
#' \code{\link[locfit]{locfit}},
#' \code{\link[splines]{bs}}
#'
#' @author
#' Fiona Cornet, Charly Favier
#'
#' @export
analyze_WLcore_proportions <- function(
    gmm_draws,
    depth_estimation,
    age_estimation,
    nn_locfit = 0.03,
    h_locfit  = 20,
    progress  = TRUE
){
  # ------------------------------------------------------------------------
  # Libraries
  # ------------------------------------------------------------------------
  required_pkgs <- c("dplyr", "splines", "progress", "locfit")
  
  suppressPackageStartupMessages(
    lapply(required_pkgs, require, character.only = TRUE)
  )
  
  # ------------------------------------------------------------------------
  # Input checks
  # ------------------------------------------------------------------------
  if (!is.list(gmm_draws))
    stop("'gmm_draws' must be a list returned by analyze_WLcore_GMM_quantiles().")
  
  required <- c(
    "charcoal_data", "depth_bins",
    "muH_series", "muW_series",
    "sigmaH_series", "sigmaW_series",
    "lambda_series", "n_bins"
  )
  
  miss <- setdiff(required, names(gmm_draws))
  if (length(miss) > 0)
    stop("gmm_draws is missing elements: ", paste(miss, collapse = ", "))
  # ---- Consistency of GMM draw matrices ----
  gmm_mats <- list(
    muH_series    = gmm_draws$muH_series,
    muW_series    = gmm_draws$muW_series,
    sigmaH_series = gmm_draws$sigmaH_series,
    sigmaW_series = gmm_draws$sigmaW_series,
    lambda_series = gmm_draws$lambda_series
  )
  
  # All must be matrices
  if (!all(vapply(gmm_mats, is.matrix, logical(1))))
    stop("All GMM parameter series must be matrices.")
  
  # All must have identical dimensions
  dims <- lapply(gmm_mats, dim)
  if (length(unique(dims)) != 1)
    stop("All GMM parameter matrices must have identical dimensions (same bins × samples).")
  rm(gmm_mats, dims)
  
  if (!is.numeric(depth_estimation) || anyNA(depth_estimation))
    stop("'depth_estimation' must be a numeric vector without NA.")
  
  if (!is.numeric(age_estimation) || anyNA(age_estimation))
    stop("'age_estimation' must be a numeric vector without NA.")
  
  # ---- depth_estimation / age_estimation consistency ----
  if (length(depth_estimation) != length(age_estimation))
    stop(
      "'depth_estimation' and 'age_estimation' must have the same length ",
      "(they define the same reconstruction axis)."
    )
  
  if (!is.numeric(nn_locfit) || nn_locfit <= 0 || nn_locfit > 1)
    stop("'nn_locfit' must be in (0, 1].")
  
  if (!is.numeric(h_locfit) || h_locfit <= 0)
    stop("'h_locfit' must be positive.")
  
  if (!is.logical(progress) || length(progress) != 1)
    stop("'progress' must be a single logical value.")
  
  # ------------------------------------------------------------------------
  # Setup
  # ------------------------------------------------------------------------
  charcoal_data <- gmm_draws$charcoal_data
  depth_bins    <- gmm_draws$depth_bins
  muH_series    <- gmm_draws$muH_series
  muW_series    <- gmm_draws$muW_series
  sigmaH_series <- gmm_draws$sigmaH_series
  sigmaW_series <- gmm_draws$sigmaW_series
  n_bins        <- gmm_draws$n_bins
  N_samples     <- ncol(muH_series)
  
  series_a <- matrix(NA, length(depth_estimation), N_samples)
  series_b <- matrix(NA, length(depth_estimation), N_samples)
  
  results1 <- data.frame(
    Depthbin = c(min(depth_estimation), depth_bins, max(depth_estimation))
  )
  
  if (progress) {
    pb <- progress::progress_bar$new(
      total = N_samples,
      format = "Reconstruction [:bar] :percent eta: :eta"
    )
  }
  
  # ------------------------------------------------------------------------
  # 4. Reconstruction per posterior draw
  # ------------------------------------------------------------------------
  for (i in seq_len(N_samples)) {
    
    results1$muH    <- c(muH_series[1, i],    muH_series[, i],    tail(muH_series[, i], 1))
    results1$muW    <- c(muW_series[1, i],    muW_series[, i],    tail(muW_series[, i], 1))
    results1$sigmaH <- c(sigmaH_series[1, i], sigmaH_series[, i], tail(sigmaH_series[, i], 1))
    results1$sigmaW <- c(sigmaW_series[1, i], sigmaW_series[, i], tail(sigmaW_series[, i], 1))
    
    df <- max(2, round(n_bins / 3))
    
    fit_muH    <- stats::lm(muH    ~ splines::bs(Depthbin, df = df), data = results1)
    fit_muW    <- stats::lm(muW    ~ splines::bs(Depthbin, df = df), data = results1)
    fit_sigmaH <- stats::lm(sigmaH ~ splines::bs(Depthbin, df = df), data = results1)
    fit_sigmaW <- stats::lm(sigmaW ~ splines::bs(Depthbin, df = df), data = results1)
    
    muH <- predict(fit_muH,    newdata = data.frame(Depthbin = depth_estimation))
    muW <- predict(fit_muW,    newdata = data.frame(Depthbin = depth_estimation))
    sH  <- predict(fit_sigmaH, newdata = data.frame(Depthbin = depth_estimation))
    sW  <- predict(fit_sigmaW, newdata = data.frame(Depthbin = depth_estimation))
    
    p1 <- stats::pnorm(0, muH, sH)
    p2 <- stats::pnorm(0, muW, sW)
    
    series_a[, i] <- 1 / (p1 - p2)
    series_b[, i] <- series_a[, i] * (p1 - 1)
    
    if (progress) pb$tick()
  }
  
  # ------------------------------------------------------------------------
  # 5. Empirical proportion curve
  # ------------------------------------------------------------------------
  fit_locfit_prop <- locfit::locfit(
    WL.gt.0.5 ~ lp(Age, nn = nn_locfit, h = h_locfit),
    data = charcoal_data
  )
  
  WL_smooth <- predict(
    fit_locfit_prop,
    newdata = data.frame(Age = age_estimation)
  )
  
  prop_rep <- series_a * WL_smooth + series_b
  q <- t(apply(prop_rep, 1, quantile, probs = c(0.05, 0.5, 0.95)))
  
  series_fuel <- data.frame(
    Age       = age_estimation,
    P_W_lower = pmax(0,  q[, 1]),
    P_W       = pmin(1, pmax(0, q[, 2])),
    P_W_upper = pmin(1, q[, 3])
  )
  
  series_fuel$P_H_lower <- 1 - series_fuel$P_W_upper
  series_fuel$P_H       <- 1 - series_fuel$P_W
  series_fuel$P_H_upper <- 1 - series_fuel$P_W_lower
  
  return(series_fuel)
}
