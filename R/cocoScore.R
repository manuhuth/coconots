#' @title Scoring Rule Based Model Assessment Procedure
#' @description The function computes log, quadratic and ranked probability scores for assessing relative performance of a fitted model.
#' @param coco An object of class coco
#' @param max_x An integer which is used as the maximum count for the computation
#'  of the score (default: `50`)
#' @param julia if TRUE, the scores are computed with \proglang{julia} (default: FALSE).
#' @return a list containing the log score, quadratic score and ranked probability score.
#' @details Scoring rules assign a numerical score based on the predictive
#'  distribution and the observed data  to measure the quality of probabilistic predictions.
#' They are provided here as a model selection tool and are computed as
#'  averages over the relevant set of (in-sample) predictions. Scoring rules are, generally, negatively oriented
#' penalties that one seeks to minimize. The literature has developed a large number of scoring
#' rules and, unless there is a unique and clearly defined underlying decision problem,
#' there is no automatic choice of a (proper) scoring rule to be used in any given situation.
#' Therefore, the use of a variety of scoring rules may be appropriate to take advantage of
#' specific emphases and strengths. Three proper scoring rules
#' (for a definition of the concept of propriety see Gneiting and Raftery, 2007),
#' which Jung, McCabe and Tremayne (2016) found to be particularly useful, are implemented.
#' For more information see the references listed below.
#' 
#' @importFrom JuliaConnectoR juliaGet
#' @references 
#' Czado, C. and Gneitling, T. and Held, L. (2009) Predictive Model Assessment for Count Data. \emph{Biometrics}, \bold{65},  1254--1261.
#'
#' Gneiting, T. and Raftery, A. E. (2007) Strictly proper scoring rules, prediction, and estimation. \emph{Journal
#' of the American Statistical Association}, 102:359-378.
#'
#' Jung, R. C., McCabe, B.P.M. and Tremayne, A.R. (2016). Model validation and diagnostics. \emph{In Handbook of Discrete
#' Valued Time Series}. Edited by Davis, R.A., Holan, S.H., Lund, R. and Ravishanker, N.. Boca Raton: Chapman and
#' Hall, pp. 189--218.
#' 
#' Jung, R. C. and Tremayne, A. R. (2011) Convolution-closed models for count time series with applications. \emph{Journal of Time Series Analysis}, \bold{32}, 268--280.
#' @author Manuel Huth
#' @examples
#' lambda <- 1
#' alpha <- 0.4
#' set.seed(12345)
#' data <- cocoSim(order = 1, type = "Poisson", par = c(lambda, alpha), length = 100)
#' fit <- cocoReg(order = 1, type = "Poisson", data = data)
#'
#' # scoring rules - R implementation
#' score_r <- cocoScore(fit)
#' @export

cocoScore <- function(coco, max_x = 50, julia=FALSE) {
  
  # the pmf is evaluated on 0..max(max_x, x_t), so max_x is bounded by the
  # largest supported count; at least one transition must be scored
  coco_check_support_limit(max_x, "max_x")
  coco_check_count_limit(coco$ts)
  coco_check_length(length(coco$ts), coco$order + 1, "cocoScore",
                    sprintf("an order-%d model scores the observations after the first %d",
                            coco$order, coco$order))
  
  if (!is.null(coco$julia_reg) & julia){
    addJuliaFunctions()
    # one round trip: the three scores (read by key) as one Float64 vector
    coco_score <- as.numeric(JuliaConnectoR::juliaLet(
      'd = compute_scores(fit, max_x); Float64[d["logarithmic_score"], d["quadratic_score"], d["ranked_probability_score"]]',
      fit = coco$julia_reg, max_x = max_x))
    return(list("log.score" = coco_score[1],
                "quad.score" = coco_score[2],
                "rps.score" = coco_score[3],
                "aic" = 2*length(coco$par) - 2 * coco$likelihood,
                "bic" = length(coco$par) * log(length(coco$ts)) - 2 * coco$likelihood
                ))
  } else {
  T <- length(coco$ts)
  par <- coco$par
  seas <- c(1, 2) #will be used as argument in future versions
  data <- coco$ts

  if (coco$type == "Poisson"){
    if (is.null(coco$cov)){
      par <- c(par, 0)
    } else {
      par <- c( par[1:(length(par)-ncol(coco$cov))], 0, par[(length(par)-ncol(coco$cov)+1):length(par)])
    }
  }
  
  if (!is.null(coco$cov)){
    xreg <- coco$cov
    betas <- par[(length(par)-ncol(xreg)+1):length(par)]
    lambda_raw <- xreg %*% betas
    lambdas <- applyLinkFunctionVec(lambda_raw, coco$link_function)

    par_original <- par

  }

  if (any((round(data) != data) | (data < 0))) {
    stop("y must be a non-negative integer")
  }

  # Without covariates the predictive pmf of x_t depends on t only through
  # the history (x_{t-1}[, x_{t-2}]): one density vector per distinct history
  if (is.null(coco$cov) && coco$order == 1) {
    validate_gp1_par(par[1], par[2], par[3])
    tt <- (seas[1] + 1):T
    y_t <- data[tt - seas[1]]
    scores <- coco_scores_by_history(
      data[tt], y_t,
      function(xmax, i) dGP1h_vec(xmax, y_t[i], par[1], par[2], par[3]),
      max_x
    )
    log.score <- mean(scores$log)
    quad.score <- mean(scores$quad)
    rps.score <- mean(scores$rps)
  }

  if (is.null(coco$cov) && coco$order == 2) {
    validate_gp2_par(par[1], par[2], par[3], par[4], par[5])
    tt <- (seas[2] + 1):T
    y_t <- data[tt - seas[1]]
    z_t <- data[tt - seas[2]]
    scores <- coco_scores_by_history(
      data[tt], paste(y_t, z_t),
      function(xmax, i) dGP2h_vec(xmax, y_t[i], z_t[i],
                                  par[1], par[2], par[3], par[4], par[5]),
      max_x
    )
    log.score <- mean(scores$log)
    quad.score <- mean(scores$quad)
    rps.score <- mean(scores$rps)
  }

  #Poisson/GP 1 (covariates)
  if  (!is.null(coco$cov) && coco$order == 1){
    calculate_scores <- function(t) {
      par <- c(lambdas[t], par_original[1:(length(par_original) - ncol(xreg))])
      validate_gp1_par(par[1], par[2], par[3])

      # all three scores share the transition densities P(x | y_t) for
      # x = 0..max_x (plus data[t] if beyond), computed in one C++ call
      dvec <- dGP1h_vec(max(max_x, data[t]), data[t - seas[1]], par[1], par[2], par[3])
      density_obs <- dvec[data[t] + 1]
      coco_check_finite_pmf(c(dvec[1:(max_x + 1)], density_obs), "predictive pmf")

      # log score
      log_score <- -log(density_obs)

      # quadratic score
      norm_p <- sum(dvec[1:(max_x + 1)]^2)
      quad_score <- -2 * density_obs + norm_p

      # ranked probability score
      cdf <- cumsum(dvec)[1:(max_x + 1)]
      rps_score <- sum((cdf - (0:max_x >= data[t]))^2)

      return(c(log_score, quad_score, rps_score))
    }

    # Use sapply for the outer loop
    scores <- sapply((seas[1] + 1):T, calculate_scores)

    # Summarize scores
    log.score <- mean(scores[1, ])
    quad.score <- mean(scores[2, ])
    rps.score <- mean(scores[3, ])

  }# Poisson/GP 1


  #Poisson/GP 2 (covariates)
  if (!is.null(coco$cov) && coco$order == 2){
    calculate_scores <- function(t) {
      par <- c(lambdas[t], par_original[1:(length(par_original) - ncol(xreg))])
      validate_gp2_par(par[1], par[2], par[3], par[4], par[5])

      dvec <- dGP2h_vec(max(max_x, data[t]), data[t - seas[1]], data[t - seas[2]],
                        par[1], par[2], par[3], par[4], par[5])
      density_obs <- dvec[data[t] + 1]
      coco_check_finite_pmf(c(dvec[1:(max_x + 1)], density_obs), "predictive pmf")

      # log score
      log_score <- -log(density_obs)

      # quadratic score
      norm_p <- sum(dvec[1:(max_x + 1)]^2)
      quad_score <- -2 * density_obs + norm_p

      # ranked probability score
      cdf <- cumsum(dvec)[1:(max_x + 1)]
      rps_score <- sum((cdf - (0:max_x >= data[t]))^2)

      return(c(log_score, quad_score, rps_score))
    }

    # Use sapply for the outer loop
    scores <- sapply((seas[2] + 1):T, calculate_scores)

    # Summarize scores
    log.score <- mean(scores[1, ])
    quad.score <- mean(scores[2, ])
    rps.score <- mean(scores[3, ])
  }# Poisson/GP 2

  list_out <- list("log.score" = log.score, "quad.score" = quad.score, "rps.score" = rps.score)

  }#end julia if
  
  list_out["bic"] <- length(coco$par) * log(length(coco$ts)) - 2 * coco$likelihood
  list_out["aic"] <- 2*length(coco$par) - 2 * coco$likelihood
  return(list_out)

}# end function

#' Scores of the no-covariate models, one density vector per history
#'
#' `x` are the observations x_t scored, `key` identifies the history of each
#' x_t (equal keys, equal predictive pmf), and `dens(xmax, i)` returns the
#' predictive pmf at 0..xmax for the history of observation i. For each
#' distinct history the pmf is computed once, up to max(max_x, largest x of
#' that history); its entries do not depend on xmax. Per observation this
#' gives the same log score -log p(x), quadratic score -2 p(x) + sum p^2 and
#' ranked probability score sum((F - 1{k >= x})^2) over k = 0..max_x as the
#' former per-t loop; the latter is evaluated as
#' sum(F^2) - 2 sum_{k >= x} F_k + #{k >= x} (rounding-level differences).
#' Returns list(log, quad, rps), each in the order of `x`.
#' @noRd
coco_scores_by_history <- function(x, key, dens, max_x) {
  n <- length(x)
  log_s <- numeric(n)
  quad_s <- numeric(n)
  rps_s <- numeric(n)
  groups <- split(seq_len(n), match(key, unique(key)))
  for (ix in groups) {
    xs <- x[ix]
    dvec <- dens(max(max_x, xs), ix[1])
    d_obs <- dvec[xs + 1]
    coco_check_finite_pmf(c(dvec[1:(max_x + 1)], d_obs), "predictive pmf")

    log_s[ix] <- -log(d_obs)

    norm_p <- sum(dvec[1:(max_x + 1)]^2)
    quad_s[ix] <- -2 * d_obs + norm_p

    cdf <- cumsum(dvec)[1:(max_x + 1)]
    nk <- length(cdf)
    s2 <- sum(cdf^2)
    suf <- rev(cumsum(rev(cdf)))  # suf[j] = sum(cdf[j:nk])
    r <- rep(s2, length(xs))
    inside <- xs < nk
    r[inside] <- s2 - 2 * suf[xs[inside] + 1] + (nk - xs[inside])
    rps_s[ix] <- r
  }
  list(log = log_s, quad = quad_s, rps = rps_s)
}
