##########################################################################
#  Function for IBMT (Intensity-based Moderated T-statistic)
#  Written by: Maureen Sartor, University of Cincinnati, 2006
##########################################################################
##
##  This function adjusts the T-statistics and p-values from a linear
##  model analysis of microarrays.  The method contains elements similar in
##  nature both to Smyth's eBayes function in limma and to the Cyber-T 
##  program (Baldi, 2001).  It is an empirical hierarchical Bayesian method. 
##  Local regression and empirical bayesian theory are used to
##  determine the prior degrees of freedom and the predicted background (prior)
##  variance for each gene dependent on average spot intensity level.  
##  The moderated T-statistic uses a weighted average of prior and likelihood
##  variances, and the posterior degrees of freedom are simply the sum of
##  prior and likelihood degrees of freedom.
##
##  Please acknowledge your use of IBMT in publications by referencing:
##  Sartor MA, Tomlinson CR, Wesselkamper SC, Sivaganesan S, Leikauf GD, and
##  Medvedovic M. Intensity-based hierarchical Bayes method improves testing for
##  differentially expressed genes in microarray experiments. BMC Bioinformatics, 
##  2006.
##
##  Inputs:
##  2 objects: mdata and testcol
##  "mdata" should be a list object from the lmFit or eBayes fcn. in  
##       limma, or at least have attributes named sigma, Amean,  
##	   df.residual, coefficients, and stdev.unscaled.
##  "testcol" is an integer or vector indicating the column(s) of
##       mdata$coefficients for which the function is to be performed.
##
##  Outputs:
##  object is augmented form of "mdata" (the input), with the additions being:
##	IBMT.t	 - posterior t-value for IBMT
##	IBMT.p	 - P-value for IBMT
##	IBMT.dfprior - prior degrees of freedom for IBMT
##	IBMT.priorvar- prior variance for IBMT
##	IBMT.postvar - posterior variance for IBMT
##
##  Example Function Call:
##      IBMT.results <- IBMT(eBayes.output,1:4)
##  For further help on implementing function, contact sartorma@ucmail.uc.edu
###########################################################################

IBMT<-function(mdata,testcol) {
   library("stats")
   library("limma")
  
logVAR<-log(mdata$sigma^2)
	df<-mdata$df.residual
	numgenes<-length(logVAR[df>0])	
	df[df==0]<-NA
	eg<-logVAR-digamma(df/2)+log(df/2)
	egpred<-loessFit(eg,mdata$Amean,iterations=1,span=0.3)$fitted
	myfct<- (eg-egpred)^2 - trigamma(df/2)
	print("Local regression fit")

	mean.myfct<-mean(myfct,na.rm=TRUE)
	priordf<-vector(); testd0<-vector()
	for (i in 1:(numgenes*10)) {
		testd0[i]<-i/10
		priordf[i]= abs(mean.myfct-trigamma(testd0[i]/2))
		if (i>2) {
			if (priordf[i-2]<priordf[i-1]) { break }
		}
	}
	d0<-testd0[match(min(priordf),priordf)]
	print("Prior degrees freedom found")

	s02<-exp(egpred + digamma(d0/2) - log(d0/2))

	post.var<- (d0*s02 + df*mdata$sigma^2)/(d0+df)
	post.df<-d0+df
	IBMTt<-mdata$coefficients[,testcol]/(mdata$stdev.unscaled[,testcol]*sqrt(post.var))
	IBMTp<-2*(1-pt(abs(IBMTt),post.df))
	print("P-values calculated")

    output<-mdata
	output$IBMT.t<-IBMTt
	output$IBMT.p<-IBMTp
	output$IBMT.postvar<-post.var
	output$IBMT.priorvar<-s02
	output$IBMT.dfprior<-d0
	output
}

###############################################################################
#pick picking

fit_gaussians_mod <- function (chromatogram, n_gaussians, min_iterations = 5, max_iterations = 10, min_R_squared = 0.5, 
  method = c("guess", "random"), filter_gaussians_center = TRUE, 
  filter_gaussians_height = 0.15, filter_gaussians_variance_min = 0.1, 
  filter_gaussians_variance_max = 50, filter_gaussians_min_dist = 1, random_seed=12345) # new parameter, filter_gaussians_min_dist (all gaussians need to be at least this distant from each other)
{
  indices <- seq_along(chromatogram)
  iter <- 0
  bestR2 <- 0
  bestCoefs <- NULL
  set.seed(random_seed)
  while ((iter < min_iterations) | (iter < max_iterations & bestR2 < min_R_squared)) { #modified so a guaranteed number of iterations are done
    iter <- iter + 1
    initial_conditions <- make_initial_conditions(chromatogram, 
      n_gaussians, method)
    A <- initial_conditions$A
    mu <- initial_conditions$mu
    sigma <- initial_conditions$sigma
    p_model <- function(x, A, mu, sigma) {
      rowSums(sapply(seq_len(n_gaussians), function(i) A[i] * 
        exp(-((x - mu[i])/sigma[i])^2)))
    }
    fit <- tryCatch({
      suppressWarnings(nls(chromatogram ~ p_model(indices, 
        A, mu, sigma), start = list(A = A, mu = mu, 
        sigma = sigma), trace = FALSE, control = list(warnOnly = TRUE, 
        minFactor = 1/2048)))
    }, error = function(e) {
      e
    }, simpleError = function(e) {
      e
    })
    if ("error" %in% class(fit)) 
      next
    coefs <- coef(fit)
    coefs <- split(coefs, rep(seq_len(3), each = n_gaussians))
    coefs <- setNames(coefs, c("A", "mu", "sigma"))
    if (filter_gaussians_variance_min > 0) {
      sigmas <- coefs[["sigma"]]
      drop <- which(sigmas < filter_gaussians_variance_min)
      if (length(drop) > 0) 
        coefs <- lapply(coefs, `[`, -drop)
    }
    if (filter_gaussians_variance_max > 0) {
      sigmas <- coefs[["sigma"]]
      drop <- which(sigmas > filter_gaussians_variance_max)
      if (length(drop) > 0) 
        coefs <- lapply(coefs, `[`, -drop)
    }
    if (filter_gaussians_center) {
      means <- coefs[["mu"]]
      drop <- which(means < 0 | means > length(chromatogram))
      if (length(drop) > 0) 
        coefs <- lapply(coefs, `[`, -drop)
    }
    if (filter_gaussians_height > 0) {
      minHeight <- max(chromatogram) * filter_gaussians_height
      heights <- coefs[["A"]]
      drop <- which(heights < minHeight)
      if (length(drop) > 0) 
        coefs <- lapply(coefs, `[`, -drop)
    }
    if (filter_gaussians_min_dist > 0){
      peak_dists <- outer(coefs[["mu"]], coefs[["mu"]], "-")
      diag(peak_dists) <- NA
      if(TRUE %in% (abs(peak_dists) < filter_gaussians_min_dist)){
        next
      }
    }
    if (length(coefs[["A"]]) == 0) 
      next
    curveFit <- fit_curve(coefs, indices)
    R2 <- cor(chromatogram, curveFit)^2
    if (R2 > bestR2 & R2 > min_R_squared) {
      bestR2 <- R2
      bestCoefs <- coefs
    }
  }
  if (!is.null(bestCoefs)) {
    curveFit <- fit_curve(bestCoefs, indices)
  }
  else {
    curveFit <- NULL
  }
  results <- list(n_gaussians = n_gaussians, R2 = bestR2, 
    iterations = iter, coefs = bestCoefs, curveFit = curveFit)
  return(results)
}

#choose_gaussians_corr_mod

choose_gaussians_corr_mod <- function (chromatogram, points = NULL, max_gaussians = 5, criterion = c("AICc", 
  "AIC", "BIC"), min_iterations=5, max_iterations = 10, min_R_squared = 0.5, 
  method = c("guess", "random"), filter_gaussians_center = TRUE, 
  filter_gaussians_height = 0.15, filter_gaussians_variance_min = 0.1, 
  filter_gaussians_variance_max = 50, filter_gaussians_min_dist=1, random_seed=12345) 
{
  criterion <- match.arg(criterion)
  if (!is.null(points)) {
    max_gaussians <- min(max_gaussians, floor(points/3))
  }
  fits <- list()
  for (n_gaussians in seq_len(max_gaussians)) fits[[n_gaussians]] <- fit_gaussians_mod(chromatogram, 
    n_gaussians, min_iterations, max_iterations, min_R_squared, method = method, 
    filter_gaussians_center, filter_gaussians_height, filter_gaussians_variance_min, 
    filter_gaussians_variance_max, filter_gaussians_min_dist, random_seed=random_seed)
  models <- map(fits, "coefs")
  drop <- map_lgl(models, is.null)
  fits <- fits[!drop]
  coefs <- map(fits, "coefs")
  if (criterion == "AICc") {
    criteria <- lapply(coefs, gaussian_aicc, chromatogram)
  }
  else if (criterion == "AIC") {
    criteria <- lapply(coefs, gaussian_aic, chromatogram) # corrected!
  }
  else if (criterion == "BIC") {
    criteria <- lapply(coefs, gaussian_bic, chromatogram) # corrected!
  }
  best <- which.min(criteria)
  if (length(best) == 0) {
    return(NULL)
  }
  else {
    return(fits[[best]])
  }
}
build_gaussians_corr_mod <- function (profile_matrix, min_points = 1, min_consecutive = 5, 
  impute_NA = TRUE, smooth = TRUE, smooth_width = 4, max_gaussians = 5, 
  criterion = c("AICc", "AIC", "BIC"), min_iterations=5, max_iterations = 50, 
  min_R_squared = 0.5, method = c("guess", "random"), filter_gaussians_center = TRUE, 
  filter_gaussians_height = 0.15, filter_gaussians_variance_min = 0.5, 
  filter_gaussians_variance_max = 50, filter_gaussians_min_dist=1,  random_seed=12345) 
{
  if (is(profile_matrix, "MSnSet")) {
    profile_matrix <- exprs(profile_matrix)
  }
  filtered <- filter_profiles(profile_matrix, min_points = min_points, 
    min_consecutive = min_consecutive)
  cleaned <- clean_profiles(filtered, impute_NA = impute_NA, 
    smooth = smooth, smooth_width = smooth_width)
  gaussians <- list()
  proteins <- rownames(cleaned)
  P <- length(proteins)
  message(".. fitting Gaussian mixture models to ", P, " profiles")
  pb <- progress_bar$new(format = "fitting :what [:bar] :percent eta: :eta", 
    clear = FALSE, total = P, width = 80)
  max_len <- max(nchar(proteins))
  for (i in seq_len(P)) {
    protein <- proteins[i]
    pb$tick(tokens = list(what = sprintf(paste0("%-", max_len, 
      "s"), protein)))
    chromatogram <- cleaned[protein, ]
    points <- sum(!is.na(profile_matrix[protein, ]))
    gaussian <- choose_gaussians_corr_mod(chromatogram, points, max_gaussians, 
      criterion, min_iterations, max_iterations, min_R_squared, method, 
      filter_gaussians_center, filter_gaussians_height, 
      filter_gaussians_variance_min, filter_gaussians_variance_max, filter_gaussians_min_dist, random_seed=random_seed) # Changed to use choose_gaussians_corr
    gaussians[[protein]] <- gaussian
  }
  return(gaussians)
}

gaussian_aicc <- function(coefs, chromatogram) {
  # first, calculate AIC
  AIC <- gaussian_aic(coefs, chromatogram)
  # second, calculate AICc
  N <- length(chromatogram)
  k <- length(unlist(coefs)) + 1
  AICc <- AIC + (2 * k * (k + 1)) / (N - k - 1)
  return(AICc)
}

gaussian_aic <- function (coefs, chromatogram) 
{
  N <- length(chromatogram)
  indices <- seq_len(N)
  fit <- fit_curve(coefs, indices)
  res <- chromatogram - fit
  w <- rep_len(1, N)
  zw <- w == 0
  loglik <- -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw)) + 
    log(sum(w * res^2)))/2
  k <- length(unlist(coefs)) + 1
  AIC <- 2 * k - 2 * loglik
  return(AIC)
}
gaussian_bic <- function (coefs, chromatogram) 
{
  N <- length(chromatogram)
  indices <- seq_len(N)
  fit <- fit_curve(coefs, indices)
  res <- chromatogram - fit
  w <- rep_len(1, N)
  zw <- w == 0
  loglik <- -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw)) + 
    log(sum(w * res^2)))/2
  k <- length(unlist(coefs)) + 1
  BIC <- log(N) * k - 2 * loglik
  return(BIC)
}


