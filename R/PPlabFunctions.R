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

############################################################################
################################ pick picking GAUSSIANS
############################################################################

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

############################################################################
####################### testDifferentialExpression_beniFix
############################################################################

testDifferentialExpression_beniFix <- function (featureVals, compare_between = "Condition", level = c("protein", 
                                                                                              "proteoform", "peptide", "complex"), measuredOnly = TRUE) 
{
  level <- match.arg(level)
  featVals <- copy(featureVals)
  if ("complex_id" %in% names(featVals)) {
    setkeyv(featVals, c("feature_id", "complex_id", "apex", 
                        "id", "fraction"))
  }
  else {
    setkeyv(featVals, c("feature_id", "apex", "id", "fraction"))
  }
  message("Excluding peptides only found in one condition...")
  if (measuredOnly) {
    featVals <- subset(featVals, imputedFraction == FALSE)
    featureValsBoth <- filterValsByFractionOverlap(featVals, 
                                                   compare_between)
    featureValsBoth[, `:=`(n_frac, .N), by = c("id", "feature_id", 
                                               "apex", compare_between)]
    featureValsBoth <- subset(featureValsBoth, n_frac > 
                                2)
  }
  else {
    featureValsBoth <- filterValsByOverlap(featVals, compare_between)
  }
  featureValsBoth <- getQuantTraces(featureValsBoth, compare_between)
  message("Testing peptide-level differential expression")
  if ("complex_id" %in% names(featureValsBoth)) {
    grpn = uniqueN(featureValsBoth[, .(id, feature_id, complex_id, 
                                       apex)])
    pb <- txtProgressBar(min = 0, max = grpn, style = 3)
    tests <- featureValsBoth[, {
      setTxtProgressBar(pb, .GRP)
      samples = unique(.SD[, get(compare_between)])
      qints = .SD[, .(s = sum(intensity)), by = .(get(compare_between), 
                                                  Replicate)]
      if (length(unique(design_matrix$Replicate)) > 1) {
        a = t.test(formula = log(qints$s) ~ qints$get, 
                   var.equal = FALSE)
      }
      else {
        a = t.test(formula = intensity ~ get(compare_between), 
                   var.equal = FALSE)
      }
      ints = .SD[imputedFraction == F, .(s = sum(intensity)), 
                 by = .(get(compare_between))]
      int1 = max(0, mean(ints[get == samples[1]]$s), na.rm = T)
      int2 = max(0, mean(ints[get == samples[2]]$s), na.rm = T)
      qint1 = mean(qints[get == samples[1]]$s)
      qint2 = mean(qints[get == samples[2]]$s)
      global_ints = .SD[, .(s = unique(global_intensity)), 
                        by = .(get(compare_between), Replicate)]
      global_ints_imp = .SD[, .(s = unique(global_intensity_imputed)), 
                            by = .(get(compare_between), Replicate)]
      global_int1 = mean(global_ints[get == samples[1]]$s)
      global_int2 = mean(global_ints[get == samples[2]]$s)
      global_int1_imp = mean(global_ints_imp[get == samples[1]]$s)
      global_int2_imp = mean(global_ints_imp[get == samples[2]]$s)
      if (length(unique(design_matrix$Replicate)) > 1) {
        b = t.test(formula = log(global_ints_imp$s) ~ 
                     global_ints_imp$get, var.equal = FALSE)
        global_pVal = b$p.value
        meanDiff = a$estimate[1] - a$estimate[2]
      }
      else {
        global_pVal = 1
        meanDiff = a$estimate
      }
      .(pVal = a$p.value, int1 = int1, int2 = int2, meanDiff = meanDiff, 
        qint1 = qint1, qint2 = qint2, log2FC = log2(qint1/qint2), 
        n_replicates = a$parameter + 1, Tstat = a$statistic, 
        testOrder = paste0(samples[1], ".vs.", samples[2]), 
        global_int1 = global_int1, global_int2 = global_int2, 
        global_log2FC = log2(global_int1/global_int2), 
        global_int1_imp = global_int1_imp, global_int2_imp = global_int2_imp, 
        global_log2FC_imp = log2(global_int1_imp/global_int2_imp), 
        global_pVal = global_pVal)
    }, by = .(id, feature_id, complex_id, apex)]
    close(pb)
  }
  else {
    grpn = uniqueN(featureValsBoth[, .(id, feature_id, apex)])
    pb <- txtProgressBar(min = 0, max = grpn, style = 3)
    tests <- featureValsBoth[, {
      setTxtProgressBar(pb, .GRP)
      samples = unique(.SD[, get(compare_between)])
      qints = .SD[, .(s = sum(intensity)), by = .(get(compare_between), 
                                                  Replicate)]
      if (length(unique(design_matrix$Replicate)) > 1) {
        a = t.test(formula = log(qints$s) ~ qints$get, 
                   var.equal = FALSE)
      }
      else {
        a = t.test(formula = intensity ~ get(compare_between), 
                   var.equal = FALSE)
      }
      ints = .SD[imputedFraction == F, .(s = sum(intensity)), 
                 by = .(get(compare_between))]
      int1 = max(0, mean(ints[get == samples[1]]$s), na.rm = T)
      int2 = max(0, mean(ints[get == samples[2]]$s), na.rm = T)
      qint1 = mean(qints[get == samples[1]]$s)
      qint2 = mean(qints[get == samples[2]]$s)
      global_ints = .SD[, .(s = unique(global_intensity)), 
                        by = .(get(compare_between), Replicate)]
      global_ints_imp = .SD[, .(s = unique(global_intensity_imputed)), 
                            by = .(get(compare_between), Replicate)]
      global_int1 = mean(global_ints[get == samples[1]]$s)
      global_int2 = mean(global_ints[get == samples[2]]$s)
      global_int1_imp = mean(global_ints_imp[get == samples[1]]$s)
      global_int2_imp = mean(global_ints_imp[get == samples[2]]$s)
      if (length(unique(design_matrix$Replicate)) > 1) {
        b = t.test(formula = log(global_ints_imp$s) ~ 
                     global_ints_imp$get, var.equal = FALSE)
        global_pVal = b$p.value
        meanDiff = a$estimate[1] - a$estimate[2]
      }
      else {
        global_pVal = 1
        meanDiff = a$estimate
      }
      .(pVal = a$p.value, int1 = int1, int2 = int2, meanDiff = meanDiff, 
        qint1 = qint1, qint2 = qint2, log2FC = log2(qint1/qint2), 
        n_replicates = a$parameter + 1, Tstat = a$statistic, 
        testOrder = paste0(samples[1], ".vs.", samples[2]), 
        global_int1 = global_int1, global_int2 = global_int2, 
        global_log2FC = log2(global_int1/global_int2), 
        global_int1_imp = global_int1_imp, global_int2_imp = global_int2_imp, 
        global_log2FC_imp = log2(global_int1_imp/global_int2_imp), 
        global_pVal = global_pVal)
    }, by = .(id, feature_id, apex)]
    close(pb)
  }
  tests[is.na(log2FC) & (int1 == 0 | int2 == 0) & (meanDiff == 
                                                     0)]$log2FC <- 0
  tests[is.na(log2FC) & (int1 == 0 | int2 == 0) & (meanDiff > 
                                                     0)]$log2FC <- Inf
  tests[is.na(log2FC) & (int1 == 0 | int2 == 0) & (meanDiff < 
                                                     0)]$log2FC <- -Inf
  if ("proteoform_id" %in% names(featVals)) {
    proteoform_ann <- unique(subset(featVals, select = c("id", 
                                                         "proteoform_id")))
    tests <- merge(tests, proteoform_ann, by = c("id"), 
                   all.x = T, all.y = F, sort = F)
  }
  if (level == "peptide") {
    tests$pBHadj <- p.adjust(tests$pVal, method = "BH")
    pQv <- qvalue::qvalue(tests$pVal, lambda = 0.4)
    tests$qVal <- pQv$qvalues
    if (length(unique(design_matrix$Replicate)) > 1) {
      tests$global_pBHadj <- p.adjust(tests$global_pVal, 
                                      method = "BH")
      global_pQv <- qvalue::qvalue(tests$global_pVal, 
                                   lambda = 0.4)
      tests$global_qVal <- global_pQv$qvalues
    }
    else {
      tests$global_pBHadj <- 1
      tests$global_qVal <- 1
    }
    return(tests)
  }
  else if (level == "proteoform") {
    message("Aggregating to proteoform-level...")
    proteoformtests <- aggregatePeptideTestsToProteoform(tests)
    return(proteoformtests)
  }
  else if (level == "protein") {
    message("Aggregating to protein-level...")
    prottests <- aggregatePeptideTests(tests)
    return(prottests)
  }
  else if (level == "complex") {
    message("Aggregating to complex-level...")
    prottests <- aggregatePeptideTests(tests)
    complextests <- aggregateProteinTests(prottests)
    return(complextests)
  }
  else {
    stop("Specified level is not valid. Please chose between peptide, protein and complex.")
  }
}


############################################################################
###################### getMassAssemblyChange_aljazfix
############################################################################
	
getMassAssemblyChange_aljazfix <- function(tracesList, design_matrix,
                                  compare_between = "Condition",
                                  quantLevel = "protein_id",
                                  plot = FALSE,
                                  PDF = FALSE,
                                  name = "beta_pvalue_histogram"){
  .tracesListTest(tracesList)
  samples <- unique(design_matrix$Sample)
  if(! all(samples %in% names(tracesList))) {
    stop("tracesList and design_matrix do not match. Pleas check sample names.")
  }
  if (! "sum_assembled_norm" %in% names(tracesList[[1]]$trace_annotation)) {
    stop("No assembled mass annotation available, please run annotateMassDistribution first.")
  }
  if (! quantLevel %in% names(tracesList[[1]]$trace_annotation)) {
    stop("quantLevel not available in provided traces.")
  }

  res <- lapply(names(tracesList), function(tr){
    vals <- subset(tracesList[[tr]]$trace_annotation,select=c(quantLevel,"sum_assembled_norm"))
    vals[,Sample := tr]
    return(vals)
  })

  res <- do.call(rbind, res)

  if(length(unique(design_matrix$Replicate)) < 2) {
    if (quantLevel == "protein_id") {
      res_cast <- dcast(res, formula = protein_id ~ Sample, value.var=c("sum_assembled_norm"))
    } else if (quantLevel == "proteoform_id") {
      res_cast <- dcast(res, formula = proteoform_id ~ Sample, value.var=c("sum_assembled_norm"))
    } else {
      stop("Functionality only available for quantLevel proetin_id or proteoform_id.")
    }
    #res_cast[, change := log2(get(samples[1])/(get(samples[2])))]
    #res_cast[change=="NaN", change := 0]
    res_cast[, meanDiff := get(samples[1])-(get(samples[2]))]
    res_cast[, betaPval := 1]
    res_cast[, betaPval_BHadj := 1]
    res_cast[, testOrder := paste0(samples[1],".vs.",samples[2])]
    res_cast <- subset(res_cast, select = c("protein_id","meanDiff","betaPval", "betaPval_BHadj","testOrder"))
    return(res_cast[])
  } else {
    res <- merge(res, design_matrix, by.x="Sample", by.y="Sample_name")
    res[,n_conditions:=length(unique(Condition)), by=c("protein_id")]
    res[,n:=.N, by=c("protein_id")]
    res[,replicates_perCondition:=.N, by=c("protein_id", "Condition")]
    #res[,sum_assembled_norm := ifelse(sum_assembled_norm>0.999,sum_assembled_norm-0.001,sum_assembled_norm)]
    #res[,sum_assembled_norm := ifelse(sum_assembled_norm<0.001,sum_assembled_norm+0.001,sum_assembled_norm)]
    res[,sum_assembled_norm_t := (sum_assembled_norm * (n - 1) + 0.5)/n, by=c("protein_id")]
    res[,unique_perCondition := length(unique(round(sum_assembled_norm, digits = 3))), by=c("protein_id","Condition")]

    diff <- res[, {
      samples = unique(.SD[,get(compare_between)])
      meanX = .SD[, .(m = mean(sum_assembled_norm)), by = .(get(compare_between))]
      meanDiff = meanX[get==samples[1]]$m - meanX[get==samples[2]]$m
      n_conditions <- unique(.SD$n_conditions)
      n_perCondition <- min(.SD$replicates_perCondition)
      n_unique_perCondition <- min(.SD$unique_perCondition)
      if( (n_conditions > 1) & (n_perCondition > 1) & (n_unique_perCondition > 1) ) {
        model = betareg(.SD$sum_assembled_norm_t ~ .SD$Condition)
        stat = lrtest(model)
        p = stat$`Pr(>Chisq)`[2]
        w = wilcox.test(formula = .SD$sum_assembled_norm ~ .SD$Condition) #CHANGED - Removed paired=F, does not work with formula objects anymore
        wilcoxPval = w$p.value
      } else {
        p = 2
        wilcoxPval = 2
      }
      .(meanDiff = meanDiff,
        meanAMF1 = meanX[get==samples[1]]$m,
        meanAMF2 = meanX[get==samples[2]]$m,
        betaPval = p,
        wilcoxPval = wilcoxPval,
        testOrder = paste0(samples[1],".vs.",samples[2]))},
      by = .(get(quantLevel))]

    diff[betaPval==2, betaPval := NA ]
    diff[wilcoxPval==2, wilcoxPval := NA ]
    if (length(unique(design_matrix$Replicate)) > 1) {
      diff[, betaPval_BHadj := p.adjust(betaPval, method = "fdr")]
      Qv <- qvalue::qvalue(diff$betaPval, lambda = 0.4)
      diff[, betaQval := Qv$qvalues]
      #diff[, wilcoxPval_BHadj := p.adjust(wilcoxPval, method = "fdr")]
      #wilcoxPvalQv <- qvalue::qvalue(diff$wilcoxPval, lambda = 0.4)
      #diff[, wilcoxQval := wilcoxPvalQv$qvalues]
    } else {
      diff[, betaPval_BHadj := 1]
      diff[, betaQval := 1]
      #diff[, wilcoxPval_BHadj := 1]
      #diff[, wilcoxQval := 1]
    }

    if(plot==TRUE){
      if(PDF){
        pdf(paste0(name,".pdf"))
      }
      hist(diff$betaPval, breaks = 100)
      hist(diff$wilcoxPval, breaks = 100)
      hist(diff$meanAMF1, breaks = 100)
      hist(diff$meanAMF2, breaks = 100)
      if(PDF){
        dev.off()
      }
    }

    setnames(diff, "get(quantLevel)", quantLevel)
    tests <- subset(diff, select = c("protein_id","meanDiff",
                                     "meanAMF1","meanAMF2",
                                     "betaPval", "betaPval_BHadj",
                                     "betaQval","testOrder",
                                     "wilcoxPval"))

    return(tests[])

  }
}

############################################################################
############################# normalizeByCyclicLoess
############################################################################

normalizeByCyclicLoess <- function(traces_list, window = 3, step = 1, plot = TRUE, PDF = TRUE, name = "normalizeByCyclicLoess") {
  .tracesListTest(traces_list, type = "peptide")
  trace_intensities_long <- lapply(traces_list, extractvaluesForNorm)
  combi_table <- rbindlist(trace_intensities_long, use.names=TRUE, fill=FALSE, idcol="sample")
  combi_table[, filename := paste0(sample,"_",fraction_number)]
  combi_table$fraction_number <- as.numeric(combi_table$fraction_number)
  combi_table <- unique(combi_table)

  combi_table_forPlot <- copy(combi_table)
  combi_table_forPlot$intensity = as.numeric(combi_table_forPlot$intensity)
  combi_table_forPlot[, total_intensity:=sum(intensity), by=c("filename","sample")]
  combi_table_forPlot <- unique(subset(combi_table_forPlot, select=c("fraction_number","total_intensity","sample")))
  pnormdata<-ggplot(combi_table_forPlot, aes(x=fraction_number, y=total_intensity, group=sample)) +
    geom_line(aes(color=sample)) +
    geom_point(aes(color=sample)) +
    theme_classic()
  ggsave(pnormdata,filename=paste0(name,"_priorNormalization.pdf"),width=7,height=3.5)

  combi_table[, intensity := log2(intensity)]
  combi_table[, intensity := ifelse(intensity == -Inf, NA, intensity)]
  #combi_table[, intensity := ifelse(intensity < 0.000001, NA, intensity)]
  combi_table_norm <- normalize_sn(combi_table, window, step)
  combi_table_norm[, intensity := 2^(intensity)]
  combi_table_norm[, intensity := ifelse(intensity == 1, 0, intensity)]

  #saveRDS(combi_table_norm,"combi_table_norm.rds")

  combi_table_toMerge <- subset(combi_table, select = c("filename","id", "sample", "fraction_number"))
  combi_table_norm_final <- merge(combi_table_toMerge, combi_table_norm, by=c("filename","id"))

  combi_table_norm_forPlot <- copy(combi_table_norm_final)
  combi_table_norm_forPlot[, total_intensity:=sum(intensity), by=c("filename","sample")]
  combi_table_norm_forPlot <- unique(subset(combi_table_norm_forPlot, select=c("fraction_number","total_intensity","sample")))
  pnormdata<-ggplot(combi_table_norm_forPlot, aes(x=fraction_number, y=total_intensity, group=sample)) +
    geom_line(aes(color=sample)) +
    geom_point(aes(color=sample)) +
    theme_classic()
  ggsave(pnormdata,filename=paste0(name,"_postNormalization.pdf"),width=7,height=3.5)

  list_norm <- split(combi_table_norm_final, by="sample")
  list_norm_wide <- lapply(list_norm, dcast_backToTraces)

  traces_list_norm <- copy(traces_list)
  sample_names <- names(traces_list)
  for(s_name in sample_names){
    traces_list_norm[[s_name]]$traces <- list_norm_wide[[s_name]]
    traces_list_norm[[s_name]]$fraction_annotation <- subset(traces_list_norm[[s_name]]$fraction_annotation, id %in% names(traces_list_norm[[s_name]]$traces))
    traces_list_norm[[s_name]]$trace_annotation <- subset(traces_list_norm[[s_name]]$trace_annotation, id %in% traces_list_norm[[s_name]]$traces$id)
  }
  .tracesListTest(traces_list_norm, type = "peptide")
  return(traces_list_norm)
}


############# from Benni

	normalize_sn <- function(X, window, step) {
  mx<-dcast(X, id~filename, value.var='intensity', sum)
  ## changes made:
  #  mx[mx<0.00001]=NA this led to the ID column to be all NA. for some reason -
  #  the logical evaluation did not work as previously intended
  #  use the code below instead
  mx <- mx %>%
    dplyr::mutate(across(where(is.numeric), ~ ifelse(. < 0.00001, NA, .)))

  mxs<-as.matrix(mx[,-1])
  rownames(mxs)<-mx$id
  id_mapping<-unique(X[,c("filename","fraction_number")])
  max_sec <- max(X$fraction_number)
  windows_sets<-SlidingWindow("data.frame",c(0:max_sec+1), window, step)
  #lmxn<-lapply(windows_sets,function(X){normalizeMedianValues(mxs[,subset(id_mapping, fraction_number %in% X)$filename])})
  lmxn<-lapply(windows_sets,function(X){limma::normalizeCyclicLoess(mxs[,subset(id_mapping, fraction_number %in% X)$filename])})
  lln<-do.call("rbind",lapply(lmxn, melt, na.rm=TRUE))
  names(lln)<-c("id", "filename", "intensity")
  lln_dt <- as.data.table(lln)
  lln_dt[,mean_intensity := mean(intensity, na.rm=T), by=c("id","filename")]
  lln_dt_sub <- unique(subset(lln_dt, select = c("id","filename","mean_intensity")))
  names(lln_dt_sub)<-c("id", "filename", "intensity")
  #lxn<-ddply(lln, .(id,filename),function(X){mean(X$intensity)})
  #names(lxn)<-c("id", "filename", "intensity")
  #return(lxn)
  return(lln_dt_sub)
}


############################################################################
#################### SlidingWindow
############################################################################
	
SlidingWindow <- function (FUN, data, window, step)
   {
     total <- length(data)
     spots <- seq(from = 1, to = (total - window), by = step)
     result <- vector(length = length(spots))
     for (i in 1:length(spots)) {
       result[i] <- match.fun(FUN)(data[spots[i]:(spots[i] +
                                                    window - 1)])
     }
     return(result)
 }
	
############################################################################
#################### testDifferentialExpression_1repfix_chatgpt
############################################################################
	
testDifferentialExpression_1repfix_chatgpt <- function(featureVals,
                                             compare_between = "Condition",
                                             level = c("protein", "proteoform", "peptide", "complex"),
                                             measuredOnly = TRUE) {
  level <- match.arg(level)
  featVals <- copy(featureVals)
   # Set key based on presence of complex_id
  if ("complex_id" %in% names(featVals)) {
    setkeyv(featVals, c("feature_id", "complex_id", "apex", "id", "fraction"))
  } else {
    setkeyv(featVals, c("feature_id", "apex", "id", "fraction"))
  }
  # Filter based on measuredOnly flag
  message("Excluding peptides only found in one condition...")
  if (measuredOnly) {
    featVals <- subset(featVals, imputedFraction == FALSE)
    featureValsBoth <- filterValsByFractionOverlap(featVals, compare_between)
    featureValsBoth[, n_frac := .N, by = c("id", "feature_id", "apex", compare_between)]
    featureValsBoth <- subset(featureValsBoth, n_frac > 2)
  } else {
    featureValsBoth <- filterValsByOverlap(featVals, compare_between)
  }
  # Get quantitative traces
  featureValsBoth <- getQuantTraces(featureValsBoth, compare_between)

  # Perform differential expression testing
  message("Testing peptide-level differential expression")
############################################################# if "complex_id" #####################################
    if ("complex_id" %in% names(featureValsBoth)) {
    grpn = uniqueN(featureValsBoth[,.(id, feature_id, complex_id, apex)])
    pb <- txtProgressBar(min = 0, max = grpn, style = 3)
    tests <- featureValsBoth[, {
      setTxtProgressBar(pb, .GRP)
      samples = unique(.SD[,get(compare_between)])
      # qints = .SD[useForQuant == T, .(s = sum(intensity)), by = .(get(compare_between))] # this disables a lot of comparisons
      qints = .SD[, .(s = sum(intensity)), by = .(get(compare_between), Replicate)]
      if (length(unique(design_matrix$Replicate)) > 1) {
        a = t.test(formula = log(qints$s) ~ qints$get, var.equal = FALSE)
      } else {
        cond1 <- .SD[get(compare_between) == samples[1], intensity]
        cond2 <- .SD[get(compare_between) == samples[2], intensity]
        a <- t.test(cond1, cond2, paired = T, var.equal = FALSE)
      }
      ints = .SD[imputedFraction == F, .(s = sum(intensity)), by = .(get(compare_between))] # this creates quantitative discrepancies depending on how many fractions are used
      int1 = max(0, mean(ints[get==samples[1]]$s), na.rm=T)
      int2 = max(0, mean(ints[get==samples[2]]$s), na.rm=T)
      qint1 = mean(qints[get==samples[1]]$s)
      qint2 = mean(qints[get==samples[2]]$s)
      global_ints = .SD[, .(s = unique(global_intensity)), by = .(get(compare_between), Replicate)]
      global_ints_imp = .SD[, .(s = unique(global_intensity_imputed)), by = .(get(compare_between), Replicate)]
      global_int1 = mean(global_ints[get==samples[1]]$s)
      global_int2 = mean(global_ints[get==samples[2]]$s)
      global_int1_imp = mean(global_ints_imp[get==samples[1]]$s)
      global_int2_imp = mean(global_ints_imp[get==samples[2]]$s)
      #local_FC_all = log2(qints[get==samples[1]]$s/qints[get==samples[2]]$s)
      #global_FC_all = log2(global_ints_imp[get==samples[1]]$s/global_ints_imp[get==samples[2]]$s)
      #local_vs_global_FC_all = data.table(fc=c(local_FC_all,global_FC_all),sam=c(rep("local",length(local_FC_all)),rep("global",length(global_FC_all))))
      if (length(unique(design_matrix$Replicate)) > 1) {
        b = t.test(formula = log(global_ints_imp$s) ~ global_ints_imp$get, var.equal = FALSE)
        global_pVal = b$p.value
        #c = t.test(formula = local_vs_global_FC_all$fc ~ local_vs_global_FC_all$sam , paired = F, var.equal = FALSE)
        #local_vs_global_pVal = c$p.value
        meanDiff=a$estimate[1]-a$estimate[2]
      } else {
        global_pVal = 1
        #local_vs_global_pVal = 1
        meanDiff=a$estimate
      }

      .(pVal = a$p.value,
        int1 = int1, int2 = int2,
        meanDiff = meanDiff,
        qint1 = qint1, qint2 = qint2, log2FC =  log2(qint1/qint2),
        n_replicates = a$parameter + 1,  Tstat = a$statistic, testOrder = paste0(samples[1],".vs.",samples[2]),
        global_int1 = global_int1, global_int2 = global_int2, global_log2FC = log2(global_int1/global_int2),
        global_int1_imp = global_int1_imp, global_int2_imp = global_int2_imp, global_log2FC_imp = log2(global_int1_imp/global_int2_imp),
        #local_vs_global_log2FC = log2(qint1/qint2)-log2(global_int1/global_int2), local_vs_global_log2FC_imp = log2(qint1/qint2)-log2(global_int1_imp/global_int2_imp),
        global_pVal = global_pVal#, local_vs_global_pVal = local_vs_global_pVal
       )},
      by = .(id, feature_id, complex_id, apex)]
    close(pb)
  } else {
############################################################# if not "complex_id" ##################################
   grpn <- uniqueN(featureValsBoth[, .(id, feature_id, apex)])
    pb <- txtProgressBar(min = 0, max = grpn, style = 3)
  
    tests <- featureValsBoth[, {
      setTxtProgressBar(pb, .GRP)
      samples <- unique(.SD[, get(compare_between)])
      qints = .SD[, .(s = sum(intensity)), by = .(get(compare_between), Replicate)] 
      if (length(unique(design_matrix$Replicate)) > 1) {
          a = t.test(formula = log(qints$s) ~ qints$get, var.equal = FALSE)
        } else {
          cond1 <- .SD[get(compare_between) == samples[1], intensity]
          cond2 <- .SD[get(compare_between) == samples[2], intensity]
          a <- t.test(cond1, cond2, paired = T, var.equal = FALSE)
        }
      
      ints <- .SD[imputedFraction == FALSE, .(s = sum(intensity)), by = .(get(compare_between))]
      int1 <- max(0, mean(ints[get == samples[1]]$s), na.rm = TRUE)
      int2 <- max(0, mean(ints[get == samples[2]]$s), na.rm = TRUE)
      qint1 <- mean(qints[get == samples[1]]$s)
      qint2 <- mean(qints[get == samples[2]]$s)
      
      global_ints = .SD[, .(s = unique(global_intensity)), by = .(get(compare_between), Replicate)]
      global_ints_imp = .SD[, .(s = unique(global_intensity_imputed)), by = .(get(compare_between), Replicate)]
      global_int1 = mean(global_ints[get==samples[1]]$s)
      global_int2 = mean(global_ints[get==samples[2]]$s)
      global_int1_imp = mean(global_ints_imp[get==samples[1]]$s)
      global_int2_imp = mean(global_ints_imp[get==samples[2]]$s)
      #local_FC_all = log2(qints[get==samples[1]]$s/qints[get==samples[2]]$s)
      #global_FC_all = log2(global_ints_imp[get==samples[1]]$s/global_ints_imp[get==samples[2]]$s)
      #local_vs_global_FC_all = data.table(fc=c(local_FC_all,global_FC_all),sam=c(rep("local",length(local_FC_all)),rep("global",length(global_FC_all))))
      if (length(unique(design_matrix$Replicate)) > 1) {
        b = t.test(formula = log(global_ints_imp$s) ~ global_ints_imp$get, var.equal = FALSE) 
        global_pVal = b$p.value
        #c = t.test(formula = local_vs_global_FC_all$fc ~ local_vs_global_FC_all$sam , paired = F, var.equal = FALSE) 
        #local_vs_global_pVal = c$p.value
        meanDiff=a$estimate[1]-a$estimate[2]
      } else {
        global_pVal = 1
        #local_vs_global_pVal = 1
        meanDiff=a$estimate
      }
      .(pVal = a$p.value, 
        int1 = int1, int2 = int2, 
        meanDiff = meanDiff,
        qint1 = qint1, qint2 = qint2, log2FC =  log2(qint1/qint2),
        n_replicates = a$parameter + 1,  Tstat = a$statistic, testOrder = paste0(samples[1],".vs.",samples[2]),
        global_int1 = global_int1, global_int2 = global_int2, global_log2FC = log2(global_int1/global_int2),
        global_int1_imp = global_int1_imp, global_int2_imp = global_int2_imp, global_log2FC_imp = log2(global_int1_imp/global_int2_imp),
        #local_vs_global_log2FC = log2(qint1/qint2)-log2(global_int1/global_int2), local_vs_global_log2FC_imp = log2(qint1/qint2)-log2(global_int1_imp/global_int2_imp),
        global_pVal = global_pVal#, local_vs_global_pVal = local_vs_global_pVal
      )},
      by = .(id, feature_id, apex)]
    close(pb)
  }
##############################################################################################################
  tests[is.na(log2FC) & (int1 == 0 | int2  == 0) & (meanDiff == 0)]$log2FC <- 0
  tests[is.na(log2FC) & (int1 == 0 | int2  == 0) & (meanDiff > 0)]$log2FC <- Inf
  tests[is.na(log2FC) & (int1 == 0 | int2  == 0) & (meanDiff < 0)]$log2FC <- -Inf

  if ("proteoform_id" %in% names(featVals)) {
    proteoform_ann <- unique(subset(featVals,select=c("id","proteoform_id")))
    tests <- merge(tests,proteoform_ann,by=c("id"),all.x=T,all.y=F,sort=F)
  }

  if(level == "peptide"){
    tests$pBHadj <- p.adjust(tests$pVal, method = "BH")
    pQv <- qvalue::qvalue(tests$pVal, lambda = 0.4)
    tests$qVal <- pQv$qvalues
    if (length(unique(design_matrix$Replicate)) > 1) {
      tests$global_pBHadj <- p.adjust(tests$global_pVal, method = "BH")
      global_pQv <- qvalue::qvalue(tests$global_pVal, lambda = 0.4)
      tests$global_qVal <- global_pQv$qvalues
      #tests$local_vs_global_pBHadj <- p.adjust(tests$local_vs_global_pVal, method = "BH")
      #local_vs_global_pQv <- try(qvalue::qvalue(tests$local_vs_global_pVal, lambda = 0.4), silent = T)
      #if (is(local_vs_global_pQv, "try-error")) {
      #  tests$local_vs_global_qVal <- NA
      #} else {
      #  tests$local_vs_global_qVal <- local_vs_global_pQv$qvalues
      #}
    } else {
      tests$global_pBHadj <- 1
      tests$global_qVal <- 1
      #tests$local_vs_global_pBHadj <- 1
      #tests$local_vs_global_qVal <- 1
    }
    return(tests)
  } else if (level == "proteoform") {
    message("Aggregating to proteoform-level...")
    proteoformtests <- aggregatePeptideTestsToProteoform(tests)
    return(proteoformtests)
  } else if (level == "protein") {
    message("Aggregating to protein-level...")
    prottests <- aggregatePeptideTests(tests)
    return(prottests)
  } else if (level == "complex") {
    message("Aggregating to complex-level...")
    prottests <- aggregatePeptideTests(tests)
    complextests <- aggregateProteinTests(prottests)
    return(complextests)
  } else {
    stop("Specified level is not valid. Please chose between peptide, protein and complex.")
  }
}

