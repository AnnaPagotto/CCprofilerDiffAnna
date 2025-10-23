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
