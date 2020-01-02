require(glmmADMB)
require(mice)
require(coda) # for mcmc post hoc fit assessment
require(scapeMCMC) #mcmc plots

############################
#
# GLMM analysis for Addition paper
# By WR Long 16 Jan 2013
# 
#
# Pre-amble:
# We use glmmADMB which is one of the only ways
# to fit a GLMM to zero-inflated data
#
# Note: glmmADMB allows post hoc MCMC estimation
#       to assess model fit.
#
#############################

# care needed
# rm(list=ls())

# load up the imputed datasets
load("final.imp.RData")

# loop over the complete datasets
# fitting a zero-inflated negative binomial
# mixed model (random intercepts for Subject)
# to each one


options(warn=2)
zi <- list(length(final.imp))
for (i in 1:length(final.imp)) {
	cat("\nimputation:",i,"\n")
	zi[[i]] <- try(glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
		( 1 | Subject), data=final.imp[[i]],family="nbinom2",zeroInflation=TRUE))
	# comment out the next line if you don't want to see the model output for each dataset
	# of course, it can be access later if needed
	print(summary(zi[[i]]))
}
options(warn=1)


use.pos <- NULL
for (i in 1:length(final.imp)) {
	foo <- summary(zi[[i]])
	if (length(foo) > 3 )
		use.pos <- rbind(use.pos,foo$coefficients[6,])
}


# Apply Rubin's rules.

# number of parameter estimates
n.est <- 6

# just used to get the names of the parameters
bar  <- summary(zi[[4]])
# note: this will fail if the chosen one happens to be one that is discarded !
# in that case, just change it to another one.

est.names <- rownames(as.data.frame(bar$coefficients))
est.Estimate <- numeric(n.est)
est.StdError <- numeric(n.est)

# loop over each estimated coefficient
for (i in 1:n.est) {
	int.pos <- NULL

	# for each estimate loop over the results from the usable datasets
	# and compile the results
	for (j in 1:length(final.imp)) {
		foo <- summary(zi[[j]])
		if (length(foo) > 3 )
			int.pos <- rbind(int.pos,foo$coefficients[i,])
	}
	int.pos <- as.data.frame(int.pos)

	theta <- mean(int.pos$Estimate)
	U.hat <- mean(int.pos$`Std. Error`^2)
	B.hat <- sum((theta-int.pos$Estimate)^2)/(nrow(use.pos)-1)
	sv <- U.hat + ( 1 + (1/nrow(use.pos)))*B.hat
	se <- sqrt(sv)

	est.Estimate[i] <- round(theta,4)
	est.StdError[i] <- round(se,4)
	est.PValue[i] <- round(2*pnorm(-abs(est.Estimate[i]/est.StdError[i])),4)

}

# store and output the result
(glmm.full.result <- as.data.frame(cbind(est.names, est.Estimate, est.StdError,est.PValue )))

cat("Used ", nrow(use.pos), " out of ", length(final.imp), " datasets\n")

###############################################
#
#  Check convergence with post hoc MCMC
#
#  Care is needed ! 50000 iterations will take a LONG time
#  Change mcmc=50000 is necessary
#
fit1.mcmc <- glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
		( 1 | Subject), data=final.imp[[4]],family="nbinom2",zeroInflation=TRUE,
		mcmc=TRUE,mcmc.opts=mcmcControl(mcmc=500000),verbose=FALSE)

m <- as.mcmc(fit1.mcmc$mcmc)

# trace plots
plotTrace(m[,2:7],layout=c(2,3))

# 95% credible intervals (highest probability density)
HPDinterval(m)[2:7,]

#
# End MCMC
#
##############################################

################
#
#  Repeat for covariate-adjusted model
#
#  Would be better to do both models in a loop
#  to avoid duplicating the code. The only thing that
#  changes is the number of parameter estimates !
#

options(warn=2)
zi <- list(length(final.imp))
for (i in 1:length(final.imp)) {
	cat("\nimputation:",i)
	zi[[i]] <- try(glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
		age_1 +
		race_1 +
		gender_1 +
		dep_abu_current +
		anxiety_current +
		mood_current +
		mq_1_transform +
		mj_mean_use_b +
		al_mean_use_b +
		ci_mean_use_b +
		al_mean_use_28day +
		ci_mean_use_28day +
		( 1 | Subject), data=final.imp[[i]],family="nbinom2",zeroInflation=TRUE))
	# comment out the next line if you don't want to see the model output for each dataset
	# of course, it can be access later if
	print(summary(zi[[i]]))
}
options(warn=1)


use.pos <- NULL
for (i in 1:length(final.imp)) {
	foo <- summary(zi[[i]])
	if (length(foo) > 3 )
		use.pos <- rbind(use.pos,foo$coefficients[6,])
}


# Apply Rubin's rules.

# number of parameter estimates
n.est <- 18

bar  <- summary(zi[[1]])

est.names <- rownames(as.data.frame(bar$coefficients))
est.Estimate <- numeric(n.est)
est.StdError <- numeric(n.est)
est.PValue <- numeric(n.est)
for (i in 1:n.est) {
	int.pos <- NULL
	for (j in 1:length(final.imp)) {
		foo <- summary(zi[[j]])
		if (length(foo) > 3 )
			int.pos <- rbind(int.pos,foo$coefficients[i,])
	}

	int.pos <- as.data.frame(int.pos)
	theta <- mean(int.pos$Estimate)
	U.hat <- mean(int.pos$`Std. Error`^2)
	B.hat <- sum((theta-int.pos$Estimate)^2)/(nrow(use.pos)-1)
	sv <- U.hat + ( 1 + (1/nrow(use.pos)))*B.hat
	se <- sqrt(sv)

	est.Estimate[i] <- round(theta,4)
	est.StdError[i] <- round(se,4)
	est.PValue[i] <- round(2*pnorm(-abs(est.Estimate[i]/est.StdError[i])),4)

}

(glmm.full.result <- as.data.frame(cbind(est.names, est.Estimate, est.StdError,est.PValue )))
cat("Used ", nrow(use.pos), " out of ", length(final.imp), " datasets\n")






