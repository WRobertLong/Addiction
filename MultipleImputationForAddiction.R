require(foreign)
require(lme4)
require(mice)
require(rrp)
require(StatMatch)  # for hot deck
require(coda) # for mcmc post hoc fit assessment
require(scapeMCMC) #mcmc plots

############################
#
# Multiple Imputation for Addiction paper
# By WR Long 16 Jan 2013
# 
#
# Pre-amble:
# This script performs MI prior to
# the main GLMM modelling.
# 
#
#############################


rm(list=ls())

setwd("C:\\Users\\rob\\Documents\\Matt.Boden")

dt <- as.data.frame(read.spss("expectancies_impute_restricted3_for glmm_transposed.sav"))

data1.1 <- dt

diff <- lmer(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects + (time | Subject), data = data1.1)
summary(diff)
plot(resid(diff) ~ fitted(diff),main="residual plot for base model")
qqnorm(resid(diff),main="Normal QQ plot for base model");qqline(resid(diff))

cftest (diff)

diff1 <- lmer(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects + 
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
(time | Subject), data = data1.1)

summary(diff1)


data1.1$anxiety_current <- as.factor(data1.1$anxiety_current)
data1.1$dep_abu_current <- as.factor(data1.1$dep_abu_current)
data1.1$gender_1 <- as.factor(data1.1$gender_1)
data1.1$mood_current <- as.factor(data1.1$mood_current)
data1.1$race_1 <- as.factor(data1.1$race_1)

md.pattern(data1.1)

#ini <- mice(data1.1, maxit = 10)
#pred <- ini$pred

pred <- quickpred(data1.1)


pred["trans_mj_mean", ] <- c(1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


data1.1$Subject <- as.factor(data1.1$Subject)

data1.1$trans_mj_mean <- data1.1$trans_mj_mean * 7
data1.1$trans_mj_mean[222] <- 1
data1.1$trans_mj_mean[365] <- 2

imp.all <- mice(data1.1, 
	me = c(
		"", 			# Subject
		"", 			# time
		"pmm", 		# trans_mj_mean
		"", 			# MEQ_negativeeffects
		"", 			# MEQ_positiveeffects
		"pmm", 		# age_1
		"logreg", 		# race_1
		"", 			# gender_1  - none missing
		"", 			# dep_abu_current  - none missing
		"",			# anxiety_current - none missing
		"", 			# mood_current - none missing
		"pmm", 		# al_mean_use_28day
		"pmm", 		# ci_mean_use_28day
		"pmm", 		# al_mean_use_b
		"pmm", 		# ci_mean_use_b
		"pmm", 		# mj_mean_use_b 
		"pmm", 		# mq_1_transform
		""			# constant for 2l.norm method in mice
	), pred = pred, m = 5, maxit = 10, seed = 8000
)


fit1 <- with(imp, lmer(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects + (time | Subject)))
summary(pool(fit1))


##################

ini <- mice(data1.1, maxit = 10)
pred <- ini$pred

pred <- quickpred(data1.1)


#pred["time", ] <- c(-2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
#pred["trans_mj_mean", ] <- c(-2, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
pred["trans_mj_mean", ] <- c(-2, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)


imp <- mice(data1.1, 
	me = c(
		"", 			# Subject
		"", 		# time
		"2l.norm", 			# trans_mj_mean
		"", 
		"", 
		"pmm", 
		"logreg", 
		"", 
		"", 
		"",
		"", 
		"pmm", 
		"pmm", 
		"pmm", 
		"pmm", 
		"pmm", 
		"pmm", 
		""
), pred = pred, m = 1, maxit = 5, seed = 8000)


complete002 <- complete(imp,1)
save(complete002,file="complete002.RData")

foo.miss <- is.na(data1.1$trans_mj_mean)

foo.miss <- which(is.na(data1.1$trans_mj_mean))

data1.1$trans_mj_mean[foo.miss]
complete001$trans_mj_mean[foo.miss]
complete002$trans_mj_mean[foo.miss]
complete003$trans_mj_mean[foo.miss]

#########################################
#  short dataset
dt1 <- data1.1[c(1,2,3,4,5,18)]
str(dt1)
smallfit <- lmer(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects + (time | Subject), data = dt)

ini <- mice(dt1, maxit = 10)
predsmall <- ini$pred

#predsmall <- quickpred(dt1)

predsmall["trans_mj_mean", ] <- c(-2, 2, 0, 1, 1, 2)

imp <- mice(dt1, 
	me = c(
		"", 			# Subject
		"", 			# time
		"2l.norm", 		# trans_mj_mean
		"", 
		"", 
		""
), pred = predsmall, m = 1, maxit = 5, seed = 8000)


mice(dt1,  m = 1, maxit = 5, seed = 8000)

###################################
##  pmm imputation
##
##

dt1 <- data1.1[c(1,2,3,4,5,18)]
str(dt1)
dt1$trans_mj_mean <- dt1$trans_mj_mean * 7
dt1$trans_mj_mean[222] <- 1
dt1$trans_mj_mean[365] <- 2
dt1$Subject <- as.factor(dt1$Subject)

predsmall <- quickpred(dt1)

predsmall["trans_mj_mean", ] <- c(1, 1, 0, 1, 1, 0)

n.imps <- 20
n.iter <- 50

imp <- mice(dt1, 
	me = c(
		"", 			# Subject
		"", 			# time
		"pmm", 		# trans_mj_mean
		"", 
		"", 
		""
), pred = predsmall, m = n.imps, maxit = n.iter, seed = 8000)


f1 <- glmmadmb(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
	(time | Subject), data=na.omit(dt1),family="nbinom2",zeroInflation=TRUE)
summary(f1)

#################################
# test without random slopes
f1 <- glmmadmb(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
	(1 | Subject), data=na.omit(dt1),family="nbinom2",zeroInflation=TRUE)
summary(f1)

for (i in 1: n.imps) {
	f2 <- glmmadmb(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
		(1 | Subject), data=complete(imp,i),family="nbinom2",zeroInflation=TRUE)
	print(summary(f2))
}





#use all vars


for (i in 1: 20) {
	f2 <- glmmadmb(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
		(1 | Subject), data=complete(imp.all,i),family="nbinom2",zeroInflation=TRUE)
	print(summary(f2))
}

f2 <- glmmadmb(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
		(1 | Subject), data=complete(imp.all,3),family="nbinom2",zeroInflation=TRUE)
print(summary(f2))

library(bbmle)
AICtab(diff)


# end test without random slopes
#################################

fit1 <- with(imp, glmmadmb(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects + 
	(time | Subject),family="nbinom2",zeroInflation=TRUE))
summary(pool(fit1))


fit1 <- with(imp, lmer(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects + (time | Subject)))
summary(pool(fit1))



############################
#
#   Transformation attempts
#
############################

foo <- data1.1$trans_mj_mean
par(mfrow=c(4,4))
hist(foo)
hist(sqrt(foo)+1)



#####################
#
#  Zero-inflation modelling
#
######################

install.packages("glmmADMB", repos="http://r-forge.r-project.org")
library(glmmADMB) 

dt.zi <- data1.1

dt.zi$out <- dt.zi$trans_mj_mean*7
dt.zi$out[222] <- 1
dt.zi$out[365] <- 2

#dt.zi$time <- as.factor(dt.zi$time)
dt.zi$Subject <- as.factor(dt.zi$Subject)

dt.zi <- na.omit(dt.zi)


zi1 <- glmmadmb(out ~ time +  
	(time | Subject), data=dt.zi,family="binomial",zeroInflation=TRUE)

zi1 <- glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
	(time | Subject), data=dt.zi,family="nbinom1",zeroInflation=TRUE)

zi1 <- glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
	(time | Subject), data=dt.zi,family="nbinom2",zeroInflation=TRUE)


#######################################################
# try out hot deck imputation with rrp

############### example
set.seed(1)
key <- 1:100
## create random values
value1 <- 10 + 2 * key + rnorm(100, 0, 10)
## make 5 values into NAs
missing <- sample( key, 5)
value1[missing] <- NA
## build a dataframe
df <- data.frame(key, value1)
## do a nearest neighbor hot deck interpolation
imputed <- rrp.impute( df )$new.data

plot( df)
points(missing, imputed$value1[missing], col="red")
#######################################################

dt.zi <- data1.1

dt.zi$out <- dt.zi$trans_mj_mean*7
dt.zi$out[222] <- 1
dt.zi$out[365] <- 2

dt.zi$Subject <- as.factor(dt.zi$Subject)



imp.hd <- rrp.impute( dt.zi )$new.data


zi1 <- glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
	(time | Subject), data=imp.hd,family="nbinom2",zeroInflation=TRUE)


#######################################################
# try out hot deck imputation with StatMatch
#

# Example of Imputation of missing values
# introducing missing vales in iris

miss <- rbinom(nrow(iris), 1, 0.3)
ir.mat <- iris
ir.mat[miss==1,"Sepal.Length"] <- NA
iris.rec <- ir.mat[miss==1,-1]
iris.don <- ir.mat[miss==0,]
#search for NND donors
imp.NND <- RANDwNND.hotdeck(data.rec=iris.rec, data.don=iris.don,
	match.vars=c("Sepal.Width","Petal.Length", "Petal.Width"),
	don.class="Species")
# imputing missing values
iris.rec.imp <- create.fused(data.rec=iris.rec, data.don=iris.don,
mtc.ids=imp.NND$mtc.ids, z.vars="Sepal.Length")
# rebuild the imputed data.frame
final <- rbind(iris.rec.imp, iris.don)
head(final)

########

#Now for the real thing.
dt.zi <- data1.1

miss <- which(is.na(dt.zi$trans_mj_mean))

dt.zi$out <- dt.zi$trans_mj_mean*7
dt.zi$out[222] <- 1
dt.zi$out[365] <- 2

dt.zi$Subject <- as.factor(dt.zi$Subject)

dt.zi.rec <- dt.zi[is.na(dt.zi$out),]
dt.zi.rec$out <- NULL

dt.zi.don <- dt.zi[!is.na(dt.zi$out),]

final.imp <- NULL

n.imps <- 50
set.seed(201)
for (i in 1:n.imps) {
	#search for NND donors
	imp.NND <- RANDwNND.hotdeck(data.rec=dt.zi.rec, data.don=dt.zi.don,
		match.vars=c(
			"MEQ_negativeeffects", 
			"MEQ_positiveeffects",
			"age_1",
			"race_1",
			"gender_1",
			"dep_abu_current",
			"anxiety_current",
			"mood_current"
		),
		don.class=NULL
	)

	# imputing missing values
	dt.zi.imp <- create.fused(data.rec=dt.zi.rec, data.don=dt.zi.don,
		mtc.ids=imp.NND$mtc.ids, z.vars="out")
	# rebuild the imputed data.frame
	final.imp[[i]] <- as.data.frame(rbind(dt.zi.imp, dt.zi.don))

#	zi1 <- glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
#		(1 | Subject), data=final.imp[[i]],family="nbinom2",zeroInflation=TRUE)
#	cat("imputation:",i)
#	print(summary(zi1))

}

options(warn=2)
zi <- NULL
for (i in 1:n.imps) {
	zi[[i]] <- try(glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
		( 1 | Subject), data=final.imp[[i]],family="nbinom2",zeroInflation=TRUE))
	cat("imputation:",i)
	print(summary(zi1[i]))
}
options(warn=1)

use.pos <- NULL
for (i in 1:n.imps) {
	foo <- summary(zi[[i]])
	if (length(foo) > 3 )
		use.pos <- rbind(use.pos,foo$coefficients[6,])
}




# Apply Rubin's rules.

# number of estimates
n.est <- 6

bar  <- summary(zi[[1]])

est.names <- rownames(as.data.frame(bar$coefficients))
est.Estimate <- numeric(n.est)
est.StdError <- numeric(n.est)
for (i in 1:n.est) {
	int.pos <- NULL
	for (j in 1:n.imps) {
		foo <- summary(zi[[j]])
		if (length(foo) > 3 )
			int.pos <- rbind(int.pos,foo$coefficients[i,])
	}

	int.pos <- as.data.frame(int.pos)
	theta <- mean(int.pos$Estimate)
	U.hat <- mean(int.pos$`Std. Error`^2)
	B.hat <- sum((theta-int.pos$Estimate)^2)/(n.imps-1)
	sv <- U.hat + ( 1 + (1/n.imps))*B.hat
	se <- sqrt(sv)

	est.Estimate[i] <- round(theta,4)
	est.StdError[i] <- round(se,4)

}

as.data.frame(cbind(est.names, est.Estimate, est.StdError))



######
#
# big model

diff1 <- lmer(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects + 
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
(time | Subject), data = data1.1)

#
#
######



zi1 <- glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
	(1 | Subject), data=final,family="nbinom2",zeroInflation=TRUE)


zi1 <- glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
	(time | Subject), data=final,family="nbinom2",zeroInflation=TRUE)
summary(zi1)



zi1 <- glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
	(1 | Subject), data=na.omit(dt.zi),family="nbinom2",zeroInflation=TRUE)


####################
#
# Random slope neeeed ?

df1 <- data1.1[c(1,2,3,4,5,18)]

df1$trans_mj_mean <- df1$trans_mj_mean * 7
df1$trans_mj_mean[222] <- 1
df1$trans_mj_mean[365] <- 2
df1$Subject <- as.factor(dt1$Subject)

ff1 <- glmmadmb(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
		(1 | Subject), data=na.omit(df1),family="nbinom2",zeroInflation=TRUE)

ff2 <- glmmadmb(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
		(time | Subject), data=na.omit(df1),family="nbinom2",zeroInflation=TRUE)


library(bbmle) 
AICtab(ff1,ff2,diff)


#####################
#
#  Ignore errors (when using random slopes

dt.zi <- data1.1

miss <- which(is.na(dt.zi$trans_mj_mean))

dt.zi$out <- dt.zi$trans_mj_mean*7
dt.zi$out[222] <- 1
dt.zi$out[365] <- 2

dt.zi$Subject <- as.factor(dt.zi$Subject)

dt.zi.rec <- dt.zi[is.na(dt.zi$out),]
dt.zi.rec$out <- NULL

dt.zi.don <- dt.zi[!is.na(dt.zi$out),]

final.imp <- NULL

n.imps <- 50

# perform RHDI
for (i in 1:n.imps) {
	#search for NND donors
	imp.NND <- RANDwNND.hotdeck(data.rec=dt.zi.rec, data.don=dt.zi.don,
		match.vars=c(
			"MEQ_negativeeffects", 
			"MEQ_positiveeffects",
			"age_1",
			"race_1",
			"gender_1",
			"dep_abu_current",
			"anxiety_current",
			"mood_current"
		),
		don.class="time")

	# imputing missing values
	dt.zi.imp <- create.fused(data.rec=dt.zi.rec, data.don=dt.zi.don,
		mtc.ids=imp.NND$mtc.ids, z.vars="out")
	# rebuild the imputed data.frame
	final.imp[[i]] <- as.data.frame(rbind(dt.zi.imp, dt.zi.don))


}

zi <- NULL
for (i in 1:n.imps) {
	zi[[i]] <- try(glmmadmb(out ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
		(time | Subject), data=final.imp[[i]],family="nbinom2",zeroInflation=TRUE))
	cat("imputation:",i)
	print(summary(zi1[i]))
}



int.pos <- NULL
for (i in 1:n.imps) {
	foo <- summary(zi[[i]])
	if ( mode(summary(zi[[i]])) == "list" )
		int.pos <- rbind(int.pos,foo$coefficients[6,])
}




# Apply Rubin's rules.

int.pos <- as.data.frame(int.pos)
theta <- mean(int.pos$Estimate)
U.hat <- mean(int.pos$`Std. Error`^2)
B.hat <- sum((theta-int.pos$Estimate)^2)/(n.imps-1)
sv <- U.hat + ( 1 + (1/n.imps))*B.hat
se <- sqrt(sv)

# number of estimates
n.est <- 6

bar  <- summary(zi[[1]])

est.names <- rownames(as.data.frame(bar$coefficients))
est.Estimate <- numeric(n.est)
est.StdError <- numeric(n.est)
for (i in 1:n.est) {
	int.pos <- NULL
	for (j in 1:n.imps) {
		if ( mode(summary(zi[[j]])) == "list" ) {
			foo <- summary(zi[[j]])
			int.pos <- rbind(int.pos,foo$coefficients[i,])
		}
	}

	int.pos <- as.data.frame(int.pos)
	theta <- mean(int.pos$Estimate)
	U.hat <- mean(int.pos$`Std. Error`^2)
	B.hat <- sum((theta-int.pos$Estimate)^2)/(n.imps-1)
	sv <- U.hat + ( 1 + (1/n.imps))*B.hat
	se <- sqrt(sv)

	est.Estimate[i] <- round(theta ,4)
	est.StdError[i] <- round(se,4)

}

res.new <- as.data.frame(cbind(est.names,est.Estimate,est.StdError))

# end ignore errors
########################






# end random slope needed
########################

####
#
#   mcmc

fit1.mcmc <- glmmadmb(trans_mj_mean ~ time*MEQ_negativeeffects + time*MEQ_positiveeffects +  
	(1 | Subject), data=na.omit(dt1),family="nbinom2",zeroInflation=TRUE,
	mcmc=TRUE,mcmc.opts=mcmcControl(mcmc=1000),verbose=TRUE)

m <- as.mcmc(fit1.mcmc$mcmc)

plotTrace(m)
HPDinterval(m)

####

