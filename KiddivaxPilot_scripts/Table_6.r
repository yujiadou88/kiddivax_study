#
# R syntax to reproduce information for Table 6 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective efficacy of seasonal influenza vaccination against seasonal and
# pandemic influenza virus infection during 2009 in Hong Kong
# CID, 2010 (in press).
#
# Last updated by Fang VJ, Ng S and Cowling BJ.
# October 12, 2010

source("../kiddivaxPilot_scripts/dataframe_cross.r")

require(Hmisc)   

tab6 <- matrix(rep(NA, 4*11), nrow=11,
               dimnames=list(c("Age<16", "Age16-45", "Age>45", "Female", "Male", "No flu A", "Flu A",
                                "no vax", "seasonal vax", "Before Oct 1", "After Oct 1"),
                              c("n",  "OR",  "CI_low", "CI_up")))

#########################################################
#Multiple Imputation for missing index+member info      #
#########################################################
#redo agegp for regressopm
all$agegp<-1*(all$age<=45)+1*(all$age<=15)

set.seed(231)
all.i <- aregImpute( ~ factor(intervention) +
                         sh1w + sh1s + sh3w + sh3s + ph1s + ariw + iliw + aris + ilis + ilip + I(day.p) +
                         flo + male + I(agegp) + vac09, data=all,n.impute=10)

all.nomiss <- list(all, all, all, all, all, all, all, all, all, all)

for(i in 1:10){
    all.nomiss[[i]]$flo[is.na(all.nomiss[[i]]$flo)] <- all.i$imputed$flo[,i]
    all.nomiss[[i]]$sh1w[is.na(all.nomiss[[i]]$sh1w)] <- all.i$imputed$sh1w[,i]
    all.nomiss[[i]]$sh1s[is.na(all.nomiss[[i]]$sh1s)] <- all.i$imputed$sh1s[,i]
    all.nomiss[[i]]$sh3w[is.na(all.nomiss[[i]]$sh3w)] <- all.i$imputed$sh3w[,i]
    all.nomiss[[i]]$sh3s[is.na(all.nomiss[[i]]$sh3s)] <- all.i$imputed$sh3s[,i]
    all.nomiss[[i]]$ilip[is.na(all.nomiss[[i]]$ilip)] <- all.i$imputed$ilip[,i]
    all.nomiss[[i]]$ph1s[is.na(all.nomiss[[i]]$ph1s)] <- all.i$imputed$ph1s[,i]
    all.nomiss[[i]]$ariw[is.na(all.nomiss[[i]]$ariw)] <- all.i$imputed$ariw[,i]
    all.nomiss[[i]]$aris[is.na(all.nomiss[[i]]$aris)] <- all.i$imputed$aris[,i]
    all.nomiss[[i]]$iliw[is.na(all.nomiss[[i]]$iliw)] <- all.i$imputed$iliw[,i]
    all.nomiss[[i]]$ilis[is.na(all.nomiss[[i]]$ilis)] <- all.i$imputed$ilis[,i]
    all.nomiss[[i]]$vac09[is.na(all.nomiss[[i]]$vac09)] <- all.i$imputed$vac09[,i]
    all.nomiss[[i]]$day.p[is.na(all.nomiss[[i]]$day.p)] <- all.i$imputed$day.p[,i]
    }

# Logistic regression

require(rms)
risk.fit <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
temp <- all.nomiss[[i]]
temp$agegp <- factor(temp$agegp)
temp$vac <- (temp$sh1w==1|temp$sh1s==1|temp$sh3w==1|temp$sh3s==1)
temp$oct1 <- cut(temp$day.p,c(0,30,60))
risk.fit[[i]] <- lrm(ph1s~agegp+male+vac+vac09+oct1,data=temp)
}

combine.mi <- function(model, n.impute){
	betas <- matrix(c(model[[1]]$coef, model[[2]]$coef, model[[3]]$coef, model[[4]]$coef, model[[5]]$coef,
	                	model[[6]]$coef, model[[7]]$coef, model[[8]]$coef, model[[9]]$coef, model[[10]]$coef), byrow=FALSE, ncol=10)
	vars <- matrix(c(diag(vcov(model[[1]])), diag(vcov(model[[2]])), diag(vcov(model[[3]])), diag(vcov(model[[4]])),
	                 diag(vcov(model[[5]])), diag(vcov(model[[6]])), diag(vcov(model[[7]])), diag(vcov(model[[8]])),
	                 diag(vcov(model[[9]])), diag(vcov(model[[10]]))), byrow=FALSE, ncol=10)
	coef.names <- names(model[[1]]$coef)
	mean.coefs <- rowMeans(betas)
	Ubar <- rowMeans(vars)
	B <- rowSums((betas - mean.coefs)*(betas-mean.coefs) /
		(n.impute - 1))
	T <- (1 + 1/n.impute) * B + Ubar
	degf <- (n.impute - 1)*(1 + Ubar / ((1 + 1/n.impute)*B))*
		(1 + Ubar / ((1 + 1/n.impute)*B))
	output <- data.frame(OR = exp(mean.coefs),
		                   lowerCI = exp(mean.coefs - qt(0.975, df=degf)*sqrt(T)),
	                     upperCI = exp(mean.coefs + qt(0.975, df=degf)*sqrt(T)),
	                     row.names=coef.names)
  round(output,2)
}
out <- combine.mi(risk.fit,10)

tab6[c(3,4,6,8,10),2] <- 1.00
tab6[,1] <- c(rev(table(all$agegp)),table(all$male),table(all$sh1w==1|all$sh1s==1|all$sh3w==1|all$sh3s==1),table(all$vac09),table(cut(all$day.p,c(0,30,60))))
tab6[c(2,1,5,7,9,11),2:4] <- as.matrix(out[-1,])
tab6


#
# End of script.
#
