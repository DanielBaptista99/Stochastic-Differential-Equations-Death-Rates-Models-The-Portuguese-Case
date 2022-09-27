## Import Dataset##

require(readxl)
datadeath <- read_excel("Death_Rates_1X1.xlsx",sheet="Folha1",col_names = FALSE)
attach(datadeath)

# Organize Data

years <- datadeath[,1] #Time horizon: 1940 - 2020
Frates <- datadeath[,2:dim(datadeath)[2]] #Death rates by Age: columns 2-101 Female gender; columns 102-201 Male gender
                                                           
Fyears <- dim(Frates)[1]
Numrates <- dim(Frates)[2]

#Define time period for forecasts : 11 and 21 years
nfore <- 11 # time period : 2010 to 2020
nfore2 <- 21 # time period : 2010 to 2030
T <- datadeath[1:(Fyears-nfore),1]
finalT <- c(t(T),c(2010:2030)) #create vector for forecasts : 21 years 
ltf <- length(finalT)

#Select data for adjustments (leave 11 years out to estimate Mean Squared Error)

Rates <- datadeath[1:(Fyears-nfore),2:201]
RowRates <- dim(Rates)[1]
ColRates <- dim(Rates)[2]
n <- RowRates-1

#Identify 11 years of data not considered for adjustment

Tfore <- years[RowRates:(RowRates+nfore),]
fore.rates <- Frates[RowRates:(RowRates+nfore),]


##GEOMETRIC BROWNIAN MOTION: ESTIMATES OF THE PARAMETERS OF THE UNIVARIATE MODEL, BY AGE , BY GENDER##

#apply transformation , using logarithmic scale, to death rates#

YF <- matrix(nrow = (RowRates + nfore), ncol = ColRates)
for(p in 1:(RowRates + nfore)){
for(q in 1:ColRates){
YF[[p,q]] <- log(Frates[[p,q]] / Frates[[1,q]])
}
}

#Select transformed data to adjust model

Y <- YF[1:(Fyears-nfore),]

#R = mu - (sigma^2) / 2

R <- vector(length=ColRates)

for (i in 1:ColRates){
R[i] <- Y[RowRates,i] / RowRates
}

#V = sigma^2

V <- vector(length = ColRates)

for(i in 1:ColRates){
V[i] <- (1 / n)*sum((Y[(2:RowRates),i] - Y[(1:(RowRates-1)),i] - R[[i]])^2)
}

#Estimate Parameter mu of the GBM

mu <- vector(length = ColRates)
for (i in 1:ColRates){
mu[i] <- R[[i]] + V[i] / 2
}

#Estimate Parameter sigma of the GBM

sig <- vector(length = ColRates)
for (i in 1:ColRates){
sig[i] <- sqrt(V[i])
}

#create table with all the parameters of the GBM

p.GBM <- cbind(R, V, mu, sig)

##GBM : ADJUSTMENT##

#Adjustment Long Term : Estimate Death rates in logarithmic scale

Y.adjust.GBM <- matrix(nrow = RowRates, ncol = ColRates)
for (p in 1:RowRates){
for (q in 1:ColRates){
Y.adjust.GBM[[1,q]] <- Y[[1,q]]
}
}
for (p in 2:RowRates){
for (q in 1:ColRates){
Y.adjust.GBM[[p,q]] <- Y[[1,q]] + R[[q]] * (T[[p,1]] - T[[1,1]])
}
}

#Invert transformation Y to adjusted values in the Long Term

Rate.adjust.GBM <- matrix(nrow = RowRates, ncol = ColRates)
for (p in 1:RowRates){
for (q in 1:ColRates){
Rate.adjust.GBM[[p,q]] <- Rates[[1,q]] * exp(Y.adjust.GBM[[p,q]])
}
}

# Estimate Squared Mean Error of the long term adjustment

MSE.adjust.GBM <- vector(length = RowRates)
error.adjust.GBM <- matrix(nrow = (RowRates + 1), ncol = ColRates)
for (p in 1:RowRates){
for (q in 1:ColRates){
error.adjust.GBM[[p,q]] <- ((Rate.adjust.GBM[[p,q]] - Rates[[p,q]])^2)
}
}

for (p in 1:RowRates){
for (q in 1:ColRates){
error.adjust.GBM[(RowRates+1),q] <- sum(error.adjust.GBM[(1:p),q]) / RowRates
}
}

MSE.adjust.GBM <- error.adjust.GBM[(RowRates+1),]

#Adjustment Step by Step : estimate death rates in logarithmic scale

Y.adjust.GBM.SS <- matrix(nrow = RowRates, ncol = ColRates)
for (p in 1:RowRates){
for (q in 1:ColRates){
Y.adjust.GBM.SS[[1,q]] <- Y[[1,q]]
}
}

for (p in 2:RowRates){
for (q in 1:ColRates){
Y.adjust.GBM.SS[[p,q]] <- Y[[p-1,q]] + R[[q]] * (T[[p,1]] - T[[p-1,1]])
}
}

#Invert transformation Y to adjusted values SS

Rate.adjust.GBM.SS <- matrix(nrow = RowRates, ncol = ColRates)
for (p in 1:RowRates){
for (q in 1:ColRates){
Rate.adjust.GBM.SS[[p,q]] <- Rates[[1,q]] * exp(Y.adjust.GBM.SS[[p,q]])
}
}

# Estimate Squared Mean Error of the Step by Step Adjustment

MSE.adjust.GBM.SS <- vector(length = RowRates)
error.adjust.GBM.SS <- matrix(nrow = (RowRates + 1), ncol = ColRates)
for (p in 1:RowRates){
for (q in 1:ColRates){
error.adjust.GBM.SS[[p,q]] <- ((Rate.adjust.GBM.SS[[p,q]] - Rates[[p,q]])^2)
}
}

for (p in 1:RowRates){
for (q in 1:ColRates){
error.adjust.GBM.SS[(RowRates+1),q] <- sum(error.adjust.GBM.SS[(1:p),q]) / RowRates
}
}
MSE.adjust.GBM.SS <- error.adjust.GBM.SS[(RowRates+1),]

##GBM: Forecasts##

#Long Term Forecasts : 11 years

Y.fore.LT.GBM <- matrix(nrow = (nfore + 1), ncol = ColRates)
for (p in 1:(nfore+1)){
for (q in 1:ColRates){
Y.fore.LT.GBM[[1,q]] <- YF[[RowRates,q]]
}
}
for (p in 2:(nfore+1)){
for (q in 1:ColRates){
Y.fore.LT.GBM[[p,q]] <- YF[[RowRates,q]] + R[[q]] * (years[[(RowRates + p - 1),1]] - years[[RowRates,1]])
}
}

#Invert transformation Y to long term forecasts : 11 years

Rate.fore.LT.GBM <- matrix(nrow = (nfore + 1), ncol = ColRates)
for (p in 1:(nfore+1)){
for (q in 1:ColRates){
Rate.fore.LT.GBM[[p,q]] <- Frates[[1,q]] * exp(Y.fore.LT.GBM[[p,q]])
}
}

#Estimate Mean Squared Error of the Long term forecasts : 11 years

MSE.fore.LT.GBM <- vector(length = (nfore + 1))
error.fore.LT.GBM <- matrix(nrow = (nfore + 2), ncol = ColRates)
for(p in 1:(nfore+1)){
for(q in 1:ColRates){
error.fore.LT.GBM[[p,q]] <- ((Rate.fore.LT.GBM[[p,q]] - fore.rates[[p,q]])^2)
}
}

for (p in 1:(nfore+1)){
for (q in 1:ColRates){
error.fore.LT.GBM[(nfore+2),q] <- sum(error.fore.LT.GBM[(1:p),q]) / (nfore + 1)
}
}

MSE.fore.LT.GBM <- error.fore.LT.GBM[(nfore + 2),]

#Long Term Forecasts : 21 years

Y.fore.LT21.GBM <- matrix(nrow=(nfore2+1),ncol = ColRates)
for (p in 1:1){
for (q in 1:ColRates){
  Y.fore.LT21.GBM[[1,q]] <- YF[[RowRates,q]]
}
}

for (p in 2:(nfore2+1)){
for (q in 1:ColRates){
Y.fore.LT21.GBM[[p,q]] <- YF[[RowRates,q]]+R[[q]]*(finalT[[RowRates+p-1]]-finalT[[RowRates]]) 
}
}

#Invert transformation Y to long term forecasts : 21 years

Rate.fore.LT21.GBM <- matrix(nrow = (nfore2 + 1), ncol = ColRates)
for (p in 1:(nfore2+1)){
for (q in 1:ColRates){
Rate.fore.LT21.GBM[p,q] <- Frates[[1,q]] * exp(Y.fore.LT21.GBM[[p,q]])
}
}

#Forecasts Step by Step : 11 years

f.n <- (dim(Frates)[1]) - 1
k <- f.n - nfore

f.YF <- YF

f.R <- matrix(nrow = nfore, ncol = ColRates)
f.V <- matrix(nrow = nfore, ncol = ColRates)
f.mu <- matrix(nrow = nfore, ncol = ColRates)
f.sig <- matrix(nrow = nfore, ncol = ColRates)

Y.fore.1S.GBM <- matrix(nrow = (nfore + 1), ncol = ColRates)
for(p in 1:(nfore+1)){
for(q in 1:ColRates){
Y.fore.1S.GBM[[1,q]] <- YF[[(k+1),q]]
}
}


for (i in 1:ColRates){
  #update annually forecasts
for (j in 1:(nfore)){
f.R[j,i] <- f.YF[(k+j),i] / (k+j)
f.V[j,i] <- (1 / (k + j)) * sum((f.YF[(2:(k+j)),i] - f.YF[(1:(k+j-1)),i] - f.R[j,i])^2)
f.mu[j,i] <- f.R[j,i] + f.V[j,i] / 2
f.sig[j,i] <- sqrt(f.V[j,i])

#forecasts SS

Y.fore.1S.GBM[(j+1),i] <- YF[(k+j),i] + f.R[j,i]

}
}

#Invert transformation Y to SS forecasts : 11 years

Rate.fore.1S.GBM <- matrix(nrow = (nfore+1),ncol = ColRates)
for (p in 1:(nfore+1)){
for (q in 1:ColRates){
Rate.fore.1S.GBM[[p,q]] <- Frates[[1,q]] * exp(Y.fore.1S.GBM[[p,q]])
}
}

#Estimate Squared Mean Error of the SS forecasts : 11 years

MSE.fore.1S.GBM <- vector(length = (nfore + 1))
error.fore.1S.GBM <- matrix(nrow = (nfore + 2), ncol = ColRates)
for (p in 1:(nfore+1)){
for (q in 1:ColRates){
error.fore.1S.GBM[[p,q]] <- ((Rate.fore.1S.GBM[[p,q]] - fore.rates[[p,q]])^2)
}
}

for (p in 1:(nfore+1)){
for (q in 1:ColRates){
error.fore.1S.GBM[[(nfore+2),q]] <- sum(error.fore.1S.GBM[(1:p),q]) / (nfore + 1)
}
}

MSE.fore.1S.GBM <- error.fore.1S.GBM[(nfore+2),]

##STOCHASTIC GOMPERTZ MODEL : ESTIMATES OF THE PARAMETERS OF THE UNIVARIATE MODEL , BY AGE , BY GENDER ##

#Apply logarithmic transformation to the death rates

YF <- log(Frates)

#Select transformed data to adjust model

Y <- YF[1:(Fyears-nfore),]

#Create Matrix with 1 year gap

Yt <- matrix(nrow = n, ncol = ColRates)

for (p in 1:n){
for (q in 1:ColRates){
Yt[[p,q]] <- Y[[(p+1),q]]
}
}

Yt_1 <- matrix(nrow = n, ncol = ColRates)

for (p in 1:n){
for (q in 1:ColRates){
Yt_1[[p,q]] <- Y[[p,q]]
}
}

#Estimate parameters : A = ln(a) , b and sigma

parameters.SGM <- vector(length = ColRates)

library(stats)

for (i in 1:ColRates){
f.A <- function(beta){
((n * (1 - exp(-beta)))^(-1)) * (sum(Yt[,i] - Yt_1[,i] * exp(-beta)))
}

f.sigma <- function(Alpha,beta){
(((2*beta) * ((n * (1 - exp(-2 * beta)))^(-1))) *
    sum((Yt[,i] - Alpha - (Yt_1[,i] - Alpha) * exp(-beta))^2))^(1/2)
}

f.LogV <- function(beta){
((n / 2) * log(2 * pi)) + ((n / 2) * log(((f.sigma(f.A(beta), beta))^2) / (2 * beta))) +
       ((n / 2) * log(1 - exp(-2 * beta))) + ((beta * ((((f.sigma(f.A(beta), beta))^2) *
             (1 - exp(-2 * beta)))^(-1))) * (sum((Yt[,i] - f.A(beta) - (Yt_1[,i] - f.A(beta)) * exp(-beta))^2)))
}

parameters.SGM[i] <- optimize(f.LogV, int = c(0.0001, 10))$min
}

A <- vector(length = ColRates)

b <- vector(length = ColRates)

sigma <- vector(length = ColRates)

a <- vector(length = ColRates)

for (i in 1:ColRates){
b[i] <- parameters.SGM[i]
A[i] <- f.A(b[i])
sigma[i] <- f.sigma(f.A(b[i]), b[i])
a[i] <- exp(A[i])
}

#Create table with estimated parameters

p.SGM <- cbind(A, b, sigma, a)

##SGM : ADJUSTMENT##

#Long term adjustment : estimate Logarithms of the Death Rates

Y.adjust.SGM <- matrix(nrow = RowRates, ncol = ColRates)
for (p in 1:RowRates){
for (q in 1:ColRates){
Y.adjust.SGM[[1,q]] <- Y[[1,q]]
}
}

for (p in 2:RowRates){
for (q in 1:ColRates){
Y.adjust.SGM[[p,q]] <- A[q] + (Y[[1,q]] - A[q]) * exp(-b[q] * (T[[p,1]] - T[[1,1]]))
}
}

#Invert transformation Y for long term adjusted values

Rate.adjust.SGM <- matrix(nrow = RowRates, ncol = ColRates)
for (p in 1:RowRates){
for (q in 1:ColRates){
Rate.adjust.SGM[[p,q]] <- exp(Y.adjust.SGM[[p,q]])
}
}


#Estimate Mean Squared Error of the Long term Adjustment

MSE.adjust.SGM <- vector(length = RowRates)
error.adjust.SGM <- matrix(nrow = (RowRates + 1), ncol = ColRates)
for(p in 1:RowRates){
for(q in 1:ColRates){
error.adjust.SGM[[p,q]] <- ((Rate.adjust.SGM[[p,q]] - Rates[[p,q]])^2)
}
}

for(p in 1:RowRates){
for(q in 1:ColRates){
error.adjust.SGM[[(RowRates+1),q]] <- sum(error.adjust.SGM[(1:p),q]) / RowRates
}
}

MSE.adjust.SGM <- error.adjust.SGM[(RowRates+1),]

#Adjustment SS: estimate Logarithms of the Death Rates

Y.adjust.SGM.SS <- matrix(nrow = RowRates, ncol = ColRates)
for(p in 1:RowRates){
for(q in 1:ColRates){
Y.adjust.SGM.SS[[1,q]] <- Y[[1,q]]
}
}

for (p in 2:RowRates){
for (q in 1:ColRates){
Y.adjust.SGM.SS[[p,q]] <- A[q] + (Y[[p-1,q]] - A[q]) * exp(-b[q] * (T[[p,1]] - T[[p-1,1]]))
}
}

#Invert transformation Y to SS adjusted values

Rate.adjust.SGM.SS <- matrix(nrow = RowRates, ncol = ColRates)
for(p in 1:RowRates){
for(q in 1:ColRates){
Rate.adjust.SGM.SS[[p,q]] <- exp(Y.adjust.SGM.SS[[p,q]])
}
}

#Estimate Mean Squared Error of the SS adjustment

MSE.adjust.SGM.SS <- vector(length = RowRates)
error.adjust.SGM.SS <- matrix(nrow = (RowRates + 1), ncol = ColRates)
for (p in 1:RowRates){
for (q in 1:ColRates){
error.adjust.SGM.SS[[p,q]] <- ((Rate.adjust.SGM.SS[[p,q]] - Rates[[p,q]])^2)
}
}

for(p in 1:RowRates){
for(q in 1:ColRates){
error.adjust.SGM.SS[[(RowRates+1),q]] <- sum(error.adjust.SGM.SS[(1:p),q]) / RowRates
}
}

MSE.adjust.SGM.SS <- error.adjust.SGM.SS[(RowRates+1),]


##SGM : FORECASTS##

#Forecasts Long Term : 11 Years

Y.fore.LT.SGM <- matrix(nrow = (nfore + 1), ncol = ColRates)
for(p in 1:(nfore+1)){
for(q in 1:ColRates){
Y.fore.LT.SGM[[1,q]] <- YF[[RowRates,q]]
}
}

for(p in 2:(nfore+1)){
for(q in 1:ColRates){
Y.fore.LT.SGM[[p,q]] <- A[q] + (YF[[RowRates,q]] - A[q]) * exp(-b[q] * (years[[(RowRates + p - 1),1]] - years[[RowRates,1]]))
}
}

#Invert transformation Y to long term forecasts : 11 years

Rate.fore.LT.SGM <- matrix(nrow = (nfore+1), ncol = ColRates)
for (p in 1:(nfore+1)){
for (q in 1:ColRates){
Rate.fore.LT.SGM[[p,q]] <- exp(Y.fore.LT.SGM[[p,q]])
}
}

#Estimate Mean Squared Error of long term forecasts : 11 years 

MSE.fore.LT.SGM <- vector(length = (nfore+1))
error.fore.LT.SGM <- matrix(nrow = (nfore + 2), ncol = ColRates)
for(p in 1:(nfore+1)){
for(q in 1:ColRates){
error.fore.LT.SGM[[p,q]] <- ((Rate.fore.LT.SGM[[p,q]] - fore.rates[[p,q]])^2)
}
}

for(p in 1:(nfore+1)){
for(q in 1:ColRates){
error.fore.LT.SGM[[(nfore+2),q]] <- sum(error.fore.LT.SGM[(1:p),q]) / (nfore)
}
}

MSE.fore.LT.SGM <- error.fore.LT.SGM[(nfore+2),]

#long term forecasts : 21 years 

Y.fore.LT21.SGM <- matrix(nrow = (nfore2+1), ncol = ColRates)
for (p in 1:1){
for (q in 1:ColRates){
Y.fore.LT21.SGM[[1,q]] <- YF[[RowRates,q]]
}
}

for(p in 2:(nfore2+1)){
for(q in 1:ColRates){
Y.fore.LT21.SGM[[p,q]] <- A[q] + (YF[[RowRates,q]] - A[q]) * exp(-b[q] * (finalT[RowRates+p-1] - finalT[RowRates])) 
}
}

#Invert transformation Y for forecasts : 21 years

Rate.fore.LT21.SGM <- matrix(nrow = (nfore2 + 1), ncol = ColRates)
for (p in 1:(nfore2+1)){
for (q in 1:ColRates){
Rate.fore.LT21.SGM[[p,q]] <- exp(Y.fore.LT21.SGM[[p,q]])
}
}

#Forecasts Step by Step : 11 years

f.n <- (dim(YF)[1]) - 1

f.Yt <- matrix(nrow = f.n, ncol = ColRates)
for(p in 1:f.n){
for(q in 1:ColRates){
f.Yt[[p,q]] <- YF[[(p+1),q]]
}
}

f.Yt_1 <- matrix(nrow = f.n, ncol = ColRates)
for(p in 1:f.n){
for(q in 1:ColRates){
f.Yt_1[[p,q]] <- YF[[p,q]]
}
}

k <- f.n - nfore

f.b <- matrix(nrow = nfore, ncol = ColRates)
f.A <- matrix(nrow = nfore, ncol = ColRates)
f.sigma <- matrix(nrow = nfore, ncol = ColRates)
f.a <- matrix(nrow = nfore, ncol = ColRates)

f.parameter.SGM <- matrix(nrow = nfore, ncol = ColRates)

f.Y.fore.1S.SGM <- matrix(nrow = (nfore + 1), ncol = ColRates)
for(p in 1:(nfore+1)){
for(q in 1:ColRates){
f.Y.fore.1S.SGM[[1,q]] <- YF[[k+1,q]]
}
}

for(i in 1:ColRates){
   #Update annually forecasts
     for (j in 1:(nfore)){
       
       SS.f.A <- function(SS.beta){
        (((k + j) * (1 - exp(-SS.beta)))^(-1)) * (sum(f.Yt[1:(k+j),i] -
           f.Yt_1[1:(k+j),i] * exp(-SS.beta)))
     }
  
       SS.f.sigma <- function(SS.Alpha,SS.beta){
      (((2*SS.beta) * (((k + j) * (1 - exp(-2 * SS.beta)))^(-1))) *
       sum((f.Yt[1:(k+j),i] - SS.Alpha - (f.Yt_1[1:(k+j),i] -
        SS.Alpha) * exp(-SS.beta))^2))^(1/2)
       }
       
       SS.fLogV <- function(SS.beta){
       (((k + j) / 2) * log(2 * pi)) + (((k + j) / 2) *
        log(((SS.f.sigma(SS.f.A(SS.beta), SS.beta))^2) / (2 * SS.beta))) +
         (((k + j) / 2) * log(1 - exp(-2 * SS.beta))) +
            ((SS.beta * ((((SS.f.sigma(SS.f.A(SS.beta), SS.beta))^2) *
           (1 - exp(-2 * SS.beta)))^(-1))) * (sum((f.Yt[1:(k+j),i] -
          SS.f.A(SS.beta) - (f.Yt_1[1:(k+j),i] - SS.f.A(SS.beta)) * exp(-SS.beta))^2)))
       }
       
       f.parameter.SGM[j,i] <- optimize(SS.fLogV, int=c(0.0001,10))$min
       
       #Estimates of parameters: A=Ln(a), b e sigma, para c=0 e Y=Ln(X)
       
       f.b[j,i] <- f.parameter.SGM[j,i]
       f.A[j,i] <- SS.f.A(f.b[j,i])
       f.sigma[j,i] <- SS.f.sigma(SS.f.A(f.b[j,i]), f.b[j,i])
       f.a[j,i] <- exp(f.A[j,i])
       
       #Forecasts SS
       
       f.Y.fore.1S.SGM[[(j+1),i]] <- f.A[j,i] + (YF[[(k+j),i]] - f.A[j,i]) * exp(-f.b[j,i])
       
}
}

#Invert transformation Y for SS forecasts : 11 years

Rate.fore.1S.SGM <- matrix(nrow = (nfore+1), ncol = ColRates)
for (p in 1:(nfore+1)){
for (q in 1:ColRates){
Rate.fore.1S.SGM[[p,q]] <- exp(f.Y.fore.1S.SGM[[p,q]])
}
}

#Estimate Mean Squared Error of the SS forecasts : 11 years

MSE.fore.1S.SGM <- vector(length = (nfore + 1))
error.fore.1S.SGM <- matrix(nrow = (nfore+2), ncol = ColRates)
for (p in 1:(nfore+1)){
for (q in 1:ColRates){
error.fore.1S.SGM[[p,q]] <- ((Rate.fore.1S.SGM[[p,q]] - fore.rates[[p,q]])^2)
}
}

for (p in 1:(nfore+1)){
for (q in 1:ColRates){
error.fore.1S.SGM[(nfore+2),q] <- sum(error.fore.1S.SGM[(1:p),q]) / (nfore)
}
}

MSE.fore.1S.SGM <- error.fore.1S.SGM[(nfore+2),]


## PLOTS ##

##PAGE 2##

par(mfrow=c(1,2))
plot(0:99,Frates[55,1:100],type="l",xlab="Age",ylab="Death Rate",cex.lab = 1.5,cex.axis = 1.5)
plot(0:99,Frates[55,101:200],type="l",xlab="Age",ylab="Death Rate",cex.lab = 1.5,cex.axis = 1.5)

fem66 <- Frates[1:81,67]
mal66 <- Frates[1:81,167]
par(mfrow=c(1,2))
plot(1940:2020,t(fem66[1:81,1]),type="l",ylab="Death Rate",xlab="Year",cex.lab = 1.5,cex.axis = 1.5)
plot(1940:2020,t(mal66[1:81,1]),type="l",ylab="Death Rate",xlab="Year",cex.lab = 1.5,cex.axis = 1.5)

##PAGE 3##

par(mfrow=c(1,1))
forcemortalilty <- function(x){-log(1-x)}
fm <- forcemortalilty(t(Frates[1:81,124]))
dat <- cbind((Frates[1:81,124]),t(fm))
plot(1940:2020,dat[,1],type="l",xlab="Year",ylab="Death rates and µ values") 
lines(1940:2020,dat[,2],col="grey")
legend("topright",c("Observed death rates","µ"),pch="_",col = c("black","grey"),cex=1.2)

##PAGE 18##

limCIinf95R <- function(R,V){R-qnorm(0.975)*sqrt(V/81)} 
limCIsup95R <- function(R,V){R+qnorm(0.975)*sqrt(V/81)} 
limCIeinf95R <- function(R,V){R-qt(0.975,80)*sqrt(81/80*V/81)} 
limCIesup95R <- function(R,V){R+qt(0.975,80)*sqrt(81/80*V/81)}

CIinf95RF <- t(limCIinf95R(t(p.GBM[1:100,1]),t(p.GBM[1:100,2])))   #F
CIsup95RF <- t(limCIsup95R(t(p.GBM[1:100,1]),t(p.GBM[1:100,2])))
CIeinf95RF <- t(limCIeinf95R(t(p.GBM[1:100,1]),t(p.GBM[1:100,2])))
CIesup95RF <- t(limCIesup95R(t(p.GBM[1:100,1]),t(p.GBM[1:100,2])))

CIinf95RM <- t(limCIinf95R(t(p.GBM[101:200,1]),t(p.GBM[101:200,2])))   #M
CIsup95RM <- t(limCIsup95R(t(p.GBM[101:200,1]),t(p.GBM[101:200,2])))
CIeinf95RM <- t(limCIeinf95R(t(p.GBM[101:200,1]),t(p.GBM[101:200,2])))
CIesup95RM <- t(limCIesup95R(t(p.GBM[101:200,1]),t(p.GBM[101:200,2])))

RF <- p.GBM[1:100,1]
RM <- p.GBM[101:200,1]
Rfemale <- cbind(CIinf95RF,CIeinf95RF,RF,CIesup95RF,CIsup95RF)
Rmale <- cbind(CIinf95RM,CIeinf95RM,RM,CIesup95RM,CIsup95RM)

par(mfrow = c(1,2))

plot(0:99,Rfemale[,1],type="l",xlab = "Age",ylab = "R",col = "grey",ylim=c(-0.15,0.15),cex.lab=1.5,cex.axis=1.5)
lines(0:99,Rfemale[,2],col="black",lty=5)
lines(0:99,Rfemale[,3],col="black")
lines(0:99,Rfemale[,4],col="black",lty=5)
lines(0:99,Rfemale[,5],col="grey")
legend("topleft",c("Estimate of R","CI 95%","Exact CI 95%"),col=c("black","black","grey"),lty = c(1,5,1),cex=1.5)

plot(0:99,Rmale[,1],type="l",xlab = "Age",ylab = "R",col = "grey",ylim = c(-0.15,0.15),cex.lab=1.5,cex.axis=1.5)
lines(0:99,Rmale[,2],col="black",lty=5)
lines(0:99,Rmale[,3],col="black")
lines(0:99,Rmale[,4],col="black",lty=5)
lines(0:99,Rmale[,5],col="grey")
legend("topleft",c("Estimate of R","CI 95%","Exact CI 95%"),col=c("black","black","grey"),lty = c(1,5,1),cex=1.5)

limCIinf95V <- function(V){V - qnorm(0.975)*sqrt((2*V^2)/81)}
limCIsup95V <- function(V){V + qnorm(0.975)*sqrt((2*V^2)/81)}
limCIeinf95V <- function(V){(81*V)/(qchisq(0.975,80))}
limCIesup95V <- function(V){(81*V)/(qchisq(0.025,80))}

CIinf95VF <-  t(limCIinf95V(t(p.GBM[1:100,2])))  
CIsup95VF <-  t(limCIsup95V(t(p.GBM[1:100,2])))  
CIeinf95VF <- t(limCIeinf95V(t(p.GBM[1:100,2])))  
CIesup95VF <- t(limCIesup95V(t(p.GBM[1:100,2]))) 

CIinf95VM <-  t(limCIinf95V(t(p.GBM[101:200,2])))  
CIsup95VM <-  t(limCIsup95V(t(p.GBM[101:200,2])))  
CIeinf95VM <- t(limCIeinf95V(t(p.GBM[101:200,2])))  
CIesup95VM <- t(limCIesup95V(t(p.GBM[101:200,2]))) 

VF <- p.GBM[1:100,2]
VM <- p.GBM[101:200,2]
Vfemale <- cbind(CIinf95VF,CIeinf95VF,VF,CIesup95VF,CIsup95VF)
Vmale <- cbind(CIinf95VM,CIeinf95VM,VM,CIesup95VM,CIsup95VM)

par(mfrow = c(1,2))

plot(0:99,Vfemale[,1],type="l",xlab = "Age",ylab = "V",col = "grey",ylim = c(0,0.2),cex.lab=1.5,cex.axis=1.5)
lines(0:99,Vfemale[,2],col="black",lty=5)
lines(0:99,Vfemale[,3],col="black")
lines(0:99,Vfemale[,4],col="black",lty=5)
lines(0:99,Vfemale[,5],col="grey")
legend("topright",c("Estimate of V","CI 95%","Exact CI 95%"),col=c("black","black","grey"),lty = c(1,5,1),cex=1.5)

plot(0:99,Vmale[,1],type="l",xlab = "Age",ylab = "V",col = "grey",ylim = c(0,0.2),cex.lab=1.5,cex.axis=1.5)
lines(0:99,Vmale[,2],col="black",lty=5)
lines(0:99,Vmale[,3],col="black")
lines(0:99,Vmale[,4],col="black",lty=5)
lines(0:99,Vmale[,5],col="grey")
legend("top",c("Estimate of V","CI 95%","Exact CI 95%"),col=c("black","black","grey"),lty = c(1,5,1),cex=1.5)

##PAGE 19##

Age15Mrate <- c(t(Frates[,116]),rep(NA,11))
Age15Mratefore21 <- c(t(Rate.adjust.GBM[,116]),t(Rate.fore.LT21.GBM[,116]))
rates15 <- cbind(Age15Mrate,Age15Mratefore21)
age15MSS <- cbind(Frates[70:81,116],Rate.fore.1S.GBM[1:12,116])

#Simulation 1940 - 2020

par(mfrow=c(1,2))

require(sde)

set.seed(123)

mu15 <- as.numeric(p.GBM[116,3])
sig15 <- as.numeric(p.GBM[116,4])
x0 <- as.numeric(Frates[1,116])
TIME <- as.numeric(Fyears)-1
r <- 2000
byn <- as.numeric(Fyears)-1

#Simulate r trajectories

delta_t <- TIME/byn
time <- seq(0,TIME,by=delta_t)

G <- matrix(rep(0,length(time)*r),nrow=r)

for (i in 1:r){G[i,] = GBM(x=x0,r=mu15,sigma = sig15,T=TIME,N=byn)}

plot(1940:2020,G[1,],t='l', col="grey", ylab="Death Rates",xlab="Year",ylim=c(0,0.1),cex.lab=1.5,cex.axis=1.5) 
for(i in 2:r){lines(1940:2020,G[i,], t='l',col="grey")}
lines(1940:2020,t(Frates[1:81,116]),t="l",col="black")
legend("topright",c("Observed death rates"),pch="_",col="black",cex=1.5)

# Simulation 2009 - 2020

set.seed(12345)
x2009 <- as.numeric(Frates[70,116])
TIME2 <- nfore
byn2 <- nfore

#Simulate r trajectories

delta_t2 <- TIME2/byn2
time2 <- seq(0,TIME2,by=delta_t2)

LT15 <- matrix(rep(0,length(time2)*r),nrow=r)

for (i in 1:r){LT15[i,] = GBM(x=x2009,r=mu15,sigma = sig15,T=TIME2,N=byn2)}

plot(2009:2020,LT15[1,],t='l', col="grey", ylab="Death Rates",xlab="Year",ylim=c(0,0.0015),cex.lab=1.5,cex.axis=1.5)
for(i in 2:r){lines(2009:2020,LT15[i,], t='l',col="grey")}
lines(2009:2020,t(Frates[70:81,116]),t="l",col="black")
legend("topright",c("Observed death rates"),pch="_",col="black",cex=1.5)

Means15 <- rep(NA,12)
sd15 <- rep(NA,12)
for(i in 1:12){
Means15[i] <- mean(LT15[,i])
sd15[i] <- sd(LT15[,i])
}

CIsupsimulation95 <- function(M,S){M + qnorm(0.975)*S}
CIinfsimulation95 <- function(M,S){M - qnorm(0.975)*S}

ValuesCIsimulationsup <- CIsupsimulation95(Means15[1:12],sd15[1:12])
ValuesCIsimulationinf <- CIinfsimulation95(Means15[1:12],sd15[1:12])

age15MLT <- cbind(Frates[70:81,116],Rate.fore.LT.GBM[1:12,116],ValuesCIsimulationinf,ValuesCIsimulationsup)

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

matplot(1940:2031,rates15,type="l",pch="_",col=c("black","black"),ylim = c(0,0.003),xlab = "Year",ylab = "Death Rates",cex.lab=1.5,cex.axis=1.5)
legend("topright",col=c("black","black"),c("Observed","Adjustment/Forecast"),lty=c(1,2),cex=1.2)

matplot(2009:2020,age15MSS,type="l",pch="_",col=c("black","black"),ylim = c(0,0.0005),xlab = "Year",ylab="Death Rates",cex.lab=1.5,cex.axis=1.5)
legend("topleft",col=c("black","black"),c("Observed","SS Forecasts"),lty=c(1,2),cex=1.2)

matplot(2009:2020,age15MLT,type="l",col=c("black","black","black"),lty=c(1,2,3,3),ylim = c(0,0.0005),xlab = "Year",ylab="Death Rates",cex.lab=1.5,cex.axis=1.5)
legend("topleft",col=c("black","black","black","black"),lty=c(1,2,3,3),c("Observed","LT forecasts","CI 95% Simulation"),cex=1.2) 


##PAGES 20 AND 21##

layout(matrix(c(1,1,1,1,1,1,2,3,4),3,3, byrow = TRUE),widths = c(c(3,6),c(3,6)))

MSEF <- MSE.adjust.GBM[1:100]
MSEM <- MSE.adjust.GBM[101:200]
mseadj <- cbind(MSEF,MSEM)

plot(0:99,mseadj[,1],type="l",col="grey",ylim = c(0,0.03),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:99,mseadj[,2],col="black")
legend("topleft",pch="_",col=c("grey","black"),c("Female","Male"),cex = 3)

plot(0:4,mseadj[1:5,1],type="l",col="grey",ylim = c(0,0.0006),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:4,mseadj[1:5,2],col="black")

plot(5:68,mseadj[6:69,1],type="l",col="grey",ylim = c(0,0.000031),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(5:68,mseadj[6:69,2],col="black")

plot(69:99,mseadj[70:100,1],type="l",col="grey",ylim = c(0,0.034),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(69:99,mseadj[70:100,2],col="black")

MSEFLT <- MSE.fore.LT.GBM[1:100]
MSEMLT <- MSE.fore.LT.GBM[101:200]
mseLT <- cbind(MSEFLT,MSEMLT)

plot(0:99,mseLT[,1],type="l",pch="_",col="grey",ylim = c(0,0.004),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:99,mseLT[,2],col="black")
legend("topleft",pch="_",col=c("grey","black"),c("Female","Male"),cex = 3)

plot(0:4,mseLT[1:5,1],type="l",pch="_",col="grey",ylim = c(0,0.0000004),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:4,mseLT[1:5,2],col="black")

plot(5:68,mseLT[6:69,1],type="l",pch="_",col="grey",ylim = c(0,0.000004),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(5:68,mseLT[6:69,2],col="black")

plot(69:99,mseLT[70:100,1],type="l",pch="_",col="grey",ylim = c(0,0.004),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(69:99,mseLT[70:100,2],col="black")


MSEFSS <- MSE.fore.1S.GBM[1:100]
MSEMSS <- MSE.fore.1S.GBM[101:200]
mseSS <- cbind(MSEFSS,MSEMSS)

plot(0:99,mseSS[,1],type="l",pch="_",col="grey",ylim = c(0,0.004),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:99,mseSS[,2],col="black")
legend("topleft",pch="_",col=c("grey","black"),c("Female","Male"),cex = 3)

plot(0:4,mseSS[1:5,1],type="l",pch="_",col="grey",ylim = c(0,0.0000004),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:4,mseSS[1:5,2],col="black")

plot(5:68,mseSS[6:69,1],type="l",pch="_",col="grey",ylim = c(0,0.000002),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(5:68,mseSS[6:69,2],col="black")

plot(69:99,mseSS[70:100,1],type="l",pch="_",ylim = c(0,0.004),col="grey",xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(69:99,mseSS[70:100,2],col="black")

par(mfrow = c(2,2))

DRM49GBMLT <- cbind(Frates[1:81,150],c(t(Rate.adjust.GBM[1:70,150]),t(Rate.fore.LT.GBM[2:12,150])))

matplot(1940:2020,DRM49GBMLT,type="l",col=1:2,ylim = c(0,0.015),xlab = "Year",ylab="Death Rates",cex.lab=1.5,cex.axis=1.5)
legend("topright",c("Observed","Adjustment/Forecast LT"),lty=1:2,col=1:2,cex = 1.5)

DRM99GBMLT <- cbind(Frates[1:81,200],c(t(Rate.adjust.GBM[1:70,200]),t(Rate.fore.LT.GBM[2:12,200])))

matplot(1940:2020,DRM99GBMLT,type="l",col=1:2,ylim = c(0,0.95),xlab = "Year",ylab="Death Rates",cex.lab=1.5,cex.axis=1.5)
legend("topright",c("Observed","Adjustment/Forecasts LT"),lty=1:2,col=1:2,cex = 1.5)

DRM49GBMSS <- cbind(Frates[1:81,150],c(t(Rate.adjust.GBM[1:70,150]),t(Rate.fore.1S.GBM[2:12,150])))

matplot(1940:2020,DRM49GBMSS,type="l",col=1:2,ylim = c(0,0.015),xlab = "Year",ylab="Death Rates",cex.lab=1.5,cex.axis=1.5)
legend("topright",c("Observed","Adjustment/Forecast SS"),lty=1:2,col=1:2,cex = 1.5)

DRM99GBMSS <- cbind(Frates[1:81,200],c(t(Rate.adjust.GBM[1:70,200]),t(Rate.fore.1S.GBM[2:12,200])))

matplot(1940:2020,DRM99GBMSS,type="l",col=1:2,ylim = c(0,0.95),xlab = "Year",ylab="Death Rates",cex.lab=1.5,cex.axis=1.5)
legend("topright",c("Observed","Adjustment/Forecasts SS"),lty=1:2,col=1:2,cex = 1.5)

#PAGE 28##

par(mfrow=c(3,2))

aplot <- cbind(a[1:100],a[101:200])

plot(0:99,aplot[,1],type="l",col="grey",ylim = c(0,0.5),xlab = "Age",ylab="a",cex.lab=1.6,cex.axis=1.6)
lines(0:99,aplot[,2],col="black")
legend("topleft",c("Female","Male"),pch="_",col=c("grey","black"),cex=2)

plot(0:89,aplot[1:90,1],type="l",col="grey",ylim = c(0,0.3),xlab = "Age",ylab="a",cex.lab=1.6,cex.axis=1.6)
lines(0:89,aplot[1:90,2],col="black")
legend("topleft",c("Female","Male"),pch="_",col=c("grey","black"),cex=2)

bplot <- cbind(b[1:100],b[101:200])

plot(0:99,bplot[,1],type="l",col="grey",ylim = c(0,10),xlab = "Age",ylab="b",cex.lab=1.6,cex.axis=1.6)
lines(0:99,bplot[,2],col="black")
legend("topleft",c("Female","Male"),pch="_",col=c("grey","black"),cex=2)

plot(0:89,bplot[1:90,1],type="l",col="grey",ylim = c(0,0.3),xlab = "Age",ylab="b",cex.lab=1.6,cex.axis=1.6)
lines(0:89,bplot[1:90,2],col="black")
legend("topleft",c("Female","Male"),pch="_",col=c("grey","black"),cex=2)

sigmaplot <- cbind(sigma[1:100],sigma[101:200])

plot(0:99,sigmaplot[,1],type="l",col="grey",ylim = c(0,0.8),xlab = "Age",ylab=expression(sigma),cex.lab=1.6,cex.axis=1.6)
lines(0:99,sigmaplot[,2],col="black")
legend("topleft",c("Female","Male"),pch="_",col=c("grey","black"),cex=2)

plot(0:89,sigmaplot[1:90,1],type="l",col="grey",ylim = c(0,0.4),xlab = "Age",ylab=expression(sigma),cex.lab=1.6,cex.axis=1.6)
lines(0:89,sigmaplot[1:90,2],col="black")
legend("topright",c("Female","Male"),pch="_",col=c("grey","black"),cex=2)

#PAGE 29##

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

Age39FSGM <- c(t(Frates[,30]),rep(NA,11))

Age39FFORESGM <- c(t(Rate.adjust.SGM[,30]),t(Rate.fore.LT21.SGM[,30]))

rates39F <- cbind(Age39FSGM,Age39FFORESGM)

matplot(1940:2031,rates39F,type="l",col=1:2,xlab = "Year",ylab="Death Rates",ylim = c(0,0.005),cex.lab=1.5,cex.axis=1.5)
legend("topright",c("Observed","Adjustament/Forecast LT"),lty=1:2,col=1:2,cex=1.5)

rates39Fss <- cbind(Frates[70:81,30],Rate.fore.1S.SGM[1:12,30])

matplot(2009:2020,rates39Fss,type="l",pch="_",col=1:2,xlab = "Year",ylab="Death Rates",ylim = c(0,0.001),cex.lab=1.5,cex.axis=1.5)
legend("topleft",c("Observed","SS forecasts"),lty=1:2,col=1:2,cex=1.5)

#Simulation

set.seed(12345)

y0 <- as.numeric(log(Frates[70,30]))
A29 <- as.numeric(p.SGM[30,1])
b29 <- as.numeric(p.SGM[30,2])
sig29 <- as.numeric(p.SGM[30,3])

simxtsgm <- matrix(rep(NA,24000),nrow=2000,ncol=12,byrow=TRUE)
t <- rep(NA,11)

simxtsgm[,1] <- rep(as.numeric(Frates[70,30]),2000)

for(t in 1:11){
simxtsgm[,(t+1)] <- exp(rnorm(2000,mean = (A29 + (y0 - A29)*exp(-b29*t)),sd = (sqrt(sig29^2*(1-exp(-2*b29*t))/(2*b29)))))
}

M29 <- rep(NA,12)
S29 <- rep(NA,12)

for(i in 1:12){
  M29[i] <- mean(simxtsgm[,i])
  S29[i] <- sd(simxtsgm[,i])
}

CIsgm29sup <- function(M,S){M + qnorm(0.975)*S}
CIsgm29inf <- function(M,S){M - qnorm(0.975)*S}

ValuesCIsgm29sup <- CIsgm29sup(M29[1:12],S29[1:12])
ValuesCIsgm29inf <- CIsgm29inf(M29[1:12],S29[1:12])

rates39FLT <- cbind(Frates[70:81,30],Rate.fore.LT.SGM[1:12,30],ValuesCIsgm29inf,ValuesCIsgm29sup)

matplot(2009:2020,rates39FLT,type="l",col=c(1:2,"black","black"),lty=c(1,2,3,3),xlab = "Year",ylab="Death Rates",ylim = c(0,0.001),cex.lab=1.5,cex.axis=1.5)
legend("topleft",c("Observed","LT forecasts","CI 95% simulation"),lty=c(1:3),col=c(1:2,"black"),cex=1.5) 

#PAGE 30#

layout(matrix(c(1,1,1,1,1,1,2,3,4),3,3, byrow = TRUE),widths = c(c(3,6),c(3,6)))

MSEFSGM <- MSE.adjust.SGM[1:100]

MSEMSGM <- MSE.adjust.SGM[101:200]

mseadjsgm <- cbind(MSEFSGM,MSEMSGM)

plot(0:99,mseadjsgm[,1],type="l",col="grey",ylim = c(0,0.03),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:99,mseadjsgm[,2],col="black")
legend("topleft",pch="_",col=c("grey","black"),c("Female","Male"),cex=3)

plot(0:4,mseadjsgm[1:5,1],type="l",col="grey",ylim = c(0,0.0006),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:4,mseadjsgm[1:5,2],col="black")

plot(5:68,mseadjsgm[6:69,1],type="l",col="grey",ylim = c(0,0.00004),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(5:68,mseadjsgm[6:69,2],col="black")

plot(69:99,mseadjsgm[70:100,1],type="l",col="grey",ylim = c(0,0.03),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(69:99,mseadjsgm[70:100,2],col="black")

MSEFLTSGM <- MSE.fore.LT.SGM[1:100]

MSEMLTSGM <- MSE.fore.LT.SGM[101:200]

mseadjLTSGM <- cbind(MSEFLTSGM,MSEMLTSGM)

plot(0:99,mseadjLTSGM[,1],type="l",col="grey",ylim = c(0,0.009),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:99,mseadjLTSGM[,2],col="black")
legend("topleft",pch="_",col=c("grey","black"),c("Female","Male"),cex=3)

plot(0:4,mseadjLTSGM[1:5,1],type="l",col="grey",ylim = c(0,0.0000004),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:4,mseadjLTSGM[1:5,2],col="black")

plot(5:68,mseadjLTSGM[6:69,1],type="l",col="grey",ylim = c(0,0.000015),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(5:68,mseadjLTSGM[6:69,2],col="black")

plot(69:99,mseadjLTSGM[70:100,1],type="l",col="grey",ylim = c(0,0.009),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(69:99,mseadjLTSGM[70:100,2],col="black")

MSEFSSSGM <- MSE.fore.1S.SGM[1:100]

MSEMSSSGM <- MSE.fore.1S.SGM[101:200]

mseadjSSSGM <- cbind(MSEFSSSGM,MSEMSSSGM)

plot(0:99,mseadjSSSGM[,1],type="l",col="grey",ylim = c(0,0.004),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:99,mseadjSSSGM[,2],col="black")
legend("topleft",pch="_",col=c("grey","black"),c("Female","Male"),cex=3)

plot(0:4,mseadjSSSGM[1:5,1],type="l",col="grey",ylim = c(0,0.0000004),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(0:4,mseadjSSSGM[1:5,2],col="black")

plot(5:68,mseadjSSSGM[6:69,1],type="l",col="grey",ylim = c(0,0.000002),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(5:68,mseadjSSSGM[6:69,2],col="black")

plot(69:99,mseadjSSSGM[70:100,1],type="l",col="grey",ylim = c(0,0.004),xlab = "Age",ylab="MSE",cex.lab=1.5,cex.axis=1.5)
lines(69:99,mseadjSSSGM[70:100,2],col="black")

#PAGE 31#

par(mfrow = c(1,2))

ObservedF23 <- c(t(Frates[,24]),rep(NA,11))

F23GBM <- c(t(Rate.adjust.GBM[1:70,24]),t(Rate.fore.LT21.GBM[1:22,24]))

F23SGM <- c(t(Rate.adjust.SGM[1:70,24]),t(Rate.fore.LT21.SGM[1:22,24]))

F23 <- cbind(ObservedF23,F23GBM,F23SGM)

matplot(1940:2031,F23,type="l",col=c("black","black","black"),ylim = c(0,0.006),xlab = "Year",ylab="Death Rates",cex.lab=1.5,cex.axis=1.5)
legend("topright",col=c("black","black","black"),lty=1:3,c("Observed","GBM","SGM"),cex=2)

ObservedM23 <- c(t(Frates[,124]),rep(NA,11))

M23GBM <- c(t(Rate.adjust.GBM[1:70,124]),t(Rate.fore.LT21.GBM[1:22,124]))

M23SGM <- c(t(Rate.adjust.SGM[1:70,124]),t(Rate.fore.LT21.SGM[1:22,124]))

M23 <- cbind(ObservedM23,M23GBM,M23SGM)

matplot(1940:2031,M23,type="l",col=c("black","black","black"),ylim = c(0,0.006),xlab = "Year",ylab="Death Rates",cex.lab=1.5,cex.axis=1.5)
legend("topright",col=c("black","black","black"),lty=1:3,c("Observed","GBM","SGM"),cex=2)

#PAGES 32 and 33#

par(mfrow=c(1,4))

DIFMSEFemale <- (MSE.adjust.GBM[1:100] - MSE.adjust.SGM[1:100])*10000

DIFMSEMale <- (MSE.adjust.GBM[101:200] - MSE.adjust.SGM[101:200])*10000

barplot(DIFMSEFemale[1:25],ylim = c(-0.25,0.02),xlab = "Age 0 to 24",ylab = "(MSE(GBM)-MSE(SGM))*10000",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEFemale[26:50],ylim = c(-0.006,0.006),xlab = "Ages 25 to 49",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEFemale[51:75],ylim = c(-0.3,0),xlab = "Ages 50 to 74",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEFemale[76:100],ylim = c(-5,100),xlab = "Ages 75 to 99",cex.lab=1.5,cex.axis = 1.5)

barplot(DIFMSEMale[1:25],ylim = c(-0.3,0.02),xlab = "Ages 0 to 24",ylab = "(MSE(GBM)-MSE(SGM))*10000",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEMale[26:50],ylim = c(-0.005,0.02),xlab = "Ages 25 to 49",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEMale[51:75],ylim = c(-0.4,0),xlab = "Ages 50 to 74",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEMale[76:100],ylim = c(-20,100),xlab = "Ages 75 to 99",cex.lab=1.5,cex.axis = 1.5)

DIFMSEFemaleSS <- (MSE.fore.1S.GBM[1:100] - MSE.fore.1S.SGM[1:100])*10000

DIFMSEMaleSS <- (MSE.fore.1S.GBM[101:200] - MSE.fore.1S.SGM[101:200])*10000

barplot(DIFMSEFemaleSS[1:25],ylim = c(-0.00015,0),xlab = "Ages 0 to 24",ylab = "(MSE(GBM)-MSE(SGM))*10000",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEFemaleSS[26:50],ylim = c(-0.0001,0.00001),xlab = "Ages 25 to 49",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEFemaleSS[51:75],ylim = c(-0.0004,0.0004),xlab = "Ages 50 to 74",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEFemaleSS[76:100],ylim = c(-20,0),xlab = "Ages 75 to 99",cex.lab=1.5,cex.axis = 1.5)

barplot(DIFMSEMaleSS[1:25],ylim = c(-0.0002,0),xlab = "Age 0 to 24",ylab = "(MSE(GBM)-MSE(SGM))*10000",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEMaleSS[26:50],ylim = c(-0.0004,0),xlab = "Ages 25 to 49",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEMaleSS[51:75],ylim = c(-0.001,0.001),xlab = "Ages 50 to 74",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEMaleSS[76:100],ylim = c(-40,20),xlab = "Ages 75 to 99",cex.lab=1.5,cex.axis = 1.5)

DIFMSEFemaleLT <- (MSE.fore.LT.GBM[1:100] - MSE.fore.LT.SGM[1:100])*10000

DIFMSEMaleLT <- (MSE.fore.LT.GBM[101:200] - MSE.fore.LT.SGM[101:200])*10000

barplot(DIFMSEFemaleLT[1:25],ylim = c(-0.0004,0.00002),xlab = "Ages 0 to 24",ylab = "(MSE(GBM)-MSE(SGM))*10000",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEFemaleLT[26:50],ylim = c(-0.004,0),xlab = "Ages 25 to 49",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEFemaleLT[51:75],ylim = c(-0.1,0),xlab = "Ages 50 to 74",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEFemaleLT[76:100],ylim = c(-70,5),xlab = "Ages 75 to 99",cex.lab=1.5,cex.axis = 1.5)

barplot(DIFMSEMaleLT[1:25],ylim = c(-0.003,0),xlab = "Ages 0 to 24",ylab = "(MSE(GBM)-MSE(SGM))*10000",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEMaleLT[26:50],ylim = c(-0.02,0),xlab = "Ages 25 to 49",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEMaleLT[51:75],ylim = c(-0.2,0),xlab = "Ages 50 to 74",cex.lab=1.5,cex.axis = 1.5)
barplot(DIFMSEMaleLT[76:100],ylim = c(-60,5),xlab = "Ages 75 to 99",cex.lab=1.5,cex.axis = 1.5)


