"clvarselsampgr" <-
function(X,G,emModels1=c("E","V"),emModels2=c("EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"),sampsize=2000,allow.EEE=TRUE,forcetwo=TRUE,itermax=100)
{
#number of rows=number of observations
n<-nrow(X)
#number of columns=number of variables
d<-ncol(X)
#Sample the subset of observations for hierarchical clustering 
sub<-sample(1:n,sampsize,replace=FALSE)
#sub<-sample(c(FALSE,TRUE),n,TRUE,c(1-sampsize/n,sampsize/n))
#First Step - selecting single variable
maxBIC<-rep(0,d)
maxdiff<-rep(0,d)
oneBIC<-rep(NA,d)
for(i in 1:d)
{
xBIC<-NULL
#Fit the cluster models from 2 to G groups
try(xBIC<-Mclust(X[,i],2:G,emModels1,initialization=list(subset=sub)),TRUE)
if(is.null(xBIC)) try(xBIC<-Mclust(X[,i],2:G,emModels1),TRUE)
#If we get all NA's from "V" starting hierarchical values use "E"
if((allow.EEE)&sum(is.finite(xBIC$BIC))==0) try(xBIC<-Mclust(X[,i],2:G,emModels1,initialization=list(hcPairs = hcE(X[,i],subset=sub))),TRUE)
if((allow.EEE)&sum(is.finite(xBIC$BIC))==0) try(xBIC<-Mclust(X[,i],2:G,emModels1,initialization=list(hcPairs = hcE(X[,i]))),TRUE)
#maxBIC is the maximum BIC over all clustering models (2 to G groups) fit
if(sum(is.finite(xBIC$BIC))==0) maxBIC[i]<-NA else maxBIC[i]<-max(xBIC$BIC[is.finite(xBIC$BIC)])
#Fit and get BIC for a single component no-cluster normal model
try(oneBIC[i]<-Mclust(X[,i],1,"V",initialization=list(subset=sub))$BIC[1],TRUE)
if(is.na(oneBIC[i])) try(oneBIC[i]<-Mclust(X[,i],1,"V")$BIC[1],TRUE)
#Difference between maximum BIC for clustering and BIC for no clustering
maxdiff[i]<-c(maxBIC[i]-oneBIC[i])
}

#Find the variable with the biggest difference between clustering and no clustering
m<-max(maxdiff[is.finite(maxdiff)])
arg<-which(maxdiff==m,arr.ind=TRUE)[1]
#This is our first selected variable/S is the matrix of currently selected clustering variables
S<-matrix(c(X[,arg]),n,1)
#BICS is the BIC value for the clustering model with the variable(s) in S
BICS<-maxBIC[arg]
colnames(S)<-colnames(X)[arg]
#NS is the matrix of currently not selected variables
NS<-as.matrix(X[,-arg])
colnames(NS)<-colnames(X)[-arg]
#mat records the proposed variable, BIC for the S matrix and difference in BIC for clustering versus no clustering on S, whether it was an addition step and if it was accepted
mat<-matrix(c(colnames(S),BICS,maxdiff[arg],"Add","Accepted"),1,5)

#Second Step - selecting second variable
regBIC<-rep(0,ncol(NS))
depBIC<-rep(0,ncol(NS))
cindepBIC<-rep(0,ncol(NS))
for(i in 1:ncol(NS))
{
 sBIC<-NULL
#Fit the regression of the proposed variable on the variable in S
 fm<-lm(NS[,i]~S)
 sigma<-(sum((summary(fm)$resid)^2)/n)^0.5
#Calculate the BIC for the regression
 regBIC[i]<--n*log(2*pi)-2*n*log(sigma)-n-log(n)*3
#Fit the cluster model on the two variables for 2 to G groups 
 try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(subset=sub)),TRUE)
 if(is.null(sBIC)) sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2)
#If we get all NA's from "VVV" starting hierarchical values use "EEE"
 if((allow.EEE)&sum(is.finite(sBIC$BIC))==0) try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(hcPairs = hcEEE(cbind(S,NS[,i])[sub,]),subset=sub)),TRUE)
  if((allow.EEE)&sum(is.finite(sBIC$BIC))==0) try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(hcPairs = hcEEE(cbind(S,NS[,i])))),TRUE)
#depBIC is the BIC for the clustering model with both variables
 if(sum(is.finite(sBIC$BIC))==0) depBIC[i]<-NA else depBIC[i]<-max(sBIC$BIC[is.finite(sBIC$BIC)])
#cindepBIC is the BIC for the clustering model on S and the regression model of the new variable on S
 cindepBIC[i]<-regBIC[i]+BICS
}
#cdiff is the difference between BIC for the models with variables' being clustering variables versus them being conditionally independent of the clustering
cdiff<-depBIC-cindepBIC
#Choose the variable with the largest difference
m<-max(cdiff[is.finite(cdiff)])
arg<-which(cdiff==m,arr.ind=TRUE)[1]

#if forcetwo is true automatically add the best second variable, otherwise only add it if its difference is positive
if(forcetwo||cdiff[arg]>0){
k<-c(colnames(S),colnames(NS)[arg])
nks<-c(colnames(NS)[-arg])
BICS<-depBIC[arg]
mat<-rbind(mat,c(colnames(NS)[arg],BICS,cdiff[arg],"Add","Accepted"))
S<-cbind(S,NS[,arg])
NS<-as.matrix(NS[,-arg])
colnames(S)<-k
colnames(NS)<-nks
} else{
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not

  mat<-rbind(mat,c(colnames(NS)[arg],BICS,cdiff[arg],"Add","Rejected"))
}

criterion<-1
iter<-0
while((criterion==1)&iter<itermax)
{
iter<-iter+1
check1<-colnames(S)

#Addition step
#For the special case where we have removed all the clustering variables/S is empty
if((ncol(S)==0||is.null(ncol(S)))){
#We simply choose the same variable as in the first step and check whether the difference between the BIC for clustering versus not clustering is positive or not
m<-max(maxdiff[is.finite(maxdiff)])
arg<-which(maxdiff==m,arr.ind=TRUE)[1]
if(maxdiff[arg]>0) 
{
#if the difference is positive this variable is selected as a clustering variable
S<-matrix(c(X[,arg]),n,1)
BICS<-maxBIC[arg]
colnames(S)<-colnames(X)[arg]
NS<-as.matrix(X[,-arg])
colnames(NS)<-colnames(X)[-arg]
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not

mat<-rbind(mat,c(colnames(S),BICS,maxdiff[arg],"Add","Accepted"))
} else{
#if the difference is not positive no clustering variables exist
BICS<-NA
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not

mat<-rbind(mat,c(colnames(X)[arg],BICS,maxdiff[arg],"Add","Rejected"))
}
} else{

#Addition Step in general (for all cases except when S is empty)
if((ncol(NS)!=0&!is.null(ncol(NS)))){
regBIC<-rep(0,ncol(NS))
depBIC<-rep(0,ncol(NS))
cindepBIC<-rep(0,ncol(NS))
#p=no of regression parameters
p<-ncol(S)+2
for(i in 1:ncol(NS))
{
 sBIC<-NULL
#Fit the regression of the proposed variable on the variable(s) in S
 fm<-lm(NS[,i]~S)
 sigma<-(sum((summary(fm)$resid)^2)/n)^0.5
#Calculate the BIC for the regression
 regBIC[i]<--n*log(2*pi)-2*n*log(sigma)-n-log(n)*p
#Fit the cluster model on the S variables with the proposed variable for 2 to G groups 
 try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(subset=sub)),TRUE)
 if(is.null(sBIC)) try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2),TRUE)
#If we get all NA's from "VVV" starting hierarchical values use "EEE"
 if((allow.EEE)&(sum(is.finite(sBIC$BIC))==0)) try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(hcPairs = hcEEE(cbind(S,NS[,i])[sub,]),subset=sub)),TRUE)
 if((allow.EEE)&(sum(is.finite(sBIC$BIC))==0)) try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(hcPairs = hcEEE(cbind(S,NS[,i])))),TRUE)
#depBIC is the BIC for the clustering model with both S and proposed variable
 if(sum(is.finite(sBIC$BIC))==0) depBIC[i]<-NA else depBIC[i]<-max(sBIC$BIC[is.finite(sBIC$BIC)])
#cindepBIC is the BIC for the clustering model on S and the regression model of the new variable on S
 cindepBIC[i]<-regBIC[i]+BICS
}
#cdiff is the difference between BIC for the models with variables' being clustering variables versus them being conditionally independent of the clustering
cdiff<-depBIC-cindepBIC
#Choose the variable with the largest difference
m<-max(cdiff[is.finite(cdiff)])
arg<-which(cdiff==m,arr.ind=TRUE)[1]

if(cdiff[arg]>0){ 
#if this difference is positive add this variable to S and update the clustering model's BICS
BICS<-depBIC[arg] 
k<-c(colnames(S),colnames(NS)[arg])
nks<-c(colnames(NS)[-arg])
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not

mat<-rbind(mat,c(colnames(NS)[arg],BICS,cdiff[arg],"Add","Accepted")) 
S<-cbind(S,NS[,arg])
NS<-as.matrix(NS[,-arg])
colnames(S)<-k
colnames(NS)<-nks
} else{ 
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not

mat<-rbind(mat,c(colnames(NS)[arg],BICS,cdiff[arg],"Add","Rejected"))
}
}
}

#Removal Step for the special case where S contains only a single variable
if(ncol(S)==1){
cdiff<-0
oneBIC<-NA
try(oneBIC<-Mclust(S,1,"V")$BIC[1],TRUE)
#Difference between maximum BIC for clustering and BIC for no clustering
cdiff<-c(BICS-oneBIC)
if(is.na(cdiff)) cdiff<-0 else cdiff<-cdiff
#Check if difference is negative
if(cdiff<=0)
{
#if negative remove the variable from S and set the BIC for the model to NA
BICS<-NA
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not

mat<-rbind(mat,c(colnames(S),BICS,cdiff,"Remove","Accepted"))
k<-c(colnames(NS),colnames(S))
NS<-cbind(NS,S)
S<-NULL
colnames(NS)<-k
} else{
#Otherwise leave S and BICS alone
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not

mat<-rbind(mat,c(colnames(S),BICS,cdiff,"Remove","Rejected"))
}
} else{
#Removal step in general (for all cases except when S is a single variable or empty)
if(ncol(S)>=2){
rdep<-rep(0,ncol(S))
regBIC<-rep(0,ncol(S))
cindepBIC<-rep(0,ncol(S))
#Check if the data is at least 3 dimensional
mult<-ncol(S)>2
ifelse(mult,name<-emModels2,name<-emModels1)
p<-ncol(S)+1
for(i in 1:ncol(S))
{
 sBIC<-NULL
#Fit the regression of the proposed variable from S on the other variable(s) in S
 fm<-lm(S[,i]~S[,-i])
 sigma<-(sum((summary(fm)$resid)^2)/n)^0.5
#Calculate the BIC for the regression
 regBIC[i]<--n*log(2*pi)-2*n*log(sigma)-n-log(n)*p
#Fit the cluster model on the S variables without the proposed variable for 2 to G groups 
 try(sBIC<-Mclust(S[,-i],2:G,name,initialization=list(subset=sub)),TRUE) 
 if(is.null(sBIC)) try(sBIC<-Mclust(S[,-i],2:G,name),TRUE) 
#If we get all NA's from "VVV" starting hierarchical values use "EEE"
 if((allow.EEE)&ncol(S)>=3&sum(is.finite(sBIC$BIC))==0){try(sBIC<-Mclust(S[,-i],2:G,name,initialization=list(hcPairs = hcEEE(S[sub,-i]),subset=sub)),TRUE)} else{if((allow.EEE)&ncol(S)==2&sum(is.finite(sBIC$BIC))==0){try(sBIC<-Mclust(S[,-i],2:G,name,initialization=list(hcPairs = hcE(S[sub,-i]),subset=sub)),TRUE)}}
 if((allow.EEE)&ncol(S)>=3&sum(is.finite(sBIC$BIC))==0){try(sBIC<-Mclust(S[,-i],2:G,name,initialization=list(hcPairs = hcEEE(S[,-i]))),TRUE)} else{if((allow.EEE)&ncol(S)==2&sum(is.finite(sBIC$BIC))==0){try(sBIC<-Mclust(S[,-i],2:G,name,initialization=list(hcPairs = hcE(S[,-i]))),TRUE)}}
 if(sum(is.finite(sBIC$BIC))==0) rdep[i]<-NA else rdep[i]<-max(sBIC$BIC[is.finite(sBIC$BIC)])
#cindepBIC is the BIC for the clustering model on the other variables in S and the regression model of the proposed variable on the other variables in S
 cindepBIC[i]<-regBIC[i]+rdep[i]
}
#depBIC is the BIC for the clustering model with all variables in S
depBIC<-BICS
#cdiff is the difference between BIC for the models with variables' being clustering variables versus them being conditionally independent of the clustering
cdiff<-depBIC-cindepBIC
#Choose the variable with the smallest difference
m<-min(cdiff[is.finite(cdiff)])
arg<-which(cdiff==m,arr.ind=TRUE)[1]

if(cdiff[arg]<=0){
#if this difference is negative remove this variable from S and update the clustering model's BICS
BICS<-rdep[arg] 
k<-c(colnames(NS),colnames(S)[arg])
nks<-c(colnames(S)[-arg])
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not

mat<-rbind(mat,c(colnames(S)[arg],BICS,cdiff[arg],"Remove","Accepted")) 
NS<-cbind(NS,S[,arg])
S<-as.matrix(S[,-arg])
colnames(S)<-nks
colnames(NS)<-k
} else{ 
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not

mat<-rbind(mat,c(colnames(S)[arg],BICS,cdiff[arg],"Remove","Rejected"))
}
}
}

#Check if the variables in S have changed or not
check2<-colnames(S)
#if they have changed (either added one or removed one or changed one) then continue the algorithm (criterion is 1) otherwise stop (criterion is 0)
if(length(check2)!=length(check1)) criterion<-1 else{ if(sum(check1==check2)!=length(check1)) criterion<-1 else {criterion<-0}}
}
#List the selected variables and the matrix of steps' information
colnames(mat)<-c("Variable proposed","BIC of new clustering variables set","BIC difference","Type of step","Decision")
if(iter==itermax) print("Warning: Algorithm stopped because maximum number of iterations was reached")
return(list(sel.var=S,steps.info=mat))
}

