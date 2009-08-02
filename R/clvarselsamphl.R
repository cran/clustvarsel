"clvarselsamphl" <-
function(X,G,emModels1=c("E","V"),emModels2=c("EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"),sampsize=2000,allow.EEE=TRUE,forcetwo=TRUE,upper=0,lower=-10,itermax=100)
{
#number of rows=number of observations
n<-nrow(X)
#number of columns=number of variables
d<-ncol(X)
#Sample the subset of observations for hierarchical clustering 
sub<-sample(1:n,sampsize,replace=FALSE)
#First Step - selecting single variable and ordering
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
if((allow.EEE)&sum(is.finite(xBIC$BIC))==0) try(xBIC<-Mclust(X[,i],2:G,emModels1,initialization=list(hcPairs = hcE(X[sub,i]),subset=sub)),TRUE) 
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
temp<-maxdiff[-arg]
temp2<-sort(temp,decreasing=TRUE,index.return=TRUE)$ix
#NS is the matrix of currently not selected variables
NS<-as.matrix(X[,-arg])
#This orders NS in terms of strongest evidence of univariate clustering versus no clustering
colnames(NS)<-colnames(X)[-arg]
NS<-as.matrix(NS[,temp2])
#mat records the proposed variable, BIC for the S matrix and difference in BIC for clustering versus no clustering on S, whether it was an addition step and if it was accepted
mat<-matrix(c(colnames(S),BICS,maxdiff[arg],"Add","Accepted"),1,5)

#Second Step - selecting second variable
regBIC<-0
depBIC<-0
DepBIC<-NULL
cindepBIC<-0
crit<--10
cdiff<-0
Cdiff<-NULL
i<-0
#We only run until we find a variable whose difference in BIC between being included in the clustering variables versus conditionally independent of the clustering is greater than upper
while(crit<=upper&i<ncol(NS))
{
 i<-i+1
 sBIC<-NULL
#Fit the regression of the proposed variable on the variable in S
 fm<-lm(NS[,i]~S)
 sigma<-(sum((summary(fm)$resid)^2)/n)^0.5
#Calculate the BIC for the regression
 regBIC<--n*log(2*pi)-2*n*log(sigma)-n-log(n)*3
#Fit the cluster model on the two variables for 2 to G groups 
 try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(subset=sub)),TRUE)
 if(is.null(sBIC)) try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2),TRUE)
#If we get all NA's from "VVV" starting hierarchical values use "EEE"
 if((allow.EEE)&sum(is.finite(sBIC$BIC))==0) try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(hcPairs = hcEEE(cbind(S,NS[,i])[sub,]),subset=sub)),TRUE)
 if((allow.EEE)&sum(is.finite(sBIC$BIC))==0) try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(hcPairs = hcEEE(cbind(S,NS[,i])))),TRUE)
#depBIC is the BIC for the clustering model with both variables
 if(sum(is.finite(sBIC$BIC))==0) depBIC<-NA else depBIC<-max(sBIC$BIC[is.finite(sBIC$BIC)])
 DepBIC<-c(DepBIC,depBIC)
 cindepBIC<-regBIC+BICS
#cindepBIC is the BIC for the clustering model on S and the regression model of the new variable on S
 cdiff<-depBIC-cindepBIC
 if(!is.finite(cdiff)) cdiff<-upper
 Cdiff<-c(Cdiff,cdiff)
 crit<-cdiff 
}
if(cdiff>upper)
{
#i.e. evidence is stronger for including variable in S
 k<-c(colnames(S),colnames(NS)[i])
 S<-cbind(S,NS[,i])
 colnames(S)<-k
 BICS<-depBIC
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not
 mat<-rbind(mat,c(colnames(NS)[i],BICS,cdiff,"Add","Accepted"))
 ns<-NULL
 s<-NULL
# i is the index of those variables not selected but whose evidence of clustering BIC did not fall below "lower" or those not looked at yet
 if(i<ncol(NS)) ns<-(i+1):ncol(NS)
 if(i>1) s<-c(1:(i-1))[which(Cdiff[-i]>lower)]

 ind<-c(s,ns)
 if(!is.null(ind)){
 nks<-c(colnames(NS)[ind])
#NS is the not selected clustering variables whose recently calculated evidence of clustering BIC was higher than "lower" or variables not yet looked at
 NS<-as.matrix(NS[,ind])
 colnames(NS)<-nks} else{
 	NS<-NULL}
} else{
if(cdiff<upper&(forcetwo)){
#if the evidence is weaker but we're forcing choice of second variable
 m<-max(Cdiff[is.finite(Cdiff)])
 i<-which(Cdiff==m,arr.ind=TRUE)[1]
 k<-c(colnames(S),colnames(NS)[i])
 S<-cbind(S,NS[,i])
 colnames(S)<-k
 BICS<-DepBIC[i]
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not
 mat<-rbind(mat,c(colnames(NS)[i],BICS,Cdiff[i],"Add","Accepted"))
 nks<-c(colnames(NS)[-i])
 NS<-as.matrix(NS[,-i])
 temp<-Cdiff[-i]
#NS is the not selected clustering variables whose recently calculated evidence of clustering BIC was higher than "lower"
 if(sum(temp>lower)!=0){
 NS<-as.matrix(NS[,c(which(temp>lower))])
 colnames(NS)<-nks[c(which(temp>lower))]
 } else{
 	NS<-NULL
 	} 
} else{
 m<-max(Cdiff[is.finite(Cdiff)])
 i<-which(Cdiff==m,arr.ind=TRUE)[1]
mat<-rbind(mat,c(colnames(NS)[i],BICS,Cdiff[i],"Add","Rejected"))
}}

criterion<-1
iter<-0
while((criterion==1)&iter<itermax)
{
iter<-iter+1
check1<-colnames(S)

#Addition step
#For the special case where we have removed all the clustering variables/S is empty

if((ncol(NS)!=0&!is.null(ncol(NS)))&(ncol(S)==0)||(is.null(ncol(S)))){
depBIC<-0
DepBIC<-NULL
crit<--10
cdiff<-0
Cdiff<-NULL
#oneBIC<-rep(NA,d)
i<-0
crit<--10
while(crit<=upper&i<ncol(NS))
{
xBIC<-NULL
oneBIC<-NULL
i<-i+1
#Fit the cluster models from 2 to G groups
try(xBIC<-Mclust(X[,i],2:G,emModels1,initialization=list(subset=sub)),TRUE)
if(is.null(xBIC)) try(xBIC<-Mclust(X[,i],2:G,emModels1),TRUE)
#If we get all NA's from "V" starting hierarchical values use "E"
if((allow.EEE)&sum(is.finite(xBIC$BIC))==0) try(xBIC<-Mclust(X[,i],2:G,emModels1,initialization=list(hcPairs = hcE(X[sub,i]),subset=sub)),TRUE)
if((allow.EEE)&sum(is.finite(xBIC$BIC))==0) try(xBIC<-Mclust(X[,i],2:G,emModels1,initialization=list(hcPairs = hcE(X[,i]))),TRUE)
#depBIC is the maximum BIC over all clustering models (2 to G groups) fit
if(sum(is.finite(xBIC$BIC))==0) depBIC<-NA else depBIC<-max(xBIC$BIC[is.finite(xBIC$BIC)])
DepBIC<-c(DepBIC,depBIC)
#Fit and get BIC for a single component no-cluster normal model
try(oneBIC<-Mclust(X[,i],1,"V",initialization=list(subset=sub))$BIC[1],TRUE)
if(is.null(oneBIC)) try(oneBIC<-Mclust(X[,i],1,"V")$BIC[1],TRUE)
#Difference between maximum BIC for clustering and BIC for no clustering
cdiff<-c(depBIC-oneBIC)
if(!is.finite(cdiff)) cdiff<-upper
Cdiff<-c(Cdiff,cdiff)
crit<-cdiff
}
if(cdiff>upper)
{
#ie. evidence is stronger for including variable in S
 k<-c(colnames(NS)[i])
 S<-as.matrix(NS[,i])
 colnames(S)<-k
 BICS<-depBIC
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not
 mat<-rbind(mat,c(colnames(NS)[i],BICS,cdiff,"Add","Accepted"))
 ns<-NULL
 s<-NULL
# i is the index of those variables not selected but whose evidence of clustering BIC did not fall below "lower" or those not looked at yet
 if(i<ncol(NS)) ns<-(i+1):ncol(NS)
 if(i>1) s<-c(1:(i-1))[which(Cdiff[-i]>lower)]
 ind<-c(s,ns)
 if(!is.null(ind)){
 	nks<-c(colnames(NS)[ind])
#NS is the not selected clustering variables whose recently calculated evidence of clustering BIC was higher than "lower" or variables not yet looked at
 NS<-as.matrix(NS[,ind])
 colnames(NS)<-nks
 } else{
 	NS<-NULL
 	}
} else{
 m<-max(Cdiff[is.finite(Cdiff)])
 i<-which(Cdiff==m,arr.ind=TRUE)[1]
 mat<-rbind(mat,c(colnames(NS)[i],BICS,Cdiff[i],"Add","Rejected"))
 ind<-c(1:ncol(NS))[which(Cdiff>lower)]
 if(!is.null(ind)){
 k<-colnames(NS)[ind]
 #Exclude variables in NS whose evidence of clustering in this step was lower than "lower"
 NS<-as.matrix(NS[,ind])
 colnames(NS)<-k
} else{
NS<-NULL
}
}
} else{

#Addition Step in general (for all cases except when S is empty)
if(ncol(NS)!=0&!is.null(ncol(NS))){
regBIC<-0
depBIC<-0
DepBIC<-NULL
cindepBIC<-0
crit<--10
cdiff<-0
Cdiff<-NULL
i<-0
#p=no of regression parameters
p<-ncol(S)+2
#We only run until we find a variable whose difference in BIC between being included in the clustering variables versus conditionally independent of the clustering is greater than upper
while(crit<=upper&i<ncol(NS))
{
 sBIC<-NULL
 i<-i+1
#Fit the regression of the proposed variable on the variable(s) in S
 fm<-lm(NS[,i]~S)
 sigma<-(sum((summary(fm)$resid)^2)/n)^0.5
#Calculate the BIC for the regression
 regBIC<--n*log(2*pi)-2*n*log(sigma)-n-log(n)*p
#Fit the cluster model on the S variables with the proposed variable for 2 to G groups 
 try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(subset=sub)),TRUE)
 if(is.null(sBIC)) try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2),TRUE)
#If we get all NA's from "VVV" starting hierarchical values use "EEE"
 if((allow.EEE)&(sum(is.finite(sBIC$BIC))==0)) try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(hcPairs = hcEEE(cbind(S,NS[,i])[sub,]),subset=sub)),TRUE) 
 if((allow.EEE)&(sum(is.finite(sBIC$BIC))==0)) try(sBIC<-Mclust(cbind(S,NS[,i]),2:G,emModels2,initialization=list(hcPairs = hcEEE(cbind(S,NS[,i])))),TRUE) 
#depBIC is the BIC for the clustering model with both S and proposed variable
 if(sum(is.finite(sBIC$BIC))==0) depBIC<-NA else depBIC<-max(sBIC$BIC[is.finite(sBIC$BIC)])
 DepBIC<-c(DepBIC,depBIC)
#cindepBIC is the BIC for the clustering model on S and the regression model of the new variable on S
 cindepBIC<-regBIC+BICS
 cdiff<-depBIC-cindepBIC
 if(!is.finite(cdiff)) cdiff<-upper
 Cdiff<-c(Cdiff,cdiff)
 crit<-cdiff
}
if(cdiff>upper)
{
#ie. evidence is stronger for including variable in S
 k<-c(colnames(S),colnames(NS)[i])
 nks<-c(colnames(NS)[-i])
 S<-cbind(S,NS[,i])
 colnames(S)<-k
 BICS<-depBIC
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not
 mat<-rbind(mat,c(colnames(NS)[i],BICS,cdiff,"Add","Accepted"))
 ns<-NULL
 s<-NULL
 #Exclude variables in NS whose evidence of clustering in this step was lower than "lower"
 if(i<ncol(NS)) ns<-(i+1):ncol(NS)
 if(i>1) s<-c(1:(i-1))[which(Cdiff[-i]>lower)]
 ind<-c(s,ns)
 if(!is.null(ind)){
 nks<-colnames(NS)[ind]
 NS<-as.matrix(NS[,ind])
 colnames(NS)<-nks
} else{
NS<-NULL
}} else{
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not
m<-max(Cdiff[is.finite(Cdiff)])
i<-which(Cdiff==m,arr.ind=TRUE)[1]
mat<-rbind(mat,c(colnames(NS)[i],BICS,Cdiff[i],"Add","Rejected"))
ind<-c(1:ncol(NS))[which(Cdiff>lower)]
if(!is.null(ind)){
k<-colnames(NS)[ind]
#Exclude variables in NS whose evidence of clustering in this step was lower than "lower"
NS<-as.matrix(NS[,ind])
colnames(NS)<-k
} else{
NS<-NULL
}}

}}

#Removal Step for the special case where S contains only a single variable

if(ncol(S)==1){
cdiff<-0
oneBIC<-NA
try(oneBIC<-Mclust(S,1,"V",initialization=list(subset=sub))$BIC[1],TRUE)
if(is.na(oneBIC)) try(oneBIC<-Mclust(S,1,"V")$BIC[1],TRUE)
#Difference between maximum BIC for clustering and BIC for no clustering
cdiff<-c(BICS-oneBIC)
if(is.na(cdiff)) cdiff<-upper
#check if difference is negative
if(cdiff<=upper)
{
#if negative remove the variable from S and set the BIC for the model to NA
BICS<-NA
mat<-rbind(mat,c(colnames(S),BICS,cdiff,"Remove","Accepted"))
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not
#Only return variable to NS if difference is greater than "lower"
if(cdiff>lower){
k<-c(colnames(NS),colnames(S))
NS<-cbind(NS,S)
colnames(NS)<-k
S<-NULL
} else{
	S<-NULL
	}} else{
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not
mat<-rbind(mat,c(colnames(S),BICS,cdiff,"Remove","Rejected"))
}} else{

#Removal step in general (for all cases except when S is a single variable or empty)
if(ncol(S)>=2){
depBIC<-BICS
rdep<-0
regBIC<-0
cindepBIC<-0
cdiff<-0
Cdiff<-NULL
crit<-10
i<-0
#Check if the data is at least 3 dimensional
mult<-ncol(S)>2
ifelse(mult,name<-emModels2,name<-emModels1)
p<-ncol(S)+1
#We only run until we find a variable whose difference in BIC between being included in the clustering variables versus conditionally independent of the clustering is lower than "upper"
while(crit>upper&i<ncol(S))
{
 i<-i+1
 sBIC<-NULL
#Fit the regression of the proposed variable from S on the other variable(s) in S
 fm<-lm(S[,i]~S[,-i])
 sigma<-(sum((summary(fm)$resid)^2)/n)^0.5
#Calculate the BIC for the regression
 regBIC<--n*log(2*pi)-2*n*log(sigma)-n-log(n)*p
#Fit the cluster model on the S variables without the proposed variable for 2 to G groups 
 try(sBIC<-Mclust(as.matrix(S[,-i]),2:G,name,initialization=list(subset=sub)),TRUE)
 if(is.null(sBIC)) try(sBIC<-Mclust(as.matrix(S[,-i]),2:G,name),TRUE)
#If we get all NA's from "VVV" starting hierarchical values use "EEE"
 if((allow.EEE)&ncol(S)>=3&sum(is.finite(sBIC$BIC))==0){try(sBIC<-Mclust(S[,-i],2:G,name,initialization=list(hcPairs = hcEEE(S[sub,-i]),subset=sub)),TRUE)} else{if((allow.EEE)&ncol(S)==2&sum(is.finite(sBIC$BIC))==0){try(sBIC<-Mclust(as.matrix(S[,-i]),2:G,name,initialization=list(hcPairs = hcE(S[sub,-i]),subset=sub)),TRUE)}} 
 if((allow.EEE)&ncol(S)>=3&sum(is.finite(sBIC$BIC))==0){try(sBIC<-Mclust(S[,-i],2:G,name,initialization=list(hcPairs = hcEEE(S[,-i]))),TRUE)} else{if((allow.EEE)&ncol(S)==2&sum(is.finite(sBIC$BIC))==0){try(sBIC<-Mclust(as.matrix(S[,-i]),2:G,name,initialization=list(hcPairs = hcE(S[,-i]))),TRUE)}} 
 if(sum(is.finite(sBIC$BIC))==0) rdep<-NA else rdep<-max(sBIC$BIC[is.finite(sBIC$BIC)])
#cindepBIC is the BIC for the clustering model on the other variables in S and the regression model of the proposed variable on the other variables in S
cindepBIC<-regBIC+rdep
cdiff<-depBIC-cindepBIC
if(!is.finite(cdiff)) cdiff<-0
Cdiff<-c(Cdiff,cdiff)
crit<-cdiff
}
if(cdiff<upper&cdiff>lower)
{
#i.e. evidence is stronger for excluding variable from S but still including it in NS
BICS<-rdep
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not
mat<-rbind(mat,c(colnames(S)[i],BICS,cdiff,"Remove","Accepted"))
k<-c(colnames(NS),colnames(S)[i])
nk<-colnames(S)[-i]
NS<-cbind(NS,S[,i])
S<-as.matrix(S[,-i])
colnames(NS)<-k
colnames(S)<-nk
} else{
if(cdiff<lower){
#exclude variable entirely
BICS<-rdep
#mat is the matrix recording the best variable proposed, the BIC value at the end of the step, the difference between clustering versus conditional independence for that variable, whether it was an addition or removal step and whether the step was accepted or not
mat<-rbind(mat,c(colnames(S)[i],BICS,cdiff,"Remove","Accepted"))
nk<-colnames(S)[-i]
S<-as.matrix(S[,-i])
colnames(S)<-nk
} else{
m<-min(Cdiff[is.finite(Cdiff)])
i<-which(Cdiff==m,arr.ind=TRUE)[1]
mat<-rbind(mat,c(colnames(S)[i],BICS,Cdiff[i],"Remove","Rejected"))
}}
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

