"clustvarsel" <-
function(X,G,emModels1=c("E","V"),emModels2=c("EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"),samp=FALSE,sampsize=2000,allow.EEE=TRUE,forcetwo=TRUE,search="greedy",upper=0,lower=-10,itermax=100)
{
require(mclust02)
#Check whether there are variable names to identify selected variables
ifelse(is.null(colnames(X)),colnames(X)<-1:ncol(X),colnames(X)<-colnames(X))
#Apply varselection function depending on whether subsampling is required, samp=T or greedy or headlong search is required
if((samp)&search=="greedy") m<-clvarselsampgr(X,G,emModels1=emModels1, emModels2=emModels2,sampsize=sampsize,allow.EEE=allow.EEE,forcetwo=forcetwo,itermax=itermax)
if((!samp)&search=="greedy") m<-clvarselnosampgr(X,G,emModels1=emModels1,emModels2=emModels2,allow.EEE=allow.EEE,forcetwo=forcetwo,itermax=itermax)
if((samp)&search=="headlong") m<-clvarselsamphl(X,G,emModels1=emModels1,emModels2=emModels2,sampsize=sampsize,allow.EEE=allow.EEE,forcetwo=forcetwo,upper=upper,lower=lower,itermax=itermax)
if((!samp)&search=="headlong") m<-clvarselnosamphl(X,G,emModels1=emModels1,emModels2=emModels2,allow.EEE=allow.EEE,forcetwo=forcetwo,upper=upper,lower=lower,itermax=itermax)

m
}

