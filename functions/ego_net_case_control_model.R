#code to estimate initial homophily coefficients from ego network data
#Inputs are: 
#formula=formula used to run logistic regression, in the same form as ergm formulas
#possible terms are: 
#nodematch, nodemix, absdiff, nodefactor, nodecov

#data=input ego network data

#var.name.degree=name of variable in data that has the degree of the respondents

#var.name.characs= vector of names of characterisictics for ego

#var.name.characs.alter=list of variable names for alter characteristics, same order as var.name.charac

#case.control.type=type of null hypothesis used to construct the control portion
#of the data frame (the 0s in the regression), one of:
#"case.case.matching", "weighted.random.matching","one_one_matching","one_one_pair_matching"
#"weighted.random.matching" if you want to randomly pair respondents together based on prob. weights
#"case.case.matching" if want to simply match all respondents with all other respondents in case part of dataset
#"one_one_matching" if you want to match each case with only one other case for the control part of the dataset
#"one_one_pair_matching" if only want each respondent in one pair so in control part once-either
#as sender or reciever but not both

#max.control.data.N=max size used when constructing baseline random comparison
#in case control model. 

#max.alter=max number of alters the respondents were allowed to name

#control.type=either "none" or "differential.degree", if "differential.degree" then control part
#is constructed by selecting people proportionally to degree

#remove.isolates.control=should isolates be removed from control part of dataset? T/F

#weight.var.name=name of variable specifying vector of weights to create representative population  (or NULL)

#weight.var.name.control=name of variable for weights of respondents (if not NULL)
#determines the probability of selection for the control part of the dataset
#if NULL, default is to use weight.var.name 

#N.type=type of N used to set size of control part of data
#if no max.control.dyad.N, either "dyads" or "respondents" "respondents.diff"
#dyads just gives you all dyads, respondents gives you the number of respondents
#and respondents.diff gives you the number of sum(respondents*degree)

#case.type="none" or "split_sample", case.type is only used if case.control.type="one_one_pair_matching"
#split sample if you want only half of cases to be present in case control dataset
#if case.control.type="one_one_pair_matching" then will get you nested cases
#two observatoins per respondent, for half the sample
#
#num.iterations=number of differnet times should reconstruct control portion and restimate model
#note here we are not taking bootstrap sampling, just capturing stochastic variation in the 
#construction of the control portion

#bootstrap.sample=T/F, should run model multiple times using different samples each time? 

#num.bootstrap.samples=if bootstrap.sample=T then how many bootstrap samples to take?

#useparallel=T/F should employ multiple CPUs when running analyses?

#num.cores=if useparallel=T, then how many cores to utilize?

#maxit=number of maximum iterations to use in bigglm function

#nodemix.reference=vector specifying the reference category if using nodemix terms

#tie.data=data frame indicating which alters should be considered as present in the analysis
#if NULL, then all in data frame assumed to be in analysis

#adjust.intercept=T/F should attempt to adjust intercept to map onto known, actual size of full population?

#true.pop.size=size of population from which sample is drawn; only neccessary if adjust.intercept=T

#output_case_control_data=T/F should output data used to run case control regression


egonet_case_control_model=function(formula,ego_data,var.name.degree,var.name.characs,
var.name.characs.alter,
case.control.type="weighted.random.matching", max.alter=NULL,
max.control.data.N=100000, control.type=NULL, remove.isolates.control=T,
weight.var.name=NULL, weight.var.name.control=NULL, 
N.type="dyads", case.type="none",
num.iterations=10, bootstrap.sample=F, num.bootstrap.samples=NULL,
useparallel=F, num.cores=2, maxit=10,
nodemix.reference=NULL, tie.data=NULL, 
adjust.intercept=F, true.pop.size=NULL){

#getting weights for respondents
if (!is.null(weight.var.name)){
resp.weights=ego_data[,weight.var.name]
}

if (is.null(weight.var.name)){
resp.weights=NULL
}

#getting weights for contstructing control portion of case-control data frame
#if specified
if (!is.null(weight.var.name.control)){
resp.control.weights=ego_data[,weight.var.name.control]
}

if (is.null(weight.var.name.control)){
resp.control.weights=NULL
}

#getting number of dyads to use when constructing case-control data frame
#set as max.control.data.N unless this number is more than total number of pairs
#possible from ego network input data?

if (N.type=="dyads" & !is.null(max.control.data.N)){
num.dyads=nrow(ego_data)*(nrow(ego_data)-1)/2
if (num.dyads<max.control.data.N){max.control.data.N=num.dyads}
}

#getting network structure to get variable names
#

net=network.initialize(n=100, directed=F)

if (length(var.name.characs)>1){
attribute_list=do.call(list, ego_data[,var.name.characs])
for (p in 1:length(var.name.characs)){

if (is.factor(attribute_list[[p]])){
attribute_list[[p]]=as.character(attribute_list[[p]])
 }
}

net=network(x=net, vertex.attr=attribute_list) 

 }

if (length(var.name.characs)==1){
set.vertex.attribute(x=net, attrname=var.name.characs, value=ego_data[,var.name.characs])
 }


if (class(formula)!="list"){
temp_flag=length(paste(formula))
form_temp=as.formula(paste("net~", paste(formula)[[temp_flag]]))

sum.vars=summary(form_temp)
name.vars=names(sum.vars) #getting name of variables in formula

homoph=unique(c(grep("mix",name.vars),grep("match",name.vars)))

#bit of code to set up formula for case control logistic regression
#putting in reference categories where appropriate

temp.name.list=c()
form.list=c()
name.list=attr(terms.formula(form_temp),"variables")

for (i in 3:length(name.list)){
form.list[[i]]=attr(terms.formula(form_temp),"variables")[[i]]

mix=length(grep("mix",paste(form.list[[i]])))>0

if (mix){

temp.stats=sum.vars[grep("mix",name.vars)[ which(grep("mix",name.vars) %in%grep(paste(form.list[[i]])[2],name.vars))]]
temp.name=names(which(temp.stats==max(temp.stats)))[1]

if (!is.null(nodemix.reference)){
temp.name3=paste("mix.", paste(form.list[[i]])[2], ".", nodemix.reference, sep="")
temp.name=temp.name3[temp.name3 %in% names(temp.stats)]
}

temp.name2=strsplit(temp.name,paste("mix.",paste(form.list[[i]][2]),".",sep=""))[[1]][2]

form.list[[i]][3]=temp.name2 #using first category, based on alphbetical order
#or using category set as reference
temp.name.list[i]=temp.name
 }
}


form.list=paste(unlist(form.list),collapse="+")
net.name=paste(form_temp[[2]])
formula.case=as.formula(paste(net.name, "~",form.list)) #updating formula to use in case control logistic regression 
}



if (class(formula)=="list"){

formula.case.list=list()

for (p.flag in 1:length(formula)){

formula_temp=formula[[p.flag]]

temp_flag=length(paste(formula_temp))
form_temp=as.formula(paste("net~", paste(formula_temp)[[temp_flag]]))

sum.vars=summary(form_temp)
name.vars=names(sum.vars) #getting name of variables in formula

homoph=unique(c(grep("mix",name.vars),grep("match",name.vars)))

#bit of code to set up formula for case control logistic regression
#putting in reference categories where appropriate

temp.name.list=c()
form.list=c()
name.list=attr(terms.formula(form_temp),"variables")

for (i in 3:length(name.list)){
form.list[[i]]=attr(terms.formula(form_temp),"variables")[[i]]

mix=length(grep("mix",paste(form.list[[i]])))>0

if (mix){

temp.stats=sum.vars[grep("mix",name.vars)[ which(grep("mix",name.vars) %in%grep(paste(form.list[[i]])[2],name.vars))]]
temp.name=names(which(temp.stats==max(temp.stats)))[1]

if (!is.null(nodemix.reference[[p.flag]])){
temp.name3=paste("mix.", paste(form.list[[i]])[2], ".", nodemix.reference[[p.flag]], sep="")
temp.name=temp.name3[temp.name3 %in% names(temp.stats)]
}

temp.name2=strsplit(temp.name,paste("mix.",paste(form.list[[i]][2]),".",sep=""))[[1]][2]

form.list[[i]][3]=temp.name2 #using first category, based on alphbetical order
#or using category set as reference
temp.name.list[i]=temp.name
 }
}


form.list=paste(unlist(form.list),collapse="+")
net.name=paste(form_temp[[2]])
formula.case=as.formula(paste(net.name, "~",form.list)) #updating formula to use in case control logistic regression 
formula.case.list[[p.flag]]=formula.case
 }

}


#here start running over iterations
if (class(formula)!="list"){
#this runs if not using parallel processing and no bootstrap sampling
if (useparallel==F & bootstrap.sample==F){

coef_list=list()
fit_list=list()
for (p in 1:num.iterations){

#creating case control data frame
case.list=case.control.data.create.vars(data=ego_data,var.name.degree=var.name.degree,
,var.name.characs=var.name.characs, case.control.type=case.control.type,
,var.name.characs.alter=var.name.characs.alter,formula=formula.case
,resp.weights=resp.weights, max.alter=max.alter
,max.control.data.N=max.control.data.N, control.type=control.type,
remove.isolates.control=remove.isolates.control,resp.control.weights=resp.control.weights,
N.type=N.type, case.type=case.type, tie.data=tie.data)


#running logistic regression
mod1<-bigglm(case.list[[2]],data=case.list[[1]],family=binomial(link="logit"),maxit=maxit)

coefs=summary(mod1)$mat[,1]
coef_list[[p]]=coefs

#getting fit statistics
N=mod1$n
aic=mod1$deviance+2*(N-mod1$df.resid)
bic=mod1$deviance+log(N)*(N-mod1$df.resid)

fit_list[[p]]=c(N.dyads=N, N.respondents=nrow(ego_data), deviance=mod1$deviance,aic=aic,bic=bic)

 }

#Putting together results
coefs=data.frame(do.call(rbind, coef_list), check.names=F)
coefs$iteration=1:nrow(coefs)
colnames(coefs)[1]="Intercept"

fit=data.frame(do.call(rbind, fit_list))
fit=apply(fit, 2, mean)
fit[1]=round(fit[1])

####
#adjusting intercept, if specified
if (adjust.intercept==T){

#getting names of vars from case control data frame
temp.dat=case.list[[1]]
temp.dat=temp.dat[temp.dat$y==0,]
cnames=colnames(temp.dat)

coefs1=coefs[,2:(ncol(coefs)-1)]
temp.terms=colnames(coefs1)

#what to do with nodefactor terms when imputing intercept? 
#right now ignore them?

#which not nodematch, nodemix or absdiff?
which.keep.temp=c(grep("nodematch",temp.terms),grep("abs.diff",temp.terms),grep("nodemix",temp.terms),
grep("mix",temp.terms))
temp.terms=temp.terms[which.keep.temp]

#getting tables of vars from data
which.keep.cnames=c(grep("nodematch",cnames),grep("abs.diff",cnames),grep("nodemix",cnames),
grep("mix",cnames))

cnames.temp=cnames[which.keep.cnames]
#here we get the statistics of interest from the case control data frame
#for the 0s, putting in the proper names to match the coefficients
vals=list()

for (p.temp in 1:length(cnames.temp)){
tab.temp=prop.table(table(temp.dat[,cnames.temp[p.temp]]))

if (names(tab.temp)[1]=="0"|names(tab.temp)[1]=="1"){
 val=tab.temp["1"]
 names(val)=cnames.temp[p.temp]
 }

if (names(tab.temp)[1]!="0"&names(tab.temp)[1]!="1"){
 val=tab.temp
 names(val)=paste(cnames.temp[p.temp], names(val), sep="")
 }
vals[[p.temp]]=val
}

#temp.terms=temp.terms[which(temp.terms %in% colnames(case.list[[1]]))]

vals=unlist(vals)
names.vals=names(vals)
vals_mat=as.matrix(vals[temp.terms])

coefs1=as.matrix(coefs1[,temp.terms])

#adjusting coefficinets
if (length(temp.terms)>1){
plus_int_1=coefs1%*%vals_mat
}

if ((length(temp.terms)==0)){
plus_int_1=0
}

mean.deg=mean(ego_data[,var.name.degree], na.rm=T)
prob.tie=true.pop.size*mean.deg/(true.pop.size*(true.pop.size-1)/2) 
log_prob.tie=log(prob.tie /(1-prob.tie))

int1=log_prob.tie-plus_int_1

coefs[,"Intercept"]=int1
 } #adjust intercept

}

##############################################
#here not using parallel processing but are doing bootstrap sampling

if (useparallel==F & bootstrap.sample==T){

coefs_list_all=list()
fit_list_all=list()

#if using bootstrap then set weights to 1 for all cases, as already 
#used weights to put them in data.

resp.weights.temp=rep(1, nrow(ego_data))
resp.control.weights.temp=resp.weights.temp

#running over different samples, otherwise, same as above
for (x in 1:num.bootstrap.samples){

#take bootstrap sample from ego network data

ids_samp=sample(1:nrow(ego_data), size=nrow(ego_data), replace=T, prob=resp.weights)

ego_data_samp=ego_data[ids_samp,]

coef_list=list()
fit_list=list()

for (p in 1:num.iterations){

case.list=case.control.data.create.vars(data=ego_data_samp,var.name.degree=var.name.degree,
,var.name.characs=var.name.characs, case.control.type=case.control.type,
,var.name.characs.alter=var.name.characs.alter,formula=formula.case
,resp.weights=resp.weights.temp, max.alter=max.alter
,max.control.data.N=max.control.data.N, control.type=control.type,
remove.isolates.control=remove.isolates.control,resp.control.weights=resp.control.weights.temp,
N.type=N.type, case.type=case.type, tie.data=tie.data)

mod1<-bigglm(case.list[[2]],data=case.list[[1]],family=binomial(link="logit"),maxit=maxit)

coefs=summary(mod1)$mat[,1]
coef_list[[p]]=coefs

N=mod1$n
aic=mod1$deviance+2*(N-mod1$df.resid)
bic=mod1$deviance+log(N)*(N-mod1$df.resid)

fit_list[[p]]=c(N.dyads=N, N.respondents=nrow(ego_data), deviance=mod1$deviance,aic=aic,bic=bic)

 }#p

coefs_temp=data.frame(do.call(rbind, coef_list), check.names=F)
coefs_temp$iteration=1:nrow(coefs_temp)
coefs_temp$sample=x
colnames(coefs_temp)[1]="Intercept"

fit=data.frame(do.call(rbind, fit_list))
fit=apply(fit, 2, mean)
fit=c(fit, sample=x)
fit[1]=round(fit[1])

fit_list_all[[x]]=fit

if (adjust.intercept==T){

temp.dat=case.list[[1]]
temp.dat=temp.dat[temp.dat$y==0,]
cnames=colnames(temp.dat)

coefs1=coefs_temp[,2:(ncol(coefs_temp)-2)]
temp.terms=colnames(coefs1)

#what to do with nodefactor terms when imputing intercept? 
#right now ignore them?

#which not nodematch, nodemix or absdiff?
which.keep.temp=c(grep("nodematch",temp.terms),grep("abs.diff",temp.terms),grep("nodemix",temp.terms),
grep("mix",temp.terms))
temp.terms=temp.terms[which.keep.temp]

#getting tables of vars from data
which.keep.cnames=c(grep("nodematch",cnames),grep("abs.diff",cnames),grep("nodemix",cnames),
grep("mix",cnames))

cnames.temp=cnames[which.keep.cnames]
vals=list()

for (p.temp in 1:length(cnames.temp)){
tab.temp=prop.table(table(temp.dat[,cnames.temp[p.temp]]))

if (names(tab.temp)[1]=="0"|names(tab.temp)[1]=="1"){
 val=tab.temp["1"]
 names(val)=cnames.temp[p.temp]
 }

if (names(tab.temp)[1]!="0"&names(tab.temp)[1]!="1"){
 val=tab.temp
 names(val)=paste(cnames.temp[p.temp], names(val), sep="")
 }
vals[[p.temp]]=val
}

#temp.terms=temp.terms[which(temp.terms %in% colnames(case.list[[1]]))]

vals=unlist(vals)
names.vals=names(vals)
vals_mat=as.matrix(vals[temp.terms])

coefs1=as.matrix(coefs1[,temp.terms])

if (length(temp.terms)>1){
plus_int_1=coefs1%*%vals_mat
}

if ((length(temp.terms)==0)){
plus_int_1=0
}

mean.deg=mean(ego_data[,var.name.degree], na.rm=T)
prob.tie=true.pop.size*mean.deg/(true.pop.size*(true.pop.size-1)/2) 
log_prob.tie=log(prob.tie /(1-prob.tie))

int1=log_prob.tie-plus_int_1

coefs_temp[,"Intercept"]=int1
 } #adjust intercept

coefs_list_all[[x]]=coefs_temp

 } #x

coefs=do.call(rbind, coefs_list_all)
fit=do.call(rbind, fit_list_all)
}


#########################################
#here doing parallel processing but not bootstrap sampling

if (useparallel==T & bootstrap.sample==F){

#setting up parallel processing
cl <- makePSOCKcluster(num.cores)
registerDoParallel(cl)
vals.temp=1:num.iterations

temp_dat <- foreach(p = vals.temp, .export=c('case.control.data.create.vars',
'case.control.data.function', 'rep.func','rep.func2')) %dopar% {

library(biglm)

case.list=case.control.data.create.vars(data=ego_data,var.name.degree=var.name.degree,
,var.name.characs=var.name.characs, case.control.type=case.control.type,
,var.name.characs.alter=var.name.characs.alter,formula=formula.case
,resp.weights=resp.weights, max.alter=max.alter
,max.control.data.N=max.control.data.N, control.type=control.type,
remove.isolates.control=remove.isolates.control,resp.control.weights=resp.control.weights,
N.type=N.type, case.type=case.type, tie.data=tie.data)

mod1<-bigglm(case.list[[2]],data=case.list[[1]],family=binomial(link="logit"),maxit=maxit)

coefs=summary(mod1)$mat[,1]

N=mod1$n
aic=mod1$deviance+2*(N-mod1$df.resid)
bic=mod1$deviance+log(N)*(N-mod1$df.resid)

fit=c(N.dyads=N, N.respondents=nrow(ego_data), deviance=mod1$deviance,aic=aic,bic=bic)


#adjusting intercept
if (adjust.intercept==T){

temp.dat=case.list[[1]]
temp.dat=temp.dat[temp.dat$y==0,]
cnames=colnames(temp.dat)

coefs1=coefs[2:(length(coefs)-1)]
temp.terms=names(coefs1)

#what to do with nodefactor terms when imputing intercept? 
#right now ignore them?

#which not nodematch, nodemix or absdiff?
which.keep.temp=c(grep("nodematch",temp.terms),grep("abs.diff",temp.terms),grep("nodemix",temp.terms),
grep("mix",temp.terms))
temp.terms=temp.terms[which.keep.temp]

#getting tables of vars from data
which.keep.cnames=c(grep("nodematch",cnames),grep("abs.diff",cnames),grep("nodemix",cnames),
grep("mix",cnames))

cnames.temp=cnames[which.keep.cnames]
vals=list()

for (p.temp in 1:length(cnames.temp)){
tab.temp=prop.table(table(temp.dat[,cnames.temp[p.temp]]))

if (names(tab.temp)[1]=="0"|names(tab.temp)[1]=="1"){
 val=tab.temp["1"]
 names(val)=cnames.temp[p.temp]
 }

if (names(tab.temp)[1]!="0"&names(tab.temp)[1]!="1"){
 val=tab.temp
 names(val)=paste(cnames.temp[p.temp], names(val), sep="")
 }
vals[[p.temp]]=val
}

#temp.terms=temp.terms[which(temp.terms %in% colnames(case.list[[1]]))]

vals=unlist(vals)
names.vals=names(vals)
vals_mat=as.matrix(vals[temp.terms])

coefs1=matrix(coefs1[temp.terms], nrow=1)

if (length(temp.terms)>1){
plus_int_1=coefs1%*%vals_mat
}

if ((length(temp.terms)==0)){
plus_int_1=0
}

mean.deg=mean(ego_data[,var.name.degree], na.rm=T)
prob.tie=true.pop.size*mean.deg/(true.pop.size*(true.pop.size-1)/2) 
log_prob.tie=log(prob.tie /(1-prob.tie))

int1=log_prob.tie-plus_int_1

coefs[1]=int1
 } #adjust intercept

###

list(coefs, fit)
}

#ending cluster and then putting things together
stopCluster(cl)

coefs=lapply(temp_dat, '[[', 1)
coefs=data.frame(do.call(rbind, coefs))
coefs$iteration=1:nrow(coefs)
colnames(coefs)[1]="Intercept"


#fit=temp_dat[[1]][[2]]
fit=lapply(temp_dat, '[[', 2)
fit=data.frame(do.call(rbind, fit))
fit=apply(fit, 2, mean)
fit[1]=round(fit[1])

}

################################################
#here doing parallel processing with bootstrap sampling

if (useparallel==T & bootstrap.sample==T){

#if using bootstrap then set weights to 1 for all cases, as already 
#used weights to put them in data

resp.weights.temp=rep(1, nrow(ego_data))
resp.control.weights.temp=resp.weights.temp

#setting up parallel processing
cl <- makePSOCKcluster(num.cores)
registerDoParallel(cl)
vals.temp=1:num.bootstrap.samples

temp_dat <- foreach(x = vals.temp, .export=c('case.control.data.create.vars',
'case.control.data.function','rep.func','rep.func2')) %dopar% {

library(biglm)

#take bootstrap sample from ego network data

ids_samp=sample(1:nrow(ego_data), size=nrow(ego_data), replace=T, prob=resp.weights)
ego_data_samp=ego_data[ids_samp,]

coef_list=list()
fit_list=list()

for (p in 1:num.iterations){

case.list=case.control.data.create.vars(data=ego_data_samp,var.name.degree=var.name.degree,
,var.name.characs=var.name.characs, case.control.type=case.control.type,
,var.name.characs.alter=var.name.characs.alter,formula=formula.case
,resp.weights=resp.weights.temp, max.alter=max.alter
,max.control.data.N=max.control.data.N, control.type=control.type,
remove.isolates.control=remove.isolates.control,resp.control.weights=resp.control.weights.temp,
N.type=N.type, case.type=case.type, tie.data=tie.data)

mod1<-bigglm(case.list[[2]],data=case.list[[1]],family=binomial(link="logit"),maxit=maxit)

coefs=summary(mod1)$mat[,1]
coef_list[[p]]=coefs

N=mod1$n
aic=mod1$deviance+2*(N-mod1$df.resid)
bic=mod1$deviance+log(N)*(N-mod1$df.resid)

fit_list[[p]]=c(N.dyads=N, N.respondents=nrow(ego_data), deviance=mod1$deviance,aic=aic,bic=bic)

 }#p

coefs_temp=data.frame(do.call(rbind, coef_list), check.names=F)
coefs_temp$iteration=1:nrow(coefs_temp)
coefs_temp$sample=x
colnames(coefs_temp)[1]="Intercept"

###adjust intercept
if (adjust.intercept==T){

temp.dat=case.list[[1]]
temp.dat=temp.dat[temp.dat$y==0,]
cnames=colnames(temp.dat)

coefs1=coefs_temp[,2:(ncol(coefs_temp)-2)]
temp.terms=colnames(coefs1)

#what to do with nodefactor terms when imputing intercept? 
#right now ignore them?

#which not nodematch, nodemix or absdiff?
which.keep.temp=c(grep("nodematch",temp.terms),grep("abs.diff",temp.terms),grep("nodemix",temp.terms),
grep("mix",temp.terms))
temp.terms=temp.terms[which.keep.temp]

#getting tables of vars from data
which.keep.cnames=c(grep("nodematch",cnames),grep("abs.diff",cnames),grep("nodemix",cnames),
grep("mix",cnames))

cnames.temp=cnames[which.keep.cnames]
vals=list()

for (p.temp in 1:length(cnames.temp)){
tab.temp=prop.table(table(temp.dat[,cnames.temp[p.temp]]))

if (names(tab.temp)[1]=="0"|names(tab.temp)[1]=="1"){
 val=tab.temp["1"]
 names(val)=cnames.temp[p.temp]
 }

if (names(tab.temp)[1]!="0"&names(tab.temp)[1]!="1"){
 val=tab.temp
 names(val)=paste(cnames.temp[p.temp], names(val), sep="")
 }
vals[[p.temp]]=val
}

#temp.terms=temp.terms[which(temp.terms %in% colnames(case.list[[1]]))]

vals=unlist(vals)
names.vals=names(vals)
vals_mat=as.matrix(vals[temp.terms])

coefs1=as.matrix(coefs1[,temp.terms])

if (length(temp.terms)>1){
plus_int_1=coefs1%*%vals_mat
}

if ((length(temp.terms)==0)){
plus_int_1=0
}

mean.deg=mean(ego_data[,var.name.degree], na.rm=T)
prob.tie=true.pop.size*mean.deg/(true.pop.size*(true.pop.size-1)/2) 
log_prob.tie=log(prob.tie /(1-prob.tie))

int1=log_prob.tie-plus_int_1

coefs_temp[,"Intercept"]=int1
 } #adjust intercept


###

fit=data.frame(do.call(rbind, fit_list))
fit=apply(fit, 2, mean)
fit=c(fit, sample=x)
fit[1]=round(fit[1])

list(coefs_temp, fit)
 } #x
 
#end cluster and putting things together
stopCluster(cl)

coefs=lapply(temp_dat, '[[', 1)
coefs=data.frame(do.call(rbind, coefs))
colnames(coefs)[1]="Intercept"

#fit=temp_dat[[1]][[2]]
fit=lapply(temp_dat, '[[', 2)
fit=data.frame(do.call(rbind, fit))


 }
} #if only one formula to run over
###
###########################################################################################
if (class(formula)=="list"){

#this runs if not using parallel processing and no bootstrap sampling
if (useparallel==F & bootstrap.sample==F){


coef_list=list()
fit_list=list()

for (mod.num in 1:length(formula.case.list)){

coef_list[[mod.num]]=list()
fit_list[[mod.num]]=list()
}

for (p in 1:num.iterations){

#creating case control data frame
#with multiple formulas, assumign that first formula is used to construct base data?
case.list=case.control.data.create.vars(data=ego_data,var.name.degree=var.name.degree,
,var.name.characs=var.name.characs, case.control.type=case.control.type,
,var.name.characs.alter=var.name.characs.alter,formula=formula.case.list[[1]]
,resp.weights=resp.weights, max.alter=max.alter
,max.control.data.N=max.control.data.N, control.type=control.type,
remove.isolates.control=remove.isolates.control,resp.control.weights=resp.control.weights,
N.type=N.type, case.type=case.type, tie.data=tie.data)

form.list2=list()
form.glm=case.list[[2]]

form.list2[[1]]=form.glm
terms.list =list()

if (length(formula.case.list) >1){

for (temp.flag in 1:length(formula.case.list) ){

# base.formula=formula.case.list[[temp.flag]]

 base.formula=formula[[temp.flag]]

#if (!check.case.control.data) {
#form.list2[[temp.flag]]=output_glm_formula(base.formula= base.formula,var.name.characs=var.name.characs,
#check.case.control.data=T,
#case.list=case.list)[[1]]
#}

case.out=output_glm_formula(base.formula= base.formula,var.name.characs=var.name.characs
,check.case.control.data=T, case.list=case.list)

form.list2[[temp.flag]]=case.out[[1]]

case.list[[1]]=case.out[[2]]
  
terms.list[[temp.flag]]=colnames(attr(terms(form.list2[[temp.flag]]),"factors"))

 } #temp.flag


terms=unique(unlist(terms.list))

 } #if form.list more than one formula


full.data=case.list[[1]]

if (length(terms)>1){
mis <- apply(is.na(full.data[,terms]), 1, any)
}

if (length(terms)==1){
mis <- is.na(full.data[,terms])
}

full.data=full.data[!mis,]


for (mod.num in 1:length(form.list2)){

current.form=form.list2[[mod.num]]


#running logistic regression
mod1<-bigglm(current.form, data=full.data, family=binomial(link="logit"), maxit=maxit)

coefs=summary(mod1)$mat[,1]
coef_list[[mod.num]][[p]]=coefs

#getting fit statistics
N=mod1$n
aic=mod1$deviance+2*(N-mod1$df.resid)
bic=mod1$deviance+log(N)*(N-mod1$df.resid)

fit_list[[mod.num]][[p]]=c(N.dyads=N, N.respondents=nrow(ego_data), deviance=mod1$deviance,aic=aic,bic=bic)

 }#mod.num
}#p

#Putting together results
coefs_temp1=list()
fit_temp1=list()

for (mod.num in 1:length(coef_list)){
coefs=data.frame(do.call(rbind, coef_list[[mod.num]]), check.names=F)
coefs$iteration=1:nrow(coefs)
colnames(coefs)[1]="Intercept"

fit=data.frame(do.call(rbind, fit_list[[mod.num]]))
fit=apply(fit, 2, mean)
fit[1]=round(fit[1])

coefs_temp1[[mod.num]]=coefs
fit_temp1[[mod.num]]=fit
}

####
#adjusting intercept, if specified
if (adjust.intercept==T){

for (mod.num in 1:length(coef_list)){

#getting names of vars from case control data frame
temp.dat=case.list[[1]]
temp.dat=temp.dat[temp.dat$y==0,]
cnames=colnames(temp.dat)

coefs1=coefs_temp1[[mod.num]][,2:(ncol(coefs_temp1[[mod.num]])-1)]
temp.terms=colnames(coefs1)

#what to do with nodefactor terms when imputing intercept? 
#right now ignore them?

#which not nodematch, nodemix or absdiff?
which.keep.temp=c(grep("nodematch",temp.terms),grep("abs.diff",temp.terms),grep("nodemix",temp.terms),
grep("mix",temp.terms))
temp.terms=temp.terms[which.keep.temp]

#getting tables of vars from data
which.keep.cnames=c(grep("nodematch",cnames),grep("abs.diff",cnames),grep("nodemix",cnames),
grep("mix",cnames))

cnames.temp=cnames[which.keep.cnames]
#here we get the statistics of interest from the case control data frame
#for the 0s, putting in the proper names to match the coefficients
vals=list()

for (p.temp in 1:length(cnames.temp)){
tab.temp=prop.table(table(temp.dat[,cnames.temp[p.temp]]))

if (names(tab.temp)[1]=="0"|names(tab.temp)[1]=="1"){
 val=tab.temp["1"]
 names(val)=cnames.temp[p.temp]
 }

if (names(tab.temp)[1]!="0"&names(tab.temp)[1]!="1"){
 val=tab.temp
 names(val)=paste(cnames.temp[p.temp], names(val), sep="")
 }
vals[[p.temp]]=val
}

#temp.terms=temp.terms[which(temp.terms %in% colnames(case.list[[1]]))]

vals=unlist(vals)
names.vals=names(vals)
vals_mat=as.matrix(vals[temp.terms])

coefs1=as.matrix(coefs1[,temp.terms])

#adjusting coefficinets
if (length(temp.terms)>1){
plus_int_1=coefs1%*%vals_mat
}

if ((length(temp.terms)==0)){
plus_int_1=0
}

mean.deg=mean(ego_data[,var.name.degree], na.rm=T)
prob.tie=true.pop.size*mean.deg/(true.pop.size*(true.pop.size-1)/2) 
log_prob.tie=log(prob.tie /(1-prob.tie))

int1=log_prob.tie-plus_int_1

#coefs[,"Intercept"]=int1
coefs_temp1[[mod.num]][,"Intercept"]=int1
  } 
 }#adjust intercept
coefs=coefs_temp1
fit=fit_temp1
}

##############################################
#here not using parallel processing but are doing bootstrap sampling

if (useparallel==F & bootstrap.sample==T){

coefs_list_all=list()
fit_list_all=list()

#if using bootstrap then set weights to 1 for all cases, as already 
#used weights to put them in data.

resp.weights.temp=rep(1, nrow(ego_data))
resp.control.weights.temp=resp.weights.temp

#running over different samples, otherwise, same as above
for (x in 1:num.bootstrap.samples){

#take bootstrap sample from ego network data

ids_samp=sample(1:nrow(ego_data), size=nrow(ego_data), replace=T, prob=resp.weights)

ego_data_samp=ego_data[ids_samp,]

coef_list=list()
fit_list=list()

for (mod.num in 1:length(formula.case.list)){

coef_list[[mod.num]]=list()
fit_list[[mod.num]]=list()
}

for (p in 1:num.iterations){

case.list=case.control.data.create.vars(data=ego_data_samp,var.name.degree=var.name.degree,
,var.name.characs=var.name.characs, case.control.type=case.control.type,
,var.name.characs.alter=var.name.characs.alter,formula=formula.case.list[[1]]
,resp.weights=resp.weights.temp, max.alter=max.alter
,max.control.data.N=max.control.data.N, control.type=control.type,
remove.isolates.control=remove.isolates.control,resp.control.weights=resp.control.weights.temp,
N.type=N.type, case.type=case.type, tie.data=tie.data)

###
form.list2=list()
form.glm=case.list[[2]]

form.list2[[1]]=form.glm
terms.list =list()

if (length(formula.case.list) >1){

for (temp.flag in 1:length(formula.case.list) ){

# base.formula=formula.case.list[[temp.flag]]

 base.formula=formula[[temp.flag]]

#if (!check.case.control.data) {
#form.list2[[temp.flag]]=output_glm_formula(base.formula= base.formula,var.name.characs=var.name.characs,
#check.case.control.data=T,
#case.list=case.list)[[1]]
#}

case.out=output_glm_formula(base.formula= base.formula,var.name.characs=var.name.characs
,check.case.control.data=T, case.list=case.list)

form.list2[[temp.flag]]=case.out[[1]]

case.list[[1]]=case.out[[2]]
  
terms.list[[temp.flag]]=colnames(attr(terms(form.list2[[temp.flag]]),"factors"))

 } #temp.flag


terms=unique(unlist(terms.list))

 } #if form.list more than one formula


full.data=case.list[[1]]

if (length(terms)>1){
mis <- apply(is.na(full.data[,terms]), 1, any)
}

if (length(terms)==1){
mis <- is.na(full.data[,terms])
}

full.data=full.data[!mis,]

###
for (mod.num in 1:length(form.list2)){

current.form=form.list2[[mod.num]]

#running logistic regression
mod1<-bigglm(current.form, data=full.data, family=binomial(link="logit"), maxit=maxit)

coefs=summary(mod1)$mat[,1]
coef_list[[mod.num]][[p]]=coefs

#getting fit statistics
N=mod1$n
aic=mod1$deviance+2*(N-mod1$df.resid)
bic=mod1$deviance+log(N)*(N-mod1$df.resid)

fit_list[[mod.num]][[p]]=c(N.dyads=N, N.respondents=nrow(ego_data), deviance=mod1$deviance,aic=aic,bic=bic)

 }#mod.num
}#p



coefs_temp1=list()
fit_temp1=list()

for (mod.num in 1:length(coef_list)){
coefs_temp=data.frame(do.call(rbind, coef_list[[mod.num]]), check.names=F)
coefs_temp$iteration=1:nrow(coefs_temp)
coefs_temp$sample=x
colnames(coefs_temp)[1]="Intercept"

fit=data.frame(do.call(rbind, fit_list[[mod.num]]))
fit=apply(fit, 2, mean)
fit=c(fit, sample=x)
fit[1]=round(fit[1])

coefs_temp1[[mod.num]]=coefs_temp
fit_temp1[[mod.num]]=fit
}

fit_list_all[[x]]=fit_temp1

if (adjust.intercept==T){

for (mod.num in 1:length(coefs_temp1)){

temp.dat=case.list[[1]]
temp.dat=temp.dat[temp.dat$y==0,]
cnames=colnames(temp.dat)

coefs1=coefs_temp1[[mod.num]][,2:(ncol(coefs_temp1[[mod.num]])-2)]
temp.terms=colnames(coefs1)

#what to do with nodefactor terms when imputing intercept? 
#right now ignore them?

#which not nodematch, nodemix or absdiff?
which.keep.temp=c(grep("nodematch",temp.terms),grep("abs.diff",temp.terms),grep("nodemix",temp.terms),
grep("mix",temp.terms))
temp.terms=temp.terms[which.keep.temp]

#getting tables of vars from data
which.keep.cnames=c(grep("nodematch",cnames),grep("abs.diff",cnames),grep("nodemix",cnames),
grep("mix",cnames))

cnames.temp=cnames[which.keep.cnames]
vals=list()

for (p.temp in 1:length(cnames.temp)){
tab.temp=prop.table(table(temp.dat[,cnames.temp[p.temp]]))

if (names(tab.temp)[1]=="0"|names(tab.temp)[1]=="1"){
 val=tab.temp["1"]
 names(val)=cnames.temp[p.temp]
 }

if (names(tab.temp)[1]!="0"&names(tab.temp)[1]!="1"){
 val=tab.temp
 names(val)=paste(cnames.temp[p.temp], names(val), sep="")
 }
vals[[p.temp]]=val
}

#temp.terms=temp.terms[which(temp.terms %in% colnames(case.list[[1]]))]

vals=unlist(vals)
names.vals=names(vals)
vals_mat=as.matrix(vals[temp.terms])

coefs1=as.matrix(coefs1[,temp.terms])

if (length(temp.terms)>1){
plus_int_1=coefs1%*%vals_mat
}

if ((length(temp.terms)==0)){
plus_int_1=0
}

mean.deg=mean(ego_data[,var.name.degree], na.rm=T)
prob.tie=true.pop.size*mean.deg/(true.pop.size*(true.pop.size-1)/2) 
log_prob.tie=log(prob.tie /(1-prob.tie))

int1=log_prob.tie-plus_int_1

coefs_temp1[[mod.num]][,"Intercept"]=int1
  } #mod.num
 } #adjust intercept

coefs_list_all[[x]]=coefs_temp1

 } #x

#coefs=do.call(rbind, coefs_list_all)
#fit=do.call(rbind, fit_list_all)

coefs=list()
fit=list()

for (mod.num in 1:length(formula.case.list)){
coefs.m=lapply(coefs_list_all, '[[', mod.num)
coefs.m=data.frame(do.call(rbind, coefs.m))
coefs[[mod.num]]=coefs.m

fit1=lapply(fit_list_all, '[[', mod.num)
fit1=data.frame(do.call(rbind, fit1))
fit[[mod.num]]=fit1

 }

}

#########################################
#here doing parallel processing but not bootstrap sampling

if (useparallel==T & bootstrap.sample==F){

#coef_list=list()
#fit_list=list()

#for (mod.num in 1:length(formula.case.list)){

#coef_list[[mod.num]]=list()
#fit_list[[mod.num]]=list()
#}

#setting up parallel processing
cl <- makePSOCKcluster(num.cores)
registerDoParallel(cl)
vals.temp=1:num.iterations

temp_dat <- foreach(p = vals.temp, .export=c('case.control.data.create.vars',
'case.control.data.function','output_glm_formula','rep.func', 'rep.func2')) %dopar% {

library(biglm)

case.list=case.control.data.create.vars(data=ego_data,var.name.degree=var.name.degree,
,var.name.characs=var.name.characs, case.control.type=case.control.type,
,var.name.characs.alter=var.name.characs.alter,formula=formula.case.list[[1]]
,resp.weights=resp.weights, max.alter=max.alter
,max.control.data.N=max.control.data.N, control.type=control.type,
remove.isolates.control=remove.isolates.control,resp.control.weights=resp.control.weights,
N.type=N.type, case.type=case.type, tie.data=tie.data)


form.list2=list()
form.glm=case.list[[2]]

form.list2[[1]]=form.glm
terms.list =list()

if (length(formula.case.list) >1){

for (temp.flag in 1:length(formula.case.list) ){

# base.formula=formula.case.list[[temp.flag]]

 base.formula=formula[[temp.flag]]

#if (!check.case.control.data) {
#form.list2[[temp.flag]]=output_glm_formula(base.formula= base.formula,var.name.characs=var.name.characs,
#check.case.control.data=T,
#case.list=case.list)[[1]]
#}

case.out=output_glm_formula(base.formula= base.formula,var.name.characs=var.name.characs
,check.case.control.data=T, case.list=case.list)

form.list2[[temp.flag]]=case.out[[1]]

case.list[[1]]=case.out[[2]]
  
terms.list[[temp.flag]]=colnames(attr(terms(form.list2[[temp.flag]]),"factors"))

 } #temp.flag


terms=unique(unlist(terms.list))

 } #if form.list more than one formula


full.data=case.list[[1]]

if (length(terms)>1){
mis <- apply(is.na(full.data[,terms]), 1, any)
}

if (length(terms)==1){
mis <- is.na(full.data[,terms])
}

full.data=full.data[!mis,]

coef_list=list()
fit_list=list()

for (mod.num in 1:length(form.list2)){

current.form=form.list2[[mod.num]]


#running logistic regression
mod1<-bigglm(current.form, data=full.data, family=binomial(link="logit"), maxit=maxit)

coefs=summary(mod1)$mat[,1]
coef_list[[mod.num]]=coefs

#getting fit statistics
N=mod1$n
aic=mod1$deviance+2*(N-mod1$df.resid)
bic=mod1$deviance+log(N)*(N-mod1$df.resid)

fit_list[[mod.num]]=c(N.dyads=N, N.respondents=nrow(ego_data), deviance=mod1$deviance,aic=aic,bic=bic)

 }#mod.num


#adjusting intercept
if (adjust.intercept==T){

for (mod.num in 1:length(coef_list)){

temp.dat=case.list[[1]]
temp.dat=temp.dat[temp.dat$y==0,]
cnames=colnames(temp.dat)

coefs1=coef_list[[mod.num]][2:(length(coef_list[[mod.num]])-1)]
temp.terms=names(coefs1)

#what to do with nodefactor terms when imputing intercept? 
#right now ignore them?

#which not nodematch, nodemix or absdiff?
which.keep.temp=c(grep("nodematch",temp.terms),grep("abs.diff",temp.terms),grep("nodemix",temp.terms),
grep("mix",temp.terms))
temp.terms=temp.terms[which.keep.temp]

#getting tables of vars from data
which.keep.cnames=c(grep("nodematch",cnames),grep("abs.diff",cnames),grep("nodemix",cnames),
grep("mix",cnames))

cnames.temp=cnames[which.keep.cnames]
vals=list()

for (p.temp in 1:length(cnames.temp)){
tab.temp=prop.table(table(temp.dat[,cnames.temp[p.temp]]))

if (names(tab.temp)[1]=="0"|names(tab.temp)[1]=="1"){
 val=tab.temp["1"]
 names(val)=cnames.temp[p.temp]
 }

if (names(tab.temp)[1]!="0"&names(tab.temp)[1]!="1"){
 val=tab.temp
 names(val)=paste(cnames.temp[p.temp], names(val), sep="")
 }
vals[[p.temp]]=val
}

#temp.terms=temp.terms[which(temp.terms %in% colnames(case.list[[1]]))]

vals=unlist(vals)
names.vals=names(vals)
vals_mat=as.matrix(vals[temp.terms])

coefs1=matrix(coefs1[temp.terms], nrow=1)

if (length(temp.terms)>1){
plus_int_1=coefs1%*%vals_mat
}

if ((length(temp.terms)==0)){
plus_int_1=0
}

mean.deg=mean(ego_data[,var.name.degree], na.rm=T)
prob.tie=true.pop.size*mean.deg/(true.pop.size*(true.pop.size-1)/2) 
log_prob.tie=log(prob.tie /(1-prob.tie))

int1=log_prob.tie-plus_int_1

coef_list[[mod.num]][1]=int1
   }
 } #adjust intercept

###

list(coef_list, fit_list)
}

#ending cluster and then putting things together
stopCluster(cl)

coefs_list2=list()
fit_list2=list()
coefs=lapply(temp_dat, '[[', 1)
fit=lapply(temp_dat, '[[', 2)


for (mod.num in 1:length(formula.case.list)){
temp=lapply(coefs, '[[', mod.num)
temp=data.frame(do.call(rbind, temp))
temp$iteration=1:nrow(temp)
colnames(temp)[1]="Intercept"
coefs_list2[[mod.num]]=temp


temp_fit=lapply(fit, '[[', mod.num)
temp_fit=data.frame(do.call(rbind, temp_fit))
temp_fit=apply(temp_fit, 2, mean)
temp_fit[1]=round(temp_fit[1])

fit_list2[[mod.num]]=temp_fit
}


coefs=coefs_list2
fit=fit_list2
}

################################################
#here doing parallel processing with bootstrap sampling

if (useparallel==T & bootstrap.sample==T){

#if using bootstrap then set weights to 1 for all cases, as already 
#used weights to put them in data

resp.weights.temp=rep(1, nrow(ego_data))
resp.control.weights.temp=resp.weights.temp

#setting up parallel processing
cl <- makePSOCKcluster(num.cores)
registerDoParallel(cl)
vals.temp=1:num.bootstrap.samples

temp_dat <- foreach(x = vals.temp, .export=c('case.control.data.create.vars',
'case.control.data.function','output_glm_formula','rep.func','rep.func2')) %dopar% {

library(biglm)

#take bootstrap sample from ego network data

ids_samp=sample(1:nrow(ego_data), size=nrow(ego_data), replace=T, prob=resp.weights)
ego_data_samp=ego_data[ids_samp,]

coef_list=list()
fit_list=list()

for (mod.num in 1:length(formula.case.list)){

coef_list[[mod.num]]=list()
fit_list[[mod.num]]=list()
}


for (p in 1:num.iterations){

case.list=case.control.data.create.vars(data=ego_data_samp,var.name.degree=var.name.degree,
,var.name.characs=var.name.characs, case.control.type=case.control.type,
,var.name.characs.alter=var.name.characs.alter,formula=formula.case.list[[1]]
,resp.weights=resp.weights.temp, max.alter=max.alter
,max.control.data.N=max.control.data.N, control.type=control.type,
remove.isolates.control=remove.isolates.control,resp.control.weights=resp.control.weights.temp,
N.type=N.type, case.type=case.type, tie.data=tie.data)

###
form.list2=list()
form.glm=case.list[[2]]

form.list2[[1]]=form.glm
terms.list =list()

if (length(formula.case.list) >1){

for (temp.flag in 1:length(formula.case.list) ){

# base.formula=formula.case.list[[temp.flag]]

 base.formula=formula[[temp.flag]]

#if (!check.case.control.data) {
#form.list2[[temp.flag]]=output_glm_formula(base.formula= base.formula,var.name.characs=var.name.characs,
#check.case.control.data=T,
#case.list=case.list)[[1]]
#}

case.out=output_glm_formula(base.formula= base.formula,var.name.characs=var.name.characs
,check.case.control.data=T, case.list=case.list)

form.list2[[temp.flag]]=case.out[[1]]

case.list[[1]]=case.out[[2]]
  
terms.list[[temp.flag]]=colnames(attr(terms(form.list2[[temp.flag]]),"factors"))

 } #temp.flag


terms=unique(unlist(terms.list))

 } #if form.list more than one formula


full.data=case.list[[1]]

if (length(terms)>1){
mis <- apply(is.na(full.data[,terms]), 1, any)
}

if (length(terms)==1){
mis <- is.na(full.data[,terms])
}

full.data=full.data[!mis,]

###
for (mod.num in 1:length(form.list2)){

current.form=form.list2[[mod.num]]

#running logistic regression
mod1<-bigglm(current.form, data=full.data, family=binomial(link="logit"), maxit=maxit)

coefs=summary(mod1)$mat[,1]
coef_list[[mod.num]][[p]]=coefs

#getting fit statistics
N=mod1$n
aic=mod1$deviance+2*(N-mod1$df.resid)
bic=mod1$deviance+log(N)*(N-mod1$df.resid)

fit_list[[mod.num]][[p]]=c(N.dyads=N, N.respondents=nrow(ego_data), deviance=mod1$deviance,aic=aic,bic=bic)

 }#mod.num
}#p


coefs_temp1=list()

for (mod.num in 1:length(coef_list)){
coefs_temp=data.frame(do.call(rbind, coef_list[[mod.num]]), check.names=F)
coefs_temp$iteration=1:nrow(coefs_temp)
coefs_temp$sample=x
colnames(coefs_temp)[1]="Intercept"
coefs_temp1[[mod.num]]=coefs_temp
}

coefs_temp=data.frame(do.call(rbind, coef_list), check.names=F)
coefs_temp$iteration=1:nrow(coefs_temp)
coefs_temp$sample=x
colnames(coefs_temp)[1]="Intercept"

###adjust intercept
if (adjust.intercept==T){

for (mod.num in 1:length(coefs_temp1)){

temp.dat=case.list[[1]]
temp.dat=temp.dat[temp.dat$y==0,]
cnames=colnames(temp.dat)

#coefs1=coefs_temp[,2:(ncol(coefs_temp)-2)]
coefs1=coefs_temp1[[mod.num]][,2:(ncol(coefs_temp1[[mod.num]])-2)]

temp.terms=colnames(coefs1)

#what to do with nodefactor terms when imputing intercept? 
#right now ignore them?

#which not nodematch, nodemix or absdiff?
which.keep.temp=c(grep("nodematch",temp.terms),grep("abs.diff",temp.terms),grep("nodemix",temp.terms),
grep("mix",temp.terms))
temp.terms=temp.terms[which.keep.temp]

#getting tables of vars from data
which.keep.cnames=c(grep("nodematch",cnames),grep("abs.diff",cnames),grep("nodemix",cnames),
grep("mix",cnames))

cnames.temp=cnames[which.keep.cnames]
vals=list()

for (p.temp in 1:length(cnames.temp)){
tab.temp=prop.table(table(temp.dat[,cnames.temp[p.temp]]))

if (names(tab.temp)[1]=="0"|names(tab.temp)[1]=="1"){
 val=tab.temp["1"]
 names(val)=cnames.temp[p.temp]
 }

if (names(tab.temp)[1]!="0"&names(tab.temp)[1]!="1"){
 val=tab.temp
 names(val)=paste(cnames.temp[p.temp], names(val), sep="")
 }
vals[[p.temp]]=val
}

#temp.terms=temp.terms[which(temp.terms %in% colnames(case.list[[1]]))]

vals=unlist(vals)
names.vals=names(vals)
vals_mat=as.matrix(vals[temp.terms])

coefs1=as.matrix(coefs1[,temp.terms])

if (length(temp.terms)>1){
plus_int_1=coefs1%*%vals_mat
}

if ((length(temp.terms)==0)){
plus_int_1=0
}

mean.deg=mean(ego_data[,var.name.degree], na.rm=T)
prob.tie=true.pop.size*mean.deg/(true.pop.size*(true.pop.size-1)/2) 
log_prob.tie=log(prob.tie /(1-prob.tie))

int1=log_prob.tie-plus_int_1

coefs_temp1[[mod.num]][,"Intercept"]=int1
  } #mod.num

 } #adjust intercept


###


#
fit_temp1=list()

for (mod.num in 1:length(fit_list)){
fit=data.frame(do.call(rbind, fit_list[[mod.num]]))
fit=apply(fit, 2, mean)
fit=c(fit, sample=x)
fit[1]=round(fit[1])
fit_temp1[[mod.num]]=fit
}

list(coefs_temp1, fit_temp1)
 } #x
 
#end cluster and putting things together
stopCluster(cl)

coefs=lapply(temp_dat, '[[', 1)
coefs=data.frame(do.call(rbind, coefs))
colnames(coefs)[1]="Intercept"

#fit=temp_dat[[1]][[2]]
fit=lapply(temp_dat, '[[', 2)
fit=data.frame(do.call(rbind, fit))

##
coefs_list2=list()
fit_list2=list()
coefs=lapply(temp_dat, '[[', 1)
fit=lapply(temp_dat, '[[', 2)


for (mod.num in 1:length(formula.case.list)){
temp=lapply(coefs, '[[', mod.num)
temp=data.frame(do.call(rbind, temp))
colnames(temp)[1]="Intercept"
coefs_list2[[mod.num]]=temp


temp_fit=lapply(fit, '[[', mod.num)
temp_fit=data.frame(do.call(rbind, temp_fit))
fit_list2[[mod.num]]=temp_fit
}


coefs=coefs_list2
fit=fit_list2


 }

} #if multiple formulas
###################################################
list(coefs=coefs, fit=fit)


} #end of function
