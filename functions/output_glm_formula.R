
output_glm_formula=function(base.formula,var.name.characs,check.case.control.data=T,
case.list=NULL,egodata=NULL){


#figuring out which terms are in secondary formulas 
var.name.characs.temp=var.name.characs
temp.TF=NA

for (x.temp2 in 1:length(var.name.characs)){

temp.TF[x.temp2]=length(grep(var.name.characs[x.temp2],attr(terms.formula(base.formula),"variables")))>0
}

var.name.characs.temp=var.name.characs.temp[temp.TF==T]

##########################
var.type=NA

for (i in 1:length(var.name.characs.temp)){

ids=grep(var.name.characs.temp[i],attr(terms.formula(base.formula),"variables"))
abs=grep("absdiff",attr(terms.formula(base.formula),"variables"))
match=grep("nodematch",attr(terms.formula(base.formula),"variables"))
mix=grep("nodemix",attr(terms.formula(base.formula),"variables"))
id.type=max(c(1*as.numeric(sum(ids %in%abs)>0)
,2*as.numeric(sum(ids %in%match)>0)
,3*as.numeric(sum(ids %in%mix)>0)
))

var.type[i]=switch(paste(id.type),"0"=NA,"1"="absdiff","2"="nodematch","3"="nodemix")
}


base.name=list()

for (i in 1:length(var.name.characs.temp)){
if (var.type[i]%in% c("absdiff", "nodematch")){base.name[[i]]=NULL}
if (var.type[i]%in% c("nodemix")){
ids=grep(var.name.characs.temp[i],attr(terms.formula(base.formula),"variables"))
mix=grep("nodemix",attr(terms.formula(base.formula),"variables"))

ids=ids[ids %in% mix]

charac.mix=as.character(attr(terms.formula(base.formula),"variables")[[ids]])
if(length(charac.mix)==3){base.name[[i]]=charac.mix[[3]]}
#if(length(charac.mix)!=3){base.name[[i]]=NULL}

 }
}

if (length(base.name)<i ){base.name[[i+1]]="hold"}


var.nodecov=NA
for ( i in 1:length(var.name.characs.temp)){
ids=grep(var.name.characs[i],attr(terms.formula(base.formula),"variables"))
node=grep("nodecov",attr(terms.formula(base.formula),"variables"))
ids=ids[ids %in% node]

var.nodecov[i]=ifelse(length(ids)>0,TRUE,FALSE)
}


var.nodefactor=NA
for ( i in 1:length(var.name.characs.temp)){
ids=grep(var.name.characs[i],attr(terms.formula(base.formula),"variables"))
node=grep("nodefactor",attr(terms.formula(base.formula),"variables"))
ids=ids[ids %in% node]

var.nodefactor[i]=ifelse(length(ids)>0,TRUE,FALSE)
}


base.name.factor=list()

for (i in 1:length(var.name.characs.temp)){
if (var.nodefactor[i] %in% FALSE){base.name.factor[[i]]=NULL}
if (var.nodefactor[i]%in% c(TRUE)){
ids=grep(var.name.characs.temp[i],attr(terms.formula(base.formula),"variables"))
node=grep("nodefactor",attr(terms.formula(base.formula),"variables"))

ids=ids[ids %in% node]

charac.mix=as.character(attr(terms.formula(base.formula),"variables")[[ids]])
if(length(charac.mix)==3){base.name.factor[[i]]=charac.mix[[3]]}

 }
}

if (length(base.name.factor)<i ){base.name.factor[[i+1]]="hold"}

#miss.ids=unlist(lapply(lapply(var.name.characs.alter,is.na),sum))
#miss.ids=which(miss.ids==1)


#getting if nodemtach is differential diff
node.diff=NA

for ( i in 1:length(var.name.characs.temp)){
ids=grep(var.name.characs[i],attr(terms.formula(base.formula),"variables"))
node=grep("diff",attr(terms.formula(base.formula),"variables"))
ids=ids[ids %in% node]
node.diff[i]=ifelse(length(ids)>0,TRUE,FALSE)
}

##
var.nodecov.factor=var.nodecov+var.nodefactor
form.list.temp=list()

for (i in 1:length(var.name.characs.temp)){

if (is.na(var.type[[i]])) {form.list.temp[[i]]=paste(var.name.characs.temp[i],1,sep="")

}

#if (!is.na(var.type[[i]]) & !(i%in%miss.ids)  ){
if (!is.na(var.type[[i]])){

if (var.type[i]=="absdiff"){

form.list.temp[[i]]=ifelse(var.nodecov.factor[i],paste(paste(var.type[i],".",var.name.characs.temp[i],sep=""),paste(var.name.characs.temp[i],1,sep=""),sep="+"),
paste(var.type[i],".",var.name.characs.temp[i],sep="")  )
}

if (var.type[i]=="nodematch"){
   if (node.diff[i]==F){
form.list.temp[[i]]=ifelse(var.nodecov.factor[i],paste(paste(var.type[i],".",var.name.characs.temp[i],sep=""),paste(var.name.characs.temp[i],1,sep=""),sep="+"),
paste(var.type[i],".",var.name.characs.temp[i],sep=""))
                 }

   if (node.diff[i]==T){
form.list.temp[[i]]=ifelse(var.nodecov.factor[i],paste(paste(var.type[i],".",var.name.characs.temp[i],".",sep=""),paste(var.name.characs.temp[i],1,sep=""),sep="+"),
paste(var.type[i],".",var.name.characs.temp[i],".",sep=""))
                 }

         }#nodematch

if (var.type[i]=="nodemix"){
form.list.temp[[i]]=ifelse(var.nodecov.factor[i],paste(paste("mix.",var.name.characs.temp[i],".",sep=""),paste(var.name.characs.temp[i],1,sep=""),sep="+"),
paste("mix.",var.name.characs.temp[i],".",sep=""))
        }
}

#if (!is.na(var.type[[i]]) & (i%in%miss.ids)  ){
#form.list.temp[[i]]=ifelse(var.nodecov.factor[i],paste(var.name.characs.temp[i],1,sep=""),"none")
#}

if (i ==1){ind.vars=form.list.temp[[i]]}
if (i>1 &form.list.temp[[i]]!="none"){ ind.vars=paste(ind.vars,form.list.temp[[i]],sep="+")}
 
} #i

glm.formula<-formula(paste("y","~",ind.vars))
##

############################################################################


if (check.case.control.data) {
#check to make sure that variables in formula on data set. if not, add it. 

case.control.data=case.list[[1]]


####
for (i in 1:length(var.name.characs.temp)){


#check to see if variable of interest is already on dataset:
#var.indata=paste(var.name.characs[i],"1",sep="") %in% colnames(case.control.data)

#if (!var.indata) {

#paste(var.name.characs[i],"1",sep="")

#person1.charac=hom.data[person1.id,var.name.characs[p]]
#person2.charac=hom.data[person2.id,var.name.characs.person2[p]]


#}



if (!(is.na(var.type[i]))){

if (var.type[i]=="absdiff") {
case.control.data[,paste("absdiff",".",var.name.characs.temp[i],sep="")]=abs(case.control.data[,paste(var.name.characs[i],"1",sep="")]-case.control.data[,paste(var.name.characs.temp[i],"2",sep="")])
     }

if (var.type[i]=="nodematch") {

 if (node.diff[i]==F){
case.control.data[,paste("nodematch.",var.name.characs.temp[i],sep="")]=as.numeric(as.character(case.control.data[,paste(var.name.characs.temp[i],"1",sep="")])==as.character(case.control.data[,paste(var.name.characs[i],"2",sep="")]))
                 }

 if (node.diff[i]==T){#creating factor appropriate for differential nodematch 

match=as.numeric(as.character(case.control.data[,paste(var.name.characs.temp[i],"1",sep="")])==as.character(case.control.data[,paste(var.name.characs[i],"2",sep="")]))

case.control.data[,paste("nodematch.",var.name.characs.temp[i],".",sep="")]=as.character(case.control.data[,paste(var.name.characs.temp[i],"1",sep="")])

case.control.data[match==0&is.na(match)==F&is.na(case.control.data[,paste("nodematch.",var.name.characs.temp[i],".",sep="")])==F,
paste("nodematch.",var.name.characs.temp[i],".",sep="")]="no match"


case.control.data[is.na(case.control.data[,paste(var.name.characs.temp[i],"1",sep="")])|is.na(case.control.data[,paste(var.name.characs.temp[i],"2",sep="")]),paste("nodematch.",var.name.characs[i],".",sep="")]=NA

case.control.data[,paste("nodematch.",var.name.characs.temp[i],".",sep="")]=as.factor(case.control.data[,paste("nodematch.",var.name.characs.temp[i],".",sep="")])
case.control.data[,paste("nodematch.",var.name.characs.temp[i],".",sep="")]=relevel(case.control.data[,paste("nodematch.",var.name.characs.temp[i],".",sep="")],"no match")
#relevling to make "no match" the base
                 }
   }#nodematch


if (var.type[i]=="nodemix"){
temp.dat1=(data.frame(case.control.data[,paste(var.name.characs.temp[i],"1",sep="")],case.control.data[,paste(var.name.characs.temp[i],"2",sep="")]))
#id=temp.dat1[,1]==temp.dat1[,2]
#id.not=!id

n.mix=unique(c(names(table(temp.dat1[,1])),names(table(temp.dat1[,2]))))
r.mix=rank(n.mix)

temp.dat=temp.dat1
temp.dat[,1]=as.numeric(as.character(factor(temp.dat1[,1],levels=n.mix,labels=r.mix)))
temp.dat[,2]=as.numeric(as.character(factor(temp.dat1[,2],levels=n.mix,labels=r.mix)))


id2=temp.dat[,1]>temp.dat[,2]
#id1=temp.dat[,1]<temp.dat[,2]
rm(temp.dat);gc()

id2[is.na(id2)]=FALSE

t2=temp.dat1[id2,2]
t1=temp.dat1[id2,1]


keep2=as.character(temp.dat1[,2])
keep2[id2]=as.character(t1)

keep1=as.character(temp.dat1[,1])
keep1[id2]=as.character(t2)

#temp.dat1[id2,2]=t1
#temp.dat1[id2,1]=t2

var=paste(keep1,keep2,sep=".")



case.control.data[,paste("mix.",var.name.characs.temp[i],".",sep="")]= var
rm(temp.dat1);gc();gc()
#taking out NA's

case.control.data[is.na(case.control.data[,paste(var.name.characs.temp[i],"1",sep="")])|is.na(case.control.data[,paste(var.name.characs[i],"2",sep="")]),paste("mix.",var.name.characs.temp[i],".",sep="")]=NA
case.control.data[,paste("mix.",var.name.characs.temp[i],".",sep="")]=as.factor(case.control.data[,paste("mix.",var.name.characs.temp[i],".",sep="")])

if (!is.null(base.name[[i]])){
case.control.data[,paste("mix.",var.name.characs.temp[i],".",sep="")]=relevel(case.control.data.temp[,paste("mix.",var.name.characs[i],".",sep="")],base.name[[i]])}
    }#mix
  }
}#i


#dealing with node factors
for (i in 1:length(var.nodefactor)){

if (var.nodefactor[i]==FALSE){}

if (var.nodefactor[i]==TRUE){

case.control.data[,paste(var.name.characs.temp[i],1,sep="")]=factor(case.control.data[,paste(var.name.characs.temp[i],1,sep="")])

if (!is.null(base.name.factor[[i]])){
case.control.data[,paste(var.name.characs.temp[i],1,sep="")]=relevel(case.control.data[,paste(var.name.characs.temp[i],1,sep="")],ref=base.name.factor[[i]])
  }
}

} #i


#if (length(miss.ids)>0){

#var.miss=paste(var.name.characs[miss.ids],2,sep="")
#for (i in 1:length(var.miss)){
#col.id=-which(colnames(case.control.data)==var.miss[i])
# if (length(col.id)>0){
# case.control.data=case.control.data[,col.id]}
 # }

#var.miss=paste("absdiff",".",var.name.characs.temp[miss.ids],sep="")
#for (i in 1:length(var.miss)){
#col.id=-which(colnames(case.control.data)==var.miss[i])
# if (length(col.id)>0){
 #case.control.data=case.control.data[,col.id]}
  #}

#var.miss=paste("nodematch",".",var.name.characs.temp[miss.ids],sep="")
#for (i in 1:length(var.miss)){
#col.id=-which(colnames(case.control.data)==var.miss[i])
# if (length(col.id)>0){
 #case.control.data=case.control.data[,col.id]}
  #}

#var.miss=paste("mix",".",var.name.characs.temp[miss.ids],".",sep="")
#for (i in 1:length(var.miss)){
#col.id=-which(colnames(case.control.data)==var.miss[i])
 #if (length(col.id)>0){
 #case.control.data=case.control.data[,col.id]}
  #}

#}



###
} #(check.case.control.data)

##################

list(glm.formula,case.control.data)


}


