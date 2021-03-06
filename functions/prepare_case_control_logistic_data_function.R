# function to take relatively standard ego network dataset 
# and get a dataset used in case control 
#logistic regression
#to get subsets of data, i.e. only kin ties
#need to alter input data so that var.name.degree=0 if no "kin" ties

#inputs:

#data=ego network data, of form chara charac1, charac2, etc.; also have ties alter info, weights and possibly other info
#var.name.degree=name of variable in data that has the degree of the respondents
#var.name.characs= vector of names of characterisictics for homophily model-case control logistic
#case.control.type=type of null hypothesis, one of:
#"case.case.matching", "weighted.random.matching","simulated.data","one_one_matching","one_one_pair_matching"
#resp.weights=vector of weights to create representative population
#if resp.control.weights is not NULL then this vector is simply passed to be used in regressions. 
#max.control.data.N=if not null, number of dyads in control part of regression 
#var.name.characs.alter=list of variable names for alter characteristics, same order as var.name.charac

#alter.tie.data=only needed if add.trans=T, a dataset indicating if alter 1 is tied to 2, 1 to 3...
#add.trans=should put transitive measure onto dataset, no for simulated data
#impute.trans.control=should put some kind of simulated transitive onto the control part of data
#control.type=either "none" or "differential.degree", if "differential.degree" then control part
#is constructed by selecting people proportionally to degree

#if control.type=differential.degree and doing one to one match the respondent is in 
#placed in as many alter-alter pairs as degree

#sim.net=simulated network from which the control part will be drawn-only neccessary if case.control.type="simulated.data"
#vector.trans.control=only needed if impute.trans.control=T; a vector of shared partner, or transtitive, values
#from which to pull values and put onto the dyads of the case part of the dataset 
#output.weights=T or F, should the weights of the respondent be outputed?

#var.name.characs.person2=name of variable of "reciever" of tie in control variable if not the same variable name as respondents
#tie.data=optional dataset which specifies if there is a tie between ego and alter for each possible alter
#remove.isolates.control=should isolates be removed from control part of dataset?
#output.ids=should the ids of the respondents be outputed in dataset?
#id.names=name of id variable if output.ids=T
#one_one_match_withreplace=if doing one_one_matching should you allows respondents to be selected more than once
#to be paired with other respondents (T), or should selectoin be done without replacement (F)
#resp.control.weights=a vector of weights if not NULL, that determines the probability of selection for the control part of the dataset
#default is to use resp.weights

#N.type=type of N used if no max.control.dyad.N, either "dyads" or "respondents" "respondents.diff"
#dyads just gives you all dyads, respondents gives you the number of respondents
#and respondents.diff gives you the number of sum(respondents*degree)

#case.type="none" or "split_sample", split sample if you want only half of cases to be present in case control dataset
# if case.control.type="one_one_pair_matching" then will get you nested cases
#two observatoins per respondent, for half the sample
#case.type is only used if case.control.type="one_one_pair_matching"




rep.func=function(x,N){rep(x,N-x)}
rep.func2=function(x,N){x:N}

case.control.data.function=function(data=NULL,var.name.degree=NULL,var.name.characs=NULL,
case.control.type=NULL,resp.weights=NULL,max.control.data.N=NULL
,var.name.characs.alter=NULL,alter.tie.data=NULL,add.trans=F,impute.trans.control=F,control.type=NULL,sim.net=NULL,vector.trans.control=NULL
,output.weights=F,var.name.characs.person2=NULL,tie.data=NULL
,remove.isolates.control=TRUE,output.ids=F,id.names=NULL,one_one_match_withreplace=T
,resp.control.weights=NULL
,N.type="dyads",case.type="none"){

#setting stuff up if inputs are null
if (is.null(var.name.characs.person2)){var.name.characs.person2=var.name.characs}
if (is.null(control.type)){control.type="none"}


#defining tie data
num.possible.alters=length(var.name.characs.alter[[1]])

if (is.null(tie.data)){
tie.data=matrix(0,nrow=nrow(data),ncol=num.possible.alters)

for (i in 1:num.possible.alters){
tie.data[data[,var.name.degree]>=i,1:i]=1
  }
} 
 

if (is.null(var.name.degree)){

data[,"num.alter"]=rowSums(tie.data,na.rm=T)
var.name.degree="num.alter"
}

#first matching all respondents to each other 
#assume everyone in dataset should be matches with eveyrone else
#after taking out the isolates-if remove.isolate.control=T


hom.data=switch(paste(remove.isolates.control),"TRUE"=data[data[,var.name.degree]>0 &!is.na(data[,var.name.degree]),],
"FALSE"=data[!is.na(data[,var.name.degree]),])


#tie.data=switch(paste(remove.isolates.control),"TRUE"=tie.data[data[,var.name.degree]>0 &!is.na(data[,var.name.degree]),],
#"FALSE"=tie.data[!is.na(data[,var.name.degree]),])

#taking out people in tie.data with no ties


tie.data=tie.data[data[,var.name.degree]>0 &!is.na(data[,var.name.degree]),]
#"FALSE"=tie.data[!is.na(data[,var.name.degree]),])

N=nrow(hom.data)
#N=nrow(tie.data)


tie.data.orig=tie.data

resp.weights=switch(paste(is.null(resp.weights)),"TRUE"=rep(1,nrow(data)),"FALSE"=resp.weights)
resp.control.weights=switch(paste(is.null(resp.control.weights)),"TRUE"=resp.weights,"FALSE"=resp.control.weights)

#resp.weights=switch(paste(remove.isolates.control),"TRUE"=resp.weights[data[,var.name.degree]>0 &!is.na(data[,var.name.degree])]
#,"FALSE"=resp.weights[!is.na(data[,var.name.degree])])

resp.weights=resp.weights[data[,var.name.degree]>0 &!is.na(data[,var.name.degree])]
#,"FALSE"=resp.weights[!is.na(data[,var.name.degree])])

resp.control.weights=switch(paste(remove.isolates.control),"TRUE"=resp.control.weights[data[,var.name.degree]>0 &!is.na(data[,var.name.degree])]
,"FALSE"=resp.control.weights[!is.na(data[,var.name.degree])])



#if want to simply match all respondents with all other respondents in case part of dataset
if (case.control.type=="case.case.matching"){

person1.id=unlist(lapply(1:(N-1),rep.func,N=N))
person2.id=unlist(lapply(2:(N),rep.func2,N=N))

for (p in 1:length(var.name.characs)){
person1.charac=hom.data[person1.id,var.name.characs[p]]
person2.charac=hom.data[person2.id,var.name.characs.person2[p]]
       
      if (p==1){
      control.data=data.frame(person1.charac,person2.charac)
rm(person1.charac,person2.charac)  
gc();gc() 
            }
if (p>1){control.data=data.frame(control.data,person1.charac,person2.charac)
 rm(person1.charac,person2.charac)  
gc();gc()
         }                
     }#p

colnames(control.data)=paste(rep(var.name.characs,each=2),1:2,sep="")
control.data$y=0
#putting degree of sender onto dataset
deg=hom.data[,var.name.degree][person1.id]
control.data[,"degree"]=deg

if (!is.null(max.control.data.N)){
 if ( max.control.data.N<nrow(control.data)) {
id.temp=sample(1:nrow(control.data),size=max.control.data.N,replace=F)
  control.data=control.data[id.temp,]
person1.id=person1.id[id.temp]
person2.id=person2.id[id.temp]

  }

if (max.control.data.N>nrow(control.data)){
#building up dataset to be size max.control.data.N 

id.temp=sample(1:nrow(control.data),size=max.control.data.N-nrow(control.data),replace=T)
  control.data=rbind(control.data,control.data[id.temp,])
person1.id=c(person1.id,person1.id[id.temp])
person2.id=c(person2.id,person2.id[id.temp])

   }
  }

 } #case.case



#if want to match respondents in data by weight, to get representative control part
if (case.control.type=="weighted.random.matching"){
#resp.weights=resp.weights[data[,var.name.degree]>0 &!is.na(data[,var.name.degree])]
if (is.null(max.control.data.N)){   N.dyads=switch(N.type,"dyads"=N*(N-1)/2,"respondents"=N,"respondents.diff"=sum(tie.data,na.rm=T) )

}

if (!is.null(max.control.data.N)){  N.dyads=max.control.data.N }

person1.id=sample(1:length(resp.control.weights),size=N.dyads,prob=resp.control.weights,replace=T)
person2.id=sample(1:length(resp.control.weights),size=N.dyads,prob=resp.control.weights,replace=T)


if (control.type=="differential.degree"){
deg=hom.data[,var.name.degree][person1.id]
person1.id=sample(person1.id,prob=deg,replace=T)
deg=hom.data[,var.name.degree][person2.id]
person2.id=sample(person2.id,prob=deg,replace=T)
   }


#making sure randomly selected pairs not same person
#ids.rm=person1.id[which((person1.id/person2.id)==1)]

ids.rm=which((person1.id/person2.id)==1)

if (length(ids.rm)>0){
#person1.id[which((person1.id/person2.id)==1)]=sample((1:length(resp.control.weights))[-ids.rm]
#,size=length(ids.rm),prob=resp.control.weights[-ids.rm],replace=T)

person1.id=person1.id[-ids.rm]
person2.id=person2.id[-ids.rm]
}


for (p in 1:length(var.name.characs)){
person1.charac=hom.data[person1.id,var.name.characs[p]]
person2.charac=hom.data[person2.id,var.name.characs.person2[p]]
       
      if (p==1){
      control.data=data.frame(person1.charac,person2.charac)
rm(person1.charac,person2.charac)  
            }
if (p>1){control.data=data.frame(control.data,person1.charac,person2.charac)
 rm(person1.charac,person2.charac)  
gc();gc()
         }                
     }#p

colnames(control.data)=paste(rep(var.name.characs,each=2),1:2,sep="")
control.data$y=0
#putting degree of sender onto dataset
deg=hom.data[,var.name.degree][person1.id]
control.data[,"degree"]=deg

   }

#control part-if want draws from simulated data
if (case.control.type=="simulated.data"){
edges=as.edgelist(sim.net)

for (p in 1:length(var.name.characs)){

var1=get.vertex.attribute(sim.net,var.name.characs[p])
var2=get.vertex.attribute(sim.net,var.name.characs.person2[p])


person1.charac=var1[edges[,1]]
person2.charac=var2[edges[,2]]

     if (p==1){
      control.data=data.frame(person1.charac,person2.charac)
rm(person1.charac,person2.charac)  
gc();gc() 
            }
if (p>1){control.data=data.frame(control.data,person1.charac,person2.charac)
 rm(person1.charac,person2.charac)  
gc();gc()
    }
}#p

colnames(control.data)=paste(rep(var.name.characs,each=2),1:2,sep="")

control.data$y=0
control.data$degree=0
   }


#if you want to match each case with only one other case for the control part of the dataset
if (case.control.type=="one_one_matching"){


if (control.type!="differential.degree"){
person1.id=1:nrow(tie.data)
person2.id=sample(person1.id,prob=resp.control.weights,replace=T)

if (one_one_match_withreplace==F){

person2.id=sample(person1.id,replace=F)
ids.rm=person1.id[which((person1.id/person2.id)==1)]
  while(length(ids.rm)==1){ #a quick fix so that we 0 or more than 1 case with id1=id2
person2.id=sample(person1.id,replace=F)
ids.rm=person1.id[which((person1.id/person2.id)==1)]
   }

}

}#no differential degree


if (control.type=="differential.degree"){

num.ties=rowSums(tie.data)
person1.id=rep(1:nrow(tie.data),times=num.ties)
person2.id=sample(1:nrow(tie.data),size=length(person1.id),prob=resp.control.weights,replace=T)

   }

#making sure we haven't selected the same person for ego and alter
ids.rm=person1.id[which((person1.id/person2.id)==1)]



if (length(ids.rm)>0){

if (one_one_match_withreplace==T){
person2.id[which((person1.id/person2.id)==1)]=sample((1:length(resp.control.weights))[-ids.rm]
,size=length(ids.rm),prob=resp.control.weights[-ids.rm],replace=T)
  }

if (one_one_match_withreplace==F){
temp.ids=person2.id[which((person1.id/person2.id)==1)]

#temp.ids=c(sample(person1.id[!person1.id%in% temp.ids] ,1),temp.ids)

person2.id[which((person1.id/person2.id)==1)]=c(temp.ids[-1],temp.ids[1])
}


}


for (p in 1:length(var.name.characs)){
person1.charac=hom.data[person1.id,var.name.characs[p]]
person2.charac=hom.data[person2.id,var.name.characs.person2[p]]
       
      if (p==1){
      control.data=data.frame(person1.charac,person2.charac)
rm(person1.charac,person2.charac)  
gc();gc() 
            }
if (p>1){control.data=data.frame(control.data,person1.charac,person2.charac)
 rm(person1.charac,person2.charac)  
gc();gc()
         }                
     }#p

colnames(control.data)=paste(rep(var.name.characs,each=2),1:2,sep="")
control.data$y=0
#putting degree of sender onto dataset
deg=hom.data[,var.name.degree][person1.id]
control.data[,"degree"]=deg

}#one to one matching


#if only want each respondent in one pair so in control part once-either
#as sender or reciever but not both
if (case.control.type=="one_one_pair_matching"){


#if (control.type!="differential.degree"){

person1.id=sample(1:nrow(tie.data),size=nrow(tie.data)/2,replace=F)
ids.temp=(1:nrow(tie.data))[!(1:nrow(tie.data)%in%person1.id)]
person2.id=sample(ids.temp,size=nrow(tie.data)/2,replace=F)

ids.rm=person1.id[which((person1.id/person2.id)==1)]
  while(length(ids.rm)==1){ #a quick fix so that we 0 or more than 1 case with id1=id2
person2.id=sample(person1.id,replace=F)
ids.rm=person1.id[which((person1.id/person2.id)==1)]
  }

#}#no differential degree


if (length(ids.rm)>0){

#if (one_one_match_withreplace==F){
temp.ids=person2.id[which((person1.id/person2.id)==1)]
person2.id[which((person1.id/person2.id)==1)]=c(temp.ids[-1],temp.ids[1])
#      }
}


for (p in 1:length(var.name.characs)){
person1.charac=hom.data[person1.id,var.name.characs[p]]
person2.charac=hom.data[person2.id,var.name.characs.person2[p]]
       
      if (p==1){
      control.data=data.frame(person1.charac,person2.charac)
rm(person1.charac,person2.charac)  
gc();gc() 
            }
if (p>1){control.data=data.frame(control.data,person1.charac,person2.charac)
 rm(person1.charac,person2.charac)  
gc();gc()
         }                
     }#p

colnames(control.data)=paste(rep(var.name.characs,each=2),1:2,sep="")
control.data$y=0
#putting degree of sender onto dataset
deg=hom.data[,var.name.degree][person1.id]
control.data[,"degree"]=deg

}#one to one pair matching



#######################################################################

#declaring hom.data as used in control part 

hom.data.control=hom.data

#now subsetting hom.data in case using different populatoin to do control

hom.data=data[data[,var.name.degree]>0 &!is.na(data[,var.name.degree]),]
#"FALSE"=data[!is.na(data[,var.name.degree]),])


#case part of data

#making 1's and 0's of ties data input as 1,2,3s etc.

tie.data=t(apply(tie.data,1,"*",1:ncol(tie.data)))


for (p in 1:length(var.name.characs)){
person2=list()
person1=list()

for (x in 1:nrow(hom.data)){
person2[[x]]=unlist(hom.data[x,var.name.characs.alter[[p]]][,tie.data[x,]])
person1[[x]]=rep(hom.data[x,var.name.characs[p]],length(person2[[x]]))
    }


if (p==1){
case.data=data.frame(unlist(person1),unlist(person2))
   }


if (p>1){
  case.data=data.frame(case.data,unlist(person1),unlist(person2))
  
  }


}#p
colnames(case.data)=paste(rep(var.name.characs,each=2),1:2,sep="")
case.data$y=1

#getting degree onto case data
for (x in 1:nrow(hom.data)){
person1[[x]]=rep(hom.data[x,var.name.degree],length(person2[[x]]))
    }

deg=unlist(person1)
case.data[,"degree"]=deg


if (add.trans){#if want transtivity onto dataset
alter.tie.data=alter.tie.data[data[,var.name.degree]>0 &!is.na(data[,var.name.degree]),]

trans.dat=tie.data.orig[,]*alter.tie.data[,]
ties.ego1=trans.dat[,1]+trans.dat[,2]
ties.ego2=trans.dat[,1]+trans.dat[,3]
ties.ego3=trans.dat[,2]+trans.dat[,3]

ties.ego=cbind(ties.ego1,ties.ego2,ties.ego3)
ties.ego.list=list()
for (x in 1:nrow(hom.data)){
ties.ego.list[[x]]=ties.ego[x,][tie.data[x,]]
    }#x
trans=unlist(ties.ego.list)
trans[is.na(trans)]=0
case.data[,"trans"]=trans
control.data[,"trans"]=0

if (impute.trans.control==T){#should imput some transtitivity level under randomness?

  control.data[,"trans"]=sample(vector.trans.control,length(control.data[,"trans"]),replace=T)
    }
}


if (output.weights==T){#should ouput weights?
control.data[,"weights.resp"]=1

p1=list()
for (x in 1:nrow(hom.data)){
p1[[x]]=rep(resp.weights[x],length(person2[[x]]))
    }

w1=unlist(p1)
case.data[,"weights.resp"]=w1

}


if (output.ids==T){#should ouput ids?

person1.ids=hom.data.control[person1.id,id.names]
person2.ids=hom.data.control[person2.id,id.names]

control.data[,"id1"]=person1.ids
control.data[,"id2"]=person2.ids


p1=list()
for (x in 1:nrow(hom.data)){
p1[[x]]=rep(hom.data[x,id.names],length(person2[[x]]))
    }

w1=unlist(p1)
case.data[,"id1"]=w1
case.data[,"id2"]=w1

}


if(case.type=="split_sample" &control.type=="one_one_pair_matching"){

case.data=case.data[person1.id,]
}

case.control.data=rbind(case.data,control.data)
rm(case.data,control.data)
gc();gc()
return(case.control.data)

}#end function

