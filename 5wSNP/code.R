#data preparation
setwd("/home/yyzhou/CUBIC/5wSNP")	 
#Load all omics data   
    load("/home/yyzhou/CUBIC/5wSNP/dataset/trainingset/snp.394.Rdata")  #50000 SNP
    load("/home/yyzhou/CUBIC/5wSNP/dataset/trainingset/Transcriptomics.394.Rdata")  #26743 Transcripts 
    load("/home/yyzhou/CUBIC/5wSNP/dataset/trainingset/Sub_Trait.394.Rdata") #19 Sub_Traits
    load("/home/yyzhou/CUBIC/5wSNP/dataset/trainingset/EW.394.Rdata")
    load("/home/yyzhou/CUBIC/5wSNP/dataset/trainingset/cv.group.Rdata")

    #creat dataframe
    CUBIC.dataframe<-data.frame(Sub_Trait.394,EW.394,snp.394,Transcriptomics.394)

#create list  
    trainingset<-list()
    validset<-list()
    
    CUBIC.G.training<-list()
    CUBIC.E.training<-list()
    CUBIC.Sub_Trait.training<-list()
    CUBIC.EW.training<-list()

    CUBIC.G.valid<-list()
    CUBIC.E.valid<-list()
    CUBIC.Sub_Trait.valid<-list()   
    CUBIC.EW.valid<-list()

for(i in 1:10){
  #create trainingset     	
      trainingset[[i]]<-CUBIC.dataframe[setdiff(rownames(snp.394),cv.group[[i]]),]
      CUBIC.G.training[[i]]<-as.matrix(trainingset[[i]][,21:50020])
      CUBIC.E.training[[i]]<-as.matrix(trainingset[[i]][,50021:76763])
      CUBIC.Sub_Trait.training[[i]]<-as.matrix(trainingset[[i]][,1:19])	
      CUBIC.EW.training[[i]]<-as.matrix(trainingset[[i]][,20])	
      colnames(CUBIC.EW.training[[i]])<-'EW'
  #create validset  
      validset[[i]]<-CUBIC.dataframe[cv.group[[i]],]  	
      CUBIC.G.valid[[i]]<-as.matrix(validset[[i]][,21:50020])
      CUBIC.E.valid[[i]]<-as.matrix(validset[[i]][,50021:76763])
      CUBIC.Sub_Trait.valid[[i]]<-as.matrix(validset[[i]][,1:19])
      CUBIC.EW.valid[[i]]<-as.matrix(validset[[i]][,20]) 
      colnames(CUBIC.EW.valid[[i]])<-'EW'		
      ifelse(i==1,CUBIC.E.valid.10merge<-CUBIC.E.valid[[i]],CUBIC.E.valid.10merge<-rbind(CUBIC.E.valid.10merge,CUBIC.E.valid[[i]]))		  	  
      ifelse(i==1,CUBIC.EW.valid.10merge<-CUBIC.EW.valid[[i]],CUBIC.EW.valid.10merge<-rbind(CUBIC.EW.valid.10merge,CUBIC.EW.valid[[i]]))
}

training.LASSO<-function(X.training,X.valid,Y.training,Y.valid){ 
    p<-ncol(Y.training)
    x<-X.training
    PredY.training<-matrix(,nrow(Y.training),ncol(Y.training))
    PredY.valid<-matrix(,nrow(Y.valid),ncol(Y.valid))
    rownames(PredY.training)<-rownames(X.training)
    colnames(PredY.training)<-colnames(Y.training)
    rownames(PredY.valid)<-rownames(X.valid)
    colnames(PredY.valid)<-colnames(Y.valid)	
    for(l in 1:p){
    	cat("\r",l)
    	y<-Y.training[,l]
    	cvfit<-cv.glmnet(x,y)
    	lambda<-cvfit$lambda.min
    	pred<-predict(cvfit,newx=x,s=lambda)
    	pred.2<-predict(cvfit,newx=X.valid,s=lambda)
    	PredY.training[,l]<-pred
    	PredY.valid[,l]<-pred.2
    	}
    print(l)
    return(list(PredY.training,PredY.valid))
  }

#calculation function of R2/PCC
evaluation<-function(Pred,True){ 
  eva.matrix<-matrix(,ncol(True),2)
  rownames(eva.matrix)<-colnames(True)
  colnames(eva.matrix)<-c("R2","PCC")
  for(l in 1:ncol(Pred)){ 
    y.Pred<-Pred[,l]
    y.True<-True[,l]
    SSE<-sum((y.Pred-y.True)^2)
    SST<-sum((y.True-mean(y.True))^2)
    eva.matrix[l,1]<-1-(SSE/SST)    
    eva.matrix[l,2]<-cor(y.Pred,y.True) 
  }
  return(eva.matrix)
}

##construct G2GPG2Sub_Trait2EW model
#construct the first layer
args=commandArgs(T)
i=as.numeric(args[1])
library("Matrix")
library("glmnet")
print('Begin')
Pred.E.training.1L<-list()
Pred.E.valid.1L<-list()
for(i in 1:10){
    Pred.E.1L<-training.LASSO(CUBIC.G.training[[i]],CUBIC.G.valid[[i]],CUBIC.E.training[[i]],CUBIC.E.valid[[i]])
    Pred.E.training.1L[[i]]<-Pred.E.1L[[1]]
    Pred.E.valid.1L[[i]]<-Pred.E.1L[[2]]
    ifelse(i==1,Pred.E.valid.1L.10merge<-Pred.E.valid.1L[[i]],Pred.E.valid.1L.10merge<-rbind(Pred.E.valid.1L.10merge,Pred.E.valid.1L[[i]]))		
}
Predictability.E.1L<-evaluation(Pred.E.valid.1L.10merge,CUBIC.E.valid.10merge)
#construct the second and third layer
Predictability.EW.3L<-matrix(,18,2) 
Pred.Sub_Trait.2L<-list()
Pred.EW.3L<-list()
for(j in 1:18){
  Pred.Sub_Trait.2L[[j]]<-list();length(Pred.Sub_Trait.2L[[j]])<-10
  Pred.EW.3L[[j]]<-list();length(Pred.EW.3L[[j]])<-10
  alpha<-(j-1)*0.05   
  ids.Pred.GPGs.1L<-rownames(Predictability.E.1L)[which(Predictability.E.1L[,1]>alpha)]
  Pred.GPGs.1L.training<-list()
  Pred.GPGs.1L.valid<-list()

  for(i in 1:10){	
      Pred.GPGs.1L.training[[i]]<-Pred.E.training.1L[[i]][,ids.Pred.GPGs.1L]
      Pred.GPGs.1L.valid[[i]]<-Pred.E.valid.1L[[i]][,ids.Pred.GPGs.1L]
      Pred.Sub_Trait.2L[[j]][[i]]<-training.LASSO(Pred.GPGs.1L.training[[i]],Pred.GPGs.1L.valid[[i]],
                                  CUBIC.Sub_Trait.training[[i]],CUBIC.Sub_Trait.valid[[i]])  
      Pred.EW.3L[[j]][[i]]<-training.LASSO(Pred.Sub_Trait.2L[[j]][[i]][[1]],Pred.Sub_Trait.2L[[j]][[i]][[2]],
                                      CUBIC.EW.training[[i]],CUBIC.EW.valid[[i]])     
      ifelse(i==1,Pred.EW.3L.valid.10merge<-Pred.EW.3L[[j]][[i]][[2]],Pred.EW.3L.valid.10merge<-rbind(Pred.EW.3L.valid.10merge,Pred.EW.3L[[j]][[i]][[2]]))
  }    
  Predictability.EW.3L[j,]<-evaluation(Pred.EW.3L.valid.10merge,CUBIC.EW.valid.10merge)
}
  rownames(Predictability.EW.3L)<-paste0("alpha",seq(0,0.85, by=0.05))
  colnames(Predictability.EW.3L)<-c("R2","PCC")
  best.alpha2<-0.05*(which(Predictability.EW.3L[,1]==max(Predictability.EW.3L[,1]))-1)  # best alpha2
  Predictability.EW.3L.394<-Predictability.EW.3L[best.alpha2,]
  Transcriptomics.GPG<-Transcriptomics.394[,rownames(Predictability.E.1L)[which(Predictability.E.1L[,1]>best.alpha2)]]
  save(Predictability.EW.3L.394,file="/home/yyzhou/CUBIC/5wSNP/result/Predictability.training/Predictability.EW.3L.394.Rdata")

##construct other models
Pred.G2EW<-list()
Pred.all_E2EW<-list()
Pred.Sub_Trait2EW<-list()
Pred.GPGs.training<-list()
Pred.GPGs.valid<-list()
Pred.GPG&Sub_Trait.training<-list()
Pred.GPG&Sub_Trait.valid<-list()
Pred.GPG&Sub_Trait2EW<-list()
#G to E to EW
Predictability.G2GPG2EW.394<-matrix(,18,2)
Pred.G2GPG2EW<-list()
  for(j in 1:18){
  alpha1<-(j-1)*0.05  
  ids.G2GPG2EW<-rownames(Predictability.E.1L)[which(Predictability.E.1L[,1]>alpha1)]
  Pred.GPGs.training<-list()
  Pred.GPGs.valid<-list()
  Pred.G2GPG2EW[[j]]<-list();length(Pred.G2GPG2EW[[j]])<-10
  for(i in 1:10){	
      Pred.GPGs.training[[i]]<-Pred.E.training.1L[[i]][,ids.G2GPG2EW]
      Pred.GPGs.valid[[i]]<-Pred.E.valid.1L[[i]][,ids.G2GPG2EW]
      Pred.G2GPG2EW[[j]][[i]]<-training.LASSO(Pred.GPGs.training[[i]],Pred.GPGs.valid[[i]],
                                  CUBIC.EW.training[[i]],CUBIC.EW.valid[[i]])  
      ifelse(i==1,Pred.G2GPG2EW.valid.10merge<-Pred.G2GPG2EW[[j]][[i]][[2]],Pred.G2GPG2EW.valid.10merge<-rbind(Pred.G2GPG2EW.valid.10merge,Pred.G2GPG2EW[[j]][[i]][[2]]))
  }    
  Predictability.G2GPG2EW.394[j,]<-evaluation(Pred.G2GPG2EW.valid.10merge,CUBIC.EW.valid.10merge)
}
rownames(Predictability.G2GPG2EW)<-paste0("alpha",seq(0,0.85, by=0.05))
colnames(Predictability.G2GPG2EW)<-c("R2","PCC")
best.alpha1<-0.05*(which(Predictability.G2GPG2EW[,1]==max(Predictability.G2GPG2EW[,1]))-1)  # best alpha1
G2E2EW.GPGs<-Transcriptomics.394[,rownames(Predictability.E.1L)[which(Predictability.E.1L[,1]>best.alpha1)]] #4256???

ids.Pred.GPGs<-rownames(Predictability.E.1L)[which(Predictability.E.1L[,1]>best.alpha2)]
for(i in 1:10){	
  #G to EW
    Pred.G2EW[[i]]<-training.LASSO(CUBIC.G.training[[i]],CUBIC.G.valid[[i]],CUBIC.EW.training[[i]],CUBIC.EW.valid[[i]]) 
    ifelse(i==1,Pred.G2EW.10merge<-Pred.G2EW[[i]][[2]],Pred.G2EW.10merge<-rbind(Pred.G2EW.10merge,Pred.G2EW[[i]][[2]]))		  
  #all Transcriptomics to EW
    Pred.all_E2EW[[i]]<-training.LASSO(CUBIC.E.training[[i]],CUBIC.E.valid[[i]],CUBIC.EW.training[[i]],CUBIC.EW.valid[[i]]) 
    ifelse(i==1,Pred.all_E2EW.10merge<-Pred.all_E2EW[[i]][[2]],Pred.all_E2EW.10merge<-rbind(Pred.all_E2EW.10merge,Pred.all_E2EW[[i]][[2]]))		  
  #Sub_Trait to EW
    Pred.Sub_Trait2EW[[i]]<-training.LASSO(CUBIC.Sub_Trait.training[[i]],CUBIC.Sub_Trait.valid[[i]],CUBIC.EW.training[[i]],CUBIC.EW.valid[[i]]) 
    ifelse(i==1,Pred.Sub_Trait2EW.10merge<-Pred.Sub_Trait2EW[[i]][[2]],Pred.Sub_Trait2EW.10merge<-rbind(Pred.Sub_Trait2EW.10merge,Pred.Sub_Trait2EW[[i]][[2]]))		  	  
  #G to GPG&Sub-Trait to EW
    Pred.GPGs.training[[i]]<-Pred.E.training.1L[[i]][,ids.Pred.GPGs]
    Pred.GPGs.valid[[i]]<-Pred.E.valid.1L[[i]][,ids.Pred.GPGs]
    Pred.GPG&Sub_Trait.training[[i]]<-cbind(Pred.GPGs.training[[i]],Pred.Sub_Trait.2L[[best.alpha2]][[i]][[1]])
    Pred.GPG&Sub_Trait.valid[[i]]<-cbind(Pred.GPGs.valid[[i]],Pred.Sub_Trait.2L[[best.alpha2]][[i]][[2]])
    Pred.GPG&Sub_Trait2EW[[i]]<-training.LASSO(Pred.GPG&Sub_Trait.training[[i]],Pred.GPG&Sub_Trait.valid[[i]],CUBIC.EW.training[[i]],CUBIC.EW.valid[[i]]) 
    ifelse(i==1,Pred.GPG&Sub_Trait2EW.10merge<-Pred.GPG&Sub_Trait2EW[[i]][[2]],Pred.GPG&Sub_Trait2EW.10merge<-rbind(Pred.GPG&Sub_Trait2EW.10merge,Pred.GPG&Sub_Trait2EW[[i]][[2]]))
}

Predictability.G2EW.394<-evaluation(Pred.G2EW.10merge,CUBIC.EW.valid.10merge)
Predictability.all_E2EW.394<-evaluation(Pred.all_E2EW.10merge,CUBIC.EW.valid.10merge)
Predictability.Sub_Trait2EW.394<-evaluation(Pred.Sub_Trait2EW.10merge,CUBIC.EW.valid.10merge)
Predictability.GPG&Sub_Trait2EW.394<-evaluation(Pred.GPG&Sub_Trait2EW.10merge,CUBIC.EW.valid.10merge)
save(Predictability.G2EW.394,file="/home/yyzhou/CUBIC/5wSNP/result/Predictability.training/Predictability.G2EW.394.Rdata")
save(Predictability.all_E2EW.394,file="/home/yyzhou/CUBIC/5wSNP/result/Predictability.training/Predictability.all_E2EW.394.Rdata")
save(Predictability.Sub_Trait2EW.394,file="/home/yyzhou/CUBIC/5wSNP/result/Predictability.training/Predictability.Sub_Trait2EW.394.Rdata")
save(Predictability.GPG&Sub_Trait2EW.394,file="/home/yyzhou/CUBIC/5wSNP/result/Predictability.training/Predictability.GPG&Sub_Trait2EW.394.Rdata")

##construct model in trainingset and test in independent testingset
load('/home/yyzhou/CUBIC20211104/datasets/test.snp.Rdata')
load('/home/yyzhou/CUBIC20211104/datasets/test.Sub_Trait.Rdata') 
load('/home/yyzhou/CUBIC20211104/datasets/test.EW.Rdata') 
testing.LASSO<-function(X.training,X.testing,Y.training,model.name){ 
    p<-ncol(Y.training)
    x<-X.training
    PredY.training<-matrix(,nrow(Y.training),ncol(Y.training))
    PredY.testing<-matrix(,nrow(X.testing),ncol(Y.training))
    rownames(PredY.training)<-rownames(Y.training)
    colnames(PredY.training)<-colnames(Y.training)
    rownames(PredY.testing)<-rownames(X.testing)
    colnames(PredY.testing)<-colnames(Y.training)		
    for(l in 1:p){
    	cat("\r",l)
    	y<-Y.training[,l]
    	cvfit<-cv.glmnet(x,y)
    	lambda<-cvfit$lambda.min
    	pred<-predict(cvfit,newx=x,s=lambda)
    	pred.2<-predict(cvfit,newx=X.testing,s=lambda)
    	PredY.training[,l]<-pred
    	PredY.testing[,l]<-pred.2
      save(cvfit,lambda,file=paste0("/home/yyzhou/CUBIC/5wSNP/model/",model.name,"/",l,".",colnames(Y.training)[l],".Rdata"))
      }
    print(l)
    return(list(PredY.training,PredY.testing))
}

#1.G to GPGs to Sub_Trait to EW
#1.1.The first layer for G to GPG
Pred.GPGs.1L<-testing.LASSO(snp.394,test.snp,Transcriptomics.GPG,'G2GPG.1L')
Pred.GPGs.1L.394<-Pred.GPGs.1L[[1]]
Pred.GPGs.1L.testing<-Pred.GPGs.1L[[2]]
#1.2.The second layer for predicted GPG to Sub_Trait 
Pred.Sub_Trait.2L<-testing.LASSO(Pred.GPGs.1L.394,Pred.GPGs.1L.testing,Sub_Trait.394,'GPG2Phe.2L')
Pred.Sub_Trait.2L.394<-Pred.Sub_Trait.2L[[1]]
Pred.Sub_Trait.2L.testing<-Pred.Sub_Trait.2L[[2]]
#1.3.The third layer for predicted Sub_Trait to EW 
Pred.EW.3L<-testing.LASSO(Pred.Sub_Trait.2L.394,Pred.Sub_Trait.2L.testing,EW.394,'Phe2EW.3L')
Predictability.Phe2EW.3L<-evaluation(Pred.EW.3L[[2]],test.EW)
save(Predictability.Phe2EW.3L,file="/home/yyzhou/CUBIC/5wSNP/result/Predictability.testing/Predictability.testing/Predictability.Phe2EW.3L.Rdata")

#2.G2GPG2EW
Pred.G2GPG<-testing.LASSO(snp.394,test.snp,G2E2EW.GPGs,'G2GPG')
Pred.G2GPG.394<-Pred.G2GPG[[1]]
Pred.G2GPG.testing<-Pred.G2GPG[[2]]
Pred.G2GPG2EW<-testing.LASSO(Pred.G2GPG.394,Pred.G2GPG.testing,EW.394,'G2GPG2EW')
Predictability.G2GPG2EW<-evaluation(Pred.G2GPG2EW[[2]],test.EW)
save(Predictability.G2GPG2EW,file="/home/yyzhou/CUBIC/5wSNP/result/Predictability.testing/Predictability.testing/Predictability.G2GPG2EW.Rdata")

#3.construct model for G to EW
Pred.G2EW<-testing.LASSO(snp.394,test.snp,EW.394,'G2EW')
Predictability.G2EW<-evaluation(Pred.G2EW[[2]],test.EW)
save(Predictability.G2EW,file="/home/yyzhou/CUBIC/5wSNP/result/Predictability.testing/Predictability.testing/Predictability.G2EW.Rdata")

#4.construct model for predicted GPG&Sub_Trait to EW  
Pred.GPG&Sub_Trait.394<-cbind(Pred.GPGs.1L.394,Pred.Sub_Trait.2L.394)
Pred.GPG&Sub_Trait.testing<-cbind(Pred.GPGs.1L.testing,Pred.Sub_Trait.2L.testing)
Pred.GPG&Sub_Trait2EW<-testing.LASSO(Pred.GPG&Sub_Trait.394,Pred.GPG&Sub_Trait.testing,EW.394,'G2GPG&Sub_Trait2EW')
Predictability.G2GPG&Sub_Trait2EW<-evaluation(Pred.GPG&Sub_Trait2EW[[2]],test.EW)
save(Predictability.G2GPG&Sub_Trait2EW,file="/home/yyzhou/CUBIC/5wSNP/result/Predictability.testing/Predictability.testing/Predictability.G2GPG&Sub_Trait2EW.Rdata")

#5.construct model for Sub_Trait to EW
id.test<-rownames(test.Sub_Trait)[complete.cases(test.Sub_Trait)]  #remove NA in the Sub_Trait,leaving 998 lines
test.snp.998<-test.snp[id.test,]
test.Sub_Trait.998<-test.Sub_Trait[id.test,]
test.EW.998<-as.matrix(test.EW[id.test,])
Pred.Sub_Trait2EW<-testing.LASSO(Sub_Trait.394,test.Sub_Trait.998,EW.394,'Phe2EW')
Predictability.Sub_Trait2EW<-evaluation(Pred.Sub_Trait2EW[[2]],test.EW.998)
save(Predictability.Sub_Trait2EW,file="/home/yyzhou/CUBIC/5wSNP/result/Predictability.testing/Predictability.testing/Predictability.Sub_Trait2EW.Rdata")
