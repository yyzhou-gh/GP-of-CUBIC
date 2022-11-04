SNP<-read.table("/home/yyzhou/CUBIC/634SNP/dataset/tmp.list.hp",head=TRUE) #634 SNPs of trainingset and testingset
load("/home/yyzhou/CUBIC/634SNP/dataset/EW.394.Rdata")
load("/home/yyzhou/CUBIC/634SNP/dataset/test.EW.Rdata")
load("/home/yyzhou/CUBIC/634SNP/dataset/cv.group.Rdata")
id.394<-rownames(as.matrix(EW.394))
snp.394<-t(SNP[,id.394])
id.testing<-rownames(test.EW)
test.snp<-t(SNP[,id.testing])

#creat dataframe
CUBIC.dataframe<-data.frame(EW.394,snp.394)

#create list  
    trainingset<-list()
    validset<-list()
    
    CUBIC.G.training<-list()
    CUBIC.EW.training<-list()

    CUBIC.G.valid<-list()  
    CUBIC.EW.valid<-list()

for(i in 1:10){
  #create trainingset     	
      trainingset[[i]]<-CUBIC.dataframe[setdiff(rownames(snp.394),cv.group[[i]]),]
      CUBIC.G.training[[i]]<-as.matrix(trainingset[[i]][,2:635])
      CUBIC.EW.training[[i]]<-as.matrix(trainingset[[i]][,1])	
      colnames(CUBIC.EW.training[[i]])<-'EW'
  #create validset  
      validset[[i]]<-CUBIC.dataframe[cv.group[[i]],]  	
      CUBIC.G.valid[[i]]<-as.matrix(validset[[i]][,2:635])
      CUBIC.EW.valid[[i]]<-as.matrix(validset[[i]][,1]) 
      colnames(CUBIC.EW.valid[[i]])<-'EW'		
      ifelse(i==1,CUBIC.EW.valid.10merge<-CUBIC.EW.valid[[i]],CUBIC.EW.valid.10merge<-rbind(CUBIC.EW.valid.10merge,CUBIC.EW.valid[[i]]))
}

library("Matrix")
library("glmnet")

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

#G to EW
Pred.G2EW<-list()
for(i in 1:10){	
    Pred.G2EW[[i]]<-training.LASSO(CUBIC.G.training[[i]],CUBIC.G.valid[[i]],CUBIC.EW.training[[i]],CUBIC.EW.valid[[i]]) 
    ifelse(i==1,Pred.G2EW.10merge<-Pred.G2EW[[i]][[2]],Pred.G2EW.10merge<-rbind(Pred.G2EW.10merge,Pred.G2EW[[i]][[2]]))		  
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
Predictability.G2EW.394<-evaluation(Pred.G2EW.10merge,CUBIC.EW.valid.10merge)
save(Predictability.G2EW.394,file="/home/yyzhou/CUBIC/634SNP/result/Predictability.G2EW.394.Rdata")

#construct model in trainingset and test in independent testingset
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
      save(cvfit,lambda,file=paste0("/home/yyzhou/CUBIC/634SNP/model/",model.name,"/",l,".",colnames(Y.training)[l],".Rdata"))
      }
    print(l)
    return(list(PredY.training,PredY.testing))
}
Pred.G2EW<-testing.LASSO(snp.394,test.snp,EW.394,'G2EW')
Predictability.G2EW.testing<-evaluation(Pred.G2EW[[2]],test.EW)
save(Predictability.G2EW.testing,file="/home/yyzhou/CUBIC/634SNP/result/Predictability.testing/Predictability.G2EW.testing.Rdata")


