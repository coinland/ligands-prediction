require(Rcpi)
require(RbioRXN)
require(rpubchem)
library("pROC")
library("caret")

active.mol<-readMolFromSDF("active.sdf")

active = cbind(extractDrugMACCSComplete(active.mol),
               extractDrugEstateComplete(active.mol),
               extractDrugExtendedComplete(active.mol),
               extractDrugGraphComplete(active.mol))

active<-as.data.frame(active)

active$outcome<-rep(c("active"),123)


inactive.mol<-readMolFromSDF("inactive.sdf")

inactive = cbind(extractDrugMACCSComplete(inactive.mol),
                 extractDrugEstateComplete(inactive.mol),
                 extractDrugExtendedComplete(inactive.mol),
                 extractDrugGraphComplete(inactive.mol)) 

inactive<-as.data.frame(inactive)

inactive$outcome<-rep(c("inactive"),412)

total<-rbind(active,inactive)

total$outcome<-as.factor(total$outcome)

library(Boruta)

boruta.train <- Boruta(outcome~., data = total, doTrace = 2)

print(boruta.train)

plot(boruta.train, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta.train$ImpHistory),function(i)
  boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])
names(lz) <- colnames(boruta.train$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta.train$ImpHistory), cex.axis = 0.7)


final.boruta <- TentativeRoughFix(boruta.train)

importance<-getSelectedAttributes(final.boruta, withTentative = F)

boruta.df <- attStats(final.boruta)



class(boruta.df)

totalnew<-as.data.frame(total$outcome)
for (i in seq_along(importance)){
  new<-get(paste("total",sep=""))[importance[i]]
  totalnew<-cbind(totalnew,new)
}
colnames(totalnew)[1]<-"outcome"  

control<-trainControl(method = "repeatedcv",
                      number = 10,
                      repeats = 3,
                      classProbs = TRUE,
                      summaryFunction = twoClassSummary)

svm.model=train(outcome~.,
                data = totalnew,
                method = "svmRadial",
                metric = "ROC",
                trControl = control)

randomforest.model=train(outcome~.,
                         data = totalnew,
                         method = 'rf',
                         metric = "ROC",
                         trControl = control)

neuralnetwork.model=train(outcome~.,
                          data = totalnew,
                          method = 'pcaNNet',
                          metric = "ROC",
                          trControl = control)
bayes.model=train(outcome~.,
                  data = totalnew,
                  method = "bartMachine",
                  metric = "ROC",
                  trControl = control)

luna.mol<-readMolFromSDF(c("prepare.sdf"))

luna = cbind(extractDrugMACCSComplete(luna.mol),
             extractDrugEstateComplete(luna.mol),
             extractDrugExtendedComplete(luna.mol),
             extractDrugGraphComplete(luna.mol))

luna<-as.data.frame(luna)

luna<-luna[importance]

luna.svm.porbs<-predict(svm.model,luna,
                   type = "prob")

luna.randomforest.porbs<-predict(randomforest.model,luna,
                            type = "prob")

luna.neuralnetwork.porbs<-predict(neuralnetwork.model,luna,
                             type = "prob")

luna.bayes.porbs<-predict(bayes.model,luna,
                                  type = "prob")

luna.result<-cbind(luna.svm.porbs,luna.neuralnetwork.porbs,
                   luna.randomforest.porbs,luna.bayes.porbs)


