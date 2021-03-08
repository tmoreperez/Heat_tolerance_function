#Thanks for checking out this code for estimating heat tolerances using variable fluorescence!
#Be sure to the read the Read.me.txt file.
#Below, the "psiiht" function can calculate heat tolerances of PSII heat tolerance. 
#For this function, define the the columns with the 1) Temperature & 2) FvFm variables, then the 3) control temperature, the 4) id/unique factor that you 
#want to estimate the heat tolerances for, 5) indicate if you want to make plots, 6) the number of bootsrap iterations for estimating 

setwd("/Users/timothyperez/Google Drive/Git_stuff/Git_projects/Heat_tolerance_function")
library(car)
htdata=read.csv("Sample_FvFm_data.csv")
psiiht=function( Temperature, FvFm, control.temp, id, plot.est, boots){
  l1=list(control.temp=control.temp, plot.est=plot.est, boots=boots)
  attach(l1)
  HTdf=data.frame(Temperature=Temperature, FvFm=FvFm, id=id)
  
 return(do.call("rbind", by(HTdf, list(HTdf$id), function(df){
   Temperature=df[,which(colnames(df)=="Temperature")]
   FvFm=df[,which(colnames(df)=="FvFm")]
   id=df[,which(colnames(df)=="id")]
   #get parameter estimates for logistic decay model
   cof=coef(lm(logit(FvFm)~Temperature)) 
   #Fit a non linear least squares model to the FvFm and Temperature data
   HT.model <- nls(FvFm ~ theta1/(1 + exp(-(theta2 + theta3*Temperature))),  start=list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
                   trace=F, control=list(maxiter=1000, tol=1e-3))
   
   #Use the parameter estimates (coef(HT.model)[#])from the HT.model to predict a new fit based on a heat treatments from 23-62 degrees celcius. Here, # = 1:3.
   y<-coef(HT.model)[1]/(1+exp(-(coef(HT.model)[2]+coef(HT.model)[3]*seq(23,62)))) 
   
   #Calculate half of the control Fv/Fm & a 95% reduction in FvFm with reference to control
   half=mean(na.omit(FvFm[which(Temperature==control.temp)]))/2  
   nine5=mean(na.omit(FvFm[which(Temperature==control.temp)]))*0.05  
   #95 Confidence Interval
   predict.boot=matrix(NA,40, boots)
   T95=T50=Tcrit=c()
   for(k in 1:boots){
     #print(k)
     srows <- sample(1:length(Temperature), length(Temperature),TRUE)
     
     if(class(try(nls(FvFm[srows] ~ theta1/(1 + exp(-(theta2 + theta3*Temperature[srows]))),  start=list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
                      trace=F, control=list(maxiter=1000, tol=1e-3)), silent=T)[[1]])=="nlsModel")
     {HT.model2 <- nls(FvFm[srows] ~ theta1/(1 + exp(-(theta2 + theta3*Temperature[srows]))),  start=list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
                       trace=F, control=list(maxiter=1000, tol=1e-6))
     predict.boot[,k]=coef(HT.model2)[1]/(1+exp(-(coef(HT.model2)[2]+coef(HT.model2)[3]*seq(23,62)))) 
     #Estimate T95
     T95[k]=(-log((coef(HT.model2)[1]/nine5)-1)-coef(HT.model2)[2])/coef(HT.model2)[3] 
     
     #Estimate T50
     T50[k]=(-log((coef(HT.model2)[[1]]/half)-1)-coef(HT.model2)[[2]])/coef(HT.model2)[[3]] #estimate ctmax
     T50k=(-log((coef(HT.model2)[[1]]/half)-1)-coef(HT.model2)[[2]])/coef(HT.model2)[[3]] #estimate ctmax
     #Use model to predict changes in FvFm & make new dataframe
     predict=data.frame(x=seq(23,62),y=coef(HT.model2)[1]/(1+exp(-(coef(HT.model2)[2]+coef(HT.model2)[3]*seq(23,62)))) ) #create the prediction data frame
     df1=cbind(predict[-1,], predict[-nrow(predict),])[,c(3,1,4,2)]
     #Use new dataframe to estimate the slope at between 1-degree intervals
     df1$slp=as.vector(apply(df1, 1, function(x) summary(lm((x[3:4]) ~ x[1:2])) [[4]][[2]] ))
     slp.at.tcrit=round(min(df1$slp), 3)*.15 #Determine where slope is 15% of max slope & round
     #Estimate the FvFm at which the slope is 15% of max slope & less than T50
     fvfv.at.tcrit=df1[which(abs(df1[which(df1[,1]<T50k),]$slp-slp.at.tcrit)==min(abs(df1[which(df1[,1]<T50k),]$slp-slp.at.tcrit))),][1,3]
     Tcrit[k]=(-log((coef(HT.model2)[[1]]/fvfv.at.tcrit)-1)-coef(HT.model2)[[2]])/coef(HT.model2)[[3]] # Estimate the temperatureat which the slope is 15% of max slope
     
     }else{(class(try(nls(FvFm ~ theta1/(1 + exp(-(theta2 + theta3*Temperature))),  start=list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
                          data=data2[srows,], trace=F, control=list(maxiter=1000, tol=1e-3)), silent=T)[[1]])=="list")
       predict.boot[,k]=NA
       T95[k]=NA
       T50[k]=NA
       Tcrit[k]=NA }}
   
   FvFm.boot=t(apply(predict.boot, 1, function(x){quantile(x,c(0.025,0.975),na.rm=T)}))
   
   #Tcrit.ci=quantile(Tcrit,c(0.025,0.975),na.rm=T)
   T50.ci=quantile(T50,c(0.025,0.975),na.rm=T)
   #T95.ci=quantile(T95,c(0.025,0.975),na.rm=T)
   
   if(plot.est==T){ 
     plot(NULL, NULL,xlab="Temperature",ylab="Fv/Fm",xlim=c(23,65),ylim=c(0,0.9), bty="l", lty=2)
     text(22,0.3, pos=4, paste(unique(paste(id))), font=4, cex=0.7)
     points(Temperature, FvFm, xlab="Temperature", ylab="Fv/Fm",xlim=c(23,65), ylim=c(0,0.85), bty="l", lty=2, pch=5,  col="black")
     lines(seq(23,62), y,lwd=1, col="purple")
     lines(seq(23,62), FvFm.boot[,1],lty=3, col="purple")
     lines(seq(23,62), FvFm.boot[,2],lty=3, col="purple")
     text(30,0, paste('Tcrit:', round(mean(na.omit(Tcrit)),1)), pos=4,col="dark blue")
     abline(v=round(mean(na.omit(Tcrit)),1), lty=2,lwd=1.5, col="dark blue",  cex=0.8)
     text(30,.1, paste('T50:',round(mean(na.omit(T50)),1)), pos=4, col="dark orange")
     abline(v=round(mean(na.omit(T50)),1), lty=2,lwd=1.5, col="dark orange",cex=0.8)
     text(30,.2, paste('T95:',round(mean(na.omit(T95)),1)),pos=4, col="dark red")
     abline(v=round(mean(na.omit(T95)),1), lty=2,lwd=1.5, col="dark red", cex=0.8)
   }
   return(data.frame(id=(unique(id)), 
                     #Tcrit.lci=round(Tcrit.ci[[1]],1),
                     Tcrit.mn=round(mean(na.omit(Tcrit)),1),  
                     #Tcrit.uci=round(Tcrit.ci[[2]],1),
                     
                     #T50.lci=round(T50.ci[[1]],1),
                     T50.mn=round(mean(na.omit(T50)),1),  
                     #T50.uci=round(T50.ci[[2]],1),
                     
                     #T95.lci=round(T95.ci[[1]],1),
                     T95.mn=round(mean(na.omit(T95)),1)  
                     #T95.uci=round(T95.ci[[2]],1)))
   ))
 })))
  detach(l1)
 }
psiiht(Temperature=htdata$Temperature, FvFm=htdata$FvFm, control.temp=23, id=htdata$id, plot.est=T, boots=100)



