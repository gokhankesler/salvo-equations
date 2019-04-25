###############################################################################
#                                                                             #
#      CHAPTER V. ANALYSIS PART 2: STOCHASTIC SALVO MODEL (CLOSED-FORM)       #
#                                                                             #
###############################################################################

#-----------------INSTALLING AND LOADING REQUAIRED LIBRARIES------------------#
if (!require("MASS")) install.packages("MASS"); library("MASS")
if (!require("visreg")) install.packages("visreg"); library("visreg")
if (!require("plot3D")) install.packages("plot3D"); library("plot3D")
if (!require("rpart.plot")) install.packages("plot3D"); library("rpart.plot")
if (!require("RcolorBrewer")) install.packages("RcolorBrewer"); library("RcolorBrewer")

#------------------------------LOADING CSV DATA-------------------------------#
sdata <- read.csv(file.choose(),header = T)

#--------------PAIR-WISE PLOTTING AND PRINTING CORRALATION TABLE--------------#
plot(sdata, pch=16, cex=.4, col="darkblue") # Page:65, Figure 35
round(cor(sdata ),2)                        # Page:66, Figure 36

#------------------------------NUMBER OF TRIALS-------------------------------#
n=50000

#---------------------------START OF THE SIMULATION---------------------------#
for (i in 1:dim(sdata)[1])
{
  if(i==1) print(Sys.time())
  print(i)
  # Predictors
  A=sdata$A[i]     # Beginning force strength for side A
  B=sdata$B[i]     # Beginning force strength for side B
  na=sdata$na[i]   # Maximum # of SSMs per ship per salvo for side A 
  nb=sdata$nb[i]   # Maximum # of SSMs per ship per salvo for side B
  pa=sdata$pa[i]   # Prob. of successful launch of single SSM for side A
  pb=sdata$pb[i]   # Prob. of successful launch of single SSM for side B
  ny=sdata$ny[i]   # Maximum # of SAMs per ship per salvo for side A
  nz=sdata$nz[i]   # Maximum # of SAMs per ship per salvo for side B
  py=sdata$py[i]   # Prob. of successful launch of single SAM for side A
  pz=sdata$pz[i]   # Prob. of successful launch of single SAM for side B
  u=sdata$u[i]     # Mean losses per hit suffered by side A
  v=sdata$v[i]     # Mean losses per hit suffered by side B
  sdu=sdata$sdu[i] # Sd of mean losses per hit suffered by side A
  sdv=sdata$sdv[i] # Sd of mean losses per hit suffered by side B
  #--------------------------------CALCULATIONS-------------------------------#
  A_offensive<-rbinom(n,A*na,pa) # Offensive power of A
  B_offensive<-rbinom(n,B*nb,pb) # Offensive power of B
  A_defensive<-rbinom(n,A*ny,py) # Defensive power of A
  B_defensive<-rbinom(n,B*nz,pz) # Defensive power of B
  
  NetBA=A_offensive-B_defensive # Nominal non-intercepted SSM's from A to B
  NetAB=B_offensive-A_defensive # Nominal non-intercepted SSM's from B to A
  
  A_damage=NULL          # Temporary variable for calculating ΔA
  B_damage=NULL          # Temporary variable for calculating ΔB
  
  for (j in 1:n) 
  {
    if(NetAB[j]>0)                           
    {                                         
      A_damage[j]=sum(rnorm(NetAB[j],u,sdu))
    }  
    else if(NetAB[j]<0)                      
    {                                                
      A_damage[j]=-sum(rnorm(-NetAB[j],u,sdu))   
    }
    else
    {
      A_damage[j]=0
    } 
    
    
    if(NetBA[j]>0)                                 
    {                                               
      B_damage[j]=sum(rnorm(NetBA[j],v,sdv))
    }  
    else if(NetBA[j]<0)                              
    {                                           
      B_damage[j]=-sum(rnorm(-NetBA[j],v,sdv)) 
    }
    else
    {
      B_damage[j]=0
    } 
    
  }
  
  
  #--------------------SAVING SUMMARIZED DATA FOR EACH DP---------------------#
  sdata$B_Loss[i]=mean(pmax(pmin(B_damage,B),0))
  sdata$B_Loss_sd[i]=sd(pmax(pmin(B_damage,B),0))
  
  sdata$FER[i]=mean((pmax(pmin(B_damage,B),0)/B)/(pmax(pmin(A_damage,A),0)/A))
  sdata$FER_sd[i]=mean((pmax(pmin(B_damage,B),0)/B)/(pmax(pmin(A_damage,A),0)/A))
  
  sdata$nolossB[i]=sum(pmax(pmin(B_damage,B),0)==0)/n  
  sdata$destroyedB[i]=sum(pmax(pmin(B_damage,B),0)==B)/n   
  sdata$Bint[i]=sum(pmax(pmin(B_damage,B),0)!=0 & pmax(pmin(B_damage,B),0)!=B)/n
  sdata$Bwin[i]=sum(pmax(pmin(B_damage,B),0)!=B & pmax(pmin(A_damage,A),0)==A)/n
  
  if(i==dim(sdata)[1]) print(Sys.time())
  
} 
#----------------------------END OF THE SIMULATION----------------------------#
write.csv(sdata, file = "RawStoResult.csv",row.names =F,quote = F)
summary(sdata)
shapiro.test(sdata$B_Loss)
shapiro.test(sdata$FER)
shapiro.test(sdata$Bwin)
#-------------------REGRESSION MODEL OF LOSS SHIP ON SIDE B-------------------#
# slm: stochastic linear model
slm1<-lm(B_Loss~(A+B+na+nb+pa+pb+ny+nz+py+pz+u+v+sdu+sdv)^2+I(A^2)+I(B^2)
         +I(na^2)+I(nb^2)+I(ny^2)+I(nz^2)+I(u^2)+I(v^2)+I(pa^2)+I(pb^2)
         +I(py^2)+I(pz^2)+I(sdu^2)+I(sdv^2)-1,sdata)
summary(slm1)
slm1<-stepAIC(slm1, direction="both")
summary(slm1)

#----------------PLOTTING PAIRWISE SCATTER PLOT THAT COMPARES-----------------#
#--------------THE SIMULATION'S OUTPUT WITH THE LM'S PREDICTION---------------#
# Page:69, Figure 38
par(oma=c(0,0,0,0))
par(mar=c(4,4,4,4))
par(mgp=c(1.4,0.4,0))

scatter2D(predict(slm1,sdata,type="response"),sdata$B_Loss
   ,colvar=sdata$B_Loss/sdata$B,cex =1,clim=c(0,1)
   ,pch=16,xlab="Predictions of the Metamodel",ylab="Armstrong’s Salvo Equation Values"
   ,clab = "ΔB/B"
   ,col=rev(c(brewer.pal(n=9,name="RdBu"))))
legend("topleft",c(paste0("R-Square:",round(summary(lm(predict(slm1,sdata
   ,type="response")~sdata$B_Loss-1))$adj.r.squared,3))),text.font=2)
abline(0,1,lwd=2)

# Page:70, Figure 39
scatter2D(predict(slm1,sdata,type="response"),sdata$B_Loss
          ,colvar=sdata$nolossB,col=c(brewer.pal(n=9,name="Blues")),xlim=c(-4,19)
          ,ylim=c(0,18),cex =1,clim=c(0,1),pch=16,colkey = F)
colkey(col=c(brewer.pal(n=9,name="Blues")),clim=c(0,1),add=TRUE,side=3
       ,length=0.9)
abline(0,1,lwd=2)
legend("topleft",c(paste0("R-Square:",round(summary(lm(predict(slm1
,sdata[which(sdata$nolossB>sdata$Bint&sdata$nolossB>sdata$destroyedB),]
,type="response")~sdata[which(sdata$nolossB>sdata$Bint&sdata$nolossB>sdata
$destroyedB),]$B_Loss-1))$r.squared,3))
,paste0("N =",sum(sdata$nolossB>sdata$Bint&sdata$nolossB>sdata$destroyedB)))
,text.font=2)
scatter2D(predict(slm1,sdata,type = "response"),sdata$B_Loss,colvar=sdata$Bint
,col=c(brewer.pal(n=9,name="Purples")),xlim=c(-4,19),ylim=c(0,18),cex=1,
clim=c(0,1),pch=16,colkey=F)
colkey(col=c(brewer.pal(n=9,name="Purples")),clim=c(0,1),add=TRUE,side=3
,length=0.9)
abline(0,1,lwd=2)
legend("topleft",c(paste0("R-Square:",round(summary(lm(predict(slm1
,sdata[which(sdata$Bint>sdata$nolossB&sdata$Bint>sdata$destroyedB),]
,type = "response")~sdata[which(sdata$Bint>sdata$nolossB&sdata$Bint>sdata
$destroyedB),]$B_Loss-1))$r.squared,3)),paste0("N = "
,sum(sdata$Bint>sdata$nolossB&sdata$Bint>sdata$destroyedB)))
,text.font=2)
scatter2D(predict(slm1,sdata,type = "response"),sdata$B_Loss
,colvar=sdata$destroyedB,col=c(brewer.pal(n=9,name="Reds")),xlim=c(-4,19)
,ylim=c(0,18),ylab ="",cex=1,clim =c(0,1),pch=16,colkey=F)
colkey(col=c(brewer.pal(n=9,name="Reds")),clim=c(0,1),add=TRUE,side=3
,length=0.9 )
abline(0,1,lwd=2)
legend("topleft",c(paste0("R-Square:",round(summary(lm(predict(slm1
,sdata[which(sdata$destroyedB>sdata$nolossB&sdata$destroyedB>sdata$Bint),]
,type = "response")~sdata[which(sdata$destroyedB>sdata$nolossB&sdata
$destroyedB>sdata$Bint),]$B_Loss-1))$r.squared,3)),paste0("N = "
,sum(sdata$destroyedB>sdata$nolossB&sdata$destroyedB>sdata$Bint)))
,text.font=2)

#--------REFITTIN LM OF LOSS SHIP WHERE EACH SIDE HAS INTERMEDIATE KILL-------#

slm2=slm2<-lm(B_Loss~(A+B+na+nb+pa+pb+ny+nz+py+pz+u+v+sdu+sdv)^2+I(A^2)+I(B^2)
   +I(na^2)+I(nb^2)+I(ny^2)+I(nz^2)+I(u^2)+I(v^2)+I(pa^2)+I(pb^2)+I(py^2) 
   +I(pz^2)+I(sdu^2)+I(sdv^2)-1,sdata[which(sdata$B_Loss/sdata$B<=0.95 & sdata$B_Loss/sdata$B>=0.02 ),])
slm2<-stepAIC(slm2, direction="both")
summary(slm2) 

#----------------PLOTTING PAIRWISE SCATTER PLOT THAT COMPARES-----------------#
#-----------THE SIMULATION'S OUTPUT WITH THE SECOND LM'S PREDICTION-----------#
# Page:71, Figure 40
par(oma=c(0,0,0,0))
par(mar=c(4,4,4,4))


scatter2D(pmax(pmin(predict(slm2,sdata),sdata$B),0),sdata$B_Loss
          ,colvar=sdata$B_Loss/sdata$B,cex =1,clim=c(0,1)
          ,pch=16,xlab="Predictions of the Metamodel",ylab="Armstrong’s Salvo Equation Values"
          ,clab = "ΔB/B"
          ,col=rev(c(brewer.pal(n=9,name="RdBu"))))
legend("topleft",c(paste0("R-Square:",round(summary(lm(pmax(pmin(predict(slm2,sdata),sdata$B),0)~sdata$B_Loss-1))$adj.r.squared,3))),text.font=2)
abline(0,1,lwd=2)

# Page:71, Figure 41
scatter2D(pmin(pmax(predict(slm2,sdata,type="response"),0),sdata$B),sdata$B_Loss,colvar=sdata$nolossB,col=c(brewer.pal(n=9,name="Blues")),xlim=c(0,18),ylim=c(0,18),cex =1,clim=c(0,1),pch=16,colkey = F)
colkey(col=c(brewer.pal(n=9,name="Blues")),clim=c(0,1),add=TRUE,side=3,length=0.9)
abline(0,1,lwd=2)
legend("topleft",c(paste0("R-Square:",round(summary(lm(pmin(pmax(predict(slm2,sdata[which(sdata$nolossB>sdata$Bint&sdata$nolossB>sdata$destroyedB),],type="response"),0),sdata[which(sdata$nolossB>sdata$Bint&sdata$nolossB>sdata$destroyedB),]$B)~sdata[which(sdata$nolossB>sdata$Bint&sdata$nolossB>sdata$destroyedB),]$B_Loss-1))$r.squared,3)),paste0("N =",sum(sdata$nolossB>sdata$Bint&sdata$nolossB>sdata$destroyedB))),text.font=2)
scatter2D(pmin(pmax(predict(slm2,sdata,type="response"),0),sdata$B),sdata$B_Loss,colvar=sdata$Bint,col=c(brewer.pal(n=9,name="Purples")),xlim=c(0,18),ylim=c(0,18),cex=1,clim=c(0,1),pch=16,colkey=F)
colkey(col=c(brewer.pal(n=9,name="Purples")),clim=c(0,1),add=TRUE,side=3,length=0.9)
abline(0,1,lwd=2)
legend("topleft",c(paste0("R-Square:",round(summary(lm(pmin(pmax(predict(slm2,sdata[which(sdata$Bint>sdata$nolossB&sdata$Bint>sdata$destroyedB),],type="response"),0),sdata[which(sdata$Bint>sdata$nolossB&sdata$Bint>sdata$destroyedB),]$B)~sdata[which(sdata$Bint>sdata$nolossB&sdata$Bint>sdata$destroyedB),]$B_Loss-1))$r.squared,3)),paste0("N = ",sum(sdata$Bint>sdata$nolossB&sdata$Bint>sdata$destroyedB))),text.font=2)
scatter2D(pmin(pmax(predict(slm2,sdata,type="response"),0),sdata$B),sdata$B_Loss,colvar=sdata$destroyedB,col=c(brewer.pal(n=9,name="Reds")),xlim=c(0,18),ylim=c(0,18),ylab ="",cex=1,clim =c(0,1),pch=16,colkey=F)
colkey(col=c(brewer.pal(n=9,name="Reds")),clim=c(0,1),add=TRUE,side=3,length=0.9 )
abline(0,1,lwd=2)
legend("topleft",c(paste0("R-Square:",round(summary(lm(pmin(pmax(predict(slm2,sdata[which(sdata$destroyedB>sdata$nolossB&sdata$destroyedB>sdata$Bint),],type="response"),0),sdata[which(sdata$destroyedB>sdata$nolossB&sdata$destroyedB>sdata$Bint),]$B)~sdata[which(sdata$destroyedB>sdata$nolossB&sdata$destroyedB>sdata$Bint),]$B_Loss-1))$r.squared,3)),paste0("N = ",sum(sdata$destroyedB>sdata$nolossB&sdata$destroyedB>sdata$Bint))),text.font=2)


#---------------------------PREDICTOR PROFILER OF ΔB--------------------------#
# Page:72, Figure 42
par(mfrow = c(4, 4))
par(oma=c(0.5, 0.5, 0.5, 0.8))
par(mar=c(3  , 3.5, 1  , 0  ))
par(mgp=c(2.1, 0.5, 0 ))
visreg(slm2,scale = "response")
par(mfrow = c(1, 1))

#------------------------------TREE DIAGRAM OF ΔB-----------------------------#
# Page:74, Figure 44
rpart.plot(rpart(lm(B_Loss~A+B+na+nb+pa+pb+ny+nz+py+pz+u+v+sdu+sdv,sdata[which(sdata$Bint>=0.42),]),
                 control=rpart.control(cp=0.012)),cex=.68)

#---------------------------REGRESSION MODEL OF FER---------------------------#

# Sensitiviy Analysis
# Page:75, Figure 45
hist(sdata$FER )
hist(log(sdata$FER),breaks = 40,main = "",col = "lightblue")

t=seq(2,100,1)
rs=NULL
lg=NULL
for (i in 1:length(t)) 
{
  print(i)
  lm5=glm(FER~(A+B+na+nb+pa+pb+ny+nz+py+pz+u+v+sdu+sdv)^2+I(A^2)+I(B^2)+
            I(na^2)+I(nb^2)+I(ny^2)+I(nz^2)+I(u^2)+I(v^2)+I(pa^2)+I(pb^2)+I(py^2) +
            I(pz^2)+I(sdu^2)+I(sdv^2)-1,sdata[which(sdata$FER<t[i]&sdata$FER>1/t[i]),],
          family=gaussian(link="log"))
  

  rs[i]=summary(lm(predict(lm5,sdata[which(is.finite(sdata$FER)),],type="response")~sdata[which(is.finite(sdata$FER)),]$FER))$adj.r.squared
  lg[i]=length(sdata[which(sdata$FER>(1/t[i])&sdata$FER<t[i]),]$FER)
}
# Page:76, Figure 46
plot(t,rs,type = "l",lwd=2,ylab = "R-squared",xlab="FER")
plot(t,lg,type = "l",lwd=2,ylab = "Sample Size",xlab="FER")#)
######


slm3=glm(FER~(A+B+na+nb+pa+pb+ny+nz+py+pz+u+v+sdu+sdv)^2+I(A^2)+I(B^2)+
          I(na^2)+I(nb^2)+I(ny^2)+I(nz^2)+I(u^2)+I(v^2)+I(pa^2)+I(pb^2)+I(py^2) +
          I(pz^2)+I(sdu^2)+I(sdv^2)-1,sdata[which(sdata$FER<100&sdata$FER>0),],
        family=gaussian(link="log"))
slm3<- stepAIC(slm3, direction="both")
1-(slm3$deviance/slm3$null.deviance)

#----------------PLOTTING PAIRWISE SCATTER PLOT THAT COMPARES-----------------#
#-----------THE SIMULATION'S OUTPUT WITH THE SECOND LM'S PREDICTION-----------#
# Page:76, Figure 47
par(mfrow=c(1,1))

scatter2D(predict(slm3,sdata[which(sdata$FER<100&sdata$FER>0),],type="response")
          ,sdata[which(sdata$FER<100&sdata$FER>0),]$FER,colvar=sdata[which(sdata$FER<100&sdata$FER>0),]$B_Loss /sdata[which(sdata$FER<100&sdata$FER>0),]$B,cex =1,clim=c(0,1)
          ,pch=16,xlab="Predictions of the Metamodel",ylab="Armstrong’s Salvo Equation Values"
          ,clab = "ΔB/B"
          ,col=rev(c(brewer.pal(n=9,name="RdBu"))))

legend("topleft",paste0("R-Square:",round((summary(lm(predict(slm3,sdata[which(sdata$FER<100&sdata$FER>0),],type="response")~sdata[which(sdata$FER<100&sdata$FER>0),]$FER))$r.squared
),3))
,cex=.8,text.font=2)
abline(0,1,lwd=2)

#---------------------------PREDICTOR PROFILER OF FER--------------------------#
# Page:77, Figure 48
par(mfrow = c(4, 4))
par(oma=c(0.5, 0.5, 0.5, 0.8))
par(mar=c(3  , 3.5, 1  , 0  ))
par(mgp=c(2.1, 0.5, 0 ))
visreg(slm3,scale = "response")

#------------------------------TREE DIAGRAM OF ΔB-----------------------------#
# Page:78, Figure 49
par(mfrow = c(1, 1))
rpart.plot(rpart(lm(FER~A+B+na+nb+pa+pb+ny+nz+py+pz+u+v+sdu+sdv,sdata[which(sdata$FER<5&sdata$FER>1/5),]),
                 control=rpart.control(cp=0.005)),cex=.85)


#---------------------------REGRESSION MODEL OF WIN---------------------------#
slm4=lm(Bwin~(A+B+na+nb+pa+pb+ny+nz+py+pz+u+v+sdu+sdv)^2+I(A^2)+I(B^2)+
           I(na^2)+I(nb^2)+I(ny^2)+I(nz^2)+I(u^2)+I(v^2)+I(pa^2)+I(pb^2)+I(py^2) +
           I(pz^2)+I(sdu^2)+I(sdv^2)-1,sdata)
slm4<- stepAIC(slm4, direction="both")
summary(slm4)

#----------------PLOTTING PAIRWISE SCATTER PLOT-----------------#
# Page:79, Figure 50

scatter2D(predict(slm4,sdata,type="response")
          ,sdata$Bwin,colvar=sdata$B_Loss/sdata$B ,cex =1,clim=c(0,1)
          ,pch=16,xlab="Predictions of the Metamodel",ylab="Armstrong’s Salvo Equation Values"
          ,clab = "ΔB/B", xlim=c(-0.2,1.2)
          ,col=rev(c(brewer.pal(n=9,name="RdBu"))))
legend("bottomright",paste0("R-Square: ",round(summary(slm4)$adj.r.squared,3)),
       cex=.8,text.font=2)

abline(0,1,lwd=2)

#-----------------SECOND LINEAR MODEL OF BWIN WITH LOGISTIC REGRESSION----#
slm5=glm(Bwin~(A+B+na+nb+pa+pb+ny+nz+py+pz+u+v+sdu+sdv)^2+I(A^2)+I(B^2)+
          I(na^2)+I(nb^2)+I(ny^2)+I(nz^2)+I(u^2)+I(v^2)+I(pa^2)+I(pb^2)+I(py^2) +
          I(pz^2)+I(sdu^2)+I(sdv^2)-1,sdata,family=binomial)
slm5<- stepAIC(slm5, direction="backward")
summary(slm5)
1-(slm5$deviance/slm5$null.deviance)

#----------------PLOTTING PAIRWISE SCATTER PLOT THAT COMPARES-----------------#
#-----------THE SIMULATION'S OUTPUT WITH THE SECOND LM'S PREDICTION-----------#
# Page:80, Figure 51
scatter2D(predict(slm5,sdata,type="response")
          ,sdata$Bwin,colvar=sdata$B_Loss/sdata$B  ,cex =1,clim=c(0,1)
          ,pch=16,xlab="Predictions of the Metamodel",ylab="Armstrong’s Salvo Equation Values"
          ,clab = "Prob. of Win", xlim=c(-0.2,1.2)
          ,col=rev(c(brewer.pal(n=9,name="RdBu"))))
legend("bottomright",paste0("R-Square: ",round(1-(slm5$deviance/slm5$null.deviance),3)),
       cex=.8,text.font=2)
abline(0,1,lwd=2)

#---------------------------PREDICTOR PROFILER OF FER--------------------------#
# Page:81, Figure 52
par(mfrow = c(4, 4))
par(oma=c(0.5, 0.5, 0.5, 0.8))
par(mar=c(3  , 3.5, 1  , 0  ))
par(mgp=c(2.1, 0.5, 0 ))
visreg(slm5,scale = "response"
       ,cond = list(A=14,B=14,  na=4, nb= 4, pa=0.9, pb=0.9,  ny=1,  nz= 1
      ,py=0.677, pz=1, u=0.4, v= 0.4, sdu=0.2, sdu=0.2))

#------------------------------TREE DIAGRAM OF ΔB-----------------------------#
# Page:82, Figure 53
par(mfrow = c(1, 1))
rpart.plot(rpart(lm(Bwin~A+B+na+nb+pa+pb+ny+nz+py+pz+u+v+sdu+sdv
                  ,sdata[which(sdata$FER<5&sdata$FER>1/5),]),
                 control=rpart.control(cp=0.0005)),cex=.75)



summary(slm5)
