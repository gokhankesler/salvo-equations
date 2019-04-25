###############################################################################
#                                                                             #
#    CHAPTER IV. ANALYSIS PART 2: DETERMINISTIC SALVO MODEL (CLOSED-FORM)     #
#                                                                             #
###############################################################################

#-----------------INSTALLING AND LOADING REQUAIRED LIBRARIES------------------#
if(!require("MASS"))install.packages("MASS");library("MASS")
if(!require("visreg"))install.packages("visreg");library("visreg")
if(!require("plot3D"))install.packages("plot3D");library("plot3D")
if(!require("rpart.plot"))install.packages("plot3D");library("rpart.plot")
if(!require("RColorBrewer"))install.packages("RColorBrewer");library("RColorBrewer")


#------------------------------LOADING CSV DATA-------------------------------#
data=read.csv(file.choose(),header=T)

#--------------PAIR-WISE PLOTTING AND PRINTING CORRALATION TABLE--------------#
plot(data,pch=16,cex=.3,col="darkblue") # Page:41, Figure 17
round(cor(data),digits=3)               # Page:42, Figure 18

#----------------------DEFINING NULL VARIABLE FOR A AND B---------------------#
#----------TO DETERMINE"Intermediate", "Over-Defense" AND "Over-Kill"---------#
data$AS=rep("intermediate",dim(data)[1])
data$BS=rep("intermediate",dim(data)[1])

#---------------------------START OF THE SIMULATION---------------------------#
for(i in 1:dim(data)[1])	#Loop for each DP
{
  A=data[i:i,1]	# Number of Ships on Side A
  B=data[i:i,2]	# Number of Ships on Side B
  a=data[i:i,3]	# Offensive power of A
  b=data[i:i,4]	# Offensive power of B
  y=data[i:i,5]	# Defensive power of A
  z=data[i:i,6]	# Defensive power of B
  w=data[i:i,7]	# Staying power of A
  x=data[i:i,8]	# Staying power of B
  
  data$deltaA[i]=(B*b-A*y)/w	# Calculation of The Nominal 
  data$deltaB[i]=(A*a-B*z)/x	# Value of ΔA and ΔB.        
  
  if(data$deltaA[i]<0)                  #######################################
  {                                     #                                     #
    data$AS[i]="Over-defense"		        #                                     #
  }                                     # Determination of whether Side A     #
  else if(data$deltaA[i]>data$A[i])     # has an Over-defense or Over-killed. #
  {                                     #                                     #
    data$AS[i]="Over-destroyed"         #                                     #
  }                                     #######################################
  
  if(data$deltaB[i]<0)                  #######################################
  {                                     #                                     #
    data$BS[i]="Over-defense"		        #                                     #
  }                                     # Determination of whether Side B     #
  else if(data$deltaB[i]>data$B[i]) 	  # has an Over-defense or Over-killed. #
  {                                     #                                     #
    data$BS[i]="Over-destroyed"		      #                                     #
  }                                     #######################################
  
                                                          #####################
  data$deltaA[i]=min(A,max(0,data$deltaA[i]))             # Calculation for   #
  data$deltaB[i]=min(B,max(0,data$deltaB[i]))             # the Actual value  #
  data$FER[i]=((data$deltaB[i])/B)/((data$deltaA[i])/A)   # of ΔA, ΔB and FER #
  data$Bratio[i]=((data$deltaB[i])/B)                     #####################
}

#----------------------------END OF THE SIMULATION----------------------------#

data$BS=factor(data$BS)        # Converting the simulation destroy #
data$AS=factor(data$AS)        # results from string to the factor #
summary(data) # Obtaining simulation statistics.

#--------------------------------NORMALITY TEST--------------------------------#
shapiro.test(sample(data$deltaB,3000))
shapiro.test(sample(data$Bratio,3000))
shapiro.test(sample(data$FER,3000))

write.csv(data, file = "RawDetResult.csv",row.names =F,quote = F)

#--THE FIRST LINEAR MODEL OF THE SIMULATION BELONGING ACTUAL LOSS ON SIDE B---#
# IV-C-1 Page:44-50
lm1=lm(deltaB~(A+B+a+b+y+z+x+w)^2+I(A^2)+I(B^2)+I(a^2)+I(b^2)+I(y^2)+I(z^2)
   +I(x^2)+I(x^2)-1,data)        
lm1=stepAIC(lm1,direction="both")    # Stepwise Regression.
summary(lm1)

#----------------PLOTTING PAIRWISE SCATTER PLOT THAT COMPARES-----------------#
#--------------THE SIMULATION'S OUTPUT WITH THE LM'S PREDICTION---------------#
# Page:45, Figure 20
{par(mfrow=c(1,1))
 par(mgp=c(1.5,0.5,0))
plot(predict(lm1,data[which(data$BS=="intermediate"),]),
   data[which(data$BS=="intermediate"),]$deltaB,pch=16,cex=.7,col="grey"
   ,xlab="Predictions of the Metamodel",ylab="Hughes’ Salvo Equation Values"
   ,xlim=c(-5,19),ylim = c(0,19))
par(new=TRUE)
plot(predict(lm1,data[which(data$BS=="Over-defense"),]),
     data[which(data$BS=="Over-defense"),]$deltaB,pch=16,cex=.7,col="#2166AC",
     xlab="",ylab="",
     xlim=c(-5,19),
     ylim = c(0,19))
par(new=TRUE)
plot(predict(lm1,data[which(data$BS=="Over-destroyed"),]),
     data[which(data$BS=="Over-destroyed"),]$deltaB,pch=16,cex=.7,col="#B2182B",
     xlab="",ylab="",xlim=c(-5,19),ylim = c(0,19))

legend("bottomright",pch=c(16,16,16),col=c("grey","#B2182B","#2166AC"),
        c("Intermediate-Kill","Over-Kill","Over-Defense")
       ,cex=.8,text.font=2)
legend("topleft",paste0("R-Square:",round(summary(lm1)$r.squared,3))
       ,cex=.8,text.font=2)
abline(0,1,lwd=3)}

#-----------------------APPLYING THE SECOND LINEAR MODEL----------------------#
#---------------FOR THE RESULTS THAT END WITH INTERMEDIATE KILL---------------#
lm2=lm(deltaB~(A+B+a+b+y+z+x+w)^2+I(A^2)+I(B^2)+I(a^2)+I(b^2)+I(y^2)+I(z^2)
       +I(x^2)+I(x^2)-1,
       data[which(data$AS=="intermediate"&data$BS=="intermediate"),])
lm2=stepAIC(lm2,direction="both")    # Stepwise Regression.
summary(lm2)

#----------SCATTER PLOT THAT COMPARES THE RESULTS OF THE SIMULATION-----------#
#---------------------WITH THE SECOND MODEL'S PREDICTION----------------------#
# Page:46, Figure 21
{par(mfrow=c(1,1))
plot(pmin(pmax(predict(lm2,data[which(data$BS=="intermediate"),]
     ,type= "response"),0),data[which(data$BS=="intermediate"),]$B),
     data[which(data$BS=="intermediate"),]$deltaB,pch=16,cex=.7,col="grey",
     xlab="Predictions of the Metamodel 2",ylab="Hughes’ Salvo Equation Values",
     xlim=c(-5,19),
     ylim = c(0,19))
par(new=TRUE)
plot(pmin(pmax(predict(lm2,data[which(data$BS=="Over-defense"),]
     ,type = "response"),0),data[which(data$BS=="Over-defense"),]$B),
     data[which(data$BS=="Over-defense"),]$deltaB,pch=16,cex=.7,col="#2166AC",
     xlab="",ylab="",
     xlim=c(-5,19),
     ylim = c(0,19))
par(new=TRUE)
plot(pmin(pmax(predict(lm2,data[which(data$BS=="Over-destroyed"),],type = "response")
     ,0),data[which(data$BS=="Over-destroyed"),]$B),
     data[which(data$BS=="Over-destroyed"),]$deltaB,pch=16,cex=.7,col="#B2182B",
     xlab="",ylab="",
     xlim=c(-5,19),
     ylim = c(0,19))
legend("bottomright",pch=c(16,16,16),col=c("grey","#B2182B","#2166AC"),
       c("Intermediate-Kill","Over-Kill","Over-Defense")
       ,text.font=2)
legend("topleft",paste0("R-Square:",round(summary(lm(pmin(pmax(predict(lm2
       ,data,type = "response"),0),data$B)~data$deltaB-1))$r.squared,3))
       ,text.font=2)
abline(0,1,lwd=3)}

#--------------------------Prediction Profiler of ΔB--------------------------#
# Page:48, Figure 22
par(mfrow=c(2,4))
visreg(lm2,scale="response")   # The code that executes the prediction profiler

#--------------------------------TREE DIAGRAM---------------------------------#
# Page:50, Figure 24
par(mfrow=c(1,1))
rpart.plot(rpart(lm(deltaB~A+B+a+b+y+z+w+x,data)
   ,control=rpart.control(cp=0.006)),cex=.7)

#--------------------BUILDING LINEAR MODEL FOR LOSS RATIO---------------------#
# IV-C-2 Page:50-54
lm3=lm(Bratio~(A+B+a+b+y+z+x+w)^2+I(A^2)+I(B^2)+I(a^2)+I(b^2)+I(y^2)+I(z^2)
       +I(x^2)+I(x^2)-1,data)        # Modeling 
lm3=stepAIC(lm3,direction="both")    # Stepwise Regression.
summary(lm3)

#----------SCATTER PLOT THAT COMPARES THE RESULTS OF THE SIMULATION-----------#
#-------------------------WITH THE MODEL'S PREDICTION-------------------------#
# Page:51 Figure 25
par(mfrow=c(1,1))
plot(predict(lm3,data[which(data$BS=="intermediate"),],type = "response"),
     data[which(data$BS=="intermediate"),]$Bratio,pch=16,cex=.5,col="grey",
     xlab="Predictions of the Metamodel",ylab="Hughes’ Salvo Equation Values",
     xlim = c(-0.4,1.5),
     ylim = c(0,1))
par(new=TRUE)
plot(predict(lm3,data[which(data$BS=="Over-defense"),]),
     data[which(data$BS=="Over-defense"),]$Bratio,pch=16,cex=.5,col="#2166AC",
     xlab="",ylab="",
     xlim = c(-0.4,1.5),
     ylim = c(0,1))
par(new=TRUE)
plot(predict(lm3,data[which(data$BS=="Over-destroyed"),]),
     data[which(data$BS=="Over-destroyed"),]$Bratio,pch=16,cex=.5,col="#B2182B",
     xlab="",ylab="",
     xlim = c(-0.4,1.5),
     ylim = c(0,1))
legend("bottomright",pch=c(16,16,16),col=c("grey","#B2182B","#2166AC"),
       c("Intermediate-Kill","Over-Kill","Over-Defense")
       ,cex=.8,text.font=2)
legend("topleft",paste0("R-Square:",round(summary(lm3)$r.squared,3))
       ,cex=.8,text.font=2)
abline(0,1,lwd=2)

#-----------------------APPLYING THE SECOND LINEAR MODEL----------------------#
#---------------FOR THE RESULTS THAT END WITH INTERMEDIATE KILL---------------#
lm4=lm(Bratio~(A+B+a+b+y+z+x+w)^2+I(A^2)+I(B^2)+I(a^2)+I(b^2)+I(y^2)+I(z^2)
       +I(x^2)+I(x^2)-1,
       data[which(data$AS=="intermediate"&data$BS=="intermediate"),])
lm4=stepAIC(lm4,direction="both")
summary(lm4)

#----------SCATTER PLOT THAT COMPARES THE RESULTS OF THE SIMULATION-----------#
#----------------------WITH THE SECOND MODEL'S PREDICTION---------------------#
# Page:52 Figure 26
par(mfrow=c(1,1))
plot(pmin(pmax(predict(lm4,data[which(data$BS=="intermediate"),]
     ,type = "response"),0),1),
     data[which(data$BS=="intermediate"),]$Bratio,pch=16,cex=.5,col="grey",
     xlab="Predictions of the Metamodel 2",ylab="Hughes’ Salvo Equation Values",
     xlim = c(-0.4,1.5),
     ylim = c(0,1))
par(new=TRUE)
plot(pmin(pmax(predict(lm4,data[which(data$BS=="Over-defense"),],type = "response"),0),1),
     data[which(data$BS=="Over-defense"),]$Bratio,pch=16,cex=.7,col="#2166AC",
     xlab="",ylab="",
     xlim = c(-0.4,1.5),
     ylim = c(0,1))
par(new=TRUE)
plot(pmin(pmax(predict(lm4,data[which(data$BS=="Over-destroyed"),],type = "response"),0),1),
     data[which(data$BS=="Over-destroyed"),]$Bratio,pch=16,cex=.7,col="#B2182B",
     xlab="",ylab="",
     xlim = c(-0.4,1.5),
     ylim = c(0,1))
legend("bottomright",pch=c(16,16,16),col=c("grey","#B2182B","#2166AC"),
       c("Intermediate-Kill","Over-Kill","Over-Defense")
       ,cex=.8,text.font=2)
legend("topleft",paste0("R-Square:",round(summary(lm(pmin(pmax(predict(lm4,data,type = "response"),0),1)~data$Bratio))$r.squared,3))
       ,cex=.8,text.font=2)
abline(0,1,lwd=3)

#----------------DRAWING THE PREDICTION PROFILER OF FIRST MODEL---------------#
#---------------THAT PREDICTS THE NUMBER OF LOSS SHIP ON SIDE B---------------#
# Page:53 Figure 27
par(mfrow=c(2,4))
par(oma=c(0.5,0.5,0.5,0.8))    # Setting the margins for best view.
par(mar=c(3,3.5,1,0))
par(mgp=c(1.5,0.5,0))
visreg(lm4,scale="response")   # The code that executes the prediction profiler

#--------------------------------TREE DIAGRAM---------------------------------#
# Page:54 Figure 28
par(mfrow=c(1,1))
rpart.plot(rpart(lm(Bratio~A+B+a+b+y+z+w+x,data),
                 control=rpart.control(cp=0.009)),cex=.7)


#---------FER SENSITIVITY ANALYSIS---------#
# Page:55 Figure 29
hist(log(data[which(data$FER>0&data$FER<4000),]$FER),breaks = 40
     ,col = "lightblue",xlab = "log(FER)",main = "")
summary(log(data[which(data$FER>0&data$FER<4000),]$FER))
shapiro.test(log(data[which(data$FER>0&data$FER<4000),]$FER))

# Page:56 Figure 30
t=seq(2,50,1)
rs=NULL
lg=NULL
for (i in 1:length(t)) 
{
  print(i)
  lm5=glm(FER~(A+B+a+b+y+z+w+x)^2+I(A^2)+I(B^2)+I(a^2)+I(b^2)+I(y^2)
          +I(z^2)+I(w^2)+I(x^2)-1,data[which(data$FER<t[i]&data$FER>1/t[i]),],
          family=gaussian(link="log"))
  
  rs[i]= summary(lm(pmax(predict(lm5,data[which(data$FER<t[i]&data$FER>1/t[i]),]
  ,type="response"),0)~data[which(data$FER<t[i]&data$FER>1/t[i]),]$FER))$adj.r.squared
  lg[i]=length(data[which(data$FER>(1/t[i])&data$FER<t[i]),]$FER)
  
}
plot(t,rs,type = "l",lwd=2,ylab = "R-squared",xlab="FER")
plot(t,lg,type = "l",lwd=2,ylab = "Sample Size",xlab="FER")


#---- Generalized Linear Model of the Simulation Belonging The FER.---#

lm5=glm(FER~+(A+B+a+b+y+z+w+x)^2+I(A^2)+I(B^2)+I(a^2)+I(b^2)+I(y^2)
        +I(z^2)+I(w^2)+I(x^2)-1,data[which(data$FER<=10&data$FER>(1/10)),],
        family=gaussian(link="log"))

lm5=stepAIC(lm5,direction="both")
summary(lm5)
1-(lm5$deviance/lm5$null.deviance)
summary(lm(predict(lm5,data[which(data$FER<10&data$FER>1/10),]
  ,type="response")~data[which(data$FER<10&data$FER>1/10),]$FER))$r.squared

#-------Scatter plot that compares the results of the simulation----------#
#---------------------with the second model's prediction.-----------------#
# Page:57 Figure 31
par(mfrow=c(1,1))
scatter2D(predict(lm5,data[which(data$FER<10&data$FER>1/10),],type="response")
          ,data[which(data$FER<10&data$FER>1/10),]$FER
          ,colvar=data[which(data$FER<10&data$FER>1/10)
          ,]$deltaB/data[which(data$FER<10&data$FER>1/10),]$B,cex =1,clim=c(0,1)
          ,pch=16,xlab="Predictions of the Metamodel",ylab="Hughes’ Salvo Equation Values"
          ,clab = "ΔB/B"
          ,col=rev(c(brewer.pal(n=9,name="RdBu"))))

legend("topleft",paste0("R-Square:",round((summary(lm(predict(lm5
  ,data[which(data$FER<10&data$FER>1/10),]
  ,type="response")~data[which(data$FER<10&data$FER>1/10)
  ,]$FER))$r.squared),3)),cex=.8,text.font=2)
abline(0,1,lwd=3)

#----------------------PREDICTION PROFILER OF FER-----------------------------#
     
par(mfrow=c(2,4))
par(oma=c(0.5,0.5,0.5,0.8))    # Setting the margins for best view.
par(mar=c(3,3.5,1,0))
par(mgp=c(1.9,0.5,0))
visreg(lm5,scale="response")   # The code that executes the prediction profiler
     
#-----------------------------TREE DIAGRAM OF FER-----------------------------#
     
par(mfrow=c(1,1))
rpart.plot(rpart(lm(FER~A+B+a+b+y+z+w+x,
data[which(data$FER<=6.7&data$FER>(1/6.7)),]),
control=rpart.control(cp=0.01)),cex=1)
 