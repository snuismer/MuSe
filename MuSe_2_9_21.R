library(ggplot2)
library(gridExtra)
library(tidyverse)
library(optimx)


filename=paste("Path\\FileName.csv",sep="")
#Ingest the Data from a csv file
Data=read.csv(filename, header = TRUE)

#Calculate number of replicates
reps=length(Data[1,])-2

#Trim the data file of trailing zeros and record the first day of "0" transgene for each replicate
TrailingTail=1
WhichRow=1
RepEnd=rep(length(Data[,1]),reps)#set the replicate end recording vector to the last row of the data and then records when each hits zero for the first time. I will calculate likelihood only over the range of non-zero values solving a couple of problems
while(TrailingTail<=2&&WhichRow<length(Data[,1]))#Note that this starts at 2 because we use the first time point to set po
{
  if((sum(Data[WhichRow,2:(reps+1)])==0))# Check to see if all replicates are zero. This is criteria for trimming the entire data file
  {
    TrailingTail=TrailingTail+1
  }
  for(k in 1:reps)#Take a look at each replicate and check to see if it has hit zero transgene frequency
  {
   if(Data[WhichRow,(k+1)]==0&&RepEnd[k]==length(Data[,1]))
   {
     RepEnd[k]=WhichRow#Note that the end is recorded as row number, not real time units...
   }
  }
  WhichRow=WhichRow+1
}
SubData=Data[1:WhichRow,1:(2+reps)]
#End of trimming

#Set intrinsic rate of increase
r=0.1
#calculate initial allele frequencies
stop=reps+1
po=c()
for(i in 2:stop)
{
  po=c(po,Data[1,i])
}
#Define the likelihood function using the vector of parameters Parms(Mu,S)
LogLike<-function(Parms)
{
  L=0
  for(k in 1:reps)#Sum likelihood over replicates
  {
    for(i in 2:(RepEnd[k]-1))#Note that this starts at 2 because we use the first time point to set po and ends at the timepoint preceeding each replicates first zero value
    {
      n=Data[i,reps+2]
      t=Data[i,1]
      p=-(exp(r*(Parms[2]*(Parms[1]-1)-Parms[1])*t)*po[k]*(Parms[2]*(Parms[1]-1)-Parms[1]))/(Parms[1]-Parms[2]*(-1+po[k]-exp(r*(Parms[2]*(Parms[1]-1)-Parms[1])*t)*po[k]+Parms[1]))
      #Nate added the next two lines to keep things from barfing on 0's. But p should never actually equal zero so unclear if the second is actually required?
      p[p==1] <- 1-(1e-12)
      p[p==0] <- 1e-12 
      x=n*Data[i,1+k]
      L=L+(n-x)*log(1-p)+x*log(p)
    }

  }
  return(-L)
}

#Use optimx to find the ML solution. Lots of testing against simulated data suggests these settings works well. 
MLSol=optimx(par=c(0.01,0.01), fn=LogLike,method = "nlminb",lower =c(0.0000001,0.0000001), upper = c(.3, .3),control=list(maxit=10000))

#Plot the Likelihood Function in the neighborhood of the solution using dynamic ranging so you can actually see things
Lo=.1*MLSol$p1
Hi=10*MLSol$p1
Width=(Hi-Lo)/20
Mu=seq(Lo,Hi,Width)
Lo=.5*MLSol$p2
Hi=2*MLSol$p2
Width=(Hi-Lo)/20
S=seq(Lo,Hi,Width)
Combos=expand.grid(Mu,S)
l=c()
NumToPlot=length(S)*length(Mu)

for(i in 1:NumToPlot){Parms=c(Combos[i,1],Combos[i,2]);l=c(l,LogLike(Parms)); print(i)}
MyData=cbind(Combos[,1],Combos[,2],l)
MyData=as.data.frame(MyData)
colnames(MyData) <- c("Mu", "S","Like")


#establish the min and max to scale figure
grandmin <- min(MyData$Like)
grandmax <- 2.0*min(MyData$Like)
contours<-20

#define the number of breaks.  In this case contours +1 
mybreaks <- round(seq(grandmin, grandmax,len = contours+1))
#Function to return the dersired number of colors
mycolors<- function(x) {
  colors<-colorRampPalette(c("darkblue", "yellow"))( contours )
  colors[1:x]
}

#Function to create labels for legend
breaklabel <- function(x){
  labels<- paste0(mybreaks[1:contours], "-", mybreaks[2:(contours+1)])
  labels[1:x]
}

#Plot the Likelihood surface and drop the estimate on top as a red point
c1=ggplot(MyData, aes(Mu, S, z=Like))+geom_contour_filled(breaks= mybreaks, show.legend = TRUE)
c2=c1+scale_fill_manual(palette=mycolors, values=breaklabel(contours), name="-Log Likelihood", drop=FALSE)
c3=c2+theme(legend.position = "right")
c4=c3+geom_point(data = data.frame(Mu = MLSol$p1, S = MLSol$p2,Like=grandmin), colour = "red", size = 3)
c5=c4+ggtitle("Likelihood surface and ML estimate") +xlab("Mutation") + ylab("Selection")+theme(plot.title = element_text(hjust = 0.5))
c5

#Next challenge is to check for lack of fit by back simulating onto the actual data and displaying this fit
#First, plot the dynamics of p using the ML estimates for mutation and selection
#Assume we start from the average initial allele frequency for simplicity of display. Should be fine unless wildly different
PredPlot=c()
pobar=mean(po)
for(i in 2:length(SubData[,1]))#Plot only up to the point where all reps have been zero for a threshold time (captured by Sub Data)
{
  t=Data[i,1]
  Plaques=Data[i,2+reps]
  p=-(exp(r*(MLSol$p2*(MLSol$p1-1)-MLSol$p1)*t)*pobar*(MLSol$p2*(MLSol$p1-1)-MLSol$p1))/(MLSol$p1-MLSol$p2*(-1+pobar-exp(r*(MLSol$p2*(MLSol$p1-1)-MLSol$p1)*t)*pobar+MLSol$p1))
  LowP=p-1.96*sqrt((p*(1-p))/Plaques)
  HiP=p+1.96*sqrt((p*(1-p))/Plaques)
  ThisData=c(t,p,LowP,HiP)
  PredPlot=rbind(PredPlot,ThisData,deparse.level = 0)
}
PredPlot=as.data.frame(PredPlot)
colnames(PredPlot) = c("Transfer","Pred","PredLo","PredHi")

#Make a plot of the re-simulated data using the ML estimate (+/- 95% CI shading) and scatterplot for the actual data
p1=ggplot() +geom_ribbon(data=PredPlot,aes(x = Transfer, ymin = PredLo, ymax = PredHi), fill = "grey70")
p2=p1+geom_line(data=PredPlot, aes(x = Transfer, y = Pred), color = "blue")
DataForPlot=SubData[,1:(reps+1)]
colnames(DataForPlot)=c("Transfer", 1:reps)
#The next line re-writes the data into a plottable format
DataForPlot <- gather(DataForPlot, "Replicate", "Freq", 2:ncol(DataForPlot))
p3=p2+geom_point(data = DataForPlot, aes(x=Transfer, y=Freq, color=Replicate), size= 2.5)
p4=p3+ ggtitle("Model fit") +xlab("Transfer time (hours)") + ylab("Frequency of transgene")+theme(plot.title = element_text(hjust = 0.5))
p4
