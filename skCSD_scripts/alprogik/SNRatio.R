#Calculating signal to noise ratio
#Let's assume the measured potential on the 
#electrode is corrupted with noise
Signal<-read.table("myLFP.txt",colClasses="numeric") #measured potential
StanDev<-3
Noise<-rnorm(prod(dim(Signal)),0,StanDev)
PowerofSignal<-1/prod(dim(Signal))*sum(Signal^2)
PowerofNoise<-1/prod(dim(Noise))*sum(Noise^2)
SignaltoNoise<-PowerofSignal/PowerofSignal

#Lets calculate the error of the reconstruction
#C.calc<-1/const*t(K.tilda)%*%ginv(K)%*%LFP 
sigma<-0.5
const<-1/(4*pi*sigma)
OriginalCurrentsSmoothed<-read.table("membcurr_line",colClasses="numeric")
OriginalCurrentsSmoothed<-read.table("membcurr_smoothed",colClasses="numeric")
TransferMatrixes<-read.table("transfermatrix_MC",colClasses="numeric")
C.calc<-1/const*TransferMatrixes%*%Signal

#Error with line current
current.error.abs<-sum((C.calc-memb.currents.vonal)^2)/sum((memb.currents.vonal)^2)