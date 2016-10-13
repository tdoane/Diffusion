 rm(list=ls())
 graphics.off()
 
 setwd('C:/Users/thdoa/Dropbox/Research/Nonlocal/Data')
 
 X_raw <- read.table('C:/Users/thdoa/Dropbox/Research/Nonlocal/Data/Sierra_Nev/B_XA.txt')
 Z_raw <- read.table('C:/Users/thdoa/Dropbox/Research/Nonlocal/Data/Sierra_Nev/B_ZA.txt')
 
 plot(X_raw[ ,1], Z_raw[ ,1])
 x11()
 
 L <- dim(X_raw)[1]
 
 interp <- approx(X_raw[ ,1],Z_raw[ ,1], method ="linear", n=L)

 
 min(X_raw)
 max(X_raw)
 
 Xtemp <- interp$x
 Ztemp <- interp$y
 ind = match(c(max(Ztemp)),Ztemp)
  
 Ls=L-ind
 Xnew=Xtemp[ind:L]
 Znew=Ztemp[ind:L]

 dx=diff(Xnew)[1]
  
 nL=501
 Lpad= ((nL-1)/2)-length(Xnew)+1
 X = seq(from= 0 , to= ((nL-1)/2)*dx, by=dx)
 X = append(seq(from= -((nL-1)/2)*dx , to=-dx, by=dx), X)
 pad = rep(0, times=Lpad)
 Ztemp2=Znew-min(Z_raw)
 Znew2=append(Ztemp2, pad)
 
 Zrev=rev(Znew2[1:((nL-1)/2)+1])
 Z=append(Zrev,Znew2)
 plot(X,Z, asp=1)
 x11()
 
 linfit=lm(Ztemp2 ~ Xnew)
 detZ=tempZ2-linfit
 plot(Xnew,detZ)
 x11()
 
 FT=fft(Z, inverse=FALSE)
 freq=(0:500)*2*pi/(500*dx)
 
 remove('Xnew', 'Ztemp', 'interp')
 
 Cost=NULL
 ztop=NULL
 
 s1=251-91
 s2=251+91
 ztop=c(29, 30, 31)
 dzT=1
 dtop=1
 k=0
 while (abs(dtop) > 0.0001){
   ztop[1]=ztop[2]-dzT
   ztop[3]=ztop[2]+dzT
   k=k+1
  for (i in 1:3){ 
  Ztemp=max(Z) + ztop[i] - 0.7*abs(X)
  Zfit=replace(Ztemp, Ztemp<0 ,0)
  DZ=Zfit-Z
  Cost[i]=sum((DZ))^2
  remove('Zfit', 'Ztemp')
  }
  con=(Cost[1] -2*Cost[2] + Cost[3])/(dzT^2)
  slp=(Cost[3] - Cost[1])/(2*dzT)
  
  dtop=-slp/con
  ztop[2]=ztop[2] + dtop;
 }
 top=ztop[2]

 
 Ztemp=max(Z) + top - 0.7*abs(X)
 Zfit=replace(Ztemp, Ztemp<0 ,0)
 Cost[4]=sum((Zfit-Z)^2)
 sum(Zfit-Z)
 
 plot(X,Zfit,type='l', asp=1, ylab=expression(zeta), xlab=expression(x))
 points(X,Z,type='l', col='purple', asp=1)
 legend(x="topright",c("Initial Condition", "Observed Evolved Moraine"), col=c('black', 'purple'), lty=1, bty="n", cex=1.0, lwd=1)
 x11()
  
 #plot(Cost)
 #x11()
 
 FTmod=fft(Zfit, inverse=FALSE)
 
 #K=c(0.01,0.02,0.03)
 delK=0.0001
 t=40000
 
 FTcost=NULL
 Kplot=NULL
 Costplot=NULL
 K=c(0.0001, 0.0002, 0.0003)
 
 plot(K[2],1)
 cnt=0
 dK=1
 
 while ( abs(dK) > 0.00001) {
   cnt=cnt+1
   K[1]=K[2] - delK
   K[3]=K[2] + delK
   FTcost=NULL
  for (i in 1:3){
    FTTime=FTmod*exp(-K[i]*t*freq^2)
    DF=Mod(FTTime[1:50])-Mod(FT[1:50])
    FTcost[i]=sum((DF)^2)
    remove('DF', 'FTTime')
  }
   conDF=(FTcost[1] - 2*FTcost[2] + FTcost[3])/(delK^2)
   slpDF=(FTcost[3] - FTcost[1])/(2*delK)
   
   dK[1]=slpDF/conDF
   Knew=K[2] - dK[1]
   K[2]=Knew
   Kplot[cnt]=K[2]
   Costplot[cnt]=FTcost[2]
   remove('FTcost')
 }
 #plot(Kplot, Costplot)
 #x11()
 Kbest=K[2]

 FTevol=FTmod*exp(-Kbest*t*freq^2)
 
 ylab.name=expression(Z(k,t))
 xlab.name=expression(k)
 
#pdf('FT_moraine.pdf', width=4, height=4)

 
 plot(freq[1:30],Mod(FTmod[1:30]),'l', ylab=ylab.name,xlab=xlab.name, cex=1.5)
 points(freq[1:30],Mod(FT[1:30]), 'l',col='purple')
 points(freq[1:30],Mod(FTevol[1:30]),col='black', pch=18)
 legend(x="topright",c("Initial Condition", "Observed Evolved Moraine"), col=c('black', 'purple'), lty=1, bty="n", cex=1.0, lwd=1)
 legend(0.09, 5900, "Analytical Solution", col='black', pch=18, bty="n", cex=1.0)

#dev.off() 
x11()
 
 
 zInv=fft(FTevol, inverse=TRUE)/251
 zInv2=Re(zInv)-min(Re(zInv))
 plot(X,zInv2,'l',asp=1, col='blue')
 points(X,Z,'l',col='red')
 points(X,Zfit,type='l')
 
 
 