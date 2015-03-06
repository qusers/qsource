#2345678901234567890123456789012345678901234567890123456789012345678901234567890
################################################################################
# Property plotter for Q runs
################################################################################
heatw <- grep("^SUM",readLines("../start.log"), value=TRUE)
nshw  <- grep("steps",readLines("../start.inp"), value=TRUE)    #numberofsteps
sshw  <- grep("stepsize",readLines("../start.inp"), value=TRUE)    #stepsize

#equiw <- grep("system  - epot     :",readLines("../water/out/equilibrate.out"), value=TRUE)
#nsew  <- grep("log_write_fq",readLines("../water/equilibrate.inp"), value=TRUE) #numberofsteps
#ssew   <- grep(" ss ",readLines("../water/equilibrate.inp"), value=TRUE)   #stepsize

#fepw  <- grep("system  - epot     :",readLines("../water/out/fep.out"), value=TRUE)
#nsfw  <- grep("log_write_fq",readLines("../water/fep.inp"), value=TRUE) #numberofsteps
#ssfw   <- grep(" ss ",readLines("../water/fep.inp"), value=TRUE)   #stepsize

#heatp <- grep("system  - epot     :",readLines("../protein/out/heat.out"), value=TRUE)
#nshp  <- grep("log_write_fq",readLines("../protein/heat.inp"), value=TRUE) #numberofsteps
#sshp   <- grep(" ss ",readLines("../protein/heat.inp"), value=TRUE)   #stepsize

#equip <- grep("system  - epot     :",readLines("../protein/out/equilibrate.out"), value=TRUE)
#nsep  <- grep("log_write_fq",readLines("../protein/equilibrate.inp"), value=TRUE) #numberofsteps
#ssep   <- grep(" ss ",readLines("../protein/equilibrate.inp"), value=TRUE)   #stepsize

#fepp  <- grep("system  - epot     :",readLines("../protein/out/fep.out"), value=TRUE)
#nsfp  <- grep("log_write_fq",readLines("../protein/fep.inp"), value=TRUE) #numberofsteps
#ssfp   <- grep(" ss ",readLines("../protein/fep.inp"), value=TRUE)   #stepsize

#     1 == water
#     2 == protein

ekinh1 <- read.table(text = heatw)[2:(length(heatw)-1), 4]
epoth1 <- read.table(text = heatw)[2:(length(heatw)-1), 3]
etoth1 <- read.table(text = heatw)[2:(length(heatw)-1), 2]
fsshw <- read.table(text=nshw)[,2]*read.table(text=sshw)[,2]

#ekine1 <- read.table(text = equiw)[2:(length(equiw)-1), 5]
#epote1 <- read.table(text = equiw)[2:(length(equiw)-1), 8]
#etote1 <- read.table(text = equiw)[2:(length(equiw)-1), 11]
#fssew <- read.table(text=nsew)[,2]*read.table(text=ssew)[,2]

## NOTE NOTE NOTE, those c%%7 and c%%6 are hard coded at the moment
## meaning that they depend on the specific log_write for this case. They
## are not general numbers for any fep file yet.
## MOLARIS outputs an average which is not differentiated from the other values 
## and then the grep command has more than it should. It also includes an additional report
## on energy zero and that is why the line c <- 1:length(ekinf1); ekinf1 <- ekinf1[!(c%%6==0)]
## is needed.
#ekinf1 <- read.table(text = fepw)[, 5]
#c <- 1:length(ekinf1); ekinf1 <- ekinf1[!(c%%7==0)]
#c <- 1:length(ekinf1); ekinf1 <- ekinf1[!(c%%6==0)]
#epotf1 <- read.table(text = fepw)[, 8]
#c <- 1:length(epotf1); epotf1 <- epotf1[!(c%%7==0)]
#c <- 1:length(epotf1); epotf1 <- epotf1[!(c%%6==0)]
#etotf1 <- read.table(text = fepw)[, 11]
#c <- 1:length(etotf1); etotf1 <- etotf1[!(c%%7==0)]
#c <- 1:length(etotf1); etotf1 <- etotf1[!(c%%6==0)]
#fssfw <- read.table(text=nsfw)[,2]*read.table(text=ssfw)[,2]


#ekinh2 <- read.table(text = heatp)[2:(length(heatp)-1), 5]
#epoth2 <- read.table(text = heatp)[2:(length(heatp)-1), 8]
#etoth2 <- read.table(text = heatp)[2:(length(heatp)-1), 11]
#fsshp <- read.table(text=nshp)[,2]*read.table(text=sshp)[,2]
#ekine2 <- read.table(text = equip)[2:(length(equip)-1), 5]
#epote2 <- read.table(text = equip)[2:(length(equip)-1), 8]
#etote2 <- read.table(text = equip)[2:(length(equip)-1), 11]
#fssep <- read.table(text=nsep)[,2]*read.table(text=ssep)[,2]

#ekinf2 <- read.table(text = fepp)[, 5]
#c <- 1:length(ekinf2); ekinf2 <- ekinf2[!(c%%7==0)]
#c <- 1:length(ekinf2); ekinf2 <- ekinf2[!(c%%6==0)]
#epotf2 <- read.table(text = fepp)[, 8]
#c <- 1:length(epotf2); epotf2 <- epotf2[!(c%%7==0)]
#c <- 1:length(epotf2); epotf2 <- epotf2[!(c%%6==0)]
#etotf2 <- read.table(text = fepp)[, 11]
#c <- 1:length(etotf2); etotf2 <- etotf2[!(c%%7==0)]
#c <- 1:length(etotf2); etotf2 <- etotf2[!(c%%6==0)]
#fssfp <- read.table(text=nsfp)[,2]*read.table(text=ssfp)[,2]


pdf(file="moldiag.pdf",family="Helvetica", width=10,height=6)

par(mfcol=c(2,3))
par(cex=0.6)
par(mar=c(1, 2, 0, 0), oma=c(4,4,4,2))
par(tcl=-0.25)
par(mgp=c(2,0.6,0))

ndim <- (nrow(as.matrix(etoth1))*fsshw)
xlim=range(0, ndim[1])
ylim=range(min(ekinh1),max(epoth1))
yrange=abs(max(epoth1)-min(ekinh1))
plot(seq(fsshw[1],ndim[1],fsshw[1]),etoth1, type="o", col="green",
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i", pch=18)
lines(seq(fsshw[1],ndim[1],fsshw[1]),ekinh1,type="b",col="red" ,xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
lines(seq(fsshw[1],ndim[1],fsshw[1]),epoth1,type="b",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
legend("topright",col=c("red","blue","green"),lty=1,legend=c("ekin","epot","etot"))
#axis(1, cex.axis=1.2, at=seq(0,ndim,(ndim/10)))
axis(2, at=seq(min(ekinh1),max(epoth1),yrange/10),las=1)
box(which="plot",col="black")

ndim <- nrow(as.matrix(etoth2))*fsshp
xlim=range(0, ndim)
ylim=range(min(ekinh2),max(epoth2))
yrange=abs(max(epoth2)-min(ekinh2))
plot(seq(fsshp,ndim,fsshp),etoth2, type="o", col="green",
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i", pch=18)
lines(seq(fsshp,ndim,fsshp),ekinh2,type="b",col="red" ,xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
lines(seq(fsshp,ndim,fsshp),epoth2,type="b",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
legend("topright",col=c("red","blue","green"),lty=1,legend=c("ekin","epot","etot"))
axis(1, cex.axis=1.0, at=seq(0,ndim,fsshp*10))
axis(2, cex.axis=1.0, at=seq(min(ekinh2),max(epoth2),yrange/10),las=1)
box(which="plot",col="black")




ndim <- nrow(as.matrix(etote1))*fssew
xlim=range(0, ndim)
ylim=range(min(ekinh1),max(epoth1))
yrange=abs(max(epoth1)-min(ekinh1))
plot(seq(fssew,ndim,fssew),etote1, type="o", col="green",
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i", pch=18)
lines(seq(fssew,ndim,fssew),ekine1,type="b",col="red" ,xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
lines(seq(fssew,ndim,fssew),epote1,type="b",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
legend("topright",col=c("red","blue","green"),lty=1,legend=c("ekin","epot","etot"), cex=1.0)
#axis(1, cex.axis=1.2, at=seq(0,ndim,(ndim/10)))
#axis(2, cex.axis=1.0, at=seq(min(ekinh1),max(epoth1),yrange/10),las=1)
box(which="plot",col="black")

ndim <- nrow(as.matrix(etote2))*fssep
xlim=range(0, ndim)
ylim=range(min(ekinh2),max(epoth2))
yrange=abs(max(epoth2)-min(ekinh2))
plot(seq(fssep,ndim,fssep),etote2, type="o", col="green",
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i", pch=18)
lines(seq(fssep,ndim,fssep),ekine2,type="b",col="red" ,xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
lines(seq(fssep,ndim,fssep),epote2,type="b",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
legend("topright",col=c("red","blue","green"),lty=1,legend=c("ekin","epot","etot"), cex=1.0)
axis(1, cex.axis=1.0, at=seq(0,ndim,fssep*10))
#axis(2, cex.axis=1.0, at=seq(min(ekinh2),max(epoth2),yrange/10),las=1)
box(which="plot",col="black")




ndim <- nrow(as.matrix(etotf1))*fssfw
xlim=range(0, ndim)
ylim=range(min(ekinh1),max(epoth1))
yrange=abs(max(epoth1)-min(ekinh1))
plot(seq(fssfw,ndim,fssfw),etotf1, type="o", col="green",
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i", pch=18)
lines(seq(fssfw,ndim,fssfw),ekinf1,type="b",col="red" ,xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
lines(seq(fssfw,ndim,fssfw),epotf1,type="b",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
legend("topright",col=c("red","blue","green"),lty=1,legend=c("ekin","epot","etot"), cex=1.0)
#axis(1, cex.axis=1.2, at=seq(0,ndim,(ndim/10)))
#axis(2, cex.axis=1.0, at=seq(min(ekinh1),max(epoth1),yrange/10),las=1)
box(which="plot",col="black")

ndim <- nrow(as.matrix(etotf2))*fssfp
xlim=range(0, ndim)
ylim=range(min(ekinh2),max(epoth2))
yrange=abs(max(epoth2)-min(ekinh2))
plot(seq(fssfp,ndim,fssfp),etotf2, type="o", col="green",
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i", pch=18)
lines(seq(fssfp,ndim,fssfp),ekinf2,type="b",col="red" ,xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
lines(seq(fssfp,ndim,fssfp),epotf2,type="b",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
legend("topright",col=c("red","blue","green"),lty=1,legend=c("ekin","epot","etot"), cex=1.0)
axis(1, cex.axis=1.0, at=seq(0,ndim,fssfp*10))
#axis(2, cex.axis=1.0, at=seq(min(ekinh2),max(epoth2),yrange/10),las=1)
box(which="plot",col="black")

par(las=0)
mtext("Water", side=3, outer=T, cex=1.0, line=1.6)
mtext("Heat", side=3, outer=T, cex=0.6, line=0.6, adj=0.16)
mtext("Equilibrate", side=3, outer=T, cex=0.6, line=0.6, adj=0.5)
mtext("FEP-US", side=3, outer=T, cex=0.6, line=0.6, adj=0.86)
mtext("Protein", side=1, outer=T, cex=1.0, line=3.0)
mtext("Simulation time in ps", side=1, outer=T, cex=0.8, line=1.0)
dev.off()
