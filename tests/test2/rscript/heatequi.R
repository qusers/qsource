#2345678901234567890123456789012345678901234567890123456789012345678901234567890
################################################################################
# Property plotter for Q runs
################################################################################
# An example of the "juice" in a Q output file is in the following summary lines:

# =======================  Energy summary at step  24500 ========================
#   el       vdW      bond     angle   torsion  improper
# solute            -7230.75   -235.91    912.03   1328.94    531.60     64.47
# solvent           -8127.62   1396.65      0.00      0.00      0.00      0.00
# solute-solvent    -2938.63   -151.95
# LRF                -231.30
# Q-atom            -1275.33     30.90   -241.13      8.34      2.86      0.00
# 
# total       fix slvnt_rad slvnt_pol     shell    solute
# restraints          631.26    113.34   -385.36     25.87     67.21    810.21
# 
# total potential   kinetic
# SUM              -10847.32 -15525.57   4678.25
# ======================= Q-atom energies at step  24500 ========================
#   type   st lambda        el       vdW      bond     angle   torsion  improper
# Q-Q     1 1.0000    329.24     15.16
# Q-Q     2 0.0000    269.55     60.38
# 
# Q-prot  1 1.0000  -1462.79      5.45
# Q-prot  2 0.0000  -1486.59      1.64
# 
# Q-wat   1 1.0000   -141.78     10.29
# Q-wat   2 0.0000   -143.81     10.35
# 
# Q-surr. 1 1.0000  -1604.57     15.75
# Q-surr. 2 0.0000  -1630.40     11.99
# 
# Q-any   1 1.0000  -1275.33     30.90   -241.13      8.34      2.86      0.00
# Q-any   2 0.0000  -1360.85     72.37    -29.70    159.91      3.82      0.00
# 
# type   st lambda     total restraint
# Q-SUM   1 1.0000  -1474.36      0.00
# Q-SUM   2 0.0000  -1154.44      0.00
# H( 1, 2) =    0.00 dist. between Q-atoms   2  31 =    1.68
# H( 1, 2) =    0.00 dist. between Q-atoms  23  31 =    2.99
# H( 1, 2) =    0.00 dist. between Q-atoms  23  26 =    1.00
# H( 1, 2) =    0.00 dist. between Q-atoms  26  33 =    2.75
# ===============================================================================
#   Temperature at step   24500:         T_tot=     295.2         T_free=     307.4

qlog <- "../start.log"
simtime         <- grep("Number of MD steps ",readLines(qlog), value=TRUE)
solute          <- grep("^solute ",readLines(qlog), value=TRUE)
solvent         <- grep("^solvent ",readLines(qlog), value=TRUE)
solutesolvent   <- grep("^solute-solvent ",readLines(qlog), value=TRUE)
lrf             <- grep("^LRF ",readLines(qlog), value=TRUE)
qatom           <- grep("^Q-atom ",readLines(qlog), value=TRUE)
restraints      <- grep("^restraints ",readLines(qlog), value=TRUE)
total           <- grep("^SUM ",readLines(qlog), value=TRUE)
temperature     <- grep("Temperature at step",readLines(qlog), value=TRUE)

simtime       <- read.table(text = simtime)[1:(length(simtime)-1), 2:10]
simtime       <- (simtime[9]*simtime[5])/1000  #simulation time in picoseconds
solute        <- read.table(text = solute)[1:(length(solute)-1), 2:7]
solvent       <- read.table(text = solvent)[1:(length(solvent)-1), 2:7]
solutesolvent <- read.table(text = solutesolvent)[1:(length(solutesolvent)-1), 2:3]
lrf           <- read.table(text = lrf)[1:(length(lrf)-1), 2]
qatom         <- read.table(text = qatom)[1:(length(qatom)-1), 2:7]
restraints    <- read.table(text = restraints)[1:(length(restraints)-1), 2:7]
total         <- read.table(text = total)[1:(length(total)-1), 2:4]
temperature   <- read.table(text = temperature)[1:(length(temperature)-1), 2:8]

qlog <- "../heat1.log"
simtime1         <- grep("Number of MD steps ",readLines(qlog), value=TRUE)
solute1          <- grep("^solute ",readLines(qlog), value=TRUE)
solvent1         <- grep("^solvent ",readLines(qlog), value=TRUE)
solutesolvent1   <- grep("^solute-solvent ",readLines(qlog), value=TRUE)
lrf1             <- grep("^LRF ",readLines(qlog), value=TRUE)
qatom1           <- grep("^Q-atom ",readLines(qlog), value=TRUE)
restraints1      <- grep("^restraints ",readLines(qlog), value=TRUE)
total1           <- grep("^SUM ",readLines(qlog), value=TRUE)
temperature1     <- grep("Temperature at step",readLines(qlog), value=TRUE)

simtime1       <- read.table(text = simtime1)[1:(length(simtime1)-1), 2:10]
simtime1       <- (simtime1[9]*simtime1[5])/1000  #simulation time in picoseconds
solute1        <- read.table(text = solute1)[1:(length(solute1)-1), 2:7]
solvent1       <- read.table(text = solvent1)[1:(length(solvent1)-1), 2:7]
solutesolvent1 <- read.table(text = solutesolvent1)[1:(length(solutesolvent1)-1), 2:3]
lrf1           <- read.table(text = lrf1)[1:(length(lrf1)-1), 2]
qatom1         <- read.table(text = qatom1)[1:(length(qatom1)-1), 2:7]
restraints1    <- read.table(text = restraints1)[1:(length(restraints1)-1), 2:7]
total1         <- read.table(text = total1)[1:(length(total1)-1), 2:4]
temperature1   <- read.table(text = temperature1)[1:(length(temperature1)-1), 2:8]

qlog <- "../heat2.log"
simtime2         <- grep("Number of MD steps ",readLines(qlog), value=TRUE)
solute2          <- grep("^solute ",readLines(qlog), value=TRUE)
solvent2         <- grep("^solvent ",readLines(qlog), value=TRUE)
solutesolvent2   <- grep("^solute-solvent ",readLines(qlog), value=TRUE)
lrf2             <- grep("^LRF ",readLines(qlog), value=TRUE)
qatom2           <- grep("^Q-atom ",readLines(qlog), value=TRUE)
restraints2      <- grep("^restraints ",readLines(qlog), value=TRUE)
total2           <- grep("^SUM ",readLines(qlog), value=TRUE)
temperature2     <- grep("Temperature at step",readLines(qlog), value=TRUE)

simtime2       <- read.table(text = simtime2)[1:(length(simtime2)-1), 2:10]
simtime2       <- (simtime2[9]*simtime2[5])/1000  #simulation time in picoseconds
solute2        <- read.table(text = solute2)[1:(length(solute2)-1), 2:7]
solvent2       <- read.table(text = solvent2)[1:(length(solvent2)-1), 2:7]
solutesolvent2 <- read.table(text = solutesolvent2)[1:(length(solutesolvent2)-1), 2:3]
lrf2           <- read.table(text = lrf2)[1:(length(lrf2)-1), 2]
qatom2         <- read.table(text = qatom2)[1:(length(qatom2)-1), 2:7]
restraints2    <- read.table(text = restraints2)[1:(length(restraints2)-1), 2:7]
total2         <- read.table(text = total2)[1:(length(total2)-1), 2:4]
temperature2   <- read.table(text = temperature2)[1:(length(temperature2)-1), 2:8]


pdf(file="qanalyze.pdf",family="Helvetica", width=10,height=6)
par(mfrow=c(2,3))
par(cex=0.6)
par(mar=c(1, 2, 0, 0), oma=c(4,4,4,2))
par(tcl=-0.25)
par(mgp=c(2,0.6,0))

xlim=range(0, simtime[1,])
ylim=range(min(solvent[,1]),max(solvent[,1]))
yrange=abs(max(solvent[,1])-min(solvent[,1]))
plot(seq(0,simtime[1,],simtime[1,]/19),solute[,1], type="o", col="green",
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i", pch=18)
lines(seq(0,simtime[1,],simtime[1,]/19),solvent[,1],type="b",col="red" ,xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
lines(seq(0,simtime[1,],simtime[1,]/19),solutesolvent[,1],type="b",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
legend("topright",col=c("red","blue","green"),lty=1,legend=c("solute-solvent","solvent","solute"))
axis(1, cex.axis=1.2, at=seq(0,simtime[1,],(simtime[1,]/10)))
axis(2, at=seq(min(solvent[,1]),max(solvent[,1]),yrange/10),las=1)
box(which="plot",col="black")


xlim=range(0, simtime[1,])
ylim=range(min(solvent[,2]),max(solvent[,2]))
yrange=abs(max(solvent[,2])-min(solvent[,2]))
plot(seq(0,simtime[1,],simtime[1,]/19),solute[,2], type="o", col="green",
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i", pch=18)
lines(seq(0,simtime[1,],simtime[1,]/19),solvent[,2],type="b",col="red" ,xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
lines(seq(0,simtime[1,],simtime[1,]/19),solutesolvent[,2],type="b",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
legend("topright",col=c("red","blue","green"),lty=1,legend=c("solute-solvent","solvent","solute"))
axis(1, cex.axis=1.2, at=seq(0,simtime[1,],(simtime[1,]/10)))
axis(2, at=seq(min(solvent[,2]),max(solvent[,2]),yrange/10),las=1)
box(which="plot",col="black")


alltemps <- rbind(temperature,temperature1,temperature2)
xlim=range(0, length(alltemps[,3]))
ylim=range(min(alltemps[,5]),max(alltemps[,5]))
yrange=abs(max(alltemps[,5])-min(alltemps[,5]))
plot(seq(1,length(alltemps[,3])),alltemps[,5], type="o", col="green",
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i", pch=18)
lines(seq(1,length(alltemps[,3])),alltemps[,7],type="b",col="red" ,xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
lines(seq(1,length(alltemps[,3])),alltemps[,5],type="b",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
legend("topright",col=c("red","blue","green"),lty=1,legend=c("solute-solvent","solvent","solute"))
axis(1, cex.axis=1.2, at=seq(0,simtime2[1,],(simtime2[1,]/10)))
axis(2, at=seq(min(alltemps[,5]),max(alltemps[,5]),yrange/10),las=1)
box(which="plot",col="black")


dev.off()



