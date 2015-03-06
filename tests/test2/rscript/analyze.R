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


simtime         <- grep("Number of MD steps ",readLines("../step1/mech1/protein/equilibration/dyn_min_wat_noshake.log"), value=TRUE)
solute          <- grep("^solute ",readLines("../step1/mech1/protein/equilibration/dyn_min_wat_noshake.log"), value=TRUE)
solvent         <- grep("^solvent ",readLines("../step1/mech1/protein/equilibration/dyn_min_wat_noshake.log"), value=TRUE)
solutesolvent   <- grep("^solute-solvent ",readLines("../step1/mech1/protein/equilibration/dyn_min_wat_noshake.log"), value=TRUE)
lrf             <- grep("^LRF ",readLines("../step1/mech1/protein/equilibration/dyn_min_wat_noshake.log"), value=TRUE)
qatom           <- grep("^Q-atom ",readLines("../step1/mech1/protein/equilibration/dyn_min_wat_noshake.log"), value=TRUE)
restraints      <- grep("^restraints ",readLines("../step1/mech1/protein/equilibration/dyn_min_wat_noshake.log"), value=TRUE)
total           <- grep("^SUM ",readLines("../step1/mech1/protein/equilibration/dyn_min_wat_noshake.log"), value=TRUE)
temperature     <- grep("Temperature at step",readLines("../step1/mech1/protein/equilibration/dyn_min_wat_noshake.log"), value=TRUE)

simtime       <- read.table(text = simtime)[1:(length(simtime)-1), 2:10]
simtime       <- (simtime[9]*simtime[5])/1000
solute        <- read.table(text = solute)[1:(length(solute)-1), 2:7]
solvent       <- read.table(text = solvent)[1:(length(solvent)-1), 2:7]
solutesolvent <- read.table(text = solutesolvent)[1:(length(solutesolvent)-1), 2:3]
lrf           <- read.table(text = lrf)[1:(length(lrf)-1), 2]
qatom         <- read.table(text = qatom)[1:(length(qatom)-1), 2:7]
restraints    <- read.table(text = restraints)[1:(length(restraints)-1), 2:7]
total         <- read.table(text = total)[1:(length(total)-1), 2:4]
temperature   <- read.table(text = temperature)[1:(length(temperature)-1), 2:8]


pdf(file="qanalyze.pdf",family="Helvetica", width=10,height=6)
par(mfcol=c(1,3))
par(cex=0.6)
par(mar=c(1, 2, 0, 0), oma=c(4,4,4,2))
par(tcl=-0.25)
par(mgp=c(2,0.6,0))

xlim=range(0, simtime[1,])
ylim=range(min(solvent[,1]),max(solvent[,1]))
yrange=abs(max(solvent[,1])-min(solvent[,1]))
plot(seq(0,simtime[1,],simtime[1,]/49),solute[,1], type="o", col="green",
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i", pch=18)
lines(seq(0,simtime[1,],simtime[1,]/49),solvent[,1],type="b",col="red" ,xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
lines(seq(0,simtime[1,],simtime[1,]/49),solutesolvent[,1],type="b",col="blue",xaxt="n",yaxt="n",xlab="",ylab="", pch=18)
legend("topright",col=c("red","blue","green"),lty=1,legend=c("solute-solvent","solvent","solute"))
axis(1, cex.axis=1.2, at=seq(0,simtime[1,],(simtime[1,]/10)))
axis(2, at=seq(min(solvent[,1]),max(solvent[,1]),yrange/10),las=1)
box(which="plot",col="black")

dev.off()

