setwd("your_path")

##source functions for Gate design
source('Functions_gate_design.R', echo=TRUE)

##load gateID training dataset (rows= single cells & columns=ClusterID followed by index data in all available channels)
data <- as.data.frame(read.table("gateID_training_dataset_test.csv", header = T, sep = "\t", stringsAsFactors = F))

##define the desired cell type 
desired.cluster <- c(1)

##define parameter for gate design
datastart <- 2  ##index data start at column 2 since colunm 1 is ClusterID
tgc <-t(combn(seq(datastart,dim(data)[2]), m = 2)) 
combs <- t(combn(seq(datastart,dim(data)[2], 1), m=2)) 
bool.desired <- data$ClusterID %in% desired.cluster
seed <- 23
thresh <- 0.3
nverts <- 4  ##define number of vertices for gates

##create dataframe for GateID results
gate_sol <- data.frame(matrix(0, ncol = 20 , nrow = nrow(tgc)))

##the for loop runs the gate design for each gate combination (runs in 5 minutes) 
##for gate design on a full dataset, run on an hpc and submit as array job
##gate_sol contains solutions
for(j in c(1:nrow(tgc))){
  gc <- tgc[j,]
  gs <- as.vector(t(combs[gc, ]))
  dimensions <- gs
  dataset <- data
  numgates <- length(gs)/2
  d2.desired <- dataset[bool.desired,dimensions]
  d2.nondesired <- dataset[!bool.desired,dimensions]
  jumpindex <- seq(1,ncol(d2.desired), 2)
  d2.desired <- lapply(jumpindex, function(x) as.matrix(d2.desired[, c(x,(x+1))]))
  d2.nondesired <- lapply(jumpindex, function(x) as.matrix(d2.nondesired[, c(x,(x+1))]))
  l = lapply(1:length(d2.desired), function(x) t(apply(d2.desired[[x]], 2, function(y) quantile(y, probs = c(0.01, 0.99)))))
  quart <- as.matrix(l[[1]])
  for (i in 2:length(l)) {quart <- rbind(quart, l[[i]])}
  initials <- lapply(l, function(y) orderCW(as.matrix(expand.grid(y[1,], y[2,])), 4))
  initials <- unlist(lapply(initials, function(x) as.numeric(x)))
  numpoints.thresh <- round(dim(d2.desired[[1]])[1] * thresh)
  nd <- length(initials)/numgates
  ji <- seq(1,length(initials),length(initials)/numgates)
  lb <-  as.numeric(unlist(lapply(quart[,1], function(y) rep(y,4))))
  ub <-  as.numeric(unlist(lapply(quart[,2], function(y) rep(y,4))))
  lb2 <-  rep(0, length(lb))
  
  ##malsChains-based optimization
  initialpop <- initials
  cmares <- malschains(li, lower = lb, upper = ub, maxEvals = 20000, seed = seed, initialpop = initialpop, verbosity = 0,
                      control = malschains.control(optimum = 0, istep = 100, ls = "sw", effort = 0.55, alpha = 0.5, popsize = 50))
  impurity <- cmares$fitness
  initialpop <- cmares$sol
  cmares2  <- malschains(my, lower = lb2, upper = ub, maxEvals = 10000, initialpop = initialpop, seed = seed, verbosity = 0,
                       control = malschains.control(istep = 300, ls = "cmaes", effort = 0.55, alpha = 0.5, popsize = 50))
  yield.cma <- abs(cmares2$fitness)
  temp = sum(bool.desired)*yield.cma
  impurity.cma <- temp/(impurity+temp)
  spec = c(gc, yield.cma, impurity.cma, cmares2$sol)
  
  gate_sol[j,] <- spec
}

##order gates based on yield (column 3) and purity (column 4)
gate_sol <- gate_sol[rev(order(gate_sol$X4, gate_sol$X3)), ]

##picks the vertex coordinates (16 values) for the best gate combination
##change the rowgate parameter to visualize different gates solutions
rowgate <- 1
bestgate <- as.numeric(gate_sol[rowgate, c(5:ncol(gate_sol))]) 

###this tells you the start of the coordinates of each gate
numgates <- 2
ji <- seq(1,length(bestgate),length(bestgate)/numgates)

###gate.ori = list of both gate coordinates
nd <- length(bestgate)/numgates
gate.ori <- list()
gate.ori <- lapply(ji, function(y) as.matrix(orderCW(bestgate[y:(y+(nd-1))],4)))
gate.ori <- lapply(gate.ori, function(y) as.matrix(rbind(y, y[1, ])))

####dims = channels in which the gates are designed (refers to the columns of the training which contain clusterID)
h.index <- as.numeric(gate_sol[rowgate, c(1,2)])
channel.comb <- t(combn(seq(datastart,dim(data)[2]), m = 2))
dims <- channel.comb[h.index, ]

# Plot gate and its parameters
cal.gate.efficacy(trainingdata = data, bool.desired = bool.desired, gates = gate.ori, dims = dims, plot = T, verbose = T)
