setwd("~/Desktop/PhD_Chloe_AvO/GateKeeper/Shareable_code")
source('Functions_GateID.R', echo=TRUE)


#_________PART 1: Gate design_________#

##load gateID training dataset (rows= single cells & columns=ClusterID followed by index data in all available channels)
data <- as.data.frame(read.table("GateID_training_dataset.csv", header = T, sep = "\t", stringsAsFactors = F))
datastart <- 2  #index data start at column 2 since colunm 1 is ClusterID

##define the desired cell type 
desired.cluster <- 1
tgc <-t(combn(seq(datastart,dim(data)[2]), m = 2)) 
bool.desired <- data$ClusterID %in% desired.cluster

##define parameters for gate design
seed <- 23
thresh <- 0.3 #yield threshold
nverts <- 4 #define number of vertices for gates

##create dataframe for GateID results
gate_sol <- data.frame(matrix(0, ncol = 20 , nrow = nrow(tgc)))

##the for loop runs the gate design for each gate combination (runs in 5 minutes) 
##for gate design on a full dataset, run on an hpc and submit as array job
##gate_sol contains solutions
for(j in c(1:nrow(tgc))){
  gc <- tgc[j,]
  gs <- as.vector(t(tgc[gc, ]))
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
rowgate <- 1 #change the rowgate parameter to visualize different gates solutions
bestgate <- as.numeric(gate_sol[rowgate, c(5:ncol(gate_sol))]) 
numgates <- 2
ji <- seq(1,length(bestgate),length(bestgate)/numgates)
nd <- length(bestgate)/numgates
gate.ori <- list()
gate.ori <- lapply(ji, function(y) as.matrix(orderCW(bestgate[y:(y+(nd-1))],4)))
gate.ori <- lapply(gate.ori, function(y) as.matrix(rbind(y, y[1, ]))) #contains predicted gate coordinates
h.index <- as.numeric(gate_sol[rowgate, c(1,2)])
channel.comb <- t(combn(seq(datastart,dim(data)[2]), m = 2))
dims <- channel.comb[h.index, ] #channels in which the gates are designed (refers to the columns of the training which contain clusterID)
cal.gate.efficacy(trainingdata = data, bool.desired = bool.desired, gates = gate.ori, dims = dims, plot = T, verbose = T) #plot gate and print gate stats


#_________PART 2: Train a machine learning model to perform desired cell type prediction_________#

##create train and test dataset from data
type.clusters.desired <- {list(
  c("Desired",desired.cluster),
  c("Other",unique(data[!data[,1] %in% desired.cluster,1]))
)}
oricl.desired <- rename.clusters(data = data, type.clusters = type.clusters.desired,name = T)

train <- (createDataPartition(oricl.desired, p = 0.8))$Resample1
traindata <- data[train, ]
trainclass <- oricl.desired[train]

testdata <- data[-train,]
testclass <- oricl.desired[-train]

set.seed(seed)

##train the random forrest model
control <- trainControl(method = "repeatedcv", number = 10, repeats = 10, sampling = "smote", classProbs = TRUE, summaryFunction = multiClassSummary)
final.model <- train(x = traindata[,-1], y = trainclass, method="rf", tuneLength = 5, trControl=control, preProcess = c("center","scale"), metric = "Kappa")
predrf <- predict(final.model, testdata[,-1]) #evaluate model on test data
p <- table(model = predrf, truth = testclass)[1,]
p[1]/(sum(as.matrix(p))) #print model prediction accuracy


#_________PART 4: Load FACS data associated with training dataset_________#

odata.fcs <- read.table("FACS_data_training_dataset.csv", header = T, sep = "\t", stringsAsFactors = F)
##for your own training dataset, this data will have to be extracted from the fcs file using the extract.presort functions as below


#_________PART 5: Load current FACS data (of new experimental dataset) and perform normalization_________#

xmlfile <- "FACS_configuration_experimental_dataset.xml" #xml of new experimental dataset
fcsfile <- "FACS_data_experimental_dataset.fcs" #fcs file of new experimental dataset (record 10000 events)
presort <- extract.presort(fcsfile = fcsfile, xmlfile = xmlfile, pregates = c("P1", "P2", "P3", "P4")) #extract non-normalized FACS parameters of new experimental dataset (pregates are gate names used to gate live cells)
presort <- presort[,c(2,3,20,5)] #selecting only FACS dimensions for the test dataset - keep all FACS parameters for real experiment
colnames(presort) <- colnames(odata.fcs) #unify colnames for all dataframes (data/odata.fcs/presort)

presort.q <- presort #create presort.q that will be the normalized new experimental dataset
for (i in 1:ncol(odata.fcs)) {
  presort.q[,i] = normalize.qspline(x = as.matrix(presort.q[,i]), target = as.matrix(odata.fcs[,i]), verbose = F, min.offset = 0.01)
}
l = unique(unlist(apply(presort.q, 2, function(x) which(x<0))))
if (length(l) > 0) {
  presort <- presort[which(!1:nrow(presort.q) %in% l),]
  presort.q <- presort.q[which(!1:nrow(presort.q) %in% l),]
}
pred <- predict(final.model, presort.q) #perform the normalization
gate.norm <- norm.class(trainingdata = data[,-1], newdata = presort, bool.desired = oricl.desired == "Desired" , bool.ml = pred == "Desired", gates = gate.ori, dims = dims-1) #normalize the gates
print(gate.norm)


#_________PART 6: Plot gates_________#

##plot new facs data before and after normalization
plot.norm.gate(prenorm = presort, postnorm = presort.q, pregates = gate.ori, postgates = gate.norm, dims = dims)
