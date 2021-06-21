source("plotFuncs.R")
source("simFuncs.R")

# Landscape
simulateFull("core.topo", paramIndex = 2573, workers = 8, noiseLevel = 0.4, scaledNoise = T,  tMax = 75) # set workers according to your system specs
dfTraj <- read.csv("core/2573_simDat.csv")
landScape(dfTraj, "2573")

# trajectories
dfTraj <- read.csv("core/2573_simDat.csv")
trajectoriesnIC(dfTraj, "2573")
simSwitching("core.topo", 2573, scaledNoise = T, tMax = 75, workers = 8)
