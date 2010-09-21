setStage("findComps", "Find Components", "cpExperiment")

setProtocol("sigma", "Overlapping sigma", fun = findComps.sigma,
            parent = "findComps")

## setProtocol("sigma_filt", "Overlapping sigma with filter", fun = findComps.sigma_filt,
##             parent = "findComps")
