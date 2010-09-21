setStage("groupComps", "Group Components", "cpExperiment")

setProtocol("density", "Density Estimation",
            representation(bw = "numeric", minfrac = "numeric",
                           minsamp = "numeric", mzwid = "numeric",
                           max = "numeric"),
            getMethodLocal("group.density", "xcmsSet"), "groupComps")

setProtocol("angle", "Spectral angle", representation(rt_window = "numeric",
                                                      dist.cutoff="numeric"),
            groupComps.angle, "groupComps")
