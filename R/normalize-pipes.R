setStage("normalize", "Normalize Quantities", "cpExperiment")

setProtocol("scale", "Scale",
            fun = normalize.scale, parent = "normalize")
