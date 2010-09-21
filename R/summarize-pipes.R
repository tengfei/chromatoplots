
setStage("summarize", "Summarize Quantities", "cpExperiment")

setProtocol("common", "Common Peaks (by m/z)",
            fun = summarize.common, parent = "summarize")
