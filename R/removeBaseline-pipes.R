setStage("removeBaseline", "Remove baseline noise", "cpSample")

setProtocol("median", "Median",
            representation(mzrad = "numeric", scanrad = "numeric"),
            removeBaseline.median, "removeBaseline")

setProtocol("rbe", "Robust Baseline Estimation",
            representation(span = "numeric", runs = "numeric", b = "numeric"),
            removeBaseline.rbe, "removeBaseline")
