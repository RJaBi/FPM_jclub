from jclub.compiled import jclub_c
import numpy as np


def spacingCalcs(flowTimeForW0, xiNumList, xig, nEval, w0PhysMean, w0PhysErr):
    ncon = np.shape(flowTimeForW0)[1]
    xiNum = len(xiNumList)
    return jclub_c.spacingcalcs_c(ncon, xiNum, flowTimeForW0, xiNumList, xig, nEval, w0PhysMean,
                                  w0PhysErr)
