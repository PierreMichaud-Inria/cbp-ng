#!/usr/bin/env python

import sys
import math

# reference core
IPCcbp0 = 6.51
CPIcbp0 = 0.0284
EPIcbp0 = 1100

# wrong-path instructions per correct-path instructions for the reference core
WPI0 = IPCcbp0 * CPIcbp0

# technology parameter (frequency vs voltage)
ALPHA = 1.625
BETA = 4*ALPHA / (ALPHA-1)**2
GAMMA = 2 / (ALPHA-1)

# ratio of CBP energy to core energy for the reference core 
cbp_energy_ratio = 0.05
EPI0 = EPIcbp0 / cbp_energy_ratio

if len(sys.argv) > 1:
    data = sys.argv[1]
else:
    data = sys.stdin.read().strip()

IPCcbp, CPIcbp, EPIcbp = map(float, data.split(','))

# wrong-path instructions per correct-path instruction
WPI0 = IPCcbp0 * CPIcbp0
WPI = IPCcbp * CPIcbp

speedup = (IPCcbp/IPCcbp0) * (1+WPI0)/(1+WPI)

LAMBDA = 1/(1+WPI0/2) - cbp_energy_ratio

normalizedEPI = ((EPIcbp/EPIcbp0) * cbp_energy_ratio + LAMBDA * speedup**GAMMA) * (1+WPI/2)

VFS = speedup * ALPHA * (1-2/(1+math.sqrt(1+BETA/(speedup*normalizedEPI))))
print(f"{VFS:.4f}")
