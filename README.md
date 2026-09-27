# SymmetryBreakingBifurcations
Code associated with the paper "Validation of Symmetry Induced Bifurcation Points in Nonlinear Diffusion Problems" by Maxime Breden, Evelyn Sander and Thomas Wanner.

You can reproduce the proofs from the paper by run script_OK.m and script_SKT.m. In each script, you can modify line 17 to change the solution considered, and line 23 if you deal with a seemingly transcritical bifurcation and want to prove it.

In order to get a fully rigorous proof, you need the Intlab toolbox (http://www.ti3.tu-harburg.de/intlab/). The code still runs without it, but then rounding errors are not controlled.
