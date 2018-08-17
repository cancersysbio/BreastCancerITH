### Tools for simulating breast tumor heterogeneity and calculating heterogeneity metrics using both simulated and clinical data

Description
---
`heterogeneityMetrics.R` is an R script that can be used to calculate relevant tumor heterogeneity metrics using mutation calls from *multiregion sequencing (MRS)* of human tumor specimens (using the output of the VAP pipeline) and using mutation calls from simulated data (using the output of `breastSimulationCode.py`). All of the metrics except for `tHFR` take as input mutation calls from multiple regional samples from within the same timepoint. `tHFR` takes as input one set of *"pre-treatment"* mutations as well sets of mutations from multiple *"post-treatment"* regions. These heterogeneity metrics can be used to measure **genetic divergence** between regions and test for clonal evolution across time (`tHFR`).

`breastSimulationCode.py` is a Python script to simulate 3D peripherally dominated tumor growth and output multi-region sequencing data (mutation calls) via an **agent-based model**. Deme subdivision is assumed in order to model cell mixing and spatial constraint. 20 samples across all octants of the deme are recorded to profile simulated **spatial tumor heterogeneity**. This version assumes at most two drivers occur on the same lineage; the fitness of Tier1 and Tier2 driver lineages is `1+s` and `(1+s)^2`, respectively.

Requirements
---
Python packages: `numpy`,`sys`,`math`,`random`,`collections`,`sets`

Usage
---
Simulation of a typical tumor (`~10^9` cells) is computationally costly. We suggest to run this script on a high performance cluster. The memory cost is also large (about `~40G` when the final_tumor_size = `10^9` and mut_rate = `0.6`).

> To run the python script
```
$ python breastSimulationCode.py s_coef repl
e.g., python breastSimulationCode.py 0.1 2
```

Since breast tumors follow patterns of heterogeneity matching growth under strong selection, we suggest simulated with selection coefficients (s_coef up to 0.4-0.5). Deme size can also be modified within the script, as can the number of public mutations and mean sequencing depth for the outputted mutations. (edited)

Contact
--
> Katherine Lee Pogrebniak: 
> kpogrebn@stanford.edu
