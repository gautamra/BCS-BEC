# BCS-BEC

A python notebook to study the BCS-BEC crossover in superconducting junctions with Bogoliubov-de Gennes (BdG) calculations on 1D tight-binding models. The workhorse for the tight-binding simulations is package `kwant`.

The base class `TBmodel` generates the Hamiltonian matrix for a nearest-neighbor 1D tight-binding model with position dependent hoppings $t_i$ and position-dependent on-site potential $u^0_i$. The method `iterate` iteratively solves the Bogoliubov-de Gennes equations self-consistently subject to the on-site Hubbard term $V_i$, corresponding to the Hamiltonian

$$\hat{H}=\sum_{i=0}^{L-1}\sum_{\sigma} (-t c^\dagger_{i,\sigma} c_{(i+1),\sigma} + h.c.) + (u^0_i - \mu)c^\dagger_{i,\sigma}c_{i,\sigma} + \frac{V_i}{2} c^\dagger_{i,\sigma}c^\dagger_{i,\bar{\sigma}}c_{i,\bar{\sigma}}c_{i,\sigma},$$

where $\sigma$ is a spin index, $i$ is a position index and meant modulo $L$, such that the last hopping term connects the first and last sites. 

## Quickstart

The class `SNS_junction` sets up a a symmetric superconductor-normal-superconductor junction. It can be thought of as three segments stitched together. The Hamiltonian is site-independent in each segment save at the edges. The superconducting segments on both sides are identical (length, attraction). The hopping is the same throughout the chain, save for the hopping connecting the first and last sites, which is zero, to implement open boundary conditions.

The class objects take the following required parameters:

| Parameter name | Meaning |
|------|----------|
|`LS`|Number of sites in both superconducting segments|
|`LN`|Number of sites in both normal segment|
|`V`|Attractive potential in the superconducitng segment|
|`t`|Hopping in all three segments|
|`U0S`|On-site potential in the superconducting segment|
|`U0N`|On-site potential in the normal segment|

The corresponding Hamiltonian is:

$$\begin{split}
\hat{H} =& \sum_{i=0}^{L-2}\sum_{\sigma} (-t c^\dagger_{i,\sigma} c_{(i+1),\sigma} + h.c.)\\
        &+ \sum_{i=0}^{LS-1}\sum_{\sigma} (U0S - \mu)c^\dagger_{i,\sigma}c_{i,\sigma} + \sum_{i=LS}^{LN+LS-1}\sum_{\sigma} (U0N - \mu)c^\dagger_{i,\sigma}c_{i,\sigma} + \sum_{i=LS_LN}^{L-1}\sum_{\sigma} (U0S - \mu)c^\dagger_{i,\sigma}c_{i,\sigma}\\
        &+ \sum_{i=0}^{LS-1} V c^\dagger_{i,\uparrow}c^\dagger_{i,\downarrow}c_{i,\downarrow}c_{i,\uparrow} + \sum_{i=LS+LN}^{L-1} V c^\dagger_{i,\uparrow}c^\dagger_{i,\downarrow}c_{i,\downarrow}c_{i,\uparrow}
\end{split}$$
