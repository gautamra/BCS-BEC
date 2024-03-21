# BCS-BEC

A python notebook to study the BCS-BEC crossover in superconducting junctions with Bogoliubov-de Gennes (BdG) calculations on 1D tight-binding models.

The base class `TBmodel` generates the Hamiltonian matrix for a nearest-neighbor 1D tight-binding model with position dependent hoppings $t_i$ and position-dependent on-site potential $u^0_i$. The method `iterate` iteratively solves the Bogoliubov-de Gennes equations self-consistently subject to the on-site Hubbard term $V_i$
