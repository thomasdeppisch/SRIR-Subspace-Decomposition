# Direct and Residual Subspace Decomposition for Spatial Room Impulse Responses

This repository contains a MATLAB implementation of the direct and residual subspace decomposition for spatial room impulse responses (SRIRs). The method decomposes SRIRs into a direct part, comprising direct sound and salient reflections, and a residual.

For more information and if you want to reference the code [please refer to the following publication](https://ieeexplore.ieee.org/document/10028731)
   
   ```
   T. Deppisch, S. Amengual Garí, P. Calamia, and J. Ahrens, 
   "Direct and Residual Subspace Decomposition of Spatial Room Impulse Responses," 
   IEEE/ACM Transactions on Audio, Speech, and Language Processing, 2023, doi:10.1109/TASLP.2023.3240657.
   ```
   
A [website with audio examples](http://www.ta.chalmers.se/srir-subspace-decomposition/) is available.   
   
## Contents
The function in `srirSubspaceDecomp.m` contains the implementation of the subspace decomposition method that is described in the paper. The script `testSrirSubspaceDecomp.m` shows how to apply the method to an SRIR. The Live Script `settingParameters.mlx` and the corresponding PDF `settingParameters.pdf` illustrate the influence of the different parameters. The SRIR that is used in the examples was measured with the Eigenmike em32 in the teaching hall at the Division of Applied Acoustics at Chalmers University of Technology and is provided as fourth-order SRIR in the spherical harmonics domain.

## Acknowledgment
We thank Meta Reality Labs Research for funding this project.

## License
This software is licensed under a Non-Commercial Software License (see [LICENSE](https://github.com/thomasdeppisch/SRIR-Subspace-Decomposition/blob/master/LICENSE) for full details).
