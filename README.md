## Model Order Reduction for Opposition Control in a Turbulent Channel Flow
### Emily Williams - 18.335 Final Project

Proper orthogonal decomposition (POD) and dynamic mode decomposition (DMD) are both used to assess the flow in a compressible turbulent channel with and without opposition control for drag. The computational setup of the turbulent channel simulation is formulated using an information-theoretic approach. Both model order reduction techniques are implemented in Matlab and contour plots are generated to visualize the 3D modes in 2D planes. Results show that the first POD mode reconstructs the mean velocity of the flow, while the second POD mode captures the velocity fluctuations. The first DMD mode provides insight into the coherent structures over time, supporting the conclusions from POD analysis. A thorough performance analysis is conducted throughout to compare the operation count, effectiveness, and usefulness of both methods for studying turbulence.

#### Files included in this repository:
_postproc.m_ - generic postprocessing script (kind of messy, but works if you don't clear the variables and iterate through the z- or x-slices manually -- better when run in Supercloud)

_compute_POD.m_ - function for computing POD modes (takes in velocity field and number of modes to calculate and returns modes, energy, and eigenvalues)

_compute_DMD.m_ - function for computing DMD modes (takes in truncated matrices for each temporal segment and returns modes and eigenvectors)

_run_POD_test.m_ and _run_POD_test_XZ.m_ - test files for computing the modes in different planes
