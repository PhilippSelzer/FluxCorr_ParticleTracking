This software consists of the programs 'GridInfo_Prisms.m', 'FVM_FluxCorr_Prisms.m', 
and 'PartTrack_Prisms.m' including functions for reading binary data files of 
HydroGeoSphere (TM), which is a propiretary software by (c) AQUANTY Inc.,
and a function for coordinate transformation.

The code 'GridInfo_Prisms.m' needs to be run once per grid. It needs to be run again,
if the extent of Dirichlet or Neumann boundary conditions are changed.

The code 'FVM_FluxCorr_Prisms.m' needs to be run once per model for all desired
time steps.

The code 'PartTrack_Prisms.m' can be run as often as liked for a specific model,
it is the actual particle tracking routine.

If some error for reading a text-file occurs, presumably the first line of the 
respective text file, which contains some meta-information needs, to be deleted
manually, and the issue will be fixed.