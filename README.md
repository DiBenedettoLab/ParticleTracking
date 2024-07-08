# ParticleTracking
Michelle's 2D particle tracking code in MATLAB which is an updated version of Nick Ouellette's PTV code. Original [source code](https://web.stanford.edu/~nto/software.shtml) here.

## Particle Finding
First use `Particle_Finder.m` to identify particles in your images. The output should be a structure with centroids of identified particles.

## Particle Tracking
Next use `Particle_Tracker.m` to link identified particles into individual trajectories.

## Example
See `example.m` for updated example use case.

## Future things to add:
- More analysis, cleaning, and plotting codes
- Particle tracking in image pair data
