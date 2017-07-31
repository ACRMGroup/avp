AVP
===

AVP ('Another Void Program') unifies the calculation of individual
voids with the calculation of packing quality. Voids are defined as
holes in the protein that are not accessible to solvent, but into
which a molecule of a given raidus (such as a water) can fit. It does
this by having two probe spheres - one that defines the solvent and
one that defines the spaces in the protein.

The program works by creating a grid around the protein. It then fills
empty grid cells with water by working from all 6 surfaces. This is
then iterated in each direction ensuring that a water is placed next
to an existing water to access grid points that might be hidden in the
first simple iteration. Typically a 1.4A probe is used to check for
absence of clashes with the protein. Initially a course grid is used
(default 1A), but close to the protein a finer grid (default 0.1A) is
used and off-grid locations are explored.

At this stage, any grid point that is not occupied by protein or water
is a void point. A second probe is then used to define these are part
of the void space. If the second probe has a zero-radius then it is a
true measure of packing quality. In practice it is recommended to use
a 0.5A probe to avoid channels running through the centre of an
alpha-helix.  Again, a finer grid is used at the boundary with the
protein.

Finally grid points assigned as void are clustered to identify
discrete voids.


