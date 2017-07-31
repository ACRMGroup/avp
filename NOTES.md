Other programs
==============

CAVENV from the CCP4 suite is grid based:

CAVENV produces a map from an input model structure, which is designed
to help visualise cavities in the protein. By default, the program
produces a "cavity" map in which each grid point is given a value
equal to the distance of the closest atom, minus its Van der Waals
radius, up to a maximum of <maxrad> (maximal probe radius to be
tested). Thus the grid values are zero within the Van der Waals
envelope of the protein. Outside this envelope, the value represents
radius of the largest probe that access that grid point. In the middle
of solvent regions, the grid value is that of the largest probe
considered <maxrad>.  1.A. Volbeda, private communication or with
reference (in french): Anne Volbeda, Speleologie des hydrogenases a
nickel et a fer. In: "Les Ecoles Physique et Chimie du Vivant, numero
1 - avril 1999, Analyse de l'organisation tridimensionnelle des
proteines", pp 47-52.  ANNE VOLBEDA, IBS/LCCP GRENOBLE CCP4 version -
Martyn Winn


The new plan..... 17.07.02

For all boundary solvent voxels, calculate actual water pseudo-atom
positions rather than assigning the whole voxel as a no-go for voxel
refinement while it doesn't prevent extra possible void points being
found.

The plan... 16.07.02

```
Iterate
{
   Iterate
   {
      Foreach protein point (a) next to a solvent point
      {
         Foreach neighbouring protein point (b)
         {
            If a solvent fits anywhere between (a) and (b)
            {
               Change (a) and (b) to solvent
            }
         }
      }
   }

   DoSolventSearch()
}
```
