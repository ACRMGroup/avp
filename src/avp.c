/* #define DEBUG 1 */
/*************************************************************************

   Program:    avp (Another Void Program)
   File:       avp.c
   
   Version:    V1.4
   Date:       31.07.17
   Function:   Find voids in proteins
   
   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 2001-17
   Author:     Dr. Andrew C. R. Martin
   Address:    School of Animal and Microbial Sciences,
               The University of Reading,
               Whiteknights,
               P.O. Box 228,
               Reading RG6 6AJ.
               England.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

   Uses a very simple algorithm as described by the following pseudo-code

   main()
   {
      ParseCmdLine()
      blOpenStdFiles()
      ReadPDB()
      SetRadii()
      [Usage()]
      FindVoids()
      {
         FindBoundaries()
         UpdateBoundaries()
         BuildGrid()
         {
            [FreeGrid()]
         }
         MarkVoidAndProteinPoints()
         {
            AtomNear()
            {
               qsort/CompareCoords()
               FindFirstAtom()
               FindLastAtom()
            }
         }
         #ifdef FLOOD_FILL_SOLVENT
            MarkSolventPoints (FLOOD-FILL VERSION)
            {
               FloodFillSolvent()
               {
                  AtomNear()
                  SetAsSolvent()
                  FloodFillSolvent()
               }
            }
         #else
            MarkSolventPoints()
            {
               FlagBoxBoundsAsSolvent()
               AdjacentIsSolvent()
               AtomNear()
               SetAsSolvent()
            }
         #endif
         PrintGrid()
         ClusterVoids()
         {
            FloodFillVoid()
            {
               FloodFillVoid()
            }
         }
         FindNearestAtoms()
         {
            FindTheNearestAtom()
         }
         RefineVoidVolumes()
         {
            RefineVoxel()
            RefineVoxelNeighbours()
            {
               RefineAVoxelNeighbour()
               {
                  RefineVoxel()
               }
            }
         }
      }
      
      PrintVoids()
      {
      }
   }
   
   i.e. 
   1. Create a grid around the protein
   2. Mark each point on the grid as protein or void. i.e. assume any
      point not occupied by protein is void
   3. Mark all void points on the surface of the box as solvent
   4. (The slightly complex bit) Walk across the grid in all 6 directions
      changing void points to solvent if they have a neighbour that is 
      solvent and no atoms within the solvent probe size distance; also
      mark any points within the solvent radius as being solvent
   5. Use flood filling to cluster void points into individual voids
   6. Run through each void, refining its volume by splitting each voxel
      into 1000 sub-voxels and reassessing these
   7. Look at protein voxels that neighbour the void voxels too see if
      parts of these are also void.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  31.10.01 Original
   V1.1  18.06.02 Fixed bug in refining void size
                  Added marginal water refinement
                  By default, when refining voxels, look at corner and 
                  edge connections as well as planar connections
   V1.2  02.07.02 Uses constructed neighbour lists in finding solvent as
                  well as in void refinement to speed things up
                  considerably!
                  Differentiates between surface and buried voids
   V1.2a xx.07.02 Added solvent point reassignment
   V1.3  18.08.08 Added OPENMP code
   V1.4  31.07.17 Cleanup for new biolip
   V1.5  27.11.20 Merged in Gil Hoben's patches for OpenMP and compilation
                  with GCC

*************************************************************************/
/* Program options
*/
/* #define FLOOD_FILL_SOLVENT */ /* Flood fill the solvent              */
/* #define NO_EXTRA_VOXEL_NEIGHBOURS */ /* Refine voxels only along planar
                                           connections                  */

/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "bioplib/MathType.h"
#include "raybox.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 256
#define MAXATOMRAD 1.9            /* Max atom radius for grid expansion */
#define FINEGRIDSTEP 10           /* Divide grid step by this to calculate
                                     refined voids                      */
#define MAXNEIGHBOURS 80          /* Max number of atom neighbours for a
                                     grid point                         */
#define RADEXPAND 2               /* Atom radius expansion to find 
                                     neighbous                          */
#define WATER_RADIUS 1.4          /* Radius of a standard water         */

#define TYPE_VOID     (1)         /* Grid point types                   */
#define TYPE_PROTEIN  (2)
#define TYPE_PROTEIN2 (4)
#define TYPE_SOLVENT  (8)
#define TYPE_SOLVENT2 (16)
#define TYPE_ASSIGNED (32768)

#define DEFAULT_GRID_SIZE  1.0    /* Default sizes                      */
#define DEFAULT_PROBE_SIZE 0.00
#define DEFAULT_WATER_SIZE 1.4

#define SOLVENT_MULTIPLIER  2     /* Assign grid points as solvent within
                                     this factor of solvent size        */
#define GRID_EXPAND 2             /* Factor by which to expand grid     */

#define SOLVENT_NO          0     /* Solvent conditions                 */
#define SOLVENT_BULK        1
#define SOLVENT_MARGIN      2

#define NSOLVFINESTEP       19    /* Number of fine steps to take in
                                     each direction when setting solvent
                                     point at the protein boundary      */
#define NUMGRIDNEIGHBOURS   3    /* Number of neighbouring points on grid
                                    to test to see if this is a surface
                                    void                                */

#define SETLOWBYTE(x, y)    (x) = ((x)&0xFF00)|(y)

/************************************************************************/
typedef struct _pointlist
{
   struct _pointlist *next;
   PDB               *nearest;
   REAL              x,
                     y,
                     z;
   int               ix,
                     iy,
                     iz;
} POINTLIST;

typedef struct _voids
{
   struct _voids *next;
   POINTLIST     *pointlist;
   REAL          volume,
                 xmin, xmax, 
                 ymin, ymax, 
                 zmin, zmax;
   int           voxelCount,
                 ix, ixmin, ixmax, 
                 iy, iymin, iymax, 
                 iz, izmin, izmax;
   BOOL          surfaceVoid;
}  VOIDS;

typedef struct
{
   USHORT ***grid;
   PDB    *****neighbours;
   REAL   xmin, 
          ymin, 
          zmin,
          xmax, 
          ymax, 
          zmax,
          gridStep,
          voxelVolume;
   int    ixmax, 
          iymax, 
          izmax;
}  GRID;

typedef struct
{  FILE *gridFile;
   FILE *gridRefinedFile;
   FILE *neighbourFile;
   BOOL printVoids;
   BOOL printSolvent;
   BOOL printAtoms;
   BOOL printRefinedVoids;
   BOOL printWaterNeighbours;
   BOOL edgeConnected;
   BOOL cornerConnected;
   BOOL printNeighbours;
   BOOL leak;
   BOOL quiet;
   BOOL reassignVoids;
   BOOL reassignSurface;
}  FLAGS;


/************************************************************************/
/* Globals
*/
FLAGS gFlags;

int  gDepth = 0,
     gMaxDepth = 0;
REAL gOffsetX = 0.0,
     gOffsetY = 0.0,
     gOffsetZ = 0.0;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
VOIDS *FindVoids(FILE *out, PDB *pdb, REAL gridStep, REAL solvSize, 
                 REAL probeSize, BOOL doRefine);
void FindBoundaries(PDB *pdb, GRID *grid);
void UpdateBoundaries(GRID *grid, REAL expandSize);
BOOL BuildGrid(GRID *grid, REAL gridStep);
void MarkVoidAndProteinPoints(PDB *pdb, GRID *grid, REAL probeSize);
BOOL MarkSolventPoints(PDB *pdb, GRID *grid, REAL solventStep,
                       REAL solventMult, BOOL quiet, BOOL doNeighbours);
int  AdjacentIsSolvent(PDB *pdb, GRID *grid, int ix, int iy, int iz,
                       REAL solvSize, REAL xoff, REAL yoff, REAL zoff);
void FlagBoxBoundsAsSolvent(GRID *grid, int ixmax, int iymax, int izmax,
                            BOOL quiet);
BOOL AtomNear(PDB *pdb, REAL x, REAL y, REAL z, REAL probeSize,
              PDB **neighbours);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *gridStep, REAL *probeSize, REAL *solvSize,
                  BOOL *doRefine);
void Usage(void);
void PrintVoids(FILE *out, VOIDS *voids);
void SetAsSolvent(GRID *grid, int ix, int iy, int iz, 
                  int numberOfSolvatedNeighbours, REAL solvRadiusSq,
                  int type);
VOIDS *ClusterVoids(GRID *grid, int *nvoids);
int FloodFillVoid(GRID *grid, int ix, int iy, int iz,
                  VOIDS *v, int voidset);
void FindNearestAtoms(PDB *pdb, VOIDS *voids);
PDB *FindTheNearestAtom(PDB *pdb, POINTLIST *point);
void SetRadii(PDB *pdb);
void PrintGrid(FILE *fp, GRID *grid, BOOL refined);
int FindFirstAtom(PDB **array, int natoms, REAL x);
int FindLastAtom(PDB **pdbArray, int natoms, REAL x);
int CompareCoords(const void *p1, const void *p2);
void FloodFillSolvent(GRID *grid, PDB *pdb, int ix, int iy, int iz, 
                      REAL solvSize, REAL solvRadiusSq, 
                      int numberOfSolvatedNeighbours);
BOOL RefineVoidVolumes(VOIDS *voidlist, GRID *grid, PDB *pdb, 
                       REAL gridStep, REAL probeSize, REAL solvSize);
void FreeGrid(GRID *grid);
int RefineVoxel(GRID *grid, PDB *pdb, POINTLIST *pt, REAL fineGridStep,
                REAL probeSize);
int RefineVoxelNeighbours(GRID *grid, PDB *pdb, POINTLIST *pt, 
                          REAL fineGridStep, REAL probeSize);
int RefineAVoxelNeighbour(GRID *grid, PDB *pdb, POINTLIST *pt,  
                          int xoff, int yoff, int zoff, 
                          REAL fineGridStep, REAL probeSize);
void PrintWaterNeighbours(FILE *out, PDB *pdb, GRID *grid);
BOOL IsSolvent(GRID *grid, int ix, int iy, int iz);
void ReassignVoidPoints(GRID *grid, REAL probeSize);
void doReassignVoidPoint(GRID *grid, int xi, int yi, int zi, 
                         PDB **neighbours, int numNeighbours,
                         REAL probeSize);
int CombineAllNeighbours(GRID *grid, int xi, int yi, int zi, 
                         PDB **neighbours);
BOOL NeighboursAreProtein(GRID *grid, int ix, int iy, int iz);
void ChangeVoxel(GRID *grid, int ix, int iy, int iz, VEC3F pt);
BOOL ReassignSurfaceProteinVoxels(PDB *pdb, GRID *grid, REAL solvSize);
BOOL CheckAndChangeNeighbours(PDB *pdb, GRID *grid, 
                              int ix, int iy, int iz, 
                              REAL solvSize,
                              int doX, int doY, int doZ);
int CombineTwoNeighbourLists(GRID *grid, int xi, int yi, int zi, 
                             int i, int j, int k,
                             PDB **neighbours);




/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Input:     
   Output:    
   Returns:   

   Main program for void finding

   17.10.01 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE  *in  = stdin,
         *out = stdout;
   char  infile[MAXBUFF],
         outfile[MAXBUFF];
   REAL  gridStep,
         probeSize,
         solvSize;
   PDB   *pdb;
   VOIDS *voidlist;
   int   natoms;
   BOOL  doRefine = FALSE;
   
   
   infile[0] = outfile[0] = '\0';
   
   if(ParseCmdLine(argc, argv, infile, outfile, &gridStep, &probeSize,
                   &solvSize, &doRefine))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=blReadPDB(in, &natoms)))
         {
            SetRadii(pdb);
            voidlist = FindVoids(out, pdb, gridStep, solvSize, probeSize,
                                 doRefine);
            if(!gFlags.quiet)
            {
               fprintf(stderr,"Printing final void list\n");
            }
            PrintVoids(out, voidlist);
#ifdef FLOOD_FILL_SOLVENT
#ifdef DEBUG
            fprintf(stderr,"Maximum recursion stack depth was %d\n", 
                    gMaxDepth);
#endif
#endif
         }
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>VOIDS *FindVoids(FILE *out, PDB *pdb, REAL gridStep, REAL solvSize, 
                    REAL probeSize, BOOL doRefine)
   -------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Main routine for finding voids. Calls routines to create a grid,
   mark protein points, identify solvent points and cluster remaining
   points into voids.

   17.10.01 Original   By: ACRM
   06.02.02 Added output file parameter and call to PrintWaterNeighbours()
*/
VOIDS *FindVoids(FILE *out, PDB *pdb, REAL gridStep, REAL solvSize, 
                 REAL probeSize, BOOL doRefine)
{
   VOIDS *voids;
   GRID  grid;
   int   nvoids,
         solventCycle = 0;
   BOOL  solventModified;
   

   if(!gFlags.quiet)
   {
      fprintf(stderr,"Finding protein box boundaries\n");
   }
   
   FindBoundaries(pdb, &grid);
   UpdateBoundaries(&grid,
                    GRID_EXPAND * (MAX(gridStep, solvSize) + MAXATOMRAD));

   if(!gFlags.quiet)
   {
      fprintf(stderr,"Creating grid\n");
   }
   BuildGrid(&grid, gridStep);

   if(!gFlags.quiet)
   {
      fprintf(stderr,"Marking protein points on grid\n");
   }
   MarkVoidAndProteinPoints(pdb, &grid, probeSize);

   if(!gFlags.quiet)
   {
      fprintf(stderr,"Marking solvent points on grid\n");
   }

   MarkSolventPoints(pdb, &grid, solvSize, SOLVENT_MULTIPLIER, 
                     FALSE, TRUE);
   if(gFlags.reassignSurface)
   {
      do
      {
         if(!gFlags.quiet)
         {
            fprintf(stderr,"Surface reassignment iteration %d\n", 
                    ++solventCycle);
         }
         
         solventModified = ReassignSurfaceProteinVoxels(pdb, &grid, 
                                                        solvSize);
         if(solventModified)
         {
            MarkSolventPoints(pdb, &grid, solvSize, SOLVENT_MULTIPLIER, 
                              FALSE, FALSE);
         }
      }  while(solventModified);
   }
   
   if(gFlags.reassignVoids)
   {
      if(!gFlags.quiet)
      {
         fprintf(stderr,"\nReassigning void points on grid\n");
      }
      ReassignVoidPoints(&grid, probeSize);
   }

   if(gFlags.printSolvent || gFlags.printAtoms)
   {
      if(!gFlags.quiet)
      {
         fprintf(stderr,"Printing grid\n");
      }
      PrintGrid(gFlags.gridFile, &grid, FALSE);
   }

   if(gFlags.printWaterNeighbours)
   {
      PrintWaterNeighbours(out, pdb, &grid);
   }

   if(!gFlags.quiet)
   {
      fprintf(stderr,"Clustering void points into distinct voids\n");
   }
   voids = ClusterVoids(&grid, &nvoids);

   if(!gFlags.quiet)
   {
      fprintf(stderr,"Found %d voids\n", nvoids);
   }
   FindNearestAtoms(pdb, voids);

   if(doRefine)
   {
      if(!gFlags.quiet)
      {
         fprintf(stderr,"Refining void volumes\n");
      }
      RefineVoidVolumes(voids, &grid, pdb, gridStep, probeSize, solvSize);
   }

   return(voids);
}


/************************************************************************/
/*>void FindBoundaries(PDB *pdb, GRID *grid)
   -----------------------------------------
   Input:     
   Output:    
   Returns:   

   Finds the boundaries of the box defined by the PDB

   17.10.01 Original   By: ACRM
*/
void FindBoundaries(PDB *pdb, GRID *grid)
{
   PDB *p;
   grid->xmin = grid->xmax = pdb->x;
   grid->ymin = grid->ymax = pdb->y;
   grid->zmin = grid->zmax = pdb->z;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->x > grid->xmax) grid->xmax = p->x;
      if(p->y > grid->ymax) grid->ymax = p->y;
      if(p->z > grid->zmax) grid->zmax = p->z;
                          
      if(p->x < grid->xmin) grid->xmin = p->x;
      if(p->y < grid->ymin) grid->ymin = p->y;
      if(p->z < grid->zmin) grid->zmin = p->z;
   }
}


/************************************************************************/
/*>void UpdateBoundaries(GRID *grid, REAL expandSize)
   --------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Update the boundaries to give space around the protein

   17.10.01 Original   By: ACRM
   19.11.01 Add offsets
*/
void UpdateBoundaries(GRID *grid, REAL expandSize)
{
   grid->xmin -= expandSize - gOffsetX;
   grid->ymin -= expandSize - gOffsetY;
   grid->zmin -= expandSize - gOffsetZ;

   grid->xmax += expandSize + gOffsetX;
   grid->ymax += expandSize + gOffsetY;
   grid->zmax += expandSize + gOffsetZ;
}


/************************************************************************/
/*>BOOL BuildGrid(GRID *grid, REAL gridStep)
   -----------------------------------------
   Input:     
   Output:    
   Returns:   

   Create the actual grid

   17.10.01 Original   By: ACRM
*/
BOOL BuildGrid(GRID *grid, REAL gridStep)
{
   int xsize, ysize, zsize, i, j, k;
   
   grid->grid = NULL;
   grid->neighbours = NULL;
   
   grid->gridStep = gridStep;
   grid->voxelVolume = gridStep * gridStep * gridStep;
   
   grid->ixmax = 1 + ((grid->xmax - grid->xmin) / grid->gridStep);
   grid->iymax = 1 + ((grid->ymax - grid->ymin) / grid->gridStep);
   grid->izmax = 1 + ((grid->zmax - grid->zmin) / grid->gridStep);
   
   xsize = 2 + ((grid->xmax-grid->xmin) / grid->gridStep);
   ysize = 2 + ((grid->ymax-grid->ymin) / grid->gridStep);
   zsize = 2 + ((grid->zmax-grid->zmin) / grid->gridStep);
   
   if((grid->grid = (USHORT ***)calloc(xsize, sizeof(USHORT **)))!=NULL)
   {
      for(i=0; i<xsize; i++)
      {
         if((grid->grid[i] = (USHORT **)calloc(ysize, sizeof(USHORT *)))
            !=NULL)
         {
            for(j=0; j<ysize; j++)
            {
               if((grid->grid[i][j] = (USHORT *)calloc(zsize, 
                                                       sizeof(USHORT)))
                  ==NULL)
                  goto DIE;
            }
         }
         else
         {
            goto DIE;
         }
      }
   }
   else
   {
      goto DIE;
   }
   
   if((grid->neighbours = (PDB *****)calloc(xsize, sizeof(PDB ****)))
      !=NULL)
   {
      for(i=0; i<xsize; i++)
      {
         if((grid->neighbours[i] = 
             (PDB ****)calloc(ysize, sizeof(PDB ***))) != NULL)
         {
            for(j=0; j<ysize; j++)
            {
               if((grid->neighbours[i][j] = 
                   (PDB ***)calloc(zsize, sizeof(PDB **)))!=NULL)
               {
                  for(k=0; k<zsize; k++)
                  {
                     if((grid->neighbours[i][j][k] = 
                         (PDB **)calloc(MAXNEIGHBOURS, sizeof(PDB *))) ==
                        NULL)
                     {
                        goto DIE;
                     }
                  }
               }
               else
               {
                  goto DIE;
               }
            }
         }
         else
         {
            goto DIE;
         }
      }
   }
   else
   {
      goto DIE;
   }
   
   return(TRUE);
   
 DIE:
   FreeGrid(grid);
   
   return(FALSE);
}


/************************************************************************/
/*>void MarkVoidAndProteinPoints(PDB *pdb, GRID *grid, REAL probeSize)
   -------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Run through all grid points assigning as protein if a protein atom
   is near enough.

   Also expands the size of each protein atom and build a list of these
   expanded protein atoms that are close to each grid point.

   17.10.01 Original   By: ACRM
*/
void MarkVoidAndProteinPoints(PDB *pdb, GRID *grid, REAL probeSize)
{
   REAL x,y,z;
   int  xi = 0;
   int  xm, ym, zm;

   xm = (int)((grid->xmax - grid->xmin)/grid->gridStep);
   if((grid->xmin + (xm*grid->gridStep)) > grid->xmax)
     xm--;

   ym = (int)((grid->ymax - grid->ymin)/grid->gridStep);
   if((grid->ymin + (ym*grid->gridStep)) > grid->ymax)
     ym--;

   zm = (int)((grid->zmax - grid->zmin)/grid->gridStep);
   if((grid->zmin + (zm*grid->gridStep)) > grid->zmax)
     zm--;

#pragma omp parallel for private(x,y,z) schedule (dynamic)
   for(xi=0; xi<=xm; xi++)
   {
      int yi;
      x=grid->xmin + xi*grid->gridStep;

      for(yi=0; yi<=ym; yi++)
      {
         int zi;
	 y=grid->ymin + yi*grid->gridStep;

         for(zi=0; zi<=zm; zi++)
         {
   	    z=grid->zmin + zi*grid->gridStep;
            if(AtomNear(pdb, x, y, z, probeSize,
                        grid->neighbours[xi][yi][zi]))
            {
               grid->grid[xi][yi][zi] = TYPE_PROTEIN;
            }
            else
            {
               grid->grid[xi][yi][zi] = TYPE_VOID;
            }
         }
      }
   }
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     REAL *gridStep, REAL *probeSize, REAL *solvSize,
                     BOOL *doRefine)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            REAL   *gridStep
            REAL   *probeSize
            REAL   *solvSize
            BOOL   *doRefine
   Returns: BOOL                Success?

   Parse the command line
   
   10.10.01 Original    By: ACRM
   06.02.02 Added -w
   09.07.02 Added -R
   17.07.02 Added -or
   20.09.02 Added -S
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *gridStep, REAL *probeSize, REAL *solvSize,
                  BOOL *doRefine)
{
   int i;
   
   argc--;
   argv++;

   gFlags.printVoids           = FALSE;
   gFlags.printAtoms           = FALSE;
   gFlags.printSolvent         = FALSE;
   gFlags.printRefinedVoids    = FALSE;
   gFlags.printWaterNeighbours = FALSE;
   gFlags.edgeConnected        = FALSE;
   gFlags.cornerConnected      = FALSE;
   gFlags.leak                 = FALSE;
   gFlags.quiet                = FALSE;
   gFlags.reassignVoids        = FALSE;
   gFlags.reassignSurface      = FALSE;
   
   *gridStep                   = DEFAULT_GRID_SIZE;
   *probeSize                  = DEFAULT_PROBE_SIZE;
   *solvSize                   = DEFAULT_WATER_SIZE;
   *doRefine                   = FALSE;

   infile[0] = outfile[0]      = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'g':
            argc--;
            argv++;
            if((!argc) || !sscanf(argv[0],"%lf", gridStep))
               return(FALSE);
            break;
         case 'p':
            argc--;
            argv++;
            if((!argc) || !sscanf(argv[0],"%lf", probeSize))
               return(FALSE);
            break;
         case 's':
            argc--;
            argv++;
            if((!argc) || !sscanf(argv[0],"%lf", solvSize))
               return(FALSE);
            break;
         case 'r':
            *doRefine = TRUE;
            break;
         case 'R':
            *doRefine = TRUE;
            gFlags.reassignVoids = TRUE;
            break;
         case 'e':
            gFlags.edgeConnected = TRUE;
            break;
         case 'l':
            gFlags.leak = TRUE;
            break;
         case 'q':
            gFlags.quiet = TRUE;
            break;
         case 'c':
            gFlags.edgeConnected   = TRUE;
            gFlags.cornerConnected = TRUE;
            break;
         case 'o':
            gFlags.printVoids = TRUE;
            for(i=1; argv[0][i]; i++)
            {
               if(argv[0][i] == 'a')
               {
                  gFlags.printAtoms = TRUE;
               }
               else if(argv[0][i] == 's')
               {
                  gFlags.printSolvent = TRUE;
               }
            }
            argc--;
            argv++;
            if((!argc) || 
               (argv[0][0] == '-') || 
               ((gFlags.gridFile = fopen(argv[0],"w"))==NULL))
            {
               if(argc)
               {
                  fprintf(stderr,"Unable to write file: %s\n", argv[0]);
               }
               return(FALSE);
            }
            break;
         case 'f':
            gFlags.printRefinedVoids = TRUE;
            argc--;
            argv++;
            if((!argc) || 
               (argv[0][0] == '-') || 
               ((gFlags.gridRefinedFile = fopen(argv[0],"w"))==NULL))
            {
               if(argc)
               {
                  fprintf(stderr,"Unable to write file: %s\n", argv[0]);
               }
               return(FALSE);
            }
            break;
         case 'n':
            gFlags.printNeighbours = TRUE;
            argc--;
            argv++;
            if((!argc) || 
               (argv[0][0] == '-') || 
               ((gFlags.neighbourFile = fopen(argv[0],"w"))==NULL))
            {
               if(argc)
               {
                  fprintf(stderr,"Unable to write file: %s\n", argv[0]);
               }
               return(FALSE);
            }
            break;
         case 'O':
            if(argv[0][2] == '\0')
            {
               return(FALSE);
            }
            else
            {
               char direction = argv[0][2];
               REAL shift;

               argc--;
               argv++;
               if((!argc) || !sscanf(argv[0],"%lf", &shift))
                  return(FALSE);
               
               switch(direction)
               {
               case 'x':
                  gOffsetX = shift;
                  break;
               case 'y':
                  gOffsetY = shift;
                  break;
               case 'z':
                  gOffsetZ = shift;
                  break;
               default:
                  return(FALSE);
               }
            }
            break;
         case 'S':
            gFlags.reassignSurface = TRUE;
            break;
         case 'w':
            gFlags.printWaterNeighbours = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void SetAsSolvent(GRID *grid, int ix, int iy, int iz, 
                     int numberOfSolvatedNeighbours, REAL solvRadiusSq,
                     int type)
   --------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Set the specified point as being a solvent point (TYPE_SOLVENT) and
   all points within the solvent radius as TYPE_SOLVENT2

   17.10.01 Original   By: ACRM
*/
void SetAsSolvent(GRID *grid, int ix, int iy, int iz, 
                  int numberOfSolvatedNeighbours, REAL solvRadiusSq,
                  int type)
{
   int jx, jy, jz;
   REAL xi,yi,zi,
        xj,yj,zj;
   
   xi = ix * grid->gridStep;
   yi = iy * grid->gridStep;
   zi = iz * grid->gridStep;
   
   SETLOWBYTE(grid->grid[ix][iy][iz], TYPE_SOLVENT);

   /* Mark other points within the solvent sphere                       */
   for(jx=ix-numberOfSolvatedNeighbours;
       jx<=ix+numberOfSolvatedNeighbours;
       jx++)
   {
      if(jx>0 && jx<grid->ixmax-1)
      {
         xj = jx * grid->gridStep;
         
         for(jy=iy-numberOfSolvatedNeighbours;
             jy<=iy+numberOfSolvatedNeighbours;
             jy++)
         {
            if(jy>0 && jy<grid->iymax-1)
            {
               yj = jy * grid->gridStep;

               for(jz=iz-numberOfSolvatedNeighbours;
                   jz<=iz+numberOfSolvatedNeighbours;
                   jz++)
               {
                  if(jz>0 && jz<grid->izmax-1)
                  {
                     zj = jz * grid->gridStep;

                     if(((xj-xi) * (xj-xi)) +
                        ((yj-yi) * (yj-yi)) +
                        ((zj-zi) * (zj-zi))
                        <= solvRadiusSq)
                     {
                        if(ISSET(grid->grid[jx][jy][jz], TYPE_VOID))
                        {
                           SETLOWBYTE(grid->grid[jx][jy][jz], TYPE_SOLVENT2);
                        }
                     }
                  }
               }
            }
         }
      }
   }
}


/************************************************************************/
/*>void PrintGrid(FILE *fp, GRID *grid, BOOL refined)
   --------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Prints non-void grid points (i.e. solvent and protein points)

   Takes an option to print the refined grid rather than the initial
   grid, but this currently does nothing

   17.10.01 Original   By: ACRM
*/
void PrintGrid(FILE *fp, GRID *grid, BOOL refined)
{
   REAL  x,y,z;
   int   atomnum=1, ix, iy, iz;

   for(x=grid->xmin,ix=0; x<=grid->xmax; x+=grid->gridStep,ix++)
   {
      for(y=grid->ymin,iy=0; y<=grid->ymax; y+=grid->gridStep,iy++)
      {
         for(z=grid->zmin,iz=0; z<=grid->zmax; z+=grid->gridStep,iz++)
         {
            switch(grid->grid[ix][iy][iz])
            {
            case TYPE_PROTEIN:
            case TYPE_PROTEIN2:
               if(!refined && gFlags.printAtoms)
               {
                  fprintf(fp, "ATOM  %5d  CA  GLY Z   \
1    %8.3f%8.3f%8.3f  1.00 20.00\n",
                          atomnum++, x,y,z);
               }
               break;
            case TYPE_SOLVENT:
               if(!refined && gFlags.printSolvent)
               {
                  fprintf(fp, "ATOM  %5d  N   GLY Z   \
1    %8.3f%8.3f%8.3f  1.00 20.00\n",
                          atomnum++, x,y,z);
               }
               break;
            default:
               break;
            }
         }
      }
   }
}


/************************************************************************/
/*>VOIDS *ClusterVoids(GRID *grid, int *nvoids)
   --------------------------------------------
   Input:     
   Output:    
   Returns:   

   First prints solvent and atom grid points if required to do so
   Place all void points into clusters by flood filling each group in
   turn

   17.10.01 Original   By: ACRM
   02.07.02 Added initialization of surfaceVoid element of VOIDS
*/
VOIDS *ClusterVoids(GRID *grid, int *nvoids)
{
   REAL  x,y,z;
   int   ix, iy, iz, 
         voidset = 0;
   VOIDS *voids  = NULL,
         *v      = NULL;
   
   for(x=grid->xmin, ix=0; x<=grid->xmax; x+=grid->gridStep, ix++)
   {
      for(y=grid->ymin, iy=0; y<=grid->ymax; y+=grid->gridStep, iy++)
      {
         for(z=grid->zmin, iz=0; z<=grid->zmax; z+=grid->gridStep, iz++)
         {
            if(grid->grid[ix][iy][iz] == TYPE_VOID)
            {
               if(!ISSET(grid->grid[ix][iy][iz], TYPE_ASSIGNED))
               {
                  if(voids==NULL)
                  {
                     INIT(voids, VOIDS);
                     v = voids;
                  }
                  else
                  {
                     ALLOCNEXT(v, VOIDS);
                  }
                  if(v==NULL)
                  {
                     fprintf(stderr,"No memory for voids list\n");
                     return(NULL);
                  }
                  v->pointlist=NULL;
                  v->ix = v->ixmin = v->ixmax = ix;
                  v->iy = v->iymin = v->iymax = iy;
                  v->iz = v->izmin = v->izmax = iz;
                  v->surfaceVoid = FALSE;
                  
                  v->voxelCount = FloodFillVoid(grid, ix, iy, iz, 
                                                v, ++voidset);
                  v->volume = v->voxelCount * grid->voxelVolume;

                  /* Calculate real bounds of void                      */
                  v->xmin = grid->xmin + grid->gridStep * v->ixmin;
                  v->ymin = grid->ymin + grid->gridStep * v->iymin;
                  v->zmin = grid->zmin + grid->gridStep * v->izmin;
                  v->xmax = grid->xmin + grid->gridStep * v->ixmax;
                  v->ymax = grid->ymin + grid->gridStep * v->iymax;
                  v->zmax = grid->zmin + grid->gridStep * v->izmax;
               }
            }
         }
      }
   }

   *nvoids = voidset;
   
   return(voids);
}


/************************************************************************/
/*>int FloodFillVoid(GRID *grid, int ix, int iy, int iz, VOIDS *v, 
                     int voidset)
   ---------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Recursive routine to flood fill a void and therefore cluster void
   points into individual voids

   Always does 6-connected filling (i.e. moving to the 6 face-adjacent
   voxels). Optionally does 18-connected filling (i.e. also do the
   12 edge-connected voxels) or 26-connected filling (i.e. also do the
   8 corner-connected voxels).

   17.10.01 Original   By: ACRM
   02.07.02 Checks for solvent boundary and setting surfaceVoid flag
*/
int FloodFillVoid(GRID *grid, int ix, int iy, int iz, VOIDS *v, 
                  int voidset)
{
   int voxelCount = 0,
       i, j, k;
   POINTLIST *p;
   
   if((ix>=0) && (ix<grid->ixmax) &&
      (iy>=0) && (iy<grid->iymax) &&
      (iz>=0) && (iz<grid->izmax))
   {
      if(grid->grid[ix][iy][iz] == TYPE_VOID)
      {
         /* Set this point as an assigned void                          */
         grid->grid[ix][iy][iz] = voidset;
         SET(grid->grid[ix][iy][iz], TYPE_ASSIGNED);
         voxelCount++;
         
         if(v->pointlist == NULL)
         {
            INIT((v->pointlist), POINTLIST);
            p=v->pointlist;
         }
         else
         {
            p=v->pointlist;
            LAST(p);
            ALLOCNEXT(p, POINTLIST);
         }
         if(p==NULL)
         {
            fprintf(stderr,"No memory for voids point list\n");
            exit(1);
         }
         
         /* Put this grid->grid point into the point list               */
         p->ix = ix;
         p->iy = iy;
         p->iz = iz;
         p->x = (ix * grid->gridStep) + grid->xmin;
         p->y = (iy * grid->gridStep) + grid->ymin;
         p->z = (iz * grid->gridStep) + grid->zmin;
         p->nearest = NULL;
         
         if(ix < v->ixmin) v->ixmin = ix;
         if(iy < v->iymin) v->iymin = iy;
         if(iz < v->izmin) v->izmin = iz;
         if(ix > v->ixmax) v->ixmax = ix;
         if(iy > v->iymax) v->iymax = iy;
         if(iz > v->izmax) v->izmax = iz;
         
         /* Recurse to neighbouring points                              */
         voxelCount += FloodFillVoid(grid,ix-1,iy,  iz,  v,voidset);
         voxelCount += FloodFillVoid(grid,ix,  iy-1,iz,  v,voidset);
         voxelCount += FloodFillVoid(grid,ix+1,iy,  iz,  v,voidset);
         voxelCount += FloodFillVoid(grid,ix,  iy,  iz-1,v,voidset);
         voxelCount += FloodFillVoid(grid,ix,  iy,  iz+1,v,voidset);
         voxelCount += FloodFillVoid(grid,ix,  iy+1,iz,  v,voidset);
         
         if(gFlags.edgeConnected)
         {
            voxelCount += FloodFillVoid(grid,ix-1,iy-1,iz,  v,voidset);
            voxelCount += FloodFillVoid(grid,ix-1,iy,  iz-1,v,voidset);
            voxelCount += FloodFillVoid(grid,ix-1,iy,  iz+1,v,voidset);
            voxelCount += FloodFillVoid(grid,ix-1,iy+1,iz,  v,voidset);
            voxelCount += FloodFillVoid(grid,ix,  iy-1,iz-1,v,voidset);
            voxelCount += FloodFillVoid(grid,ix,  iy-1,iz+1,v,voidset);
            voxelCount += FloodFillVoid(grid,ix,  iy+1,iz-1,v,voidset);
            voxelCount += FloodFillVoid(grid,ix,  iy+1,iz+1,v,voidset);
            voxelCount += FloodFillVoid(grid,ix+1,iy-1,iz,  v,voidset);
            voxelCount += FloodFillVoid(grid,ix+1,iy,  iz-1,v,voidset);
            voxelCount += FloodFillVoid(grid,ix+1,iy,  iz+1,v,voidset);
            voxelCount += FloodFillVoid(grid,ix+1,iy+1,iz,  v,voidset);
            
            if(gFlags.cornerConnected)
            {
               voxelCount += FloodFillVoid(grid,ix-1,iy-1,iz-1,v,voidset);
               voxelCount += FloodFillVoid(grid,ix-1,iy-1,iz+1,v,voidset);
               voxelCount += FloodFillVoid(grid,ix-1,iy+1,iz-1,v,voidset);
               voxelCount += FloodFillVoid(grid,ix-1,iy+1,iz+1,v,voidset);
               voxelCount += FloodFillVoid(grid,ix+1,iy-1,iz-1,v,voidset);
               voxelCount += FloodFillVoid(grid,ix+1,iy-1,iz+1,v,voidset);
               voxelCount += FloodFillVoid(grid,ix+1,iy+1,iz-1,v,voidset);
               voxelCount += FloodFillVoid(grid,ix+1,iy+1,iz+1,v,voidset);
            }
         }

         /* If not already assigned this void as being a surface void   */
         if(!v->surfaceVoid)
         {
            /* If any of the fully-connected neighbouring points is 
               solvent then set this as a surface void
            */
            for(i=-NUMGRIDNEIGHBOURS; i<=NUMGRIDNEIGHBOURS; i++)
            {
               for(j=-NUMGRIDNEIGHBOURS; j<=NUMGRIDNEIGHBOURS; j++)
               {
                  for(k=-NUMGRIDNEIGHBOURS; k<=NUMGRIDNEIGHBOURS; k++)
                  {
                     if(IsSolvent(grid, ix+i, iy+j, iz+k))
                     {
                        v->surfaceVoid = TRUE;
                        /* Force break of all nested loops              */
                        i=j=k=NUMGRIDNEIGHBOURS+1;
                        break;
                     }
                  }
               }
            }
         }
      }
   }

   return(voxelCount);
}


/************************************************************************/
/*>void FindNearestAtoms(PDB *pdb, VOIDS *voids)
   ---------------------------------------------
   Input:     
   Output:    
   Returns:   

   Identify which atom is nearest to each void point

   17.10.01 Original   By: ACRM
*/
void FindNearestAtoms(PDB *pdb, VOIDS *voids)
{
   VOIDS *v;
   POINTLIST *p;

#ifdef _OPENMP
   VOIDS **voidarray;
   int   nvoids, i;

   for(v=voids, nvoids=0; v!=NULL; NEXT(v))
   {
      nvoids++;
   }
   if((voidarray=(VOIDS **)malloc(nvoids*sizeof(VOIDS *)))==NULL)
   {
      fprintf(stderr,"Error!: No memory for voids array\n");
      exit(1);
   }
   for(v=voids, i=0; v!=NULL; NEXT(v))
   {
      voidarray[i++] = v;
   }

#pragma omp parallel for private(v,p)
   for(i=0; i<nvoids; i++)
   {
      v=voidarray[i];

      for(p=v->pointlist; p!=NULL; NEXT(p))
      {
#pragma omp atomic write
         p->nearest = FindTheNearestAtom(pdb, p);
      }
   }
#else
   for(v=voids; v!=NULL; NEXT(v))
   {
      for(p=v->pointlist; p!=NULL; NEXT(p))
      {
         p->nearest = FindTheNearestAtom(pdb, p);
      }
   }
#endif
}


/************************************************************************/
/*>PDB *FindTheNearestAtom(PDB *pdb, POINTLIST *point)
   ---------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Given a void point, find which atom is nearest

   17.10.01 Original   By: ACRM
*/
PDB *FindTheNearestAtom(PDB *pdb, POINTLIST *point)
{
   PDB *p, *nearest = NULL;
   REAL mindist = (REAL)100.0, dist;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      dist = DISTSQ(point, p);
      if(dist < mindist)
      {
         mindist = dist;
         nearest = p;
      }
   }

   return(nearest);
}


/************************************************************************/
/*>void SetRadii(PDB *pdb)
   -----------------------
   Input:     
   Output:    
   Returns:   

   Set radii of atoms in the PDB linked list. These values taken from
   ACCESS.

*** TODO add ability to read values from a file                        ***

   17.10.01 Original   By: ACRM
*/
void SetRadii(PDB *pdb)
{
   PDB *p;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"C   ",4))
      {
         p->occ = 1.76;
      }
      else if((!strncmp(p->resnam, "ARG", 3) && 
               !strncmp(p->atnam, "CZ  ", 4)) ||
              (!strncmp(p->resnam, "AS",  2) &&
               !strncmp(p->atnam, "CG  ", 4)) || /* ASP,ASN,ASX */
              (!strncmp(p->resnam, "GL",  2) &&
               !strncmp(p->atnam, "CD  ", 4)) || /* GLU,GLN,GLX */
              (!strncmp(p->resnam, "PCA", 3) &&
               !strncmp(p->atnam, "CD  ", 4)) ||
              (!strncmp(p->resnam, "HIS", 3) &&
               !strncmp(p->atnam, "CG  ", 4)) ||
              (!strncmp(p->resnam, "HIS", 3) &&
               !strncmp(p->atnam, "CE1 ", 4)) ||
              (!strncmp(p->resnam, "HIS", 3) &&
               !strncmp(p->atnam, "CD2 ", 4)) ||
              (!strncmp(p->resnam, "HIS", 3) &&
               !strncmp(p->atnam, "AD1 ", 4)) ||
              (!strncmp(p->resnam, "HIS", 3) &&
               !strncmp(p->atnam, "AD2 ", 4)) ||
              (!strncmp(p->resnam, "HIS", 3) &&
               !strncmp(p->atnam, "AE1 ", 4)) ||
              (!strncmp(p->resnam, "HIS", 3) &&
               !strncmp(p->atnam, "AE2 ", 4)) ||
              (!strncmp(p->resnam, "PHE", 3) &&
               !strncmp(p->atnam, "CG  ", 4)) ||
              (!strncmp(p->resnam, "PHE", 3) &&
               !strncmp(p->atnam, "CD1 ", 4)) ||
              (!strncmp(p->resnam, "PHE", 3) &&
               !strncmp(p->atnam, "CE1 ", 4)) ||
              (!strncmp(p->resnam, "PHE", 3) &&
               !strncmp(p->atnam, "CZ  ", 4)) ||
              (!strncmp(p->resnam, "PHE", 3) &&
               !strncmp(p->atnam, "CE2 ", 4)) ||
              (!strncmp(p->resnam, "PHE", 3) &&
               !strncmp(p->atnam, "CD2 ", 4)) ||
              (!strncmp(p->resnam, "TYR", 3) &&
               !strncmp(p->atnam, "CG  ", 4)) ||
              (!strncmp(p->resnam, "TYR", 3) &&
               !strncmp(p->atnam, "CD1 ", 4)) ||
              (!strncmp(p->resnam, "TYR", 3) &&
               !strncmp(p->atnam, "CE1 ", 4)) ||
              (!strncmp(p->resnam, "TYR", 3) &&
               !strncmp(p->atnam, "CZ  ", 4)) ||
              (!strncmp(p->resnam, "TYR", 3) &&
               !strncmp(p->atnam, "CE2 ", 4)) ||
              (!strncmp(p->resnam, "TYR", 3) &&
               !strncmp(p->atnam, "CD2 ", 4)) ||
              (!strncmp(p->resnam, "TRP", 3) &&
               !strncmp(p->atnam, "CG  ", 4)) ||
              (!strncmp(p->resnam, "TRP", 3) &&
               !strncmp(p->atnam, "CD1 ", 4)) ||
              (!strncmp(p->resnam, "TRP", 3) &&
               !strncmp(p->atnam, "CE2 ", 4)) ||
              (!strncmp(p->resnam, "TRP", 3) &&
               !strncmp(p->atnam, "CZ2 ", 4)) ||
              (!strncmp(p->resnam, "TRP", 3) &&
               !strncmp(p->atnam, "CH2 ", 4)) ||
              (!strncmp(p->resnam, "TRP", 3) &&
               !strncmp(p->atnam, "CZ3 ", 4)) ||
              (!strncmp(p->resnam, "TRP", 3) &&
               !strncmp(p->atnam, "CE3 ", 4)) ||
              (!strncmp(p->resnam, "TRP", 3) &&
               !strncmp(p->atnam, "CD2 ", 4)) ||
              (!strncmp(p->resnam, "ASN", 3) &&
               !strncmp(p->atnam, "CG  ", 4)) ||
              (!strncmp(p->resnam, "ASN", 3) &&
               !strncmp(p->atnam, "CG  ", 4)) ||
              (!strncmp(p->resnam, "ASN", 3) &&
               !strncmp(p->atnam, "CG  ", 4)) ||
              (!strncmp(p->resnam, "ASN", 3) &&
               !strncmp(p->atnam, "CG  ", 4)))
      {
         p->occ = 1.76;
      }
      else if(p->atnam[0] == 'N')
      {
         p->occ = 1.65;
      }
      else if(p->atnam[0] == 'O')
      {
         p->occ = 1.40;
      }
      else if(p->atnam[0] == 'S')
      {
         p->occ = 1.85;
      }
      else if(p->atnam[0] == 'P')
      {
         p->occ = 1.9;
      }
      else
      {
         p->occ = 1.87;
      }
   }
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Input:     
   Output:    
   Returns:   

   Prints a usage message

   17.10.01 Original   By: ACRM
   17.06.02 V1.1
   02.07.02 V1.2
   09.07.02 added -R
   18.08.08 V1.3
   31.07.17 V1.4
   27.11.20 V1.5
*/
void Usage(void)
{
   fprintf(stderr,"\navp V1.5 (c) 2001-20, Prof. Andrew C.R. Martin, \
University of Reading, UCL\n");

   fprintf(stderr,"\nUsage: avp [-q] [-g gridspacing] [-p probesize] \
[-s solventsize]\n");
   fprintf(stderr,"           [-r] [-e] [-c] [-o[a][s] file] \n");
   fprintf(stderr,"           [-f file] [-l] [-n file] [-O(xyz) value] \
[-w] [-S] [file.pdb [file.out]]\n");

   fprintf(stderr,"   -q Quiet - do not report progress\n");
   fprintf(stderr,"   -g Specify the grid spacing (Default: %f)\n",
           DEFAULT_GRID_SIZE);
   fprintf(stderr,"   -p Specify the probe size (Default: %f)\n",
           DEFAULT_PROBE_SIZE);
   fprintf(stderr,"   -s Specify the solvent size (Default: %f)\n",
           DEFAULT_WATER_SIZE);
   fprintf(stderr,"   -r Refine void sizes\n");
   fprintf(stderr,"   -R Redefine void points and refine void sizes\n");
   fprintf(stderr,"   -e Assign points to voids with edge connections\n");
   fprintf(stderr,"   -c Assign points to voids with corner connections \
(implies -e)\n");
   fprintf(stderr,"   -o Output the void grid points to file. With 'a', \
also output atom\n");
   fprintf(stderr,"      grid points; with 's', also output solvent grid \
points\n");
   fprintf(stderr,"   -f Output the refined void points to a file. Used \
with -r\n");
   fprintf(stderr,"   -l Include protein voxels next to void and solvent \
in volume refinement\n");
   fprintf(stderr,"   -n Output atom records for atoms nearest to each \
void to file.\n");
   fprintf(stderr,"   -Ox -Oy -Oz Specify an offset for the grid\n");
   fprintf(stderr,"   -w Print a list of waters that neighbour voids\n");
   fprintf(stderr,"   -S Reassign surface protein to solvent if can \
fit a solvent near\n");

   fprintf(stderr,"\navp (Another Void Program) is a program to \
calculate void volumes in\n");
   fprintf(stderr,"proteins. It uses a simple grid-based method but \
separates the probe\n");
   fprintf(stderr,"size used to define void points as being voids (-p) \
from the probe\n");
   fprintf(stderr,"size used to define channels to the surface (-s). \
Thus one can find\n");
   fprintf(stderr,"very small voids without these being connected via \
very small diameter\n");
   fprintf(stderr,"passages to the surface.\n\n");
}


/************************************************************************/
/*>void PrintVoids(FILE *out, VOIDS *voids)
   ----------------------------------------
   Input:     
   Output:    
   Returns:   

   Print summary details of each void. If required, also print the actual
   void grid points in PDB format and/or the neighbouring ATOMs from the
   original PDB file.

   17.10.01 Original   By: ACRM
   17.06.02 Added printing of largest individual void
   18.06.02 Added printing of void mid-point
   02.07.02 Added printing of surfaceVoid flag
*/
void PrintVoids(FILE *out, VOIDS *voids)
{
   VOIDS *v;
   POINTLIST *p;
   int   count  = 1,
         resnum = 2,
         atomnum = 2000;
   REAL  totalVoidVolume          = (REAL)0.0,
         totalBuriedVoidVolume    = (REAL)0.0,
         totalSurfaceVoidVolume   = (REAL)0.0,
         largestVoidVolume        = (REAL)0.0,
         largestBuriedVoidVolume  = (REAL)0.0,
         largestSurfaceVoidVolume = (REAL)0.0;
   
   
   for(v=voids, resnum=2; v!=NULL; NEXT(v), resnum++)
   {
      fprintf(out,"Void: %4d voxelCount: %3d volume: %8.3f \
midpoint: %8.3f%8.3f%8.3f %s\n",
              count++, v->voxelCount, v->volume,
              (v->xmin + v->xmax)/2.0,
              (v->ymin + v->ymax)/2.0,
              (v->zmin + v->zmax)/2.0,
              ((v->surfaceVoid)?"Surface-void":" Buried-void"));
      totalVoidVolume += v->volume;
      if(v->volume > largestVoidVolume)
      {
         largestVoidVolume = v->volume;
      }

      if(v->surfaceVoid)
      {
         totalSurfaceVoidVolume += v->volume;
         if(v->volume > largestSurfaceVoidVolume)
         {
            largestSurfaceVoidVolume = v->volume;
         }
      }
      else
      {
         totalBuriedVoidVolume += v->volume;
         if(v->volume > largestBuriedVoidVolume)
         {
            largestBuriedVoidVolume = v->volume;
         }
      }

      if(gFlags.printVoids)
      {
         for(p=v->pointlist; p!=NULL; NEXT(p))
         {
            if(v->surfaceVoid)
            {
               fprintf(gFlags.gridFile, "ATOM  %5d  O   GLY X%4d    \
%8.3f%8.3f%8.3f  1.00 20.00\n",  atomnum++, resnum, p->x, p->y, p->z);
            }
            else
            {
               fprintf(gFlags.gridFile, "ATOM  %5d  O   GLY Y%4d    \
%8.3f%8.3f%8.3f  1.00 20.00\n",  atomnum++, resnum, p->x, p->y, p->z);
            }
         }
      }

      if(gFlags.printNeighbours)
      {
         /* Set occupancies to zero for this void. We'll use this as a
            flag to say whether this atom has been printed
         */
         for(p=v->pointlist; p!=NULL; NEXT(p))
         {
            p->nearest->occ = 0.0;
         }
         
         for(p=v->pointlist; p!=NULL; NEXT(p))
         {
            if(p->nearest->occ < 0.01)
            {
               p->nearest->occ = (REAL)1.0;
               blWritePDBRecord(gFlags.neighbourFile, p->nearest);
            }
         }
      }
   }

   fprintf(out, "\nTotal buried void volume: %f\n", 
           totalBuriedVoidVolume);
   fprintf(out, "Largest buried void volume: %f\n", 
           largestBuriedVoidVolume);

   fprintf(out, "\nTotal surface void volume: %f\n", 
           totalSurfaceVoidVolume);
   fprintf(out, "Largest surface void volume: %f\n", 
           largestSurfaceVoidVolume);

   fprintf(out, "\nTotal void volume: %f\n", totalVoidVolume);
   fprintf(out, "Largest void volume: %f\n", largestVoidVolume);
}


/************************************************************************/
/*>BOOL AtomNear(PDB *pdb, REAL x, REAL y, REAL z, REAL probeSize,
                 PDB **neighbours)
   ---------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Tests whether any atom is near to the specified position
   If the neighbours array is non-NULL, then tests whether it is already
   filled in with neighbours. If not, fills it in with all neighbouring 
   atoms.

   17.10.01 Original   By: ACRM
   31.01.01 Added code to build neighbours list
   02.07.02 Modified so that the neighbour list is used internally
            if it exists rather than having to step through the
            neighbouring atoms in a separate routine
*/
BOOL AtomNear(PDB *pdb, REAL x, REAL y, REAL z, REAL probeSize,
              PDB **neighbours)

{
   static PDB **pdbArray = NULL;
#pragma omp threadprivate(pdbArray)
   static int natoms     = 0;
   VEC3F      point;
   REAL       cutoffSq,
              cutoff2Sq;
   REAL       xmin = 0.0,
              xmax = 0.0,
              ymin = 0.0,
              ymax = 0.0,
              zmin = 0.0,
              zmax = 0.0;
   int        start = 0,
              stop = 0,
              neighbCount = 0,
              i;
   BOOL       retValue = FALSE;
   static BOOL warned = FALSE;
   
   
   /* Create an index array into the PDB linked list which is sorted by
      x coordinate - done on the first call only
   */
   if((pdbArray==NULL) && (pdb != NULL))
   {
      /* Malloc and populate pdbArray[] with PDB pointers               */
      if((pdbArray = blIndexPDB(pdb, &natoms)) == NULL)
      {
         fprintf(stderr,"No memory for PDB index array\n");
         exit(1);
      }

      /* Now sort the atoms on the basis of coordinates                 */
#pragma omp critical
      qsort((void *)pdbArray, (size_t)natoms, sizeof(PDB *), 
            CompareCoords);
   }
   
   /* Atoms are only interesting if their coordinates lie between these
      limits
   */
   if((neighbours==NULL) || (neighbours[0] == NULL))
   {
      xmin = x-(probeSize + RADEXPAND*MAXATOMRAD);
      xmax = x+(probeSize + RADEXPAND*MAXATOMRAD);
      ymin = y-(probeSize + RADEXPAND*MAXATOMRAD);
      ymax = y+(probeSize + RADEXPAND*MAXATOMRAD);
      zmin = z-(probeSize + RADEXPAND*MAXATOMRAD);
      zmax = z+(probeSize + RADEXPAND*MAXATOMRAD);

      start = FindFirstAtom(pdbArray, natoms, xmin);
      stop  = FindLastAtom(pdbArray, natoms, xmax);

      if((start==(-1)) || (stop==(-1)))
         return(FALSE);
   }
   
   point.x = x;
   point.y = y;
   point.z = z;
   
   if(neighbours==NULL) /* Not using a neighbour list                   */
   {
      for(i=start; i<=stop; i++)
      {
         if((pdbArray[i]->y >= ymin) &&
            (pdbArray[i]->y <= ymax) &&
            (pdbArray[i]->z >= zmin) &&
            (pdbArray[i]->z <= zmax))
         {
            cutoffSq = (probeSize+pdbArray[i]->occ) * 
                       (probeSize+pdbArray[i]->occ);
            if(DISTSQ(pdbArray[i], &point) < cutoffSq)
            {
               return(TRUE);
            }
         }
      }
   }
   else /* neighbours pointer is non-null                               */
   {
      if(neighbours[0] == NULL) /* Neighbours list not yet created      */
      {
         for(i=start; i<=stop; i++)
         {
            if((pdbArray[i]->y >= ymin) &&
               (pdbArray[i]->y <= ymax) &&
               (pdbArray[i]->z >= zmin) &&
               (pdbArray[i]->z <= zmax))
            {
               /* Calculate cutoff for this being a clash               */
               cutoffSq  = (probeSize + pdbArray[i]->occ) * 
                           (probeSize + pdbArray[i]->occ);

               /* Calculate cutoff for neighbour list                   */
               cutoff2Sq = (probeSize + RADEXPAND*pdbArray[i]->occ) * 
                           (probeSize + RADEXPAND*pdbArray[i]->occ);

               /* Add to neighbour list if near enough                  */
               if(DISTSQ(pdbArray[i], &point) < cutoff2Sq)
               {
#pragma omp critical
                  if(neighbCount < MAXNEIGHBOURS)
                  {
                     neighbours[neighbCount++] = pdbArray[i];
                  }
                  else
                  {
                     if(!warned)
                     {
                        fprintf(stderr,"Max number of atoms neighbouring \
a grid point exceeded. Increase MAXNEIGHBOURS\n");
                        warned = TRUE;
                     }
                  }

                  /* If it was actually close enough to clash, then set
                     the return value
                  */
                  if(DISTSQ(pdbArray[i], &point) < cutoffSq)
                  {
                     retValue = TRUE;
                  }
               }
            }
         }
      }
      else /* We already have a neighbour list                          */
      {
         /* Run through the neighbour list                              */
         for(i=0; i<MAXNEIGHBOURS; i++)
         {
            PDB *p = neighbours[i];
            
            if(p==NULL) /* Exit if we have run out of entries           */
            {
               break;
            }

            /* Calculate the cutoff and see if we are in range          */
            cutoffSq = (probeSize+p->occ) * 
                       (probeSize+p->occ);
            if(DISTSQ(p, &point) < cutoffSq)
            {
               retValue = TRUE;
               break;
            }
         }
      }
   }

   return(retValue);
}

/************************************************************************/
/*>int FindFirstAtom(PDB **array, int natoms, REAL x)
   --------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Find the first atom with x coordinate >= x

   31.10.01 Original   By: ACRM
*/
int FindFirstAtom(PDB **array, int natoms, REAL x)
{
   int bottom = 0,
       top = natoms,
       cutpoint;
   
   /* Check that at least one value will satisfy                        */
   if((array[natoms-1])->x < x)
   {
      return(-1);
   }
 
   for(;;)
   {
      /* Define the new cut point                                       */
      cutpoint=(int)((top-bottom)/2)+bottom;
      
      /* If this point satisfies....                                    */
      if((array[cutpoint])->x >= x)
      {
         /* Check if we're at the end of the array                      */
         if((cutpoint-1) < 0)
         {
            return(cutpoint);
         }
         
         /* If the next point doesn't satisy we have found what we're 
            after 
         */
         if((array[cutpoint-1])->x < x)
         {
            return(cutpoint);
         }
         
         /* Otherwise make this the new top                             */
         top=cutpoint;
      }
      else
      {
         /* Make this the new bottom                                    */
         bottom=cutpoint;
      }
   }
   return(-1);
}


/************************************************************************/
/*>int FindLastAtom(PDB **array, int natoms, REAL x)
   -------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Find the first atom with x coordinate <= x 

   31.10.01 Original   By: ACRM
*/
int FindLastAtom(PDB **array, int natoms, REAL x)
{
   int bottom = 0,
       top = natoms,
       cutpoint;

   /* Check that at least one value will satisfy */
   if((array[0])->x > x)
   {
      return(-1);
   }

   for(;;)
   {
      /* Define the new cut point */
      cutpoint=(int)((top-bottom)/2)+bottom;
      
      /* If this point satisfies.... */
      if((array[cutpoint])->x <= x)
      {
         
         /* Check if we're at the end of the array */
         if((cutpoint+1) >= natoms)
         {
            return(cutpoint);
         }
         
         /* If the next point doesn't satisy we have found what we're 
            after 
         */
         if((array[cutpoint+1])->x > x)
         {
            return(cutpoint);
         }
      
         /*        Otherwise make this the new bottom */
         bottom=cutpoint;
      }
      else
      {
         /*        Make this the new top */
         top=cutpoint;
      }
      
   }

   return(-1);
}


/************************************************************************/
/*>int CompareCoords(const void *p1, const void *p2)
   -------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Compare the coordinates of two atoms for the qsort() routine

   31.10.01 Original   By: ACRM
*/
int CompareCoords(const void *p1, const void *p2)
{
   PDB *atom1 = *(PDB **)p1;
   PDB *atom2 = *(PDB **)p2;

   if(atom1->x < atom2->x)
      return(-1);
   if(atom1->x > atom2->x)
      return(1);
   
   return(0);
}


#ifdef FLOOD_FILL_SOLVENT
/************************************************************************/
/*>BOOL MarkSolventPoints(PDB *pdb, GRID *grid, REAL solventStep,
                          REAL solventMult, BOOL quiet, BOOL doNeighbours)
   -----------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Identify solvent points by flood filling. This is faster than the
   alternative 6-directional search, but runs out of recursion stack
   space on anything much bigger than crambin...

   17.10.01 Original   By: ACRM
*/
BOOL MarkSolventPoints(PDB *pdb, GRID *grid, REAL solventStep, 
                       REAL solventMult, BOOL quiet, BOOL doNeighbours)
{
   int numberOfSolvatedNeighbours;
   REAL solvRadiusSq = solventMult * solventStep * 
                       solventMult * solventStep;

   numberOfSolvatedNeighbours = 
      (int)(solventMult * solventStep / grid->gridStep);

   /* Flood fill the solvent starting from point 0,0,0                  */
   FloodFillSolvent(grid, pdb, 0, 0, 0, solventStep, solvRadiusSq,
                    (doNeighbours?numberOfSolvatedNeighbours:0));
   return(TRUE);
}


/************************************************************************/
/*>void FloodFillSolvent(GRID *grid, PDB *pdb, int ix, int iy, int iz, 
                         REAL solvSize, REAL solvRadiusSq,
                         int numberOfSolvatedNeighbours)
   --------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Recursive routine for solvent flood filling

   31.10.01 Original   By: ACRM
   02.07.02 Uses neighbour list in calls to AtomNear()
*/
void FloodFillSolvent(GRID *grid, PDB *pdb, int ix, int iy, int iz, 
                      REAL solvSize, REAL solvRadiusSq,
                      int numberOfSolvatedNeighbours)
{
   gDepth++;
   if(gDepth > gMaxDepth)
   {
      gMaxDepth = gDepth;
   }
#ifdef DEBUG
      fprintf(stderr,"Stack depth: %d maxdepth: %d\n", gDepth, gMaxDepth);
#endif

   /* if this is a valid grid point currently VOID or SOLVENT2          */
   if((ix>=0) && (ix<grid->ixmax) &&
      (iy>=0) && (iy<grid->iymax) &&
      (iz>=0) && (iz<grid->izmax) &&
      ((grid->grid[ix][iy][iz] == TYPE_VOID) ||
       (grid->grid[ix][iy][iz] == TYPE_SOLVENT2)))
   {
      /* If no atom clashes                                             */
      if(!AtomNear(pdb,
                   grid->xmin + ix*grid->gridStep,
                   grid->ymin + iy*grid->gridStep,
                   grid->zmin + iz*grid->gridStep,
                   solvSize, grid->neighbours[ix][iy][iz]))
      {
         SetAsSolvent(grid,ix,iy,iz,
                      numberOfSolvatedNeighbours, 
                      solvRadiusSq,
                      TYPE_SOLVENT);
         
         FloodFillSolvent(grid, pdb, ix-1,iy,  iz,  solvSize, 
                          solvRadiusSq, numberOfSolvatedNeighbours);
         FloodFillSolvent(grid, pdb, ix,  iy-1,iz,  solvSize, 
                          solvRadiusSq, numberOfSolvatedNeighbours);
         FloodFillSolvent(grid, pdb, ix+1,iy,  iz,  solvSize, 
                          solvRadiusSq, numberOfSolvatedNeighbours);
         FloodFillSolvent(grid, pdb, ix,  iy,  iz-1,solvSize, 
                          solvRadiusSq, numberOfSolvatedNeighbours);
         FloodFillSolvent(grid, pdb, ix,  iy,  iz+1,solvSize, 
                          solvRadiusSq, numberOfSolvatedNeighbours);
         FloodFillSolvent(grid, pdb, ix,  iy+1,iz,  solvSize, 
                          solvRadiusSq, numberOfSolvatedNeighbours);
      }
   }
   gDepth--;
}

#else

/************************************************************************/
/*>void FlagBoxBoundsAsSolvent(GRID *grid, int ixmax, int iymax, 
                               int izmax, BOOL quiet)
   -------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Runs around the outside of the grid setting all points as solvent

   17.10.01 Original   By: ACRM
*/
void FlagBoxBoundsAsSolvent(GRID *grid, int ixmax, int iymax, int izmax,
                            BOOL quiet)
{
   int ix, iy, iz;
   BOOL warned = FALSE;

#pragma omp parallel for private(iy,iz)
   for(ix=0; ix<ixmax; ix++)
   {
      for(iy=0; iy<iymax; iy++)
      {
         if((grid->grid[ix][iy][0] == TYPE_PROTEIN) || 
            (grid->grid[ix][iy][izmax-1] == TYPE_PROTEIN))
         {
            if(!warned && !quiet)
            {
               fprintf(stderr,"Protein intersects edge of box!\n");
               warned = TRUE;
            }
         }
         if(grid->grid[ix][iy][0] != TYPE_PROTEIN)
         {
            grid->grid[ix][iy][0] = TYPE_SOLVENT2;
         }
         if(grid->grid[ix][iy][izmax-1] != TYPE_PROTEIN)
         {
            grid->grid[ix][iy][izmax-1] = TYPE_SOLVENT2;
         }
      }
      for(iz=0; iz<izmax; iz++)
      {
         if((grid->grid[ix][0][iz] == TYPE_PROTEIN) || 
            (grid->grid[ix][iymax-1][iz] == TYPE_PROTEIN))
         {
            if(!warned && !quiet)
            {
               fprintf(stderr,"Protein intersects edge of box!\n");
               warned = TRUE;
            }
         }
         if(grid->grid[ix][0][iz] != TYPE_PROTEIN)
         {
            grid->grid[ix][0][iz] = TYPE_SOLVENT2;
         }
         if(grid->grid[ix][iymax-1][iz] != TYPE_PROTEIN)
         {
            grid->grid[ix][iymax-1][iz] = TYPE_SOLVENT2;
         }
      }
   }
   
#pragma omp parallel for private(iz)
   for(iy=0; iy<iymax; iy++)
   {
      for(iz=0; iz<izmax; iz++)
      {
         if((grid->grid[0][iy][iz] == TYPE_PROTEIN) || 
            (grid->grid[ixmax-1][iy][iz] == TYPE_PROTEIN))
         {
            if(!warned && !quiet)
            {
               fprintf(stderr,"Protein intersects edge of box!\n");
               warned = TRUE;
            }
         }

         if(grid->grid[0][iy][iz] != TYPE_PROTEIN)
         {
            grid->grid[0][iy][iz] = TYPE_SOLVENT2;
         }
         if(grid->grid[ixmax-1][iy][iz] != TYPE_PROTEIN)
         {
            grid->grid[ixmax-1][iy][iz] = TYPE_SOLVENT2;
         }
      }
   }
}


/************************************************************************/
/*>int AdjacentIsSolvent(PDB *pdb, GRID *grid, int ix, int iy, int iz,
                         REAL solvSize, REAL xoff, REAL yoff, REAL zoff)
   ---------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Look at all 26 neighbouring points to see if any one of them is solvent
   If so, check this point is > solvSize away from protein and if OK, 
   then return TRUE

   17.10.01 Original   By: ACRM
   18.06.02 Now distinguishes 3 cases: SOLVENT_BULK (previously 'TRUE'),
            SOLVENT_NO (previously 'FALSE') and the new case 
            SOLVENT_MARGIN where this point has solvent neighbours, but
            is too close to an atom (previously was just 'FALSE')
            Also now accepts xoff,yoff,zoff offsets 
   02.07.02 Uses neighbour list in calls to AtomNear()
*/
int AdjacentIsSolvent(PDB *pdb, GRID *grid, int ix, int iy, int iz,
                      REAL solvSize, REAL xoff, REAL yoff, REAL zoff)
{
   /* If we are not on the boundary of the box                          */
   if((ix>0) && (ix<grid->ixmax) &&
      (iy>0) && (iy<grid->iymax) &&
      (iz>0) && (iz<grid->izmax))
   {
      /* If any neighbour is solvent and there is no atom too close,
         then return SOLVENT_BULK
      */
      if((grid->grid[ix-1][iy-1][iz-1] == TYPE_SOLVENT) ||
         (grid->grid[ix-1][iy-1][iz] == TYPE_SOLVENT)   ||
         (grid->grid[ix-1][iy-1][iz+1] == TYPE_SOLVENT) ||
         (grid->grid[ix-1][iy][iz-1] == TYPE_SOLVENT)   ||
         (grid->grid[ix-1][iy][iz] == TYPE_SOLVENT)     ||
         (grid->grid[ix-1][iy][iz+1] == TYPE_SOLVENT)   ||
         (grid->grid[ix-1][iy+1][iz-1] == TYPE_SOLVENT) ||
         (grid->grid[ix-1][iy+1][iz] == TYPE_SOLVENT)   ||
         (grid->grid[ix-1][iy+1][iz+1] == TYPE_SOLVENT) ||
         (grid->grid[ix][iy-1][iz-1] == TYPE_SOLVENT)   ||
         (grid->grid[ix][iy-1][iz] == TYPE_SOLVENT)     ||
         (grid->grid[ix][iy-1][iz+1] == TYPE_SOLVENT)   ||
         (grid->grid[ix][iy][iz-1] == TYPE_SOLVENT)     ||
         (grid->grid[ix][iy][iz+1] == TYPE_SOLVENT)     ||
         (grid->grid[ix][iy+1][iz-1] == TYPE_SOLVENT)   ||
         (grid->grid[ix][iy+1][iz] == TYPE_SOLVENT)     ||
         (grid->grid[ix][iy+1][iz+1] == TYPE_SOLVENT)   ||
         (grid->grid[ix+1][iy-1][iz-1] == TYPE_SOLVENT) ||
         (grid->grid[ix+1][iy-1][iz] == TYPE_SOLVENT)   ||
         (grid->grid[ix+1][iy-1][iz+1] == TYPE_SOLVENT) ||
         (grid->grid[ix+1][iy][iz-1] == TYPE_SOLVENT)   ||
         (grid->grid[ix+1][iy][iz] == TYPE_SOLVENT)     ||
         (grid->grid[ix+1][iy][iz+1] == TYPE_SOLVENT)   ||
         (grid->grid[ix+1][iy+1][iz-1] == TYPE_SOLVENT) ||
         (grid->grid[ix+1][iy+1][iz] == TYPE_SOLVENT)   ||
         (grid->grid[ix+1][iy+1][iz+1] == TYPE_SOLVENT))
      {
         if(!AtomNear(pdb, 
                      grid->xmin + ix*grid->gridStep + xoff, 
                      grid->ymin + iy*grid->gridStep + yoff, 
                      grid->zmin + iz*grid->gridStep + zoff, 
                      solvSize, grid->neighbours[ix][iy][iz]))
         {
            return(SOLVENT_BULK);
         }
         else
         {
            /* This point has solvent and protein neighbours            */
            return(SOLVENT_MARGIN);
         }
      }
      
      /* Not in solvent at all                                          */
      return(SOLVENT_NO);
   }
   else
   {
      /* It's on the boundary of the box so must be solvent             */
      return(SOLVENT_BULK);
   }
}


/************************************************************************/
/*>BOOL MarkSolventPoints(PDB *pdb, GRID *grid, REAL solventStep,
                          REAL solventMult, BOOL quiet, BOOL doNeighbours)
   -----------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Assign all 6 surfaces of the grid as solvent then work from the 6 
   surface reassigning points as being solvent if they are next to a 
   solvent point and are not too close to a protein atom. Also set any
   grid point within the atom radius to solvent type 2.

*** TODO: Could maybe speed this up by looking at the grid points      ***
*** rather than the actual atom coordinates to detect being too close  ***

   17.10.01 Original   By: ACRM
   18.06.02 Handles SOLVENT_MARGIN points
   17.07.02 Returns a flag to indicate whether any solvent points were
            modified
            Added doNeighbours parameter
*/
BOOL MarkSolventPoints(PDB *pdb, GRID *grid, REAL solventStep, 
                       REAL solventMult, BOOL quiet, BOOL doNeighbours)
{
   BOOL modified = TRUE,
        doneSomething = FALSE;
   int  ix, iy, iz, iteration = 1,solvState;
   int  numberOfSolvatedNeighbours;
   REAL solvRadiusSq = solventMult * solventStep * 
                       solventMult * solventStep,
        smallStep    = grid->gridStep / (2*FINEGRIDSTEP),
        xoff,
        yoff,
        zoff;
   
   
   numberOfSolvatedNeighbours = 
      ((doNeighbours)?
       (int)(solventMult * solventStep / grid->gridStep):
       0);
   
   /* Assume all points on the surface of the box are solvent           */
   FlagBoxBoundsAsSolvent(grid, grid->ixmax, grid->iymax, grid->izmax,
                          quiet); 
   FlagBoxBoundsAsSolvent(grid, 
                          grid->ixmax-1, grid->iymax-1, grid->izmax-1,
                          quiet);
   
   /* Then work from each of the 6 faces marking any point as solvent if 
      it has an adjacent solvent point. Somehow need to account for being 
      able to fit solvent molecules in through any passage to the outside 
      world
   */
   
   while(modified)
   {
      if(!gFlags.quiet)
      {
         fprintf(stderr,"Iteration %d...\n", iteration++);
      }
      
      modified = FALSE;
#pragma omp parallel for private(iy,iz,xoff,yoff,zoff) schedule (dynamic)
      for(ix=0; ix<grid->ixmax; ix++)
      {
         for(iy=0; iy<grid->iymax; iy++)
         {
            /* Forwards...                                              */
            for(iz=0; iz<grid->izmax; iz++)
            {
               if((grid->grid[ix][iy][iz] == TYPE_VOID) ||
                  (grid->grid[ix][iy][iz] == TYPE_SOLVENT2))
               {
                  if((solvState=AdjacentIsSolvent(pdb, grid, ix, iy, iz, 
                                                  solventStep,
                                                  0.0, 0.0, 0.0))
                     == SOLVENT_BULK)
                  {
                     SetAsSolvent(grid,ix,iy,iz,
                                  numberOfSolvatedNeighbours, 
                                  solvRadiusSq,
                                  TYPE_SOLVENT);
                     modified = TRUE;
                     doneSomething = TRUE;
                  }
                  else if(solvState == SOLVENT_MARGIN)
                  {
                     for(xoff = (-NSOLVFINESTEP*smallStep);
                         xoff <= (NSOLVFINESTEP*smallStep);
                         xoff += smallStep)
                     {
                        for(yoff = (-NSOLVFINESTEP*smallStep);
                            yoff <= (NSOLVFINESTEP*smallStep);
                            yoff += smallStep)
                        {
                           if(AdjacentIsSolvent(pdb, grid, ix, iy, iz,
                                                solventStep,
                                                xoff, yoff, 0.0)
                              == SOLVENT_BULK)
                           {
                              SetAsSolvent(grid,ix,iy,iz,
                                           numberOfSolvatedNeighbours, 
                                           solvRadiusSq,
                                           TYPE_SOLVENT);
                              modified = TRUE;
                              doneSomething = TRUE;
                           }
                        }
                     }
                  }
               }
            }
            /* Backwards...                                             */
            for(iz=grid->izmax-1; iz>=0; iz--)
            {
               if((grid->grid[ix][iy][iz] == TYPE_VOID) ||
                  (grid->grid[ix][iy][iz] == TYPE_SOLVENT2))
               {
                  if((solvState=AdjacentIsSolvent(pdb, grid, ix, iy, iz, 
                                                  solventStep,
                                                  0.0, 0.0, 0.0))
                     == SOLVENT_BULK)
                  {
                     SetAsSolvent(grid,ix,iy,iz,
                                  numberOfSolvatedNeighbours, 
                                  solvRadiusSq,
                                  TYPE_SOLVENT);
                     modified = TRUE;
                     doneSomething = TRUE;
                  }
                  else if(solvState == SOLVENT_MARGIN)
                  {
                     for(xoff = (-NSOLVFINESTEP*smallStep);
                         xoff <= (NSOLVFINESTEP*smallStep);
                         xoff += smallStep)
                     {
                        for(yoff = (-NSOLVFINESTEP*smallStep);
                            yoff <= (NSOLVFINESTEP*smallStep);
                            yoff += smallStep)
                        {
                           if(AdjacentIsSolvent(pdb, grid, ix, iy, iz,
                                                solventStep,
                                                xoff, yoff, 0.0)
                              == SOLVENT_BULK)
                           {
                              SetAsSolvent(grid,ix,iy,iz,
                                           numberOfSolvatedNeighbours, 
                                           solvRadiusSq,
                                           TYPE_SOLVENT);
                              modified = TRUE;
                              doneSomething = TRUE;
                           }
                        }
                     }
                  }
               }
            }
         }
         
         for(iz=0; iz<grid->izmax; iz++)
         {
            /* Up...                                                    */
            for(iy=0; iy<grid->iymax; iy++)
            {
               if((grid->grid[ix][iy][iz] == TYPE_VOID) ||
                  (grid->grid[ix][iy][iz] == TYPE_SOLVENT2))
               {
                  if((solvState=AdjacentIsSolvent(pdb, grid, ix, iy, iz, 
                                                  solventStep,
                                                  0.0, 0.0, 0.0))
                     == SOLVENT_BULK)
                  {
                     SetAsSolvent(grid,ix,iy,iz,
                                  numberOfSolvatedNeighbours, 
                                  solvRadiusSq,
                                  TYPE_SOLVENT);
                     modified = TRUE;
                     doneSomething = TRUE;
                  }
                  else if(solvState == SOLVENT_MARGIN)
                  {
                     for(xoff = (-NSOLVFINESTEP*smallStep);
                         xoff <= (NSOLVFINESTEP*smallStep);
                         xoff += smallStep)
                     {
                        for(zoff = (-NSOLVFINESTEP*smallStep);
                            zoff <= (NSOLVFINESTEP*smallStep);
                            zoff += smallStep)
                        {
                           if(AdjacentIsSolvent(pdb, grid, ix, iy, iz,
                                                solventStep,
                                                xoff, 0.0, zoff)
                              == SOLVENT_BULK)
                           {
                              SetAsSolvent(grid,ix,iy,iz,
                                           numberOfSolvatedNeighbours, 
                                           solvRadiusSq,
                                           TYPE_SOLVENT);
                              modified = TRUE;
                              doneSomething = TRUE;
                           }
                        }
                     }
                  }
               }
            }
            /* Down...                                                  */
            for(iy=grid->iymax-1; iy>=0; iy--)
            {
               if((grid->grid[ix][iy][iz] == TYPE_VOID) ||
                  (grid->grid[ix][iy][iz] == TYPE_SOLVENT2))
               {
                  if((solvState=AdjacentIsSolvent(pdb, grid, ix, iy, iz, 
                                                  solventStep,
                                                  0.0, 0.0, 0.0))
                     == SOLVENT_BULK)
                  {
                     SetAsSolvent(grid,ix,iy,iz,
                                  numberOfSolvatedNeighbours, 
                                  solvRadiusSq,
                                  TYPE_SOLVENT);
                     modified = TRUE;
                     doneSomething = TRUE;
                  }
                  else if(solvState == SOLVENT_MARGIN)
                  {
                     for(xoff = (-NSOLVFINESTEP*smallStep);
                         xoff <= (NSOLVFINESTEP*smallStep);
                         xoff += smallStep)
                     {
                        for(zoff = (-NSOLVFINESTEP*smallStep);
                            zoff <= (NSOLVFINESTEP*smallStep);
                            zoff += smallStep)
                        {
                           if(AdjacentIsSolvent(pdb, grid, ix, iy, iz,
                                                solventStep,
                                                xoff, 0.0, zoff)
                              == SOLVENT_BULK)
                           {
                              SetAsSolvent(grid,ix,iy,iz,
                                           numberOfSolvatedNeighbours, 
                                           solvRadiusSq,
                                           TYPE_SOLVENT);
                              modified = TRUE;
                              doneSomething = TRUE;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      
#pragma omp parallel for private(ix,iz,xoff,yoff,zoff) schedule (dynamic)
      for(iy=0; iy<grid->iymax; iy++)
      {
         for(iz=0; iz<grid->izmax; iz++)
         {
            /* Across...                                                */
            for(ix=0; ix<grid->ixmax; ix++)
            {
               if((grid->grid[ix][iy][iz] == TYPE_VOID) ||
                  (grid->grid[ix][iy][iz] == TYPE_SOLVENT2))
               {
                  if((solvState=AdjacentIsSolvent(pdb, grid, ix, iy, iz, 
                                                  solventStep,
                                                  0.0, 0.0, 0.0))
                     == SOLVENT_BULK)
                  {
                     SetAsSolvent(grid,ix,iy,iz,
                                  numberOfSolvatedNeighbours, 
                                  solvRadiusSq,
                                  TYPE_SOLVENT);
                     modified = TRUE;
                     doneSomething = TRUE;
                  }
                  else if(solvState == SOLVENT_MARGIN)
                  {
                     for(yoff = (-NSOLVFINESTEP*smallStep);
                         yoff <= (NSOLVFINESTEP*smallStep);
                         yoff += smallStep)
                     {
                        for(zoff = (-NSOLVFINESTEP*smallStep);
                            zoff <= (NSOLVFINESTEP*smallStep);
                            zoff += smallStep)
                        {
                           if(AdjacentIsSolvent(pdb, grid, ix, iy, iz,
                                                solventStep,
                                                0.0, yoff, zoff)
                              == SOLVENT_BULK)
                           {
                              SetAsSolvent(grid,ix,iy,iz,
                                           numberOfSolvatedNeighbours, 
                                           solvRadiusSq,
                                           TYPE_SOLVENT);
                              modified = TRUE;
                              doneSomething = TRUE;
                           }
                        }
                     }
                  }
               }
            }
            /* Backwards...                                             */
            for(ix=grid->ixmax-1; ix>=0; ix--)
            {
               if((grid->grid[ix][iy][iz] == TYPE_VOID) ||
                  (grid->grid[ix][iy][iz] == TYPE_SOLVENT2))
               {
                  if((solvState=AdjacentIsSolvent(pdb, grid, ix, iy, iz, 
                                                  solventStep,
                                                  0.0, 0.0, 0.0))
                     == SOLVENT_BULK)
                  {
                     SetAsSolvent(grid,ix,iy,iz,
                                  numberOfSolvatedNeighbours, 
                                  solvRadiusSq,
                                  TYPE_SOLVENT);
                     modified = TRUE;
                     doneSomething = TRUE;
                  }
                  else if(solvState == SOLVENT_MARGIN)
                  {
                     for(yoff = (-NSOLVFINESTEP*smallStep);
                         yoff <= (NSOLVFINESTEP*smallStep);
                         yoff += smallStep)
                     {
                        for(zoff = (-NSOLVFINESTEP*smallStep);
                            zoff <= (NSOLVFINESTEP*smallStep);
                            zoff += smallStep)
                        {
                           if(AdjacentIsSolvent(pdb, grid, ix, iy, iz,
                                                solventStep,
                                                0.0, yoff, zoff)
                              == SOLVENT_BULK)
                           {
                              SetAsSolvent(grid,ix,iy,iz,
                                           numberOfSolvatedNeighbours, 
                                           solvRadiusSq,
                                           TYPE_SOLVENT);
                              modified = TRUE;
                              doneSomething = TRUE;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }  /* End iterations                                                 */

   return(doneSomething);
}

#endif


/************************************************************************/
/*>void FreeGrid(GRID *grid)
   -------------------------
   Input:     
   Output:    
   Returns:   

   Frees malloc'd memory within the grid structure

   31.10.01 Original   By: ACRM
*/
void FreeGrid(GRID *grid)
{
   int i, j, k;
   
   if(grid->grid!=NULL)
   {
      for(i=0; i<grid->ixmax; i++)
      {
         if(grid->grid[i] != NULL)
         {
            for(j=0; j<grid->iymax; j++)
            {
               if(grid->grid[i][j] != NULL)
               {
                  free(grid->grid[i][j]);
                  grid->grid[i][j] = NULL;
               }
            }
            free(grid->grid[i]);
            grid->grid[i] = NULL;
         }
      }
      free(grid->grid);
      grid->grid = NULL;
   }
   
   if(grid->neighbours!=NULL)
   {
      for(i=0; i<grid->ixmax; i++)
      {
         if(grid->neighbours[i] != NULL)
         {
            for(j=0; j<grid->iymax; j++)
            {
               if(grid->neighbours[i][j] != NULL)
               {
                  for(k=0; k<grid->izmax; k++)
                  {
                     if(grid->neighbours[i][j][k] != NULL)
                     {
                        free(grid->neighbours[i][j][k]);
                        grid->neighbours[i][j][k] = NULL;
                     }
                  }
                  free(grid->neighbours[i][j]);
                  grid->neighbours[i][j] = NULL;
               }
            }
            free(grid->neighbours[i]);
            grid->neighbours[i] = NULL;
         }
      }
      free(grid->neighbours);
      grid->neighbours = NULL;
   }
}


/************************************************************************/
/*>BOOL RefineVoidVolumes(VOIDS *voidlist, GRID *grid, PDB *pdb, 
                          REAL gridStep, REAL probeSize, REAL solvSize)
   --------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Refines the void volumes. Calls RefineVoxel() on each voxel of a void; 
   this splits it into FINEGRIDSTEP^3 (1000) sub-voxels and calculates
   which are really void.

   Then calls RefineVoxelNeighbours() on each voxel. This identifies
   any void sub-voxels within neighbouring voxels assigned as protein

   31.10.01 Original   By: ACRM
*/
BOOL RefineVoidVolumes(VOIDS *voidlist, GRID *grid, PDB *pdb, 
                       REAL gridStep, REAL probeSize, REAL solvSize)
{
   VOIDS     *v;
   int       vcount  = 0,
             nvoids  = 0;
   REAL      voxelVolume,
             fineGridStep,
             fineProbeSize;
   POINTLIST *pt;
   
   fineGridStep = grid->gridStep / FINEGRIDSTEP;
   voxelVolume = fineGridStep * fineGridStep * fineGridStep;

   /* Calculate a probe size for using with the sub-voxels. This will
      be fineGridStep/courseGridStep of the normal probe size
   */
   fineProbeSize = probeSize * fineGridStep / grid->gridStep;
   
   for(v=voidlist; v!=NULL; NEXT(v))
   {
      vcount = 0;
      if(!gFlags.quiet)
      {
         fprintf(stderr, "   refining void %d...\n", ++nvoids);
      }
      
      /* Run through the points that make up this void                  */
      for(pt=v->pointlist; pt!=NULL; NEXT(pt))
      {
         vcount += RefineVoxel(grid, pdb, pt, fineGridStep,
                               fineProbeSize);
      }
#ifdef DEBUG
      fprintf(stderr,"      Stage 1 refinement volume: %f\n",
              vcount * voxelVolume);
#endif
      
      for(pt=v->pointlist; pt!=NULL; NEXT(pt))
      {
         vcount += RefineVoxelNeighbours(grid, pdb, pt, fineGridStep,
                                         fineProbeSize);
      }

      /* Calculate the new volume for this void                         */
      v->volume = vcount * voxelVolume;
   }
   return(TRUE);
}


/************************************************************************/
/*>int RefineVoxel(GRID *grid, PDB *pdb, POINTLIST *pt, REAL fineGridStep,
                   REAL fineProbeSize)
   -----------------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Splits a voxel into FINGRIDSTEP^3 (1000) subvoxels and checks each to
   see if it is protein or void

   31.10.01 Original   By: ACRM
   18.06.02 Changed probeSize to fineProbeSize
*/
int RefineVoxel(GRID *grid, PDB *pdb, POINTLIST *pt, REAL fineGridStep,
                REAL fineProbeSize)
{
   REAL  xmin, xmax,
         ymin, ymax,
         zmin, zmax,
         cutoffSq;
   int   vcount = 0,
         i;
   static int atomnum = 1;
   BOOL  isVoid;
   VEC3F gridCoor;

   xmin = pt->x - grid->gridStep / 2.0;
   ymin = pt->y - grid->gridStep / 2.0;
   zmin = pt->z - grid->gridStep / 2.0;
   xmax = pt->x + grid->gridStep / 2.0;
   ymax = pt->y + grid->gridStep / 2.0;
   zmax = pt->z + grid->gridStep / 2.0;
   
   for(gridCoor.x=xmin+(fineGridStep/2); 
       gridCoor.x<xmax; 
       gridCoor.x+=fineGridStep)
   {
      for(gridCoor.y=ymin+(fineGridStep/2); 
          gridCoor.y<ymax; 
          gridCoor.y+=fineGridStep)
      {
         for(gridCoor.z=zmin+(fineGridStep/2); 
             gridCoor.z<zmax; 
             gridCoor.z+=fineGridStep)
         {
            /* Simply work through the list of atoms associated with
               this current (course) grid point to look whether this
               grid point is void
            */
            isVoid = TRUE;
            for(i=0; i<MAXNEIGHBOURS; i++)
            {
               PDB *p = grid->neighbours[pt->ix][pt->iy][pt->iz][i];
               
               if(p==NULL)
               {
                  break;
               }
               
               cutoffSq = (fineProbeSize+p->occ) * 
                          (fineProbeSize+p->occ); 
               if(DISTSQ(&gridCoor, p) < cutoffSq)
               {
                  isVoid = FALSE;
                  break;
               }
            }
            if(isVoid)
            {
               vcount++;

               if(gFlags.printRefinedVoids)
               {
                  fprintf(gFlags.gridRefinedFile, "ATOM  %5d  O   \
GLY Z%4d    %8.3f%8.3f%8.3f  1.00 20.00\n",  
                          atomnum++, 1, 
                          gridCoor.x, gridCoor.y, gridCoor.z);
               }
            }
         }
      }
   }
   
   return(vcount);
}

/************************************************************************/
/*>int RefineVoxelNeighbours(GRID *grid, PDB *pdb, POINTLIST *pt, 
                             REAL fineGridStep, REAL fineProbeSize)
   ----------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Looks at each of the 6 planar-neighbouring voxels to add any void
   sub-voxels they might contain

   31.10.01 Original   By: ACRM
   17.06.02 Added edge and corner neighbours
   18.06.02 Changed probeSize to fineProbeSize
*/
int RefineVoxelNeighbours(GRID *grid, PDB *pdb, POINTLIST *pt, 
                          REAL fineGridStep, REAL fineProbeSize)
{
   int vcount = 0;
   
   /* Look at the 6 (plane) neighbours of this point                    */
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  1,  0,  0, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt, -1,  0,  0, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  0,  1,  0, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  0, -1,  0, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  0,  0,  1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  0,  0, -1, 
                                   fineGridStep, fineProbeSize);
   
#ifndef NO_EXTRA_VOXEL_NEIGHBOURS
   /* Look at the 12 edge-neighbours                                    */
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  -1,  -1,  0, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  -1,  0,  -1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  -1,  1,  0, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  -1,  0,  1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  0, -1, -1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  0, -1,  1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  0,  1, -1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  0,  1,  1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  1, -1,  0, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  1,  0, -1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  1,  0,  1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  1,  1,  0, 
                                   fineGridStep, fineProbeSize);

   /* Look at the 8 corner neighbours                                   */
   vcount += RefineAVoxelNeighbour(grid, pdb, pt, -1, -1, -1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt, -1, -1,  1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt, -1,  1, -1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt, -1,  1,  1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  1, -1, -1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  1, -1,  1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  1,  1, -1, 
                                   fineGridStep, fineProbeSize);
   vcount += RefineAVoxelNeighbour(grid, pdb, pt,  1,  1,  1, 
                                   fineGridStep, fineProbeSize);
#endif

   return(vcount);
}

/************************************************************************/
/*>int RefineAVoxelNeighbour(GRID *grid, PDB *pdb, POINTLIST *pt,  
                             int xoff, int yoff, int zoff, 
                             REAL fineGridStep, REAL fineProbeSize)
   ----------------------------------------------------------------
   Input:     
   Output:    
   Returns:   

   Takes a voxel which is a neighbour of a void voxel and tests whether
   it is protein. If so, calls RefineVoxel() on this to see if some of it
   should be void. The voxel type is changed from TYPE_PROTEIN to 
   TYPE_PROTEIN2 to avoid counting the same neighbour more than once.

   If gFlags.leak is not set, then none of the neighbours of this voxel
   is allowed to be a solvent voxel.

   31.10.01 Original   By: ACRM
   18.06.02 Changed probeSize to fineProbeSize
*/
int RefineAVoxelNeighbour(GRID *grid, PDB *pdb, POINTLIST *pt,  
                          int xoff, int yoff, int zoff, 
                          REAL fineGridStep, REAL fineProbeSize)
{
   int ix, iy, iz, vcount=0;
   POINTLIST npt;
   
   ix = pt->ix + xoff;
   iy = pt->iy + yoff;
   iz = pt->iz + zoff;

   if((ix >= 0) && (ix < grid->ixmax) &&
      (iy >= 0) && (iy < grid->iymax) &&
      (iz >= 0) && (iz < grid->izmax) &&
      (grid->grid[ix][iy][iz] == TYPE_PROTEIN))
   {
      /* This is a neighbouring protein point.                          */

      /* Stop the current voxel from getting counted twice              */
      grid->grid[ix][iy][iz] = TYPE_PROTEIN2;
      
      /* Check that none of its planar neighbours is solvent            
         Note that we do no bounds checking as no protein point should be
         at the edge of the box
       */
      if(gFlags.leak ||
         ((grid->grid[ix-1][iy][iz] != TYPE_SOLVENT)  &&
          (grid->grid[ix-1][iy][iz] != TYPE_SOLVENT2) &&
          (grid->grid[ix+1][iy][iz] != TYPE_SOLVENT)  &&
          (grid->grid[ix+1][iy][iz] != TYPE_SOLVENT2) &&
          (grid->grid[ix][iy-1][iz] != TYPE_SOLVENT)  &&
          (grid->grid[ix][iy-1][iz] != TYPE_SOLVENT2) &&
          (grid->grid[ix][iy+1][iz] != TYPE_SOLVENT)  &&
          (grid->grid[ix][iy+1][iz] != TYPE_SOLVENT2) &&
          (grid->grid[ix][iy][iz-1] != TYPE_SOLVENT)  &&
          (grid->grid[ix][iy][iz-1] != TYPE_SOLVENT2) &&
          (grid->grid[ix][iy][iz+1] != TYPE_SOLVENT)  &&
          (grid->grid[ix][iy][iz+1] != TYPE_SOLVENT2)))
      {
         /* Fill in a POINTLIST entry for this neighbouring voxel       */
         npt.next    = NULL;
         npt.nearest = NULL;
         npt.x       = grid->xmin + ix * grid->gridStep;
         npt.y       = grid->ymin + iy * grid->gridStep;
         npt.z       = grid->zmin + iz * grid->gridStep;
         npt.ix      = ix;
         npt.iy      = iy;
         npt.iz      = iz;
         
         vcount = RefineVoxel(grid, pdb, &npt, fineGridStep, 
                              fineProbeSize);
      }
   }
   
   return(vcount);
}

/************************************************************************/
void PrintWaterNeighbours(FILE *out, PDB *pdb, GRID *grid)
{
   PDB  *p;
   int  ix, iy, iz,
        grid_x, grid_y, grid_z,
        nSteps;
   REAL maxDist, maxDistSq;

   
   maxDist   = WATER_RADIUS + grid->gridStep;
   maxDistSq = maxDist * maxDist;
   nSteps    = 2 * ((maxDist / grid->gridStep) + 1);

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->record_type, "HETATM", 6))
      {
         if(ISWATER(p))
         {
            /* Find the grid point nearest to this water atom           */
            grid_x = (int)(0.5 + ((p->x - grid->xmin)/grid->gridStep));
            grid_y = (int)(0.5 + ((p->y - grid->ymin)/grid->gridStep));
            grid_z = (int)(0.5 + ((p->z - grid->zmin)/grid->gridStep));

            for(ix = grid_x - nSteps; 
                ix <= grid_x + nSteps; 
                ix++)
            {
               if((ix >= 0) && (ix < grid->ixmax))
               {
                  for(iy = grid_y - nSteps; 
                      iy <= grid_y + nSteps; 
                      iy++)
                  {
                     if((iy >= 0) && (iy < grid->iymax))
                     {
                        for(iz = grid_z - nSteps; 
                            iz <= grid_z + nSteps; 
                            iz++)
                        {
                           if((iz >= 0) && (iz < grid->izmax))
                           {
                              VEC3F point;

                              if((grid->grid[ix][iy][iz] == TYPE_VOID) ||
                                 ISSET(grid->grid[ix][iy][iz], TYPE_ASSIGNED))
                              {
                                 point.x = grid->xmin+(ix*grid->gridStep);
                                 point.y = grid->ymin+(iy*grid->gridStep);
                                 point.z = grid->zmin+(iz*grid->gridStep);
                              
                                 if(DISTSQ(p,(&point)) < maxDistSq)
                                 {
                                    fprintf(out, "Water %d (atom %d) is \
touching a void\n", p->resnum, p->atnum);
                                    goto nextatom;
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
nextatom: 
      ;
   }
}


/************************************************************************/
/*>BOOL IsSolvent(GRID *grid, int ix, int iy, int iz)
   --------------------------------------------------
   Input:     GRID   *grid       Grid structure
              int    ix,iy,iz    Position on grid
   Returns:   BOOL               Is this a solvent point?

   Tests whether a specified point is on the grid and if so is it
   solvent?

   02.07.02 Original   By: ACRM
*/
BOOL IsSolvent(GRID *grid, int ix, int iy, int iz)
{
   if((ix>=0) && (ix<grid->ixmax) &&
      (iy>=0) && (iy<grid->iymax) &&
      (iz>=0) && (iz<grid->izmax))
   {
      if((grid->grid[ix][iy][iz] == TYPE_SOLVENT) ||
         (grid->grid[ix][iy][iz] == TYPE_SOLVENT2))
      {
         return(TRUE);
      }
   }
   
   return(FALSE);
}



/************************************************************************/
/*>void ReassignVoidPoints(GRID *grid, REAL probeSize)
   ---------------------------------------------------
   Checks all protein voxels to see if they contain a void that isn't 
   found by placing the test probe at the centre of the voxel
*/
void ReassignVoidPoints(GRID *grid, REAL probeSize)
{
   REAL x,  y,  z;
   int  xi, yi, zi,
        numNeighbours;
   PDB  *neighbours[MAXNEIGHBOURS];
   
   
   for(x=grid->xmin, xi=0; x<=grid->xmax; x+=grid->gridStep, xi++)
   {
      for(y=grid->ymin, yi=0; y<=grid->ymax; y+=grid->gridStep, yi++)
      {
         for(z=grid->zmin, zi=0; z<=grid->zmax; z+=grid->gridStep, zi++)
         {
            /* If this point has been assigned as protein               */
            if(grid->grid[xi][yi][zi] == TYPE_PROTEIN)
            {
               if(NeighboursAreProtein(grid,xi,yi,zi))
               {
                  /* Create a combined list of all the atoms neighbouring 
                     this protein point and all adjacent protein points
                  */
                  numNeighbours = CombineAllNeighbours(grid, xi, yi, zi, 
                                                       neighbours);
                  doReassignVoidPoint(grid, xi, yi, zi, neighbours,
                                      numNeighbours, probeSize);
               }
            }
         }
      }
   }
   
}


/************************************************************************/
/*
  the occupancy field stores the radii
 */
void doReassignVoidPoint(GRID *grid, int xi, int yi, int zi, 
                         PDB **neighbours, int numNeighbours,
                         REAL probeSize)
{
   int   i, j, k;
   PDB   *p, *q, *r;
   REAL  minDistSq,
         distSq,
         dist,
         d,
         ps;
   VEC3F pq,
         pt,
         pV,
         corner1,
         corner2;
   BOOL  isVoid;
#ifdef DEBUG
   REAL  nearestAtomDist;
   PDB   *nearestAtom = NULL;
#endif
   
   /* Enforce a minimum probe size of 0.1A - otherwise virtually every 
      protein point will get converted back to void
   */
   ps = MAX(probeSize, (REAL)0.1);

   /* Calculate the bounding corners of this voxel                      */
   corner1.x = (grid->xmin + (((REAL)xi-0.5) * grid->gridStep));
   corner1.y = (grid->ymin + (((REAL)yi-0.5) * grid->gridStep));
   corner1.z = (grid->zmin + (((REAL)zi-0.5) * grid->gridStep));
   corner2.x = (grid->xmin + (((REAL)xi+0.5) * grid->gridStep));
   corner2.y = (grid->ymin + (((REAL)yi+0.5) * grid->gridStep));
   corner2.z = (grid->zmin + (((REAL)zi+0.5) * grid->gridStep));

   /* Run through all the neighbouring atoms                            */
   for(i=0; i<numNeighbours; i++)
   {
      p = neighbours[i];
      pV.x = p->x;
      pV.y = p->y;
      pV.z = p->z;
         
      for(j=i+1; j<numNeighbours; j++)
      {
         q = neighbours[j];

         /* Calculate the sum of the two atom radii plus the diameter of
            the probe. An atom pair must be separated by at least this
            amount for us to be interested
         */
         minDistSq = (p->occ + q->occ + (ps * 2.0)) *
                     (p->occ + q->occ + (ps * 2.0));
         /* Calculate the actual separation of these two atoms          */
         distSq = DISTSQ(p, q);

         /* If these two atoms are separated by more than the sum of their
            radii plus the diameter of the probe, then walk along the 
            vector between then and see if any point is at least a probe 
            radius away from all other atoms
         */
         if(distSq > minDistSq)
         {
            /* Calculate the distance between the two points (i.e. the
               length of the vector between them)
            */
            dist = (REAL)sqrt((double)distSq);
            /* Now define the unit vector                               */
            pq.x = (q->x - p->x) / dist;
            pq.y = (q->y - p->y) / dist;
            pq.z = (q->z - p->z) / dist;

            /* Now see if the vector intersects the voxel of interest   */
            if(DoesLineCrossBox(pV, pq, dist, corner1, corner2))
            {
               /* Now make steps along the line in 0.05A units          */
               for(d = p->occ+ps; d <= (dist - q->occ - ps ); d += 0.05)
               {
                  pt.x = p->x + (d*pq.x);
                  pt.y = p->y + (d*pq.y);
                  pt.z = p->z + (d*pq.z);
                  
                  /* If this point is within the current voxel          */
/*
//                  if((pt.x >= (grid->xmin + 
//                               (((REAL)xi-0.5) * grid->gridStep))) &&
//                     (pt.x <  (grid->xmin + 
//                               (((REAL)xi+0.5) * grid->gridStep))) &&
//                     (pt.y >= (grid->ymin + 
//                               (((REAL)yi-0.5) * grid->gridStep))) &&
//                     (pt.y <  (grid->ymin + 
//                               (((REAL)yi+0.5) * grid->gridStep))) &&
//                     (pt.z >= (grid->zmin + 
//                               (((REAL)zi-0.5) * grid->gridStep))) &&
//                     (pt.z <  (grid->zmin + 
//                               (((REAL)zi+0.5) * grid->gridStep))))
*/
                  {
                     /* Check the distances to all neighbouring atoms. 
                        Assume that if all atoms are sufficiently distant
                        that there is actually a void point along this 
                        vector, once that is shown not to be true, set 
                        isVoid to FALSE and break out of the loop
                     */
                     isVoid = TRUE;
#ifdef DEBUG
                     nearestAtomDist = 100000.0;
#endif
                     for(k=0; k<numNeighbours; k++)
                     {
                        r = neighbours[k];
                        minDistSq = (r->occ + ps) * (r->occ + ps);
                        distSq = DISTSQ(r, &pt);
                        if(distSq < minDistSq)
                        {
                           isVoid = FALSE;
                           break;
                        }
#ifdef DEBUG
                        if(distSq < nearestAtomDist)
                        {
                           nearestAtomDist = distSq;
                           nearestAtom = r;
                        }
#endif
                     }
                     
                     /* We have found a point along this line which is a 
                        void point. Reassign the current grid point from 
                        protein to void and return.
                     */
                     if(isVoid)
                     {
/*  This gives worse results than just changing the current voxel
//                        ChangeVoxel(grid, xi, yi, zi, pt);
*/                      

                      grid->grid[xi][yi][zi] = TYPE_VOID;

#ifdef DEBUG
                        fprintf(stderr,"\nAtoms %d and %d are separated \
by %.1fA (required separation %.1fA)\n", 
                                p->atnum, q->atnum, dist, 
                                sqrt(minDistSq));
                        fprintf(stderr,"At point %8.3f%8.3f%8.3f\n",
                                pt.x, pt.y, pt.z);
                        fprintf(stderr, "Nearest atom distance was %f to \
atom %d\n",
                                sqrt((double)nearestAtomDist), 
                                nearestAtom->atnum);
                  fprintf(stdout, "ATOM  %5d  CA  GLY Z   \
1    %8.3f%8.3f%8.3f  1.00 20.00\n",
                          1, pt.x,pt.y,pt.z);
#endif
                        return; 
                     }
                  }
               }
            }
         }
      }
   }
}




/************************************************************************/
int CombineAllNeighbours(GRID *grid, int xi, int yi, int zi, 
                         PDB **neighbours)
{
   int  numNeighbours = 0,
        i, j, k, m, n;
   BOOL found;

   for(i=xi-1; i<=xi+1; i++)
   {
      if((i>=0) && (i<grid->ixmax))
      {
         for(j=yi-1; j<=yi+1; j++)
         {
            if((j>=0) && (j<grid->iymax))
            {
               for(k=zi-1; k<=zi+1; k++)
               {
                  if((k>=0) && (k<grid->izmax))
                  {
                     /* Step through the neighbours listed for this 
                        grid point
                     */
                     for(n=0; 
                         (n < MAXNEIGHBOURS) && 
                            (grid->neighbours[i][j][k][n] != NULL);
                         n++)
                     {
                        /* See if this neighbour is already in our
                           combined list
                        */
                        found = FALSE;
                        for(m=0; m<numNeighbours; m++)
                        {
                           if(grid->neighbours[i][j][k][n] ==
                              neighbours[m])
                           {
                              found = TRUE;
                              break;
                           }
                        }
                        /* If not in the combined list, add it          */
                        if(!found)
                        {
                           if(numNeighbours >= MAXNEIGHBOURS)
                           {
                              fprintf(stderr,"Error!: Maximum number of \
neighbours exceeded in combining neighbour lists.\n\
Increase MAXNEIGHBOURS!\n");
                              exit(1);
                           }
                           neighbours[numNeighbours++] =
                              grid->neighbours[i][j][k][n];
                        }
                     }  /* Loop through neighbours                      */
                  }
               }  /* Loop through z                                     */
            }
         }  /* Loop through y                                           */
      }
   }  /* Loop through x                                                 */

   return(numNeighbours);
}


#define NUMGRIDNEIGHBOURS2 3
BOOL NeighboursAreProtein(GRID *grid, int ix, int iy, int iz)
{
   int i, j, k;

   for(i=ix-NUMGRIDNEIGHBOURS2; i<=ix+NUMGRIDNEIGHBOURS2; i++)
   {
      if((i>=0) && (i<grid->ixmax))
      {
         for(j=iy-NUMGRIDNEIGHBOURS2; j<=iy+NUMGRIDNEIGHBOURS2; j++)
         {
            if((j>=0) && (j<grid->iymax))
            {
               for(k=iz-NUMGRIDNEIGHBOURS2; k<=iz+NUMGRIDNEIGHBOURS2; k++)
               {
                  if((k>=0) && (k<grid->izmax))
                  {
                     if(grid->grid[i][j][k] != TYPE_PROTEIN)
                     {
                        return(FALSE);
                     }
                  }
                  else
                  {
                     return(FALSE);
                  }
               }
            }
            else
            {
               return(FALSE);
            }
         }
      }
      else
      {
         return(FALSE);
      }
   }

   return(TRUE);
}


void ChangeVoxel(GRID *grid, int ix, int iy, int iz, VEC3F pt)
{
   int   i, j, k,
         iNear = ix, 
         jNear = iy, 
         kNear = iz;
   REAL  minDistSq = 1000,
         distSq = 0.0;
   VEC3F voxel;
   
   for(i=ix-1; i<=ix+1; i++)
   {
      if((i>=0) && (i<grid->ixmax))
      {
         for(j=iy-1; j<=iy+1; j++)
         {
            if((j>=0) && (j<grid->iymax))
            {
               for(k=iz-1; k<=iz+1; k++)
               {
                  if((k>=0) && (k<grid->izmax))
                  {
                     voxel.x = grid->xmin + ix*grid->gridStep;
                     voxel.y = grid->ymin + iy*grid->gridStep;
                     voxel.z = grid->zmin + iz*grid->gridStep;
                     
                     distSq = DISTSQ(&voxel, &pt);
                     if(distSq < minDistSq)
                     {
                        minDistSq = distSq;
                        iNear = i;
                        jNear = j;
                        kNear = k;
                     }
                  }
               }
            }
         }
      }
   }

   grid->grid[iNear][jNear][kNear] = TYPE_VOID;
}



/************************************************************************/
/*>BOOL ReassignSurfaceProteinVoxels(PDB *pdb, GRID *grid, REAL solvSize)
   ----------------------------------------------------------------------

   Looks at surface protein voxels and reassigns them to solvent if
   a solvent molecule can be placed at any point along the vector
   between the centre of this voxel and the centre of an adjacent
   voxel in the same plane perpendicular to the direction of the
   walk.
*/
BOOL ReassignSurfaceProteinVoxels(PDB *pdb, GRID *grid, REAL solvSize)
{
   BOOL modified = FALSE;
   int  ix, iy, iz;
   
   
   /* Walk across the grid from each of the six sides                   */
   for(ix=0; ix<grid->ixmax; ix++)
   {
      for(iy=0; iy<grid->iymax; iy++)
      {
         /* Forwards...                                                 */
         for(iz=1; iz<grid->izmax; iz++)
         {
            /* If this is a protein point                               */
            if(grid->grid[ix][iy][iz] == TYPE_PROTEIN)
            {
               /* If previous point was solvent                         */
               if((grid->grid[ix][iy][iz-1] == TYPE_SOLVENT) ||
                  (grid->grid[ix][iy][iz-1] == TYPE_SOLVENT2))
               {
                  if(CheckAndChangeNeighbours(pdb, grid, 
                                              ix, iy, iz, 
                                              solvSize,
                                              1, 1, 0))
                  {
                     modified = TRUE;
                  }
               }
               break;
            }
         }
         /* Backwards...                                                */
         for(iz=grid->izmax-2; iz>=0; iz--)
         {
            /* If this is a protein point                               */
            if(grid->grid[ix][iy][iz] == TYPE_PROTEIN)
            {
               /* If previous point was solvent                         */
               if((grid->grid[ix][iy][iz+1] == TYPE_SOLVENT) ||
                  (grid->grid[ix][iy][iz+1] == TYPE_SOLVENT2))
               {
                  if(CheckAndChangeNeighbours(pdb, grid, 
                                              ix, iy, iz, 
                                              solvSize,
                                              1, 1, 0))
                  {
                     modified = TRUE;
                  }
               }
               break;
            }
         }
      }
      
      for(iz=0; iz<grid->izmax; iz++)
      {
         /* Up...                                                       */
         for(iy=1; iy<grid->iymax; iy++)
         {
            /* If this is a protein point                               */
            if(grid->grid[ix][iy][iz] == TYPE_PROTEIN)
            {
               /* If previous point was solvent                         */
               if((grid->grid[ix][iy-1][iz] == TYPE_SOLVENT) ||
                  (grid->grid[ix][iy-1][iz] == TYPE_SOLVENT2))
               {
                  if(CheckAndChangeNeighbours(pdb, grid, 
                                              ix, iy, iz, 
                                              solvSize,
                                              1, 0, 1))
                  {
                     modified = TRUE;
                  }
               }
               break;
            }
         }
         
         /* Down...                                                     */
         for(iy=grid->iymax-2; iy>=0; iy--)
         {
            /* If this is a protein point                               */
            if(grid->grid[ix][iy][iz] == TYPE_PROTEIN)
            {
               /* If previous point was solvent                         */
               if((grid->grid[ix][iy+1][iz] == TYPE_SOLVENT) ||
                  (grid->grid[ix][iy+1][iz] == TYPE_SOLVENT2))
               {
                  if(CheckAndChangeNeighbours(pdb, grid, 
                                              ix, iy, iz, 
                                              solvSize,
                                              1, 0, 1))
                  {
                     modified = TRUE;
                  }
               }
               break;
            }
         }
      }
   }

   for(iy=0; iy<grid->iymax; iy++)
   {
      for(iz=0; iz<grid->izmax; iz++)
      {
         /* Across...                                                   */
         for(ix=1; ix<grid->ixmax; ix++)
         {
            /* If this is a protein point                               */
            if(grid->grid[ix][iy][iz] == TYPE_PROTEIN)
            {
               /* If previous point was solvent                         */
               if((grid->grid[ix-1][iy][iz] == TYPE_SOLVENT) ||
                  (grid->grid[ix-1][iy][iz] == TYPE_SOLVENT2))
               {
                  if(CheckAndChangeNeighbours(pdb, grid, 
                                              ix, iy, iz, 
                                              solvSize,
                                              0, 1, 1))
                  {
                     modified = TRUE;
                  }
               }
               break;
            }
         }
         /* Backwards...                                                */
         for(ix=grid->ixmax-2; ix>=0; ix--)
         {
            /* If this is a protein point                               */
            if(grid->grid[ix][iy][iz] == TYPE_PROTEIN)
            {
               /* If previous point was solvent                         */
               if((grid->grid[ix+1][iy][iz] == TYPE_SOLVENT) ||
                  (grid->grid[ix+1][iy][iz] == TYPE_SOLVENT2))
               {
                  if(CheckAndChangeNeighbours(pdb, grid, 
                                              ix, iy, iz, 
                                              solvSize,
                                              0, 1, 1))
                  {
                     modified = TRUE;
                  }
               }
               break;
            }
         }
      }
   }
   
   return(modified);
}

BOOL CheckAndChangeNeighbours(PDB *pdb, GRID *grid, 
                              int ix, int iy, int iz, 
                              REAL solvSize,
                              int doX, int doY, int doZ)
{
   int   i, j, k, n,
         numNeighbours;
   PDB   *neighbours[MAXNEIGHBOURS],
         *r;
   VEC3F centre,
         neighb,
         vec,
         pt;
   REAL  len, d, minDistSq = 0.0, distSq;
   BOOL  isVoid;

   doX = doY = doZ = 1;

   centre.x = grid->xmin + (ix * grid->gridStep);
   centre.y = grid->ymin + (iy * grid->gridStep);
   centre.z = grid->zmin + (iz * grid->gridStep);
   
   for(i=ix-doX; i<=ix+doX; i++)
   {
      for(j=iy-doY; j<=iy+doY; j++)
      {
         for(k=iz-doZ; k<=iz+doZ; k++)
         {
            /* If this is not our central voxel                         */
            if((i!=ix) || (j!=iy) || (k!=iz))
            {
               /* If it is also a protein voxel                         */
               if(grid->grid[i][j][k] == TYPE_PROTEIN)
               {
                  /* Build a list of neighbours to both voxels          */
                  numNeighbours = CombineTwoNeighbourLists(grid,
                                                           ix, iy, iz,
                                                           i, j, k,
                                                           neighbours);

                  /* Describe the vector between the two voxel centres  */
                  neighb.x = grid->xmin + (i * grid->gridStep);
                  neighb.y = grid->ymin + (j * grid->gridStep);
                  neighb.z = grid->zmin + (k * grid->gridStep);
                  
                  len = DIST(&centre, &neighb);
                  
                  /* Define the unit vector                             */
                  vec.x = (neighb.x - centre.x) / len;
                  vec.y = (neighb.y - centre.y) / len;
                  vec.z = (neighb.z - centre.z) / len;

                  /* Now walk from the centre of one voxel to the centre 
                     of the other and see if there are any points where 
                     a water can fit
                  */
                  for(d=0; d<len; d+=0.05)
                  {
                     pt.x = centre.x + (d*vec.x);
                     pt.y = centre.y + (d*vec.y);
                     pt.z = centre.z + (d*vec.z);

                     /* Check the distances to all neighbouring atoms. 
                        Assume that if all atoms are sufficiently distant
                        that there is actually a void point along this 
                        vector, once that is shown not to be true, set 
                        isVoid to FALSE and break out of the loop
                     */
                     isVoid = TRUE;
                     for(n=0; n<numNeighbours; n++)
                     {
                        r = neighbours[n];
                        minDistSq = (r->occ + solvSize) * 
                                    (r->occ + solvSize);
                        distSq = DISTSQ(r, &pt);
                        if(distSq < minDistSq)
                        {
                           isVoid = FALSE;
                           break;
                        }
                     }
                     
                     /* We have found a point along this line which is a 
                        void point. Reassign the two grids point from 
                        protein to solvent and return.
                     */
                     if(isVoid)
                     {
                        grid->grid[ix][iy][iz] = TYPE_SOLVENT;
                        grid->grid[i][j][k]    = TYPE_SOLVENT;
                        return(TRUE);
                     }
                  }
               }
            }
         }
      }
   }
   return(FALSE);
}


/************************************************************************/
int CombineTwoNeighbourLists(GRID *grid, int xi, int yi, int zi, 
                             int i, int j, int k,
                             PDB **neighbours)
{
   int  numNeighbours = 0,
        m, n;
   BOOL found;

   /* Copy the neighbours from the first point into our new list        */
   for(n=0; 
       (n < MAXNEIGHBOURS) && 
          (grid->neighbours[xi][yi][zi][n] != NULL);
       n++)
   {
      if(numNeighbours >= MAXNEIGHBOURS)
      {
         fprintf(stderr,"Error!: Maximum number of \
neighbours exceeded in combining neighbour lists.\n\
Increase MAXNEIGHBOURS!\n");
         exit(1);
      }
      neighbours[numNeighbours++] = grid->neighbours[xi][yi][zi][n];
   }
   
   /* Now step through the neighbours list for the second grid point    */
   for(n=0; 
       (n < MAXNEIGHBOURS) && 
          (grid->neighbours[i][j][k][n] != NULL);
       n++)
   {
      /* See if this neighbour is already in our
         combined list
      */
      found = FALSE;
      for(m=0; m<numNeighbours; m++)
      {
         if(grid->neighbours[i][j][k][n] ==
            neighbours[m])
         {
            found = TRUE;
            break;
         }
      }
      /* If not in the combined list, add it          */
      if(!found)
      {
         if(numNeighbours >= MAXNEIGHBOURS)
         {
            fprintf(stderr,"Error!: Maximum number of \
neighbours exceeded in combining neighbour lists.\n\
Increase MAXNEIGHBOURS!\n");
            exit(1);
         }
         neighbours[numNeighbours++] =
            grid->neighbours[i][j][k][n];
      }
   }  /* Loop through neighbours                      */

   return(numNeighbours);
}
