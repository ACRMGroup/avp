#define TWOWATERS
/*************************************************************************

   Program:    
   File:       
   
   Version:    
   Date:       
   Function:   
   
   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 2001
   Author:     Dr. Andrew C. R. Martin
   Address:    School of Animal and Microbial Sciences,
               The University of Reading,
               Whiteknights,
               P.O. Box 228,
               Reading RG6 6AJ.
               England.
   Phone:      +44 (0)118 987 5123 Extn. 7022
   Fax:        +44 (0)118 931 0180
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

For each triplet of atoms, where for each pair, the sum of the radii plus
the diameter of the solvent probe is > separation between the centres,
(i.e. the solvent can't just fit between the atoms) calculate the plane 
described by the 3 centres. Then find the point within the triangle
that is maximally spaced from the 3 vertices. Find the normal to that  
plane that passes through this point. Then find the points on that normal
on either side of the plane that are closest to the plane, but do not
lead to a clash with any other atom.

Alternatively, for a triplet of atoms, ignoring any other atoms, the
closest point to place a probe (on either side of the plane) can be
calculated as follows. The probe can be visualized as rolling around
the surface of any one atom. The sphere described by the centre of
the probe is thus a sphere of radius r_atom + r_probe centred about
the centre of the atom. To find where a probe sits amongst a set of
three atoms, one simply needs to find the intersections of the three
enlarged spheres. As before, draw a vector through this point normal 
to the plane described by the 3 atoms and work out the closest point
to the plane that doesn't clash with any other atom

Code URLS:
ftp://lassp-ftp.ccmr.cornell.edu/pub/
http://ggt.sourceforge.net/html/MatrixOps_8h-source.html
http://www.netlib.org/sfmm/solve.f
http://www.netlib.org/sfmm/decomp.f
http://www.nauticom.net/www/jdtaft/JavaMatrix.htm
http://calcium.emc.uq.edu.au/apr/articles/pre2cut/node10.html
http://mathforum.org/library/drmath/sets/college_3d.html
http://www.acm.org/tog/GraphicsGems/


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/MathUtil.h"
#include "bioplib/matrix.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define WAT_RADIUS 1.4
#define CLASHSQ 7.84

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/long SolveNonlinearSystem(REAL (**func)(REAL*),
                          REAL *x,
                          register long	n,
                          REAL numsig,
                          REAL *maxit);
PDB *Solvate(PDB *pdb);
void SetWater(PDB *w, VEC3F watPos);
int PlaceWater(PDB *pdb, PDB *p, PDB *q, PDB *r, VEC3F CofG, 
                VEC3F *watPos1, VEC3F *watPos2);
REAL eq1(REAL *x);
REAL eq2(REAL *x);
REAL eq3(REAL *x);
VEC3F FindVectorTangentToPlane(REAL x1, REAL y1, REAL z1,
                               REAL x2, REAL y2, REAL z2,
                               REAL x3, REAL y3, REAL z3);
void PrintAsPDB(REAL x,REAL y,REAL z, REAL r,char *resnam);
BOOL Find3SphereIntersection(REAL x1, REAL y1, REAL z1, REAL r1,
                             REAL x2, REAL y2, REAL z2, REAL r2,
                             REAL x3, REAL y3, REAL z3, REAL r3,
                             VEC3F *isect1, VEC3F *isect2);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);
void SetRadii(PDB *pdb);


/************************************************************************/

/************************************************************************/
int main(int argc, char **argv)
{
   FILE *in  = stdin, 
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   PDB  *pdb  = NULL, 
        *pdbw = NULL;
   int  natoms;

   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=ReadPDB(in, &natoms))!=NULL)
         {
            SetRadii(pdb);
            if((pdbw = Solvate(pdb))==NULL)
            {
               fprintf(stderr,"No memory for waters\n");
               return(1);
            }
            
            WritePDB(out, pdb);
            WritePDB(out, pdbw);
         }
         else
         {
            fprintf(stderr,"Unable to read input PDB file\n");
            return(1);            
         }
      }
      else
      {
         fprintf(stderr,"Unable to open I/O files\n");
         return(1);            
      }
   }
   else
   {
      Usage();
   }
   return(0);
}

/************************************************************************/
PDB *Solvate(PDB *pdb)
{
   PDB   *p, *q, *r,
         *water = NULL, 
         *w;
   VEC3F CofG,
         watPos1,
         watPos2;
   int   nwater;
   
   GetCofGPDB(pdb, &CofG);
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      for(q=p->next; q!=NULL; NEXT(q))
      {
         for(r=q->next; r!=NULL; NEXT(r))
         {
            nwater = PlaceWater(pdb, p, q, r, CofG, &watPos1, &watPos2);
            if(nwater & 1)
            {
               if(water==NULL)
               {
                  INIT(water, PDB);
                  w=water;
               }
               else
               {
                  ALLOCNEXT(w, PDB);
               }
               if(w==NULL)
               {
                  FREELIST(water, PDB);
                  return(NULL);
               }
               SetWater(w, watPos1);

#ifdef TWOWATER
               if(nwater & 2)
               {
                  ALLOCNEXT(w, PDB);
                  if(w==NULL)
                  {
                     FREELIST(water, PDB);
                     return(NULL);
                  }
                  SetWater(w, watPos2);
               }
#endif
            }
         }
      }
   }
   return(water);
}

/************************************************************************/
void SetWater(PDB *w, VEC3F watPos)
{
   static int atnum = 9000;
   
   strcpy(w->record_type, "ATOM  ");
   w->atnum = atnum % 5000;
   w->resnum = atnum % 5000;
   atnum++;
   
   strcpy(w->atnam, "O   ");
   strcpy(w->atnam_raw, "O   ");
   strcpy(w->resnam, "WAT ");
   strcpy(w->insert, " ");
   strcpy(w->chain, " ");
   w->x = watPos.x;
   w->y = watPos.y;
   w->z = watPos.z;
   w->occ = 1.0;
   w->bval = 20.0;
}

/************************************************************************/
int PlaceWater(PDB *pdb, PDB *p, PDB *q, PDB *r, VEC3F CofG, 
                VEC3F *watPos1, VEC3F *watPos2)
{
   REAL  r1, r2, r3, d1, d2, d3;
   VEC3F isect1, isect2;
   PDB   *s;
   int   nwater = 0;
   BOOL  hit1 = FALSE,
         hit2 = FALSE;
   

   d1 = DISTSQ(p, q);
   d2 = DISTSQ(p, r);
   d3 = DISTSQ(r, q);
   
   r1 = p->occ + WAT_RADIUS;
   r2 = q->occ + WAT_RADIUS;
   r3 = r->occ + WAT_RADIUS;

   /* Only look at atom triplets which are close enough to support the
      water resting on them
   */
   if((d1 <= (r1+r2)*(r1+r2)) &&
      (d2 <= (r1+r3)*(r1+r3)) &&
      (d3 <= (r2+r3)*(r2+r3)))
   {
      if(Find3SphereIntersection(p->x, p->y, p->z, r1,
                                 q->x, q->y, q->z, r2,
                                 r->x, r->y, r->z, r3,
                                 &isect1, &isect2))
      {
         /* Ensure that isect1 is the water farther from the CofG */
         d1 = DISTSQ(&isect1, &CofG);
         d2 = DISTSQ(&isect2, &CofG);
         if(d1 > d2)
         {
            *watPos1 = isect1;
            *watPos2 = isect2;
         }
         else
         {
            *watPos1 = isect2;
            *watPos2 = isect1;
         }
      
         for(s=pdb; s!=NULL; NEXT(s))
         {
            if((s!=p) && (s!=q) && (s!=r))
            {
               if(DISTSQ(s, watPos1) < 
                  ((s->occ+WAT_RADIUS)*(s->occ+WAT_RADIUS)))
               {
                  hit1=TRUE;
                  break;
               }
            }
         }
         
         for(s=pdb; s!=NULL; NEXT(s))
         {
            if((s!=p) && (s!=q) && (s!=r))
            {
               if(DISTSQ(s, watPos2) < CLASHSQ)
               {
                  hit2=TRUE;
                  break;
               }
            }
         }

         nwater=0;
         if(!hit1) nwater |= 1;
         if(!hit2) nwater |= 2;
      }
   }

   return(nwater);
}



/************************************************************************/
static REAL eq1x, eq1y, eq1z, eq1r,
            eq2x, eq2y, eq2z, eq2r,
            eq3x, eq3y, eq3z, eq3r;

/************************************************************************/
REAL eq1(REAL *x)
{
   return(((x[0] - eq1x)*(x[0] - eq1x)) +
          ((x[1] - eq1y)*(x[1] - eq1y)) +
          ((x[2] - eq1z)*(x[2] - eq1z)) -
          (eq1r * eq1r));
}
/************************************************************************/
REAL eq2(REAL *x)
{
   return(((x[0] - eq2x)*(x[0] - eq2x)) +
          ((x[1] - eq2y)*(x[1] - eq2y)) +
          ((x[2] - eq2z)*(x[2] - eq2z)) -
          (eq2r * eq2r));
}
/************************************************************************/
REAL eq3(REAL *x)
{
   return(((x[0] - eq3x)*(x[0] - eq3x)) +
          ((x[1] - eq3y)*(x[1] - eq3y)) +
          ((x[2] - eq3z)*(x[2] - eq3z)) -
          (eq3r * eq3r));
}


/************************************************************************/
/* If you know three points on a plane and you want an equation for the
   plane, you can take the difference between point 1 and point 2 and 
   call that the first vector. Then take the difference between point 1 
   and point 3 and call that the second vector. Now if you take the 
   cross product of the first vector and the second vector, you will have 
   a vector that points normal to the plane. (This is because it is a 
   fundamental property of the vector cross product that the result is 
   always perpendicular to the two vectors you are multiplying.)  Once 
   you have the normal and any one point, you can proceed as above.    
*/
VEC3F FindVectorTangentToPlane(REAL x1, REAL y1, REAL z1,
                               REAL x2, REAL y2, REAL z2,
                               REAL x3, REAL y3, REAL z3)
{
   VEC3F v1, v2, nv;
   REAL  vlen;

   
   v1.x = x2 - x1;
   v1.y = y2 - y1;
   v1.z = z2 - z1;

   v2.x = x3 - x1;
   v2.y = y3 - y1;
   v2.z = z3 - z1;
   
   CrossProd3(&nv, v1, v2);

   /* normalize this vector */
   vlen = sqrt((nv.x*nv.x) + (nv.y*nv.y) + (nv.z*nv.z));
   nv.x /= vlen;
   nv.y /= vlen;
   nv.z /= vlen;
   
   return(nv);
}


/************************************************************************/
void PrintAsPDB(REAL x,REAL y,REAL z, REAL r,char *resnam)
{
   static int atnum = 0;
   
   printf("ATOM  %5d  CA  %3s %5d    %8.3f%8.3f%8.3f  1.00%6.2f              \n",
          atnum,resnam,atnum,x,y,z,r);
   atnum++;
}


/************************************************************************/
typedef REAL (*FUNCPOINTER)(REAL *);

BOOL Find3SphereIntersection(REAL x1, REAL y1, REAL z1, REAL r1,
                             REAL x2, REAL y2, REAL z2, REAL r2,
                             REAL x3, REAL y3, REAL z3, REAL r3,
                             VEC3F *isect1, VEC3F *isect2)
{
   FUNCPOINTER funcs[3];
   REAL x[3];
   REAL niter = 1000.0;
   VEC3F nv;
   BOOL retval = TRUE;

   eq1x = x1;
   eq1y = y1;
   eq1z = z1;
   eq1r = r1;
   eq2x = x2;
   eq2y = y2;
   eq2z = z2;
   eq2r = r2;
   eq3x = x3;
   eq3y = y3;
   eq3z = z3;
   eq3r = r3;

   funcs[0] = eq1;
   funcs[1] = eq2;
   funcs[2] = eq3;

   /* Make a guess at the solution - take the average coordinates, then move along
      the normal vector by the smallest radius
   */
   nv = FindVectorTangentToPlane(x1,y1,z1,x2,y2,z2,x3,y3,z3);
   
   x[0] = ((x1+x2+x3) / 3.0);
   x[1] = ((y1+y2+y3) / 3.0);
   x[2] = ((z1+z2+z3) / 3.0);
#ifdef DEBUG
   PrintAsPDB(x[0],x[1],x[2],0.5,"CYS");
#endif
   x[0] -= r1 * nv.x;
   x[1] -= r1 * nv.y;
   x[2] -= r1 * nv.z;
   


#ifdef DEBUG
   fprintf(stderr,"Guess 1: %.3f %.3f %.3f\n", x[0],x[1],x[2]);
   PrintAsPDB(x1,y1,z1,r1,"GLY");
   PrintAsPDB(x2,y2,z2,r2,"GLY");
   PrintAsPDB(x3,y3,z3,r3,"GLY");
   PrintAsPDB(x[0],x[1],x[2],1.0,"ALA");
#endif
   
   if(!SolveNonlinearSystem(funcs, x, 3, 3.0, &niter))
   {
      retval = FALSE;
   }
   isect1->x = x[0];
   isect1->y = x[1];
   isect1->z = x[2];

   x[0] = ((x1+x2+x3) / 3.0) + (r1 * nv.x);
   x[1] = ((y1+y2+y3) / 3.0) + (r1 * nv.y);
   x[2] = ((z1+z2+z3) / 3.0) + (r1 * nv.z);

#ifdef DEBUG
   fprintf(stderr,"Guess 2: %.3f %.3f %.3f\n", x[0],x[1],x[2]);
   PrintAsPDB(x[0],x[1],x[2],1.0,"ALA");
#endif
   
   if(!SolveNonlinearSystem(funcs, x, 3, 3.0, &niter))
   {
      retval = FALSE;
   }
   isect2->x = x[0];
   isect2->y = x[1];
   isect2->z = x[2];
   
   return(retval);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
   Input:   int    argc        Argument count
            char   **argv      Argument array
   Output:  char   *infile     Input filename (or blank string)
            char   *outfile    Output filename (or blank string)
   Returns: BOOL               Success

   Parse the command line

   29.06.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
{
   argc--;
   argv++;
   
   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
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

void Usage(void)
{
   
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


