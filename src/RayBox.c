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
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define RIGHT	0
#define LEFT	1
#define MIDDLE	2

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
BOOL DoesLineCrossBox(VEC3F start, VEC3F stop, REAL length, 
                      VEC3F corner1, VEC3F corner2);
BOOL HitBoundingBox(REAL minB[3], REAL maxB[3],
                    REAL origin[3], REAL dir[3],
                    REAL coord[3]);


/************************************************************************/
/*>BOOL DoesLineCrossBox(VEC3F start, VEC3F stop, REAL length, 
                         VEC3F corner1, VEC3F corner2)
   -----------------------------------------------------------
   Tests whether a line specified by its two endpoints (start and stop)
   or its start (start), direction (stop) and length )length) intersects 
   a box. If given the two endpoints, 'length' should be set to zero.

   16.07.02 Original   By: ACRM
*/
BOOL DoesLineCrossBox(VEC3F start, VEC3F stop, REAL length, 
                      VEC3F corner1, VEC3F corner2)
{
   REAL bl[3],
        tr[3],
        pt[3],
        dir[3],
        is[3],
        lengthsq;
   VEC3F isect;
   
   
   if(length==(REAL)0.0)
   {
      lengthsq = DISTSQ(&start, &stop);
      length = (REAL)sqrt((double)lengthsq);
      
      stop.x = (stop.x - start.x)/length;
      stop.y = (stop.y - start.y)/length;
      stop.z = (stop.z - start.z)/length;
   }
   else
   {
      lengthsq = length * length;
   }
   

   bl[0]  = corner1.x;   bl[1]  = corner1.y;   bl[2]  = corner1.z;
   tr[0]  = corner2.x;   tr[1]  = corner2.y;   tr[2]  = corner2.z;
   pt[0]  = start.x;     pt[1]  = start.y;     pt[2]  = start.z;
   dir[0] = stop.x;      dir[1] = stop.y;      dir[2] = stop.z;
   if(!HitBoundingBox(bl, tr, pt, dir, is))
   {
      return(FALSE);
   }

   isect.x = is[0];   isect.y = is[1];   isect.z = is[2];
   if(DISTSQ(&start, &isect) > lengthsq)
   {
      return(FALSE);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL HitBoundingBox(REAL minB[3], REAL maxB[3],
                       REAL origin[3], REAL dir[3],
                       REAL coord[3])
   ----------------------------------------------------------

   Given the two corners of a box, the origin of a line and its direction,
   determine whether the line intersects the box and, if so, give the
   coordinates of the first intersection point.

   Tidied up from 'Fast Ray-Box Intersection' by Andrew Woo
   from "Graphics Gems", Academic Press, 1990
   http://www.acm.org/pubs/tog/GraphicsGems/gems/

   01.01.90  Original By: AW
   16.07.02  Tidied up By: ACRM
*/
BOOL HitBoundingBox(REAL minB[3], REAL maxB[3],
                    REAL origin[3], REAL dir[3],
                    REAL coord[3])
{
   register int i;
   int          whichPlane;
   char         inside = TRUE;
   char         quadrant[3];
   REAL         maxT[3];
   REAL         candidatePlane[3];
   
   /* Find candidate planes; this loop can be avoided if rays cast all 
      from the eye(assume perpsective view) 
   */
   for (i=0; i<3; i++)
   {
      if(origin[i] < minB[i]) 
      {
         quadrant[i]       = LEFT;
         candidatePlane[i] = minB[i];
         inside            = FALSE;
      }
      else if (origin[i] > maxB[i]) 
      {
         quadrant[i]       = RIGHT;
         candidatePlane[i] = maxB[i];
         inside            = FALSE;
      }
      else	
      {
         quadrant[i] = MIDDLE;
      }
   }
   
   
   /* Ray origin inside bounding box                                    */
   if(inside)
   {
      coord = origin;
      return(TRUE);
   }
   
   
   /* Calculate T distances to candidate planes                         */
   for(i=0; i<3; i++)
   {
      if((quadrant[i] != MIDDLE) && (dir[i] != 0.0))
      {
         maxT[i] = (candidatePlane[i]-origin[i]) / dir[i];
      }
      else
      {
         maxT[i] = -1.0;
      }
   }
   
   
   /* Get largest of the maxT's for final choice of intersection        */
   whichPlane = 0;
   for(i=1; i<3; i++)
   {
      if (maxT[whichPlane] < maxT[i])
      {
         whichPlane = i;
      }
   }
   
   
   /* Check final candidate actually inside box                         */
   if(maxT[whichPlane] < 0.0)
   {
      return(FALSE);
   }
   
   for(i=0; i<3; i++)
   {
      if(whichPlane != i) 
      {
         coord[i] = origin[i] + maxT[whichPlane] * dir[i];
         if((coord[i] < minB[i]) || (coord[i] > maxB[i]))
         {
            return(FALSE);
         }
      } 
      else 
      {
         coord[i] = candidatePlane[i];
      }
   }
   
   return(TRUE);				/* ray hits box         */
}	


#ifdef DEMO
#include <stdio.h>
int main(int argc, char **argv)
{
   VEC3F start, stop, corner1, corner2;
   REAL  length = 0.0;
   
   corner1.x = 0.0;  corner1.y = 0.0;  corner1.z = 0.0;
   corner2.x = 10.0; corner2.y = 10.0; corner2.z = 10.0;

   start.x = -1.0;   start.y =  5.0;   start.z =  5.0;   
   stop.x  = -0.1;  stop.y  =  5.0;   stop.z  =  5.0;   

   if(DoesLineCrossBox(start, stop, length, corner1, corner2))
   {
      printf("Intersection!\n");
   }
   else
   {
      printf("No intersection\n");
   }


   start.x = -1.0;   start.y =  5.0;   start.z =  5.0;   
   stop.x  = 1.0;    stop.y  =  0.0;   stop.z  =  0.0;   
   length = 0.1;

   if(DoesLineCrossBox(start, stop, length, corner1, corner2))
   {
      printf("Intersection!\n");
   }
   else
   {
      printf("No intersection\n");
   }

   return(0);
}

#endif
