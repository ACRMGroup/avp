/*************************************************************************

   Program:    
   File:       CSClip.c
   
   Version:    V1.2
   Date:       01.02.96
   Function:   Clips a line to fit within a box. Doesn't draw it.
   
   Copyright:  (c) SciTech Software 1991-6
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   CSClip() implements Cohen-Sutherland line clipping and returns 
   coordinates to fit within the clipping limits.

   This version simply returns the clipped coordinates (in the input
   pointers), so may be used for PostScript clipping, etc.

   Adapted from Clip() routine by Mike Nelson in AUI Jan 1991, p38--40

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.1  28.05.92 ANSIed and static'd
   V1.2  01.02.96 Cleaned up for incorporation in bioplib

*************************************************************************/
/* Includes
*/
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define  SPECIAL 255

/************************************************************************/
/* Globals
*/
static int     sCLP_LineX0,
               sCLP_LineX1,
               sCLP_LineY0,
               sCLP_LineY1;
static UBYTE   sCLP_CodeP1,
               sCLP_CodeP2,
               sCLP_NewCodeP1,
               sCLP_NewCodeP2;


/************************************************************************/
/* Prototypes
*/
static UBYTE Quadrant(int x, int y, int left, int up, int right, 
                      int down);
static void ClipEdge(int left, int up, int right, int down);
BOOL CSClip(int *x0in, int *y0in, int *x1in, int *y1in,
            int left, int up, int right, int down);


/************************************************************************/
/*>BOOL CSClip(int *x0in, int *y0in, int *x1in, int *y1in,
               int left,  int up,    int right, int down)
   -------------------------------------------------------
   Input:   int               left,up     Clip box limits
                              right,down
   I/O:     int               *x1,*y1     Line to draw. Returns clipped
                              *x2,*y2     values.
   Returns: BOOL              FALSE       Nothing to draw
                              TRUE        Draw something

   CSClip() implements Cohen-Sutherland line clipping and returns 
   coordinates to fit within the clipping limits.

   00.00.00 Original   By: ACRM
   28.05.92 ANSIed and static'd
   01.02.96 Cleaned up for incorporation in bioplib
*/
BOOL CSClip(int *x0in, int *y0in, int *x1in, int *y1in,
            int left,  int up,    int right, int down)
{
   int   x1temp,
         x0temp,
         y1temp,
         y0temp;
   int   x0,y0,x1,y1;

   /* We copy the input values from the pointers for ease of use        */
   x0 = *x0in;
   y0 = *y0in;
   x1 = *x1in;
   y1 = *y1in;
         
   /* Make copies of the input line coordinates                         */
   sCLP_LineX1 = x1temp = x1;
   sCLP_LineY1 = y1temp = y1;
   sCLP_LineX0 = x0temp = x0;
   sCLP_LineY0 = y0temp = y0;

   /* Find the quadrants for the ends of the lines                      */
   sCLP_CodeP1 = Quadrant(sCLP_LineX0,sCLP_LineY0,left,up,right,down);
   sCLP_CodeP2 = Quadrant(sCLP_LineX1,sCLP_LineY1,left,up,right,down);

   /* If one of these was set, then we need to do some clipping         */
   if(sCLP_CodeP1||sCLP_CodeP2)
   {
      if(sCLP_CodeP1 & sCLP_CodeP2)
      {  /* If both in same quadrant, we can ignore it                  */
         return(FALSE);
      }
      else
      {  /* They're not the same, we need to do some work               */
         if(!sCLP_CodeP2)      /* See if this end is in the window       */
         {
            x1 = sCLP_LineX1;
            y1 = sCLP_LineY1;
         }
         else                 /* It wasn't, so we must clip it          */
         {
            sCLP_NewCodeP1 = sCLP_CodeP1;
            sCLP_NewCodeP2 = sCLP_CodeP2;
            ClipEdge(left,up,right,down);
            if(sCLP_NewCodeP1 == SPECIAL)  /* Special case error         */
            {
               return(FALSE);
            }
            else              /* Keep these clipped points              */
            {
               x1 = sCLP_LineX0;
               y1 = sCLP_LineY0;
            }
         }
         
         /* Sort out the other point, swapping the coordinates round    */
         sCLP_LineX1     = x0temp;
         sCLP_LineX0     = x1temp;
         sCLP_LineY1     = y0temp;
         sCLP_LineY0     = y1temp;
         sCLP_NewCodeP1  = sCLP_CodeP2;
         sCLP_NewCodeP2  = sCLP_CodeP1;
         
         if(!sCLP_NewCodeP2)   /* See if this end is in the window       */
         {
            x0 = sCLP_LineX1;
            y0 = sCLP_LineY1;
         }
         else                 /* It wasn't, so we must clip it          */
         {
            ClipEdge(left,up,right,down);
            x0 = sCLP_LineX0;
            y0 = sCLP_LineY0;
         }
      }
   }
   
   /* O.K, we should now have our required coordinates,
      so we can return the fixed values
   */
   *x0in = x0;
   *y0in = y0;
   *x1in = x1;
   *y1in = y1;
   
   return(TRUE);
}
         
         
/************************************************************************/
/*>static UBYTE Quadrant(int x, int y, int left, int up, int right, 
                         int down)
   ----------------------------------------------------------------
   Return a code for the quadrant

   00.00.00 Original   By: ACRM
*/
static UBYTE Quadrant(int x, int y, int left, int up, int right, int down)
{
   UBYTE code = 0;

   if(x<left)  code |= 1;
   if(x>right) code |= 2;
   if(y>down)  code |= 4;
   if(y<up)    code |= 8;

   return(code);
}


/************************************************************************/
/*>static void ClipEdge(int left, int up, int right, int down)
   -----------------------------------------------------------
   Calculate the clip point

   00.00.00 Original   By: ACRM
*/
static void ClipEdge(int left, int up, int right, int down)
{
   int   LineCentreX,
         LineCentreY;
   UBYTE code,
         flag = 0;                                /* Used to break loop */
         
   LineCentreX = (sCLP_LineX0 + sCLP_LineX1) >> 1;  /* divide by 2      */
   LineCentreY = (sCLP_LineY0 + sCLP_LineY1) >> 1;  /* divide by 2      */
   
   while(!flag)
   {
      /* Check quadrant for middle of line                              */
      code = Quadrant(LineCentreX,LineCentreY,left,up,right,down);
      
      if(!(code & sCLP_NewCodeP2))   /* Are these in different quads?   */
      {
         sCLP_LineX0 = LineCentreX;  /* Move point 0 to the midpoint    */
         sCLP_LineY0 = LineCentreY;
      }
      else
      {
         sCLP_LineX1 = LineCentreX;  /* Move point 1 to the midpoint    */
         sCLP_LineY1 = LineCentreY;
      }
      
      /* Recalc the midpoint                                            */
      LineCentreX = (sCLP_LineX0 + sCLP_LineX1) >> 1; /* divide by 2    */
      LineCentreY = (sCLP_LineY0 + sCLP_LineY1) >> 1; /* divide by 2    */
   
      /* Check if we've finished                                        */
      if((LineCentreX == sCLP_LineX0)&&(LineCentreY == sCLP_LineY0)) 
         flag = 1;
      if((LineCentreX == sCLP_LineX1)&&(LineCentreY == sCLP_LineY1)) 
         flag = 1;
      
      /* Check for the special case                                     */
      if((code & sCLP_NewCodeP1)&&(code & sCLP_NewCodeP2))
      {
         flag = 1;
         sCLP_NewCodeP1 = SPECIAL;
      }
   }
}

