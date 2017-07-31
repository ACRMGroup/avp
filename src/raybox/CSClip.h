/*************************************************************************

   Program:    
   File:       CSClip.c
   
   Version:    V1.2
   Date:       01.02.96
   Function:   Header file for Cohen/Sutherland clipping
   
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

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
BOOL CSClip(int *x0in, int *y0in, int *x1in, int *y1in,
            int left, int up, int right, int down);

