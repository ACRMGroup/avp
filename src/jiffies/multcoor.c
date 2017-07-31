#include <stdio.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"

int main(int argc, char **argv)
{
   PDB *pdb,
      *p;
   int natoms;
   
   pdb = ReadPDB(stdin, &natoms);
   for(p=pdb; p!=NULL; NEXT(p))
   {
      p->x *= 10.0;
      p->y *= 10.0;
      p->z *= 10.0;
   }
   
   WritePDB(stdout, pdb);

   return(0);
}
