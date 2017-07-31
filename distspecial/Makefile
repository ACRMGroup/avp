CC    = cc
# To allow decompression of gzipped PDB files on the fly, uncomment
# this line
# GUNZIP = -DGUNZIP_SUPPORT
COPT = -O2 -Wall -ansi -pedantic -I./include $(GUNZIP)
OFILES = avp.o RayBox.o
IFILES = raybox.h

LFILES = bioplib/ReadPDB.o \
bioplib/WritePDB.o \
bioplib/fsscanf.o \
bioplib/OpenStdFiles.o \
bioplib/FindNextResidue.o \
bioplib/IndexPDB.o \
bioplib/chindex.o \
bioplib/padterm.o 

avp : $(OFILES) $(LFILES)
	$(CC) -o $@ $(OFILES) $(LFILES) -lm

.c.o : $(IFILES)
	$(CC) $(COPT) -o $@ -c $<

clean : 
	\rm $(OFILES) $(LFILES)

