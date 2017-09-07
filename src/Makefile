COPT = -O3 -Wall -ansi -pedantic -I$(HOME)/include  -Wno-unknown-pragmas
#COPT = -g -Wall -ansi -pedantic -I$(HOME)/include -Wno-unknown-pragmas
LOPT = -L$(HOME)/lib
LIBS = -lbiop -lgen -lm -lxml2
OFILES = avp.o RayBox.o

.c.o :
	$(CC) $(COPT) -c -o $@ $<

avp : $(OFILES)
	$(CC) $(LOPT) -o $@ $(OFILES) $(LIBS)

clean : 
	\rm $(OFILES)
