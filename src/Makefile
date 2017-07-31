COPT = -O3 -Wall -ansi -pedantic -I/home/bsm/martin/include
#COPT = -g -Wall -ansi -pedantic -I/home/bsm/martin/include
LOPT = -L/home/bsm/martin/lib
LIBS = -lbiop -lgen -lm
OFILES = avp.o RayBox.o

.c.o :
	$(CC) $(COPT) -c -o $@ $<

avp : $(OFILES)
	$(CC) $(LOPT) -o $@ $(OFILES) $(LIBS)

clean : 
	\rm $(OFILES)
