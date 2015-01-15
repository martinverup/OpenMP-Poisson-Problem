TARGET	= poisson
OBJS	= poisson.o

OPT		= -g -fast -xrestrict
PARA	= -xopenmp -xloopinfo -xvpara -xinstrument=datarace

CCC	= cc

CFLAGS	= $(OPT) $(PARA)


all: $(TARGET)

$(TARGET): $(OBJS) 
	$(CCC) $(CFLAGS) -o $@ $(OBJS)

clean:
	@/bin/rm -f *.o core

realclean: clean
	@/bin/rm -f $(TARGET)

poisson.o  : poisson.c
