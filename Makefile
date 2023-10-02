#
# Written by Jonas Juselius <jonas@iki.fi> Tue Dec  4 11:41:04 EET 2001
#
#DEBUG=1
#PROF=1
#RANGE=1

topdir:=..
builddir:=..
 
include $(builddir)/make.config

INST_PROGS:=xx2caogrid
INST_LIBS:=
INST_INCLUDES:=
INST_DATA:=
# module names are listed in smallcaps and wihtout any .mod
INST_MODULES:=

FCFLAGS+=-I../libx2c

include $(builddir)/make.rules

all: $(INST_PROGS) $(INST_LIBS)

libs:=-lx2cint -llibr

$(bindir)/xx2caogrid: $(objs) $(libs)
	$(FC) $(FLDFLAGS) -o $@ $(objs) $(libs) $(BLAS_LIBS) $(LAPACK_LIBS)



