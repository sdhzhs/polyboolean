vpath %.o obj

# compile macro
CC = gcc
CCFLAGS = -W -O3 -fPIC
LDFLAGS = -O3 -lm
INCLUDES = -I.

src:=$(sort $(wildcard *.c))
mainsrc:=main.c
libsrc:=$(filter-out $(mainsrc),$(src))
libobjects:=$(libsrc:.c=.o)
mainobject:=$(mainsrc:.c=.o)
exec:=polyop
libname:=libpolyboolean.so
bprefix:=bin/
oprefix:=obj/
bdir:=bin
odir:=obj
md:=mkdir -p
rd:=rm -r -f
del:=rm -f
cp:=cp -f
cd:=cp -r -f

vpath $(exec) bin
vpath $(libname) bin

# dependence
%.o: %.c|$(odir)
	${CC} ${CCFLAGS} ${INCLUDES} -c $< -o $(addprefix $(oprefix),$@)

$(exec): $(mainobject) $(libname)|$(bdir)
	$(CC) $(addprefix $(oprefix),$(mainobject)) -L$(bdir) -lpolyboolean -Wl,-rpath,'$$ORIGIN' $(LDFLAGS) -o $(addprefix $(bprefix),$(exec))

$(libname): $(libobjects)|$(bdir)
	$(CC) -shared $(addprefix $(oprefix),$(libobjects)) $(LDFLAGS) -o $(addprefix $(bprefix),$(libname))

$(odir):
	-$(md) $(odir)

$(bdir):
	-$(md) $(bdir)

.PHONY:clean
clean:
	-$(del) *.o $(exec)
	-$(rd) $(odir) $(bdir)
