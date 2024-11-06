vpath %.o obj

# compile macro
CC = gcc
CCFLAGS = -W -O3
LDFLAGS = -O3 -lm
INCLUDES = -I.

src:=$(sort $(wildcard *.c))
objects:=$(src:.c=.o)
exec:=polyop
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

# dependence
%.o: %.c|$(odir)
	${CC} ${CCFLAGS} ${INCLUDES} -c $< -o $(addprefix $(oprefix),$@)
$(exec):$(objects)|$(bdir)
	$(CC) $(addprefix $(oprefix),$(objects)) $(LDFLAGS) -o $(addprefix $(bprefix),$(exec))
$(odir):
	-$(md) $(odir)
$(bdir):
	-$(md) $(bdir)

.PHONY:clean
clean:
	-$(del) *.o $(exec)
	-$(rd) $(odir) $(bdir)
