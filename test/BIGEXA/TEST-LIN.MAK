PROJ	=TEST-LIN
DEBUG	=1
CC	=qcl
CFLAGS_G	= /AL /W4 /Ze /DFLOAT
CFLAGS_D	= /Zi /Zr /Gi$(PROJ).mdt /Od
CFLAGS_R	= /O /Ot /DNDEBUG
CFLAGS	=$(CFLAGS_G) $(CFLAGS_D)
LFLAGS_G	= /CP:0xffff /NOI /SE:0x80 /ST:0x2710
LFLAGS_D	= /CO /INCR
LFLAGS_R	=
LFLAGS	=$(LFLAGS_G) $(LFLAGS_D)
RUNFLAGS	=
OBJS_EXT =
LIBS_EXT =

.asm.obj: ; $(AS) $(AFLAGS) -c $*.asm

all:	$(PROJ).EXE

spaceli3.obj:	spaceli3.c $(H)

spaceli4.obj:	spaceli4.c $(H)

spacelib.obj:	spacelib.c $(H)

test-lin.obj:	test-lin.c $(H)

linear.obj:	linear.c $(H)

linear2.obj:	linear2.c $(H)

$(PROJ).EXE:	spaceli3.obj spaceli4.obj spacelib.obj test-lin.obj linear.obj linear2.obj $(OBJS_EXT)
	echo >NUL @<<$(PROJ).crf
spaceli3.obj +
spaceli4.obj +
spacelib.obj +
test-lin.obj +
linear.obj +
linear2.obj +
$(OBJS_EXT)
$(PROJ).EXE

$(LIBS_EXT);
<<
	ilink -a -e "qlink $(LFLAGS) @$(PROJ).crf" $(PROJ)

run: $(PROJ).EXE
	$(PROJ) $(RUNFLAGS)