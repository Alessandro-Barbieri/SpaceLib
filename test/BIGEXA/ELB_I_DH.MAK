PROJ	=ELB_I_DH
DEBUG	=1
CC	=qcl
CFLAGS_G	= /AM /W1 /Ze /DFLOAT
CFLAGS_D	= /Zi /Zr /Gi$(PROJ).mdt /Od
CFLAGS_R	= /O /Ot /DNDEBUG
CFLAGS	=$(CFLAGS_G) $(CFLAGS_D)
LFLAGS_G	= /CP:0xfff /NOI /SE:0x80 /ST:0x2710
LFLAGS_D	= /CO /INCR
LFLAGS_R	=
LFLAGS	=$(LFLAGS_G) $(LFLAGS_D)
RUNFLAGS	=
OBJS_EXT =
LIBS_EXT =

.asm.obj: ; $(AS) $(AFLAGS) -c $*.asm

all:	$(PROJ).EXE

elb_i_dh.obj:	elb_i_dh.c $(H)

spacelib.obj:	spacelib.c $(H)

spaceli4.obj:	spaceli4.c $(H)

spaceli3.obj:	spaceli3.c $(H)

linear.obj:	linear.c $(H)

linear2.obj:	linear2.c $(H)

$(PROJ).EXE:	elb_i_dh.obj spacelib.obj spaceli4.obj spaceli3.obj linear.obj linear2.obj $(OBJS_EXT)
	echo >NUL @<<$(PROJ).crf
elb_i_dh.obj +
spacelib.obj +
spaceli4.obj +
spaceli3.obj +
linear.obj +
linear2.obj +
$(OBJS_EXT)
$(PROJ).EXE

$(LIBS_EXT);
<<
	ilink -a -e "qlink $(LFLAGS) @$(PROJ).crf" $(PROJ)

run: $(PROJ).EXE
	$(PROJ) $(RUNFLAGS)