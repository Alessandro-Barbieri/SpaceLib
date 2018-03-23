PROJ	=ELB_D_PA
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
RUNFLAGS	=elbow.dat tmp.mot t2.t
OBJS_EXT = 	
LIBS_EXT = 	

.asm.obj: ; $(AS) $(AFLAGS) -c $*.asm

all:	$(PROJ).EXE

elb_d_pa.obj:	elb_d_pa.c $(H)

spacelib.obj:	spacelib.c $(H)

spaceli4.obj:	spaceli4.c $(H)

spaceli3.obj:	spaceli3.c $(H)

$(PROJ).EXE:	elb_d_pa.obj spacelib.obj spaceli4.obj spaceli3.obj $(OBJS_EXT)
	echo >NUL @<<$(PROJ).crf
elb_d_pa.obj +
spacelib.obj +
spaceli4.obj +
spaceli3.obj +
$(OBJS_EXT)
$(PROJ).EXE

$(LIBS_EXT);
<<
	ilink -a -e "qlink $(LFLAGS) @$(PROJ).crf" $(PROJ)

run: $(PROJ).EXE
	$(PROJ) $(RUNFLAGS)

