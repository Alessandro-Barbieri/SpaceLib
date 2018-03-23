PROJ	=SAT
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

spacelib.obj:	spacelib.c $(H)

spaceli4.obj:	spaceli4.c $(H)

spaceli3.obj:	spaceli3.c $(H)

sat.obj:	sat.c $(H)

$(PROJ).EXE:	spacelib.obj spaceli4.obj spaceli3.obj sat.obj $(OBJS_EXT)
	echo >NUL @<<$(PROJ).crf
spacelib.obj +
spaceli4.obj +
spaceli3.obj +
sat.obj +
$(OBJS_EXT)
$(PROJ).EXE

$(LIBS_EXT);
<<
	ilink -a -e "qlink $(LFLAGS) @$(PROJ).crf" $(PROJ)

run: $(PROJ).EXE
	$(PROJ) $(RUNFLAGS)

