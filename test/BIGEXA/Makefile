CC ?= clang
CFLAGS ?= -O3 -fstack-protector-all -D_FORTIFY_SOURCE=2
CFLAGS += -std=c11 -fPIE -g -I../../src -L../.. -Wall
LDFLAGS += -lm -lspacelib

OBJECTS := ELB_D_DH ELB_D_PA ELB_I_DH ROB-MAT ROBSCARA SAT TEST TEST-LIN TEST_NEW TRIANGLE

.PHONY: all clean distclean

all : $(OBJECTS)

$(OBJECTS): %: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $<

clean:
	@- $(RM) $(OBJECTS)

distclean: clean
