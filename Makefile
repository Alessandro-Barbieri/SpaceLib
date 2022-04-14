SRC_DIR := ./src
SHORT_TEST_DIR := ./test/SHORTEXA

CC ?= clang

CFLAGS ?= -O3 -fstack-protector-all -D_FORTIFY_SOURCE=2
CFLAGS += -std=c11 -g -Wall -I${SRC_DIR} -L./
LDFLAGS += -lm

TARGET = libspacelib.so
SOURCES = $(wildcard $(SRC_DIR)/*.c)
HEADERS = $(wildcard $(SRC_DIR)/*.H)
LIB_OBJECTS = $(SOURCES:$(SRC_DIR)/%.c=$(SRC_DIR)/%.o)

SHORT_TEST_SRC := $(wildcard $(SHORT_TEST_DIR)/*.c)
SHORT_TEST = $(SHORT_TEST_SRC:$(SHORT_TEST_DIR)/%.c=$(SHORT_TEST_DIR)/%)

.PHONY: all clean distclean

all: $(TARGET)

$(LIB_OBJECTS) : %.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -fPIC -c -o $@ $<

$(TARGET) : $(LIB_OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIB_OBJECTS) -fPIC -shared -Wl,-soname,${TARGET} -o $(TARGET)

$(SHORT_TEST) : %: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -fPIE -lspacelib -o $@ $<

test: $(TARGET) $(SHORT_TEST)

clean:
	@- $(RM) $(TARGET)
	@- $(RM) $(LIB_OBJECTS)
	@- $(RM) $(SHORT_TEST)

distclean: clean
