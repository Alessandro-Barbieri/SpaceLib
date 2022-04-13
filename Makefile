SRC_DIR := ./src
LIB_DIR := ./lib

CC ?= clang

CFLAGS ?= -O3 -fstack-protector-all -D_FORTIFY_SOURCE=2
CFLAGS += -std=c11 -fPIC -g -Wall
LDFLAGS += -lm

TARGET = libspacelib.so
SOURCES = $(wildcard $(SRC_DIR)/*.c)
HEADERS = $(wildcard $(SRC_DIR)/*.H)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.c=$(SRC_DIR)/%.o)

.PHONY: all clean distclean

all: $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJECTS) -shared -Wl,-soname,${TARGET} -o $(TARGET)

clean:
	@- $(RM) $(TARGET)
	@- $(RM) $(OBJECTS)

distclean: clean
