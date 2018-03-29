SRC_DIR := ./src
LIB_DIR := ./lib

CFLAGS += -std=c11 -fPIC -g -O3 -fstack-protector-all -D_FORTIFY_SOURCE=2 -DREENTRANT
LDFLAGS += -shared -lm

TARGET = $(LIB_DIR)/libspacelib.so
SOURCES = $(wildcard $(SRC_DIR)/*.c)
HEADERS = $(wildcard $(SRC_DIR)/*.H)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.c=$(SRC_DIR)/%.o)

.PHONY: all clean distclean

all: $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH) $(OBJECTS) -o $(TARGET)

clean:
	@- $(RM) $(TARGET)
	@- $(RM) $(OBJECTS)

distclean: clean