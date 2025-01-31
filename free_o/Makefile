CC := gcc
CFLAGS := -g -Wall -Wextra
LIBS := -lm 

TARGET := free
SRC := free.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(TARGET)

.PHONY: all clean