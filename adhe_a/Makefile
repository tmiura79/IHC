CC := gcc
CFLAGS := -g -Wall -Wextra
LIBS := -lm 

TARGET := adhe
SRC := adhe.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(TARGET)

.PHONY: all clean
