PROG = rbin

SRC_DIR = src
SRC = $(SRC_DIR)/main.c \
      $(SRC_DIR)/formula.c \
      $(SRC_DIR)/sturges.c \
      $(SRC_DIR)/util.c

OBJ = $(SRC:.c=.o)

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra
LDLIBS  ?= -lm

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

all: $(PROG)

$(PROG): $(OBJ)
	$(CC) -o $@ $(OBJ) $(LDLIBS)

$(OBJ): $(SRC_DIR)/config.h $(SRC_DIR)/util.h

install: $(PROG)
	install -d $(DESTDIR)$(BINDIR)
	install -m 755 $(PROG) $(DESTDIR)$(BINDIR)/$(PROG)

clean:
	rm -f $(PROG) $(OBJ)

.PHONY: all install clean
