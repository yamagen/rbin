PROG = rbin

SRC_DIR = src
SRC = $(SRC_DIR)/main.c \
      $(SRC_DIR)/formula.c \
      $(SRC_DIR)/util.c

OBJ = $(SRC:.c=.o)

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra
LDLIBS  ?= -lm

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin


all: $(PROG)

test: $(PROG)
	@echo "Testing skew..."
	@./$(PROG) -s test/skew_test0.034.dat | grep -q 'skew=0.034'
	@echo "Testing kurtosis..."
	@./$(PROG) -s test/kurt_test2.17.dat | grep -q 'kurtosis=2.170'

$(PROG): $(OBJ)
	$(CC) -o $@ $(OBJ) $(LDLIBS)

$(OBJ): $(SRC_DIR)/config.h $(SRC_DIR)/util.h $(SRC_DIR)/formula.h

install: $(PROG)
	install -d $(DESTDIR)$(BINDIR)
	install -m 755 $(PROG) $(DESTDIR)$(BINDIR)/$(PROG)

clean:
	rm -f $(PROG) $(OBJ)

.PHONY: all install clean test
