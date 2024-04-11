CC = gcc
CFLAGS = -Wall -Werror
INCLUDES = -Iincludes

SRCDIR = src
OBJDIR = obj
BINDIR = bin

SRCFILES = $(wildcard $(SRCDIR)/*.c)
OBJFILES = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCFILES))

LIBS = -lm

TARGET = $(BINDIR)/main

$(TARGET): $(OBJFILES)
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $^ -o $@ $(INCLUDES) $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean

clean:
	rm -rf $(OBJDIR) $(BINDIR)
