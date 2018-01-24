CC = gcc
CFLAGS = -O2 -Wall -ansi -pedantic
DEBUG = -ggdb -pg

LIBS = libgbfp.a
SHAREDLIB = libgbfp.so
SHAREDLIBV = libgbfp.so.1.0.1
SHAREDLIBM = libgbfp.so.1

AR = ar rcs

prefix = /usr/local
libdir = $(prefix)/lib
includedir = $(prefix)/include
mandir = $(prefix)/share/man
man3dir = $(mandir)/man3

OBJS = gbfp.o

all: $(SHAREDLIB) $(LIBS) seqext

$(OBJS): gbfp.c
	$(CC) $(CFLAGS) -D_GNU_SOURCE -fPIC -g -c gbfp.c

$(SHAREDLIB): $(OBJS)
	$(CC) -shared -Wl,-h,$(SHAREDLIBM) -o $(SHAREDLIBV) $(OBJS)

$(LIBS): $(OBJS)
	$(AR) $(LIBS) $(OBJS)

install: $(SHAREDLIB) $(LIBS)
	if [ ! -d $(libdir) ]; then mkdir -p $(libdir); fi
	if [ ! -d $(ncludedir) ]; then mkdir -p $(includedir); fi
	cp -f gbfp.h $(includedir)
	cp -f $(LIBS) $(libdir)
	cd $(libdir); chmod 755 $(LIBS)
	cp -f $(SHAREDLIBV) $(libdir)
	cd $(libdir); \
	if [ -f $(SHAREDLIBV) ]; then \
	 chmod 755 $(SHAREDLIBV); \
	 rm -f $(SHAREDLIB) $(SHAREDLIBM); \
	 ln -s $(SHAREDLIBV) $(SHAREDLIB); \
	 ln -s $(SHAREDLIBV) $(SHAREDLIBM); \
	 (ldconfig || true) > /dev/null 2>&1; \
	fi

uninstall:
	cd $(libdir); \
	rm -f $(SHAREDLIB) $(SHAREDLIBV) $(SHAREDLIBM) $(LIBS)
	cd $(includedir); \
	rm -f gbfp.h

seqext: seqext.c gbfp.c
	$(CC) $(CFLAGS) -D_GNU_SOURCE -o seqext seqext.c gbfp.c

debug: seqext.c gbfp.c
	$(CC) $(CFLAGS) $(DEBUG) -D_GNU_SOURCE -o seqext seqext.c gbfp.c

clean:
	rm -f $(OBJS) $(LIBS) $(SHAREDLIBV) $(SHAREDLIBM) seqext
