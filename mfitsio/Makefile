# File:        Makefile, Makefile.header
# Author:      Damian Ryan Eads <eads@lanl.gov>
# Description: The makefile for compiling mfitsio within a Linux environment.
# Created:     December 9, 2002
#
# This software and ancillary information (herein called ``Software'')
# called MFITSIO is made available under the terms described here. The
# SOFTWARE has been approved for release with associated LA-CC number:
# LA-CC-02-085.
#
# This SOFTWARE has been authored by an employee or employees of the
# University of California, operator of the Los Alamos National Laboratory
# under Contract No. W-7405-ENG-36 with the U.S. Department of Energy. 
# The U.S. Government has rights to use, reproduce, and distribute this
# SOFTWARE.  The public may copy, distribute, prepare derivative works and
# publicly display this SOFTWARE without charge, provided that this Notice
# and any statement of authorship are reproduced on all copies.  Neither
# the Government nor the University makes any warranty, express or
# implied, or assumes any liability or responsibility for the use of this
# SOFTWARE.  If SOFTWARE is modified to produce derivative works, such
# modified SOFTWARE should be clearly marked, so as not to confuse it with
# the version available from LANL.
#
# Consult the COPYING file for more details on licensing.


CC = mex
DOCGEN = doxygen
CFLAGS = mfitsio.c -lcfitsio -I/usr/local/include
HEADERS = mfitsio.h
C_FILES = mfitsio.c
VERSION = 1.2.4

all:
	./bin/adetect binaries

tar:
	./bin/adetect tar_binaries

tar_src:
	./bin/adetect tar_src

clean:
	./bin/adetect clean

doc:
	./bin/adetect doc
linux: linux_binaries

linux_binaries: fits_read_image.mexglx \
                fits_write_image.mexglx \
                fits_read_header.mexglx \
                fits_write_header.mexglx \
                fits_delete_keyword.mexglx \
		fits_read_image_subset.mexglx \
		fits_write_image_subset.mexglx

fits_read_image.mexglx: fits_read_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image.c $(CFLAGS)

fits_read_image_subset.mexglx: fits_read_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image_subset.c $(CFLAGS)

fits_write_image_subset.mexglx: fits_write_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image_subset.c $(CFLAGS)

fits_write_image.mexglx: fits_write_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image.c $(CFLAGS)

fits_read_header.mexglx: fits_read_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_header.c $(CFLAGS)

fits_write_header.mexglx: fits_write_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_header.c $(CFLAGS)

fits_delete_keyword.mexglx: fits_delete_keyword.c $(HEADERS) $(C_FILES)
	$(CC) fits_delete_keyword.c $(CFLAGS)

linux_tar_binaries:
	mkdir mfitsio-$(VERSION)-linux-bin
	mkdir mfitsio-$(VERSION)-linux-bin/mfitsio
	cd mfitsio-$(VERSION)-linux-bin/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../*.m ../../*.mexglx \
	   ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL \
	   ../../Makefile .
	tar cvf - mfitsio-$(VERSION)-linux-bin/* \
	        | gzip -c - > mfitsio-$(VERSION)-linux-bin.tar.gz
	rm -rf mfitsio-$(VERSION)-linux-bin

linux_tar_src:
	mkdir mfitsio-$(VERSION)-src
	mkdir mfitsio-$(VERSION)-src/mfitsio
	cd mfitsio-$(VERSION)-src/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../Makefile ../../Makefile.root ../../Makefile.header .; \
	cp ../../*.[chm] ../../ChangeLog ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL .
	tar cvf - mfitsio-$(VERSION)-src/* | gzip -c - > mfitsio-$(VERSION)-src.tar.gz
	rm -rf mfitsio-$(VERSION)-src

linux_doc:
	doxygen Doxyfile

linux_clean:
	rm *.mexglx
solaris: solaris_binaries

solaris_binaries: fits_read_image.mexsol \
                fits_write_image.mexsol \
                fits_read_header.mexsol \
                fits_write_header.mexsol \
                fits_delete_keyword.mexsol \
		fits_read_image_subset.mexsol \
		fits_write_image_subset.mexsol

fits_read_image.mexsol: fits_read_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image.c $(CFLAGS)

fits_read_image_subset.mexsol: fits_read_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image_subset.c $(CFLAGS)

fits_write_image_subset.mexsol: fits_write_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image_subset.c $(CFLAGS)

fits_write_image.mexsol: fits_write_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image.c $(CFLAGS)

fits_read_header.mexsol: fits_read_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_header.c $(CFLAGS)

fits_write_header.mexsol: fits_write_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_header.c $(CFLAGS)

fits_delete_keyword.mexsol: fits_delete_keyword.c $(HEADERS) $(C_FILES)
	$(CC) fits_delete_keyword.c $(CFLAGS)

solaris_tar_binaries:
	mkdir mfitsio-$(VERSION)-solaris-bin
	mkdir mfitsio-$(VERSION)-solaris-bin/mfitsio
	cd mfitsio-$(VERSION)-solaris-bin/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../*.m ../../*.mexsol \
	   ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL \
	   ../../Makefile .
	tar cvf - mfitsio-$(VERSION)-solaris-bin/* \
	        | gzip -c - > mfitsio-$(VERSION)-solaris-bin.tar.gz
	rm -rf mfitsio-$(VERSION)-solaris-bin

solaris_tar_src:
	mkdir mfitsio-$(VERSION)-src
	mkdir mfitsio-$(VERSION)-src/mfitsio
	cd mfitsio-$(VERSION)-src/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../Makefile ../../Makefile.root ../../Makefile.header .; \
	cp ../../*.[chm] ../../ChangeLog ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL .
	tar cvf - mfitsio-$(VERSION)-src/* | gzip -c - > mfitsio-$(VERSION)-src.tar.gz
	rm -rf mfitsio-$(VERSION)-src

solaris_doc:
	doxygen Doxyfile

solaris_clean:
	rm *.mexsol
macosx: macosx_binaries

macosx_binaries: fits_read_image.mexmac \
                fits_write_image.mexmac \
                fits_read_header.mexmac \
                fits_write_header.mexmac \
                fits_delete_keyword.mexmac \
		fits_read_image_subset.mexmac \
		fits_write_image_subset.mexmac

fits_read_image.mexmac: fits_read_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image.c $(CFLAGS)

fits_read_image_subset.mexmac: fits_read_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image_subset.c $(CFLAGS)

fits_write_image_subset.mexmac: fits_write_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image_subset.c $(CFLAGS)

fits_write_image.mexmac: fits_write_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image.c $(CFLAGS)

fits_read_header.mexmac: fits_read_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_header.c $(CFLAGS)

fits_write_header.mexmac: fits_write_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_header.c $(CFLAGS)

fits_delete_keyword.mexmac: fits_delete_keyword.c $(HEADERS) $(C_FILES)
	$(CC) fits_delete_keyword.c $(CFLAGS)

macosx_tar_binaries:
	mkdir mfitsio-$(VERSION)-macosx-bin
	mkdir mfitsio-$(VERSION)-macosx-bin/mfitsio
	cd mfitsio-$(VERSION)-macosx-bin/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../*.m ../../*.mexmac \
	   ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL \
	   ../../Makefile .
	tar cvf - mfitsio-$(VERSION)-macosx-bin/* \
	        | gzip -c - > mfitsio-$(VERSION)-macosx-bin.tar.gz
	rm -rf mfitsio-$(VERSION)-macosx-bin

macosx_tar_src:
	mkdir mfitsio-$(VERSION)-src
	mkdir mfitsio-$(VERSION)-src/mfitsio
	cd mfitsio-$(VERSION)-src/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../Makefile ../../Makefile.root ../../Makefile.header .; \
	cp ../../*.[chm] ../../ChangeLog ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL .
	tar cvf - mfitsio-$(VERSION)-src/* | gzip -c - > mfitsio-$(VERSION)-src.tar.gz
	rm -rf mfitsio-$(VERSION)-src

macosx_doc:
	doxygen Doxyfile

macosx_clean:
	rm *.mexmac
hp700: hp700_binaries

hp700_binaries: fits_read_image.mexhp7 \
                fits_write_image.mexhp7 \
                fits_read_header.mexhp7 \
                fits_write_header.mexhp7 \
                fits_delete_keyword.mexhp7 \
		fits_read_image_subset.mexhp7 \
		fits_write_image_subset.mexhp7

fits_read_image.mexhp7: fits_read_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image.c $(CFLAGS)

fits_read_image_subset.mexhp7: fits_read_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image_subset.c $(CFLAGS)

fits_write_image_subset.mexhp7: fits_write_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image_subset.c $(CFLAGS)

fits_write_image.mexhp7: fits_write_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image.c $(CFLAGS)

fits_read_header.mexhp7: fits_read_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_header.c $(CFLAGS)

fits_write_header.mexhp7: fits_write_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_header.c $(CFLAGS)

fits_delete_keyword.mexhp7: fits_delete_keyword.c $(HEADERS) $(C_FILES)
	$(CC) fits_delete_keyword.c $(CFLAGS)

hp700_tar_binaries:
	mkdir mfitsio-$(VERSION)-hp700-bin
	mkdir mfitsio-$(VERSION)-hp700-bin/mfitsio
	cd mfitsio-$(VERSION)-hp700-bin/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../*.m ../../*.mexhp7 \
	   ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL \
	   ../../Makefile .
	tar cvf - mfitsio-$(VERSION)-hp700-bin/* \
	        | gzip -c - > mfitsio-$(VERSION)-hp700-bin.tar.gz
	rm -rf mfitsio-$(VERSION)-hp700-bin

hp700_tar_src:
	mkdir mfitsio-$(VERSION)-src
	mkdir mfitsio-$(VERSION)-src/mfitsio
	cd mfitsio-$(VERSION)-src/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../Makefile ../../Makefile.root ../../Makefile.header .; \
	cp ../../*.[chm] ../../ChangeLog ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL .
	tar cvf - mfitsio-$(VERSION)-src/* | gzip -c - > mfitsio-$(VERSION)-src.tar.gz
	rm -rf mfitsio-$(VERSION)-src

hp700_doc:
	doxygen Doxyfile

hp700_clean:
	rm *.mexhp7
hpunix: hpunix_binaries

hpunix_binaries: fits_read_image.mexhpux \
                fits_write_image.mexhpux \
                fits_read_header.mexhpux \
                fits_write_header.mexhpux \
                fits_delete_keyword.mexhpux \
		fits_read_image_subset.mexhpux \
		fits_write_image_subset.mexhpux

fits_read_image.mexhpux: fits_read_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image.c $(CFLAGS)

fits_read_image_subset.mexhpux: fits_read_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image_subset.c $(CFLAGS)

fits_write_image_subset.mexhpux: fits_write_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image_subset.c $(CFLAGS)

fits_write_image.mexhpux: fits_write_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image.c $(CFLAGS)

fits_read_header.mexhpux: fits_read_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_header.c $(CFLAGS)

fits_write_header.mexhpux: fits_write_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_header.c $(CFLAGS)

fits_delete_keyword.mexhpux: fits_delete_keyword.c $(HEADERS) $(C_FILES)
	$(CC) fits_delete_keyword.c $(CFLAGS)

hpunix_tar_binaries:
	mkdir mfitsio-$(VERSION)-hpunix-bin
	mkdir mfitsio-$(VERSION)-hpunix-bin/mfitsio
	cd mfitsio-$(VERSION)-hpunix-bin/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../*.m ../../*.mexhpux \
	   ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL \
	   ../../Makefile .
	tar cvf - mfitsio-$(VERSION)-hpunix-bin/* \
	        | gzip -c - > mfitsio-$(VERSION)-hpunix-bin.tar.gz
	rm -rf mfitsio-$(VERSION)-hpunix-bin

hpunix_tar_src:
	mkdir mfitsio-$(VERSION)-src
	mkdir mfitsio-$(VERSION)-src/mfitsio
	cd mfitsio-$(VERSION)-src/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../Makefile ../../Makefile.root ../../Makefile.header .; \
	cp ../../*.[chm] ../../ChangeLog ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL .
	tar cvf - mfitsio-$(VERSION)-src/* | gzip -c - > mfitsio-$(VERSION)-src.tar.gz
	rm -rf mfitsio-$(VERSION)-src

hpunix_doc:
	doxygen Doxyfile

hpunix_clean:
	rm *.mexhpux
irix: irix_binaries

irix_binaries: fits_read_image.mexsgi \
                fits_write_image.mexsgi \
                fits_read_header.mexsgi \
                fits_write_header.mexsgi \
                fits_delete_keyword.mexsgi \
		fits_read_image_subset.mexsgi \
		fits_write_image_subset.mexsgi

fits_read_image.mexsgi: fits_read_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image.c $(CFLAGS)

fits_read_image_subset.mexsgi: fits_read_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image_subset.c $(CFLAGS)

fits_write_image_subset.mexsgi: fits_write_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image_subset.c $(CFLAGS)

fits_write_image.mexsgi: fits_write_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image.c $(CFLAGS)

fits_read_header.mexsgi: fits_read_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_header.c $(CFLAGS)

fits_write_header.mexsgi: fits_write_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_header.c $(CFLAGS)

fits_delete_keyword.mexsgi: fits_delete_keyword.c $(HEADERS) $(C_FILES)
	$(CC) fits_delete_keyword.c $(CFLAGS)

irix_tar_binaries:
	mkdir mfitsio-$(VERSION)-irix-bin
	mkdir mfitsio-$(VERSION)-irix-bin/mfitsio
	cd mfitsio-$(VERSION)-irix-bin/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../*.m ../../*.mexsgi \
	   ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL \
	   ../../Makefile .
	tar cvf - mfitsio-$(VERSION)-irix-bin/* \
	        | gzip -c - > mfitsio-$(VERSION)-irix-bin.tar.gz
	rm -rf mfitsio-$(VERSION)-irix-bin

irix_tar_src:
	mkdir mfitsio-$(VERSION)-src
	mkdir mfitsio-$(VERSION)-src/mfitsio
	cd mfitsio-$(VERSION)-src/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../Makefile ../../Makefile.root ../../Makefile.header .; \
	cp ../../*.[chm] ../../ChangeLog ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL .
	tar cvf - mfitsio-$(VERSION)-src/* | gzip -c - > mfitsio-$(VERSION)-src.tar.gz
	rm -rf mfitsio-$(VERSION)-src

irix_doc:
	doxygen Doxyfile

irix_clean:
	rm *.mexsgi
alpha: alpha_binaries

alpha_binaries: fits_read_image.mexalpha \
                fits_write_image.mexalpha \
                fits_read_header.mexalpha \
                fits_write_header.mexalpha \
                fits_delete_keyword.mexalpha \
		fits_read_image_subset.mexalpha \
		fits_write_image_subset.mexalpha

fits_read_image.mexalpha: fits_read_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image.c $(CFLAGS)

fits_read_image_subset.mexalpha: fits_read_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image_subset.c $(CFLAGS)

fits_write_image_subset.mexalpha: fits_write_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image_subset.c $(CFLAGS)

fits_write_image.mexalpha: fits_write_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image.c $(CFLAGS)

fits_read_header.mexalpha: fits_read_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_header.c $(CFLAGS)

fits_write_header.mexalpha: fits_write_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_header.c $(CFLAGS)

fits_delete_keyword.mexalpha: fits_delete_keyword.c $(HEADERS) $(C_FILES)
	$(CC) fits_delete_keyword.c $(CFLAGS)

alpha_tar_binaries:
	mkdir mfitsio-$(VERSION)-alpha-bin
	mkdir mfitsio-$(VERSION)-alpha-bin/mfitsio
	cd mfitsio-$(VERSION)-alpha-bin/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../*.m ../../*.mexalpha \
	   ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL \
	   ../../Makefile .
	tar cvf - mfitsio-$(VERSION)-alpha-bin/* \
	        | gzip -c - > mfitsio-$(VERSION)-alpha-bin.tar.gz
	rm -rf mfitsio-$(VERSION)-alpha-bin

alpha_tar_src:
	mkdir mfitsio-$(VERSION)-src
	mkdir mfitsio-$(VERSION)-src/mfitsio
	cd mfitsio-$(VERSION)-src/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../Makefile ../../Makefile.root ../../Makefile.header .; \
	cp ../../*.[chm] ../../ChangeLog ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL .
	tar cvf - mfitsio-$(VERSION)-src/* | gzip -c - > mfitsio-$(VERSION)-src.tar.gz
	rm -rf mfitsio-$(VERSION)-src

alpha_doc:
	doxygen Doxyfile

alpha_clean:
	rm *.mexalpha
rs6000: rs6000_binaries

rs6000_binaries: fits_read_image.mexrs6 \
                fits_write_image.mexrs6 \
                fits_read_header.mexrs6 \
                fits_write_header.mexrs6 \
                fits_delete_keyword.mexrs6 \
		fits_read_image_subset.mexrs6 \
		fits_write_image_subset.mexrs6

fits_read_image.mexrs6: fits_read_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image.c $(CFLAGS)

fits_read_image_subset.mexrs6: fits_read_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image_subset.c $(CFLAGS)

fits_write_image_subset.mexrs6: fits_write_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image_subset.c $(CFLAGS)

fits_write_image.mexrs6: fits_write_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image.c $(CFLAGS)

fits_read_header.mexrs6: fits_read_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_header.c $(CFLAGS)

fits_write_header.mexrs6: fits_write_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_header.c $(CFLAGS)

fits_delete_keyword.mexrs6: fits_delete_keyword.c $(HEADERS) $(C_FILES)
	$(CC) fits_delete_keyword.c $(CFLAGS)

rs6000_tar_binaries:
	mkdir mfitsio-$(VERSION)-rs6000-bin
	mkdir mfitsio-$(VERSION)-rs6000-bin/mfitsio
	cd mfitsio-$(VERSION)-rs6000-bin/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../*.m ../../*.mexrs6 \
	   ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL \
	   ../../Makefile .
	tar cvf - mfitsio-$(VERSION)-rs6000-bin/* \
	        | gzip -c - > mfitsio-$(VERSION)-rs6000-bin.tar.gz
	rm -rf mfitsio-$(VERSION)-rs6000-bin

rs6000_tar_src:
	mkdir mfitsio-$(VERSION)-src
	mkdir mfitsio-$(VERSION)-src/mfitsio
	cd mfitsio-$(VERSION)-src/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../Makefile ../../Makefile.root ../../Makefile.header .; \
	cp ../../*.[chm] ../../ChangeLog ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL .
	tar cvf - mfitsio-$(VERSION)-src/* | gzip -c - > mfitsio-$(VERSION)-src.tar.gz
	rm -rf mfitsio-$(VERSION)-src

rs6000_doc:
	doxygen Doxyfile

rs6000_clean:
	rm *.mexrs6
win32: win32_binaries

win32_binaries: fits_read_image.dll \
                fits_write_image.dll \
                fits_read_header.dll \
                fits_write_header.dll \
                fits_delete_keyword.dll \
		fits_read_image_subset.dll \
		fits_write_image_subset.dll

fits_read_image.dll: fits_read_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image.c $(CFLAGS)

fits_read_image_subset.dll: fits_read_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_image_subset.c $(CFLAGS)

fits_write_image_subset.dll: fits_write_image_subset.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image_subset.c $(CFLAGS)

fits_write_image.dll: fits_write_image.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_image.c $(CFLAGS)

fits_read_header.dll: fits_read_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_read_header.c $(CFLAGS)

fits_write_header.dll: fits_write_header.c $(HEADERS) $(C_FILES)
	$(CC) fits_write_header.c $(CFLAGS)

fits_delete_keyword.dll: fits_delete_keyword.c $(HEADERS) $(C_FILES)
	$(CC) fits_delete_keyword.c $(CFLAGS)

win32_tar_binaries:
	mkdir mfitsio-$(VERSION)-win32-bin
	mkdir mfitsio-$(VERSION)-win32-bin/mfitsio
	cd mfitsio-$(VERSION)-win32-bin/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../*.m ../../*.dll \
	   ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL \
	   ../../Makefile .
	tar cvf - mfitsio-$(VERSION)-win32-bin/* \
	        | gzip -c - > mfitsio-$(VERSION)-win32-bin.tar.gz
	rm -rf mfitsio-$(VERSION)-win32-bin

win32_tar_src:
	mkdir mfitsio-$(VERSION)-src
	mkdir mfitsio-$(VERSION)-src/mfitsio
	cd mfitsio-$(VERSION)-src/mfitsio; \
	mkdir bin doc; \
	cp ../../bin/adetect ../../bin/mkmf ./bin; \
	cp ../../Makefile ../../Makefile.root ../../Makefile.header .; \
	cp ../../*.[chm] ../../ChangeLog ../../README ../../NEWS ../../COPYING ../../AUTHORS ../../INSTALL .
	tar cvf - mfitsio-$(VERSION)-src/* | gzip -c - > mfitsio-$(VERSION)-src.tar.gz
	rm -rf mfitsio-$(VERSION)-src

win32_doc:
	doxygen Doxyfile

win32_clean:
	rm *.dll
