#-----------------------------------------------------------------------------
# GNU makefile for building WCSLIB 4.8
#
# Summary of the main targets
# ---------------------------
#   all:       Do 'make all' in each subdirectory (excluding ./doxygen).
#   check:     Do 'make check' in each subdirectory (compile and run tests).
#   tests:     Do 'make tests' in each subdirectory (compile test programs but
#              don't run them).
#   install:   Do 'make install' in each subdirectory.
#   clean:     Recursively delete intermediate files produced as part of the
#              build, e.g. object modules, core dumps, etc.
#   cleaner:   Recursively clean, and also delete test executables, test
#              input and output, and intermediates produced in compiling the
#              programmers' manual.
#   distclean (or realclean): Recursively delete all platform-dependent files
#              generated during the build, preserving only the programmers'
#              manual and man pages (which are normally provided pre-built).
#              It is the one to use between builds for multiple platforms.
#   cleanest:  Like distclean, but deletes everything that can be regenerated
#              from the source files, including the programmers' manual and
#              man pages, but excluding 'configure'.
#   show:      Print the values of important variables used in this and the
#              other makefiles.
#   writable:  Run chmod recursively to make all sources writable.
#
# Notes:
#   1) If you need to make changes then preferably modify makedefs.in instead.
#
#   2) Refer also to the makefiles in subdirectories, particularly
#      C/GNUmakefile.
#
# Author: Mark Calabretta, Australia Telescope National Facility
# http://www.atnf.csiro.au/~mcalabre/index.html
# $Id: GNUmakefile,v 4.8.1.2 2011/10/04 07:49:19 cal103 Exp cal103 $
#-----------------------------------------------------------------------------
# Get configure settings.
include makedefs

.PHONY : build check chmod clean cleaner cleanest distclean install \
         realclean show tests writable

build :
	-@ for DIR in $(SUBDIRS) ; do \
	     echo '' ; \
	     $(TIMER) ; \
	     $(MAKE) -k -C $$DIR build ; \
	   done

check tests : show
	-@ echo ''
	-@ $(TIMER)
	 @ for DIR in $(SUBDIRS) ; do \
	     echo '' ; \
	     $(MAKE) -i -C $$DIR cleaner ; \
	   done
	-@ echo ''
	 @ for DIR in $(TSTDIRS) ; do \
	     echo '' ; \
	     $(TIMER) ; \
	     $(MAKE) -k -C $$DIR $@ ; \
	   done

install :
	 @ for DIR in $(INSTDIR) ; do \
	     $(MAKE) -k -C $$DIR $@ ; \
	   done
	   $(INSTALL) -m 444 wcsconfig.h wcsconfig_f77.h $(INCDIR)
	-  if [ ! -d "$(DOCDIR)" ] ; then \
	     $(INSTALL) -d -m 2775 $(DOCDIR) ; \
	   fi
	   $(INSTALL) -m 444 CHANGES COPYING* README $(DOCDIR)
	-  if [ ! -d "$(PDFDIR)" ] ; then \
	     $(INSTALL) -d -m 2775 $(PDFDIR) ; \
	   fi
	   $(INSTALL) -m 444 wcslib.pdf $(PDFDIR)
	-  if [ ! -d "$(HTMLDIR)/html" ] ; then \
	     $(INSTALL) -d -m 2775 $(HTMLDIR)/html ; \
	   fi
	   $(INSTALL) -m 444 html/* $(HTMLDIR)/html
	   if [ ! -d "$(LIBDIR)/pkgconfig" ] ; then \
	     $(INSTALL) -d -m 2775 $(LIBDIR)/pkgconfig ; \
	   fi
	   $(INSTALL) -m 444 wcslib.pc $(LIBDIR)/pkgconfig/wcslib.pc

clean cleaner :
	   for DIR in $(SUBDIRS) doxygen ; do \
	     $(MAKE) -C $$DIR $@ ; \
	   done

cleanest distclean realclean :
	   for DIR in $(SUBDIRS) doxygen ; do \
	     $(MAKE) -C $$DIR $@ ; \
	   done
	-  $(RM) *.log
	-  $(RM) -r autom4te.cache autoscan.log
	-  $(RM) confdefs.h conftest.*
	-  $(RM) config.log config.status configure.lineno
	-  $(RM) makedefs wcslib.pc
	-  $(RM) wcsconfig.h wcsconfig_*.h
	-  $(RM) wcslib-*.tar.gz

show ::
	-@ echo 'Subdirectories to be built...'
	-@ echo '  SUBDIRS     := $(SUBDIRS)'
	-@ echo '  TSTDIRS     := $(TSTDIRS)'
	-@ echo ''

writable :
	  chmod -R u+w .

GNUmakefile : makedefs ;

makedefs : makedefs.in config.status
	-@ echo ''
	-@ $(TIMER)
	   ./config.status

config.status : configure
	-@ echo ''
	-@ $(TIMER)
	-@ echo ''
	-@ echo "Environment variables that affect 'configure':"
	-@ echo "  CC       = $${CC-(undefined)}"
	-@ echo "  CFLAGS   = $${CFLAGS-(undefined)}"
	-@ echo "  CPP      = $${CPP-(undefined)}"
	-@ echo "  CPPFLAGS = $${CPPFLAGS-(undefined)}"
	-@ echo "  F77      = $${F77-(undefined)}"
	-@ echo "  FFLAGS   = $${FFLAGS-(undefined)}"
	-@ echo "  LDFLAGS  = $${LDFLAGS-(undefined)}"
	-@ echo ''
	   ./configure --no-create


#-----------------------------------------------------------------------------
# These are for code management.

.PHONY : dist

dist :
	   $(MAKE) -C doxygen cleanest build
	   $(MAKE) -C utils man
	   $(MAKE) distclean
	-@ echo $(WCSLIBPKG)/C/RCS        >  wcslib.X
	-@ echo $(WCSLIBPKG)/C/flexed/RCS >> wcslib.X
	-@ echo $(WCSLIBPKG)/C/test/RCS   >> wcslib.X
	-@ echo $(WCSLIBPKG)/doxygen/RCS  >> wcslib.X
	-@ echo $(WCSLIBPKG)/Fortran/RCS  >> wcslib.X
	-@ echo $(WCSLIBPKG)/Fortran/test/RCS >> wcslib.X
	-@ echo $(WCSLIBPKG)/makedefs     >> wcslib.X
	-@ echo $(WCSLIBPKG)/other        >> wcslib.X
	-@ echo $(WCSLIBPKG)/pgsbox/RCS   >> wcslib.X
	-@ echo $(WCSLIBPKG)/RCS          >> wcslib.X
	-@ echo $(WCSLIBPKG)/TODO         >> wcslib.X
	-@ echo $(WCSLIBPKG)/utils/RCS    >> wcslib.X
	-@ echo $(WCSLIBPKG)/wcslib.T     >> wcslib.X
	-@ echo $(WCSLIBPKG)/wcslib.X     >> wcslib.X
	   rm -f $(WCSLIBPKG).tar.bz2
	   tar cf - -C .. -X wcslib.X $(WCSLIBPKG) | \
	     tar t | \
	     grep -v '/$$' | \
	     sort > wcslib.T
	   rm -f wcslib.X
	   tar cvf $(WCSLIBPKG).tar -C .. -T wcslib.T
	   rm -f wcslib.T
	   bzip2 $(WCSLIBPKG).tar
	   chmod 444 $(WCSLIBPKG).tar.bz2

install_dist :
	   cp -fp $(WCSLIBPKG).tar.bz2 /nfs/ftp/software/wcslib/
	   mv -f  $(WCSLIBPKG).tar.bz2 ../wcslib-releases/
	   (cd /nfs/ftp/software/wcslib/ && \
	     rm -f wcslib.tar.bz2 && \
	     ln -s $(WCSLIBPKG).tar.bz2 wcslib.tar.bz2)
	   cp -fp CHANGES wcslib.pdf ~/public_html/WCS/
	   rsync --archive --delete html/ ~/public_html/WCS/wcslib/

configure : configure.ac
	-@ echo ''
	-@ $(TIMER)
	   autoconf

# Code development overrides must be included specifically before 'configure'
# generates makedefs.
-include flavours
