#===============================================================================
# basiclu build configuration
#===============================================================================

# The file has been copied and adapted from SuiteSparse_config.mk contained
# in SuiteSparse version 4.5.3 and available from http://www.suitesparse.com.
# As with the original file, no licensing restrictions apply to this file.

BASICLU_VERSION = 2.1.0

#===============================================================================
# Defaults for any system
#===============================================================================

    #---------------------------------------------------------------------------
    # basiclu root directory
    #---------------------------------------------------------------------------

    BASICLUROOT = $(realpath $(CURDIR))

    #---------------------------------------------------------------------------
    # optimization level
    #---------------------------------------------------------------------------

    OPTIMIZATION ?= -O3 -DNDEBUG

    #---------------------------------------------------------------------------
    # compiler flags for the C compiler
    #---------------------------------------------------------------------------

    # The CF macro is used as a combination of
    # CFLAGS, CPPFLAGS, TARGET_ARCH, and system-dependent settings.
    CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(OPTIMIZATION) -fPIC

    #---------------------------------------------------------------------------
    # compiler
    #---------------------------------------------------------------------------

    # basiclu is written in the C99 standard. The CC99 macro is used to invoke a
    # standard conformant compiler. The default assumes that the system provides
    # the c99 command.

    CC99 ?= c99
    #CC99 ?= gcc -std=c99

    #---------------------------------------------------------------------------
    # required libraries
    #---------------------------------------------------------------------------

    # basiclu requires the -lm (Math) libraries.
    LDLIBS = -lm

    #---------------------------------------------------------------------------
    # shell commands
    #---------------------------------------------------------------------------

    # ranlib, and ar, for generating libraries.  If you don't need ranlib,
    # just change it to RANLIB = echo
    RANLIB = ranlib
    ARCHIVE = $(AR) $(ARFLAGS)

#===============================================================================
# System-dependent configurations
#===============================================================================

    #---------------------------------------------------------------------------
    # determine what system we are on
    #---------------------------------------------------------------------------

    # To disable these auto configurations, use 'make UNAME=custom'

    ifndef UNAME
        ifeq ($(OS),Windows_NT)
            # Cygwin Make on Windows has an $(OS) variable, but not uname.
            # Note that this option is untested.
            UNAME = Windows
        else
            # Linux and Darwin (Mac OSX) have been tested.
            UNAME := $(shell uname)
        endif
    endif

#===============================================================================
# Building the shared and static libraries
#===============================================================================

LIBRARY = libbasiclu
VERSION = 2.1.0
SO_VERSION = 2

ifeq ($(UNAME),Windows)
    # Cygwin Make on Windows (untested)
    AR_TARGET = $(LIBRARY).lib
    SO_PLAIN  = $(LIBRARY).dll
    SO_MAIN   = $(LIBRARY).$(SO_VERSION).dll
    SO_TARGET = $(LIBRARY).$(VERSION).dll
    SO_INSTALL_NAME = echo
else
    # Mac or Linux/Unix
    AR_TARGET = $(LIBRARY).a
    ifeq ($(UNAME),Darwin)
        # Mac
        SO_PLAIN  = $(LIBRARY).dylib
        SO_MAIN   = $(LIBRARY).$(SO_VERSION).dylib
        SO_TARGET = $(LIBRARY).$(VERSION).dylib
        SO_OPTS   = -dynamiclib -compatibility_version $(SO_VERSION) \
                    -current_version $(VERSION) \
                    -shared -undefined dynamic_lookup
        # When a Mac *.dylib file is moved, this command is required
        # to change its internal name to match its location in the filesystem:
        SO_INSTALL_NAME = install_name_tool -id
    else
        # Linux and other variants of Unix
        SO_PLAIN  = $(LIBRARY).so
        SO_MAIN   = $(LIBRARY).so.$(SO_VERSION)
        SO_TARGET = $(LIBRARY).so.$(VERSION)
        SO_OPTS   = -shared -Wl,-soname -Wl,$(SO_MAIN) -Wl,--no-undefined
        # Linux/Unix *.so files can be moved without modification:
        SO_INSTALL_NAME = echo
    endif
endif
