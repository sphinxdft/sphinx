# ----------------------------------------------------------------------------
#
#           The general purpose cross platform C/C++ framework
#
#                       S x A c c e l e r a t e
#
#           Home:       https://www.sxlib.de
#           License:    Apache 2
#           Authors:    see src/AUTHORS
#
# ----------------------------------------------------------------------------

NAME      = "sxdemo6"
MODE      = "libOnly"
VERSION   = "1:0:0"
NAME_HASH = 3076441707 # sxloghash



SX_DEPLIBS  = $(SXPARSERKIT) $(SXGRAPH) $(SXFS) $(SXUTIL) $(SXNONSTD)
NUM_DEPLIBS = $(REGEX) $(PTHREADS)

SUBDIRS =


# ----------------------------------------------------------------------------

include $(top_srcdir)/system/make.am.filelist

AM_CXXFLAGS += @CXX_PEDANTIC@

# --- LEX and YACC support
AUX_SRC       = SxDemo6Parser.lex.cpp
AUX_SRC      += SxDemo6Parser.tab.cpp SxDemo6Parser.tab.hpp SxDemo6Parser.grammar
BUILT_SOURCES = $(AUX_SRC)
CLEANFILES += $(AUX_SRC) SxDemo6Parser.output

%.lex.lo: AM_CXXFLAGS += @LEX_CXXFLAGS@ -Wno-unused
%.tab.lo: AM_CXXFLAGS += @LEX_CXXFLAGS@

%.lex.cpp: %.lpp %.tab.cpp %.grammar
	export LANG="" ; $(LEX) -o $@ $<

%.tab.cpp %.tab.hpp: %.ypp
	export LANG="" ; $(YACC) --report=itemset $< -o $@

%.grammar: %.tab.cpp
	@head -1 $(@:.grammar=.output) 2>/dev/null | grep '^Grammar' > $@ \
      || ( echo "Conflicts in $(@:.grammar=.ypp). Inspect $(@:.grammar=.output)."; rm $@; exit 1)
