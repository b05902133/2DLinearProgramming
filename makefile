CXX      := g++
CXXFLAGS := -std=c++11
LD       := g++
PROG     := 2dlp
CTAGS    := ctags

srcDir  := src
testDir := test

export PROG CXX CXXFLAGS LD

.PHONY: all release debug ${PROG} test tags clean

all: debug

${PROG}:
	${MAKE} -C ${srcDir}

release: CXXFLAGS += -O2 -DNDEBUG
release: ${PROG}

debug: CXXFLAGS += -g -Wall
debug: ${PROG}

test: CXXFLAGS += -DNDEBUG
test: debug
	${MAKE} -C ${testDir}

tags:
	${CTAGS} -R ${srcDir}

clean:
	-${MAKE} -C ${srcDir} clean
	-${MAKE} -C ${testDir} clean
	rm ${PROG}
