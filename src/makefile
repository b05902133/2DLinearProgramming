# PROG CXX CXXFLAGS LD from ../makefile
OBJS := $(subst .cpp,.o,$(wildcard *.cpp))

program := ../${PROG}

.PHONY: all clean

all: ${program}

-include $(subst .o,.d,${OBJS})

${program}: ${OBJS}
	${LD} -o $@ $^

%.o: %.d

%.d: %.cpp
	${CXX} ${CXXFLAGS} -MM -MT '$@ $(subst .d,.o,$@)' $< > $@

clean:
	rm *.{o,d}
