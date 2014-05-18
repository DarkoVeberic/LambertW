CPPFLAGS := -I.
#CXXFLAGS := -Wall -Wextra -ggdb3 -O0 -fno-inline -pipe
CXXFLAGS := -O2
#CXXFLAGS := -Wall -Wextra -ggdb3 -O2 -pipe
#CXXFLAGS := -Wall -Wextra -ggdb3 -march=native -Ofast -pipe
LDFLAGS := $(CXXFLAGS)

SUFFIXES := .o .cc .cxx

EXES := $(basename $(wildcard *.cxx))
OBJS := $(patsubst %.cc, %.o, $(wildcard *.cc))
DEPS := $(patsubst %.o, %.P, $(OBJS)) $(addsuffix .P, $(EXES))

define cxx_compile_with_dependency_creation
  $(COMPILE.cc) -MD -o $@ $<
  @sed -e 's|.*:|$*.o:|' <$*.d >$*.P
  @sed -e 's/.*://' -e 's/\\$$//' <$*.d | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >>$*.P
  @rm -f $*.d
endef

define cxx_link_rule
  $(CXX) $^ $(LDFLAGS) -o $@
endef

%.o: %.cc
	$(call cxx_compile_with_dependency_creation)

%.o: %.cxx
	$(call cxx_compile_with_dependency_creation)

%: %.o
	$(call cxx_link_rule)

.PHONY: all
all: $(EXES)

lambertw: lambertw.o $(OBJS)

.PHONY: check
check: $(basename $(wildcard test_*.cxx))
	for t in $^ ; do echo $$t ; ./$$t || exit $$? ; done

.PHONY: clean
clean:
	- $(RM) $(OBJS) $(addsuffix .o, $(EXES)) $(EXES) *.P

-include $(DEPS)
