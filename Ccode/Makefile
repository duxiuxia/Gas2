PROG_NAME = Janus

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
DEPENDS = $(SOURCES:.cpp=.d)

MPICXX ?= mpic++

X_CPPFLAGS = $(CPPFLAGS)
X_CXXFLAGS = $(CXXFLAGS) -O2 -g
X_LDFLAGS  = $(LDFLAGS)
X_LIBS     = $(LIBS) -lnetcdf_c++ -lsqlite3

$(PROG_NAME): $(OBJECTS)
	$(MPICXX) $(X_CXXFLAGS) $(X_LDFLAGS) $(X_LIBS) $(OBJECTS) -o $(PROG_NAME)

%.o: %.cpp
	$(MPICXX) $(X_CPPFLAGS) $(X_CXXFLAGS) -MP -MMD -c -o $@ $<

-include $(DEPENDS)

.PHONY: clean
clean:
	rm -f $(OBJECTS) $(DEPENDS) $(PROG_NAME)
