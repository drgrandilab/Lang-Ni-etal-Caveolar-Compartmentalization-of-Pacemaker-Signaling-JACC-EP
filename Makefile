DEF= -DCPU
# CXX=g++
CXX=icpx
# CXX=g++ #use g++ if icpc is not available 
CXXFLAGS=-O3 -qopenmp   -std=c++11 -Ilib# -lz  # for icpc
# CXXFLAGS=-O3 -fopenmp   -std=c++11 -Ilib -lz # for g++

# CXX = g++
srcs = $(wildcard lib/*.cpp)
BUILD_DIR=lib/build

objs := $(srcs:%=$(BUILD_DIR)/%.o)
deps := $(objs:.o=.d)

SAN_pace: SAN_pace.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

.PHONY: clean

# $(RM) is rm -f by default
clean:
	$(RM) $(objs) $(deps) SAN_pace

$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(INC_PARAMS)   $(CXXFLAGS)  -MMD -MP -c $< -o $@ 

-include $(deps)
MKDIR_P ?= mkdir -p
