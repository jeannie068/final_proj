CXX      := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -MMD
LIBS     := -lm
DEBUG_FLAGS = -g -DDEBUG
EXEC     := Staple_Insertion
SRC_DIRS := .\
            parser\
            data_structure\
            alg
SRCS     := $(wildcard $(SRC_DIRS:=/*.cpp))
OBJS     := $(SRCS:.cpp=.o)
DEPS     := $(OBJS:.o=.d)

all: $(EXEC)
	@mkdir -p input output

$(EXEC): $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(EXEC) $(OBJS) $(DEPS)

ifeq (test, $(firstword $(MAKECMDGOALS)))
  TESTCASE := $(word 2, $(MAKECMDGOALS))
  $(eval $(TESTCASE):;@:)
endif

test: $(EXEC)
	@echo Testing with $(TESTCASE)
	./$(EXEC) ../input/$(TESTCASE).txt ../output/$(TESTCASE).txt

.PHONY: all clean test
-include $(DEPS)