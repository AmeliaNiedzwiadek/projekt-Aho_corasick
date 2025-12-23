# Makefile dla projektu automaty DNA:
#  - aho_gapped.cpp
#  - aho_corasick.cpp
#  - patterns_generator.cpp
#  - mutations.cpp


CXX = g++
CXXFLAGS = -O3 -std=c++17 -march=native -Wall -Wextra -Wshadow
LDFLAGS = 

# Źródła (każdy plik .cpp kompilowany osobno)
SRCS = aho_gapped.cpp aho_corasick.cpp patterns_generator.cpp mutations.cpp

# Obiekty utworzone z powyższych plików
OBJS = $(SRCS:.cpp=.o)

# Nazwy binarek
TARGETS = aho_gapped aho_corasick patterns_generator mutations

# skompiluj wszystkie programy
all: $(TARGETS)

# Każda z binarek ma jedno źródło:
aho_gapped: aho_gapped.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

aho_corasick: aho_corasick.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

patterns_generator: patterns_generator.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

mutations: mutations.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)


# Automatyczne generowanie .o z .cpp
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

# Czyszczenie projektu
clean:
	rm -f $(OBJS) $(TARGETS) *.dot *.png *.log


# Debug build (wolniejsze, czytelniejsze)
debug: CXXFLAGS = -O0 -g -std=c++17 -Wall -Wextra
debug: clean all