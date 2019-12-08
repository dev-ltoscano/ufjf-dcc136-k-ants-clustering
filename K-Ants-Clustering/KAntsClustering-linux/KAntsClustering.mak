EXE = KAntsClustering

$(EXE) : $(EXE).cpp
	g++ -std=c++11 -o $@ $< $(CFLAGS)
