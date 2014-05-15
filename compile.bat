g++ -c -O2 -fwrapv -o matrix.o matrix.cpp
g++ -c -O2 -fwrapv -o datapoint.o datapoint.cpp
g++ -c -O2 -fwrapv -o kmeans.o kmeans.cpp
g++ -c -O2 -fwrapv -o trie.o trie.cpp
g++ -c -O2 -fwrapv -o tfidf.o tfidf.cpp
g++ -O2 -fwrapv -o test.exe test.cpp kmeans.o datapoint.o matrix.o tfidf.o trie.o
del *.o
