all: perbase pbt

pbt: process_bedtools.cpp
	g++ -O3 -Wno-write-strings -std=c++11 process_bedtools.cpp -o pbt
	#g++ -g -Wno-write-strings -std=c++11 process_bedtools.cpp -o pbt

perbase: perbase.cpp
	g++ -O3 -Wno-write-strings -std=c++11 perbase.cpp -o perbase
	#g++ -g -Wno-write-strings -std=c++11 perbase.cpp -o perbase

test: tests/test.bed perbase
	cat tests/test.bed | ./perbase -c tests/chr1.size -f tests/small_sm.tsv -m -s 5  > test.bed.test
	diff tests/test.bed.expected test.bed.test 

clean: 
	rm perbase
