#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <unistd.h>
#include <cstdlib>
#include <assert.h>
#include <vector>
#include <cstdint>
#include <fcntl.h>

static const int CHRM_COL=0;
static const int CHRM_SIZE_COL=1;
static const int START_COL=1;
static const int END_COL=2;
static const int VALUE_COL=3;
static const int NUM_CHRM=25;

template<typename T>
T* build_array(long size)
{
	//initialize to all 0's
	T* a = new T[size]();
	return a;
}

//multiple sources for this kind of tokenization, one which was useful was:
//https://yunmingzhang.wordpress.com/2015/07/14/how-to-read-file-line-by-lien-and-split-a-string-in-c/
void split_string(std::string line, char delim, std::vector<std::string>* tokens)
{
	tokens->clear();
	std::stringstream ss(line);
	std::string token;
	while(getline(ss, token, delim))
	{
		tokens->push_back(token);
	}
}

int parse_chrm_idx(const char* chrm_id)
{
	assert(strncmp(chrm_id,"chr",3)==0);
	int len = strlen(chrm_id);
	assert(len <= 5);
	//extract chrm ID as an integer
	//hardcoded to Human at this point
	int idx = -1;
	switch(chrm_id[3])
	{
		case 'X': idx = 23; break;
		case 'Y': idx = 24; break;
		case 'M': idx = 25; break;
		default:
			  if(len < 5)
				  //subtract ascii offset
				  idx = chrm_id[3]-48;
			  else
				  idx = atoi(&(chrm_id[3]));
	}
	//printf("chr idx %s %d\n",chrm_id,idx);
	return idx;
}

template<typename T>
T** build_chromosome_array(std::string chrm_file)
{
	std::ifstream infile(chrm_file.c_str());
	assert(infile);
	std::string line;
	T** chrm_array = new T*[NUM_CHRM+1];
	std::vector<std::string> tokens;
	while(getline(infile, line))
	{
		split_string(line, '\t', &tokens);
		const char* chrm_id = tokens.at(CHRM_COL).c_str();
		int idx = parse_chrm_idx(chrm_id);
		long chrm_size = atol(tokens.at(CHRM_SIZE_COL).c_str());
		//printf("chr idx %s %d %u\n",chrm_id,idx,chrm_size);
		chrm_array[idx] = build_array<T>(chrm_size);
	}
	infile.close();
	return chrm_array;
}


template <typename T>
T extract_val(const char* str) { }

template <>
uint8_t extract_val <uint8_t> (const char* str) { return atoi(str); }

template <>
double extract_val <double> (const char* str) { return atof(str); }

template <typename T>
void set_value(int cidx, long start, long end, T value, T** chrm_array)
{
	//assume BED format start-1 to end
	for(int i = start; i < end; i++)
	{
		chrm_array[cidx][i] = value;
	}
}
/*
template <typename T>
void output_line(std::string line, T summary) { } 

template <>
void output_line<uint8_t>(std::string line, uint8_t summary) { printf("%s\t%d\n",line, summary); } 

template <>
void output_line<double>(std::string line, double summary) { printf("%s\t%.3f\n",line, summary); } 
*/

//about 5x faster than the sstring/string::getline version
template <typename T>
T process_line(char* line, char* delim, int* cidx, long* start, long* end)
{
	char* line_copy = strdup(line);
	char* tok = strtok(line_copy, delim);
	int i = 0;
	char* chrm;
	T value;
	while(tok != NULL)
	{
		if(i == CHRM_COL)
		{
			chrm = strdup(tok);
			*cidx = parse_chrm_idx(chrm);
		}
		if(i == START_COL)
			*start = atol(tok);
		if(i == END_COL)
			*end = atol(tok);
		if(i == VALUE_COL)
			value = extract_val<T>(tok);
		i++;
		tok = strtok(NULL,delim);
	}
	if(line_copy)
		free(line_copy);
	if(chrm)
		free(chrm);
	if(line)
		free(line);
	return value;
}

template <typename T>
void go(std::string chrm_file, std::string perbase_file)
{
	T** chrm_array = build_chromosome_array<T>(chrm_file);
	//std::ifstream infile(perbase_file.c_str());
	//assert(infile);
	int cidx;
	int pcidx = 0;
	long start, end;
	T value;
	
	char* line = NULL;
	size_t length = 0;
	FILE* fin = fopen(perbase_file.c_str(), "r");
	assert(fin);
	ssize_t bytes_read = getline(&line, &length, fin);
	while(bytes_read != -1)
	{
		//assumes no header
		T value = process_line<T>(strdup(line), "\t", &cidx, &start, &end);
		if(cidx != pcidx)
			printf("chr idx %d done\n",pcidx);
		pcidx = cidx;
		set_value<T>(cidx, start, end, value, chrm_array);
		bytes_read = getline(&line, &length, fin);
		//printf("%d %d\n",cidx,chrm_array[cidx][end-1]);
	}
	std::cout << "building done" << " " << perbase_file << "\n";
	/*std::string line;
	//now read main file from STDIN line-by-line
	while(getline(std::cin, line))
	{
		std::vector<std::string> tokens;
		split_string(line, '\t', &tokens);
		cidx = parse_chrm_idx(tokens.at(CHRM_COL).c_str());
		start = atol(tokens.at(START_COL).c_str());
		end = atol(tokens.at(END_COL).c_str());
		T summary = summarize_region<T>(cidx, start, end, chrm_array);
		output_line<T>(line, summary);
	}
	std::cout << "matching done" << " " << perbase_file << "\n";*/
}

int main(int argc, char* argv[])
{
	int o;
	std::string perbase_file;
	std::string chrm_file;
	std::string perbase_type = "int";
	while((o  = getopt(argc, argv, "f:t:c:")) != -1) 
	{
		switch(o) 
		{
			case 'c': chrm_file = optarg; break;
			case 'f': perbase_file = optarg; break;
			case 't': perbase_type = optarg; break;
		}
	}
	std::cout << "hello" << " " << perbase_file << " " << perbase_type << "\n";
	go<uint8_t>(chrm_file, perbase_file);
}
