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

static const int ASCII_OFFSET=48;
static const int CHRM_COL=0;
static const int CHRM_SIZE_COL=1;
static const int START_COL=1;
static const int END_COL=2;
static const int VALUE_COL=3;
static const int NUM_CHRM=25;
static const int BOTH_OPPOSITE_VAL=10;
//1MB per line should be more than enough
static const int LINE_BUFFER_LENGTH=1048576;
//For future use
static int TEST_MODE = -1;
static int SPLICE_MOTIF_MODE = 0;

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
				  idx = chrm_id[3]-ASCII_OFFSET;
			  else
				  idx = atoi(&(chrm_id[3]));
	}
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
		chrm_array[idx] = build_array<T>(chrm_size);
	}
	infile.close();
	return chrm_array;
}

template<typename T>
void delete_nested_array(T** array, long length)
{
	for(long i = 0; i < length; i++)
		delete array[i];
	delete array;
}

template <typename T>
T extract_val(const char* str) { }

template <>
uint8_t extract_val <uint8_t> (const char* str) { return atoi(str); }

template <>
double extract_val <double> (const char* str) { return atof(str); }

template <typename T>
void set_value(int cidx, long start, long end, char strand, T value, T** chrm_array, int strand_col)
{
	//assume BED format start-1 to end
	//splice motifs
	//if(strand_col != -1 && end-start == 2)
	if(SPLICE_MOTIF_MODE)
	{
		value = strand == '+'?1:3;
		chrm_array[cidx][start] = value;
		return;
	}
	for(int i = start; i < end; i++)
	{
		T current_val = chrm_array[cidx][i];
		chrm_array[cidx][i] = value;
	}
}

void output_line(char* line, double summary) { fprintf(stdout,"%s\t%.3f\n", line, summary); } 
//void output_line(char* line, double summary) { fprintf(stdout,"%s\t%.0f\n", line, summary); } 

//about 3x faster than the sstring/string::getline version
template <typename T>
T process_line(char* line, char* delim, int* cidx, long* start, long* end, char* strand, int value_col, int strand_col)
{
	char* line_copy = strdup(line);
	char* tok = strtok(line_copy, delim);
	int i = 0;
	char* chrm;
	T value;
	int last_col = END_COL;
	if(strand_col != -1 || value_col != -1)
		last_col = value_col>strand_col?value_col:strand_col;
	while(tok != NULL)
	{
		if(i > last_col)
			break;
		if(i == CHRM_COL)
		{
			chrm = strdup(tok);
			*cidx = parse_chrm_idx(chrm);
		}
		if(i == START_COL)
			*start = atol(tok);
		if(i == END_COL)
			*end = atol(tok);
		if(i == value_col)
			value = extract_val<T>(tok);
		if(i == strand_col)
			*strand = tok[0];
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
double summarize_region(int* cidx, long* start, long* end, char* strand, T** chrm_array, int strand_col)
{
	double summary = 0.0;
	long length = *end - *start;
	int strand_val = *strand=='+'?1:3;
	//fprintf(stderr,"value: %d line: %c\n",strand_val,*strand);
	for(int i=*start; i < *end; i++)
	{
		if(SPLICE_MOTIF_MODE)
		{
			//only count if we're on the same strand AND the splice motif (2 nucs) is fully within the boundaries
			T value = chrm_array[*cidx][i];
			//if(strand_val == value && !((i == *start && *strand=='-') || (i+1 == *end && *strand=='+')))
			if(strand_val == value && i+1 != *end)
			{
				summary += 1;
			}
		}
		else
			summary	+= chrm_array[*cidx][i];
	}
	double final_val = summary;
	if(strand_col == -1)
		final_val /= length;
	return final_val;
}

template <typename T>
void go(std::string chrm_file, std::string perbase_file, int strand_col)
{
	T** chrm_array = build_chromosome_array<T>(chrm_file);
	int cidx;
	int pcidx = 0;
	long start, end;
	T value;
	char strand;
	
	char* line = new char[LINE_BUFFER_LENGTH];
	size_t length = LINE_BUFFER_LENGTH;
	FILE* fin = fopen(perbase_file.c_str(), "r");
	assert(fin);
	ssize_t bytes_read = getline(&line, &length, fin);
	while(bytes_read != -1)
	{
		//assumes no header
		T value = process_line<T>(strdup(line), "\t", &cidx, &start, &end, &strand, VALUE_COL, strand_col);
		if(cidx != pcidx && pcidx != 0)
		{
			fprintf(stderr,"BUILDING: chr idx %d done\n",pcidx);
			//only do chr1 and 10 for testing
			/*if(pcidx == 10)
				break;*/
		}
		pcidx = cidx;
		set_value<T>(cidx, start, end, strand, value, chrm_array, strand_col);
		bytes_read = getline(&line, &length, fin);
	}
	std::cerr << "building genome wide value array done\n";
	pcidx = 0;
	//now read main file from STDIN line-by-line
	//FILE* fin1 = fopen("q1", "r");
	//bytes_read = getline(&line, &length, fin1);
	bytes_read = getline(&line, &length, stdin);
	char* line_wo_nl = new char[LINE_BUFFER_LENGTH];
	while(bytes_read != -1)
	{
		//get rid of newline
		memcpy(line_wo_nl, line, bytes_read-1);
		line_wo_nl[bytes_read-1]='\0';
		//assumes no header
		T value = process_line<T>(strdup(line), "\t", &cidx, &start, &end, &strand, -1, strand_col);
		if(cidx != pcidx && pcidx != 0)
			fprintf(stderr,"MATCHING: chr idx %d done\n",pcidx);
		pcidx = cidx;
		double summary = summarize_region<T>(&cidx, &start, &end, &strand, chrm_array, strand_col);
		output_line(line_wo_nl, summary);
		//bytes_read = getline(&line, &length, fin1);
		bytes_read = getline(&line, &length, stdin);
	}
	std::cerr << "matching done\n";
	delete line;
	delete line_wo_nl;
	delete_nested_array(chrm_array, NUM_CHRM+1);
}

int main(int argc, char* argv[])
{
	int o;
	std::string perbase_file;
	std::string chrm_file;
	std::string perbase_type = "int";
	int strand_col = -1;
	while((o  = getopt(argc, argv, "f:t:c:s:m")) != -1) 
	{
		switch(o) 
		{
			case 'c': chrm_file = optarg; break;
			case 'f': perbase_file = optarg; break;
			case 't': perbase_type = optarg; break;
			case 's': strand_col = atoi(optarg); break;
			case 'm': SPLICE_MOTIF_MODE = 1; break;
		}
	}
	if(SPLICE_MOTIF_MODE && strand_col == -1)
	{
		std::cerr << "Missing strand column parameter (e.g. -s 5) for splice motif mode" << "\n";
		exit(-1);
	}
	std::cerr << "hello" << " " << perbase_file << " " << perbase_type << "\n";
	switch(perbase_type.c_str()[0])
	{
		case 'd': go<double>(chrm_file, perbase_file, strand_col); break;
		case 'l': go<long>(chrm_file, perbase_file, strand_col); break;
		default:
			go<uint8_t>(chrm_file, perbase_file, strand_col);
	}
}
