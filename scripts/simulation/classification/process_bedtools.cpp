#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <string.h>
#include <unistd.h>
#include <cstdlib>
#include <assert.h>
#include <vector>
#include <cstdint>
#include <fcntl.h>
//#include <stdio.h>
#include <map>

//Parses output of `bedtools intersect -storted -wao -a A.bed -b B.bed`
//where we want to quickly collapse multiple overlaps from B to one entry from A
//as well as allowing for a "wiggle" coordinate match within some window (default 20)

static const int ASCII_OFFSET=48;
static const int CHRM_COL=0;
static const int CHRM_SIZE_COL=1;
static const int START_COL=1;
static const int END_COL=2;
//1MB per line should be more than enough
static const int LINE_BUFFER_LENGTH=20048576;
//For future use
static int TEST_MODE = -1;

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

//about 3x faster than the sstring/string::getline version
//bool matched = process_line<T>(strdup(line), "\t", start2_col, last_val_col, &npline, &num_chars, &err);
//template <typename T>
int process_line(char* line, char* delim, int start2_col, int last_val_col, char** npline, int* num_chars, char** cline, int* num_chars2, int window)
{
	char* line_copy = strdup(line);
	char* tok = strtok(line_copy, delim);
	int i = 0;
	int last_col = last_val_col;
    uint64_t bps = 0;
    uint64_t start=0;
    uint64_t end=0;
    uint64_t start2=0;
    uint64_t end2=0;
    (*num_chars) = 0;
    (*num_chars2) = 0;
    int num_toks = 0;
    int num_toks2 = 0;
	while(tok != NULL)
	{
		if(i > last_val_col)
			break;
        //track length of the line upto the 2nd set of coordinates
        //overlapping from bedtools -wao
        if(i < (start2_col-1))
        {
            num_toks++;
            (*num_chars) += strlen(tok);
        }
        //track length of line upto the end of the first set
        //of coordinates, this is for the unique key
        if(i <= END_COL)
        {
            num_toks2++;
            (*num_chars2) += strlen(tok);
        }
		if(i == START_COL)
			start = atol(tok);
		if(i == END_COL)
			end = atol(tok);
		if(i == start2_col)
			start2 = atol(tok);
		if(i == (start2_col+1))
			end2 = atol(tok);
		if(i == last_val_col)
			//value = extract_val<T>(tok);
            bps = atol(tok);
		i++;
		tok = strtok(NULL,delim);
	}
    //add for number of delimiters
    (*num_chars) += (num_toks-1);
    (*num_chars2) += (num_toks2-1);
    memcpy(*npline, line, *num_chars);
    memcpy(*cline, line, *num_chars2);
    (*npline)[*num_chars]='\0';
    (*cline)[*num_chars2]='\0';
	if(line_copy)
		free(line_copy);
	if(line)
		free(line);
    uint64_t diff1 = abs(start2 - start);
    uint64_t diff2 = abs(end - end2);
    if(bps > 0 && diff1 < window && diff2 < window)
        return 1;
	return 0;
}


typedef std::unordered_map<std::string, bool> coords2matching;
void go(int start2_col, int last_col, int window)
{
	char* line_wo_nl = new char[LINE_BUFFER_LENGTH];
    uint64_t large_count = 0;
    uint64_t small_count = 0;
    uint64_t matching = 0;
    //uint64_t small_matches = 0;
    uint32_t small_matches = 0;
    char* pline = new char[LINE_BUFFER_LENGTH];
    char* npline = new char[LINE_BUFFER_LENGTH];
    char* pcline = new char[LINE_BUFFER_LENGTH];
    char* cline = new char[LINE_BUFFER_LENGTH];
    int num_chars = 0;
    int num_chars2 = 0;
    coords2matching matches; 
    int at_least_one_match = 0;

	//now read main file from STDIN line-by-line
	char* line = new char[LINE_BUFFER_LENGTH];
	size_t length = LINE_BUFFER_LENGTH;
    //FILE* fin1 = fopen("q1", "r");
	//ssize_t bytes_read = getline(&line, &length, fin1);
	ssize_t bytes_read = getline(&line, &length, stdin);
    bool new_line = false;
	while(bytes_read != -1)
	{
		//get rid of newline
		memcpy(line_wo_nl, line, bytes_read-1);
		line_wo_nl[bytes_read-1]='\0';
		//assumes no header
		bool matched = process_line(strdup(line), "\t", start2_col, last_col, &npline, &num_chars, &cline, &num_chars2, window);
        new_line = false;
        //match on just coordinates
        if(large_count > 0 && strcmp(pcline, cline) != 0)
        {
            //multiple line (same set of coordinates) match
            small_count++;
            fprintf(stdout, "%s\t%d\n", pline, at_least_one_match);
            //reset for next group
            at_least_one_match = 0;
            new_line = true;
        }
        //setup for next (or first) group of same coordinates
        if(large_count == 0 || new_line)
        {
            memcpy(pline, npline, num_chars);
            pline[num_chars]='\0';
            memcpy(pcline, cline, num_chars2);
            pcline[num_chars2]='\0';
        }
        //track every individual line
        large_count++;
        if(matched)
        {
            at_least_one_match = 1;
            //if(matches.find(cline) == matches.end())
            //    small_matches++;
            //count every individual line if it matches
            matching++;
            matches[cline] = true;
        }
		//bytes_read = getline(&line, &length, fin1);
		bytes_read = getline(&line, &length, stdin);
	}
    //pick up final group count
    if(large_count > 0)
    {
        small_count++;
        fprintf(stdout, "%s\t%d\n", pline, at_least_one_match);
    }
    //unique coordinate (left) matches
    small_matches = matches.size();
    //group matches = small_matches/small_count, includes duplicate (within group) matches=matching/large_count
    fprintf(stderr,"%.3f\t%.3f\t%u\t%lu\t%lu\t%lu\n",((double)small_matches/small_count),((double)matching/large_count),small_matches,matching,small_count,large_count);
    delete pline;
    delete npline;
    delete pcline;
    delete cline;
	delete line;
	delete line_wo_nl;
}

int main(int argc, char* argv[])
{
	int o;
	int window = 20;
	int start2_col = 7;
    int last_val_col = 9;
	while((o  = getopt(argc, argv, "w:s:l:")) != -1) 
	{
		switch(o) 
		{
			case 'w': window = atoi(optarg); break;
			case 's': start2_col = atoi(optarg); break;
			case 'l': last_val_col = atoi(optarg); break;
		}
    }
	go(start2_col, last_val_col, window);
}
