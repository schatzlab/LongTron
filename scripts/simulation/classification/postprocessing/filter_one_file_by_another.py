#!/usr/bin/python
#replaces the egrep -f patter_file file command with this one, should be faster on larger datasets

import sys
import re
import getopt
from collections import defaultdict

#make sure the user doesnt unexpected characters
def purge(str):
	return re.sub(r'\'"\n\r',r'',str)

#load the filters from the file
def read_filters(filter_file,case,target_column):
	fin=open(filter_file,"r")
	filters={}
	#filter_ranges=defaultdict(int)
	filter_ranges=defaultdict(int)
	counter=0
	for line in fin:
		counter=counter+1
		#if counter % 1000 == 0:
		#	sys.stdout.write("counter at %s\n" % str(counter))
		#	sys.stdout.flush()
		line=line.rstrip()
		if case:
			line=line.upper()	
		#sys.stdout.write("b4 filter1\n")
		filters[line]=1
		if target_column is None:
			p=re.compile("("+line+")")
			filters[line]=p
		#sys.stdout.write("b4 filter2 key check\n")
		if line not in filter_ranges:
			filter_ranges[line] = 0
		#sys.stdout.write("b4 filter3 add\n")
		filter_ranges[line] += 1
	fin.close()
	return [filters,filter_ranges]

#default split function is just str.split for speed
global split_func
split_func = lambda pattern, line: line.split(pattern)

#filter on a word delimited by a user submitted pattern or the default of \t
def filter_on_words(line,case,word,pattern,complement,filters,filter_ranges,count,target_column):
	#fields=re.split(pattern,line)
	fields = split_func(pattern,line)
	matched=False

	if target_column and int(target_column) < len(fields):
		if fields[int(target_column)] in filters:
			if count and not complement:
				for i in range(filter_ranges[fields[int(target_column)]]):
					print line
			elif not complement:
				print line
			matched=True
	else:
		for field in fields:
			if field in filters:
				if count and not complement:
					for i in range(filter_ranges[field]):
						print line
				elif not complement:
					print line
				matched=True
				break
	if not matched and complement:
		print line			

#filter on the whole line, more like normal egrep -f which also means it's slower but easier to perform and more sensitive
def filter_on_whole_line(line,case,complement,filters,filter_ranges,count):
	matched=False
	for filter in filters.values():
		f=filter.search(line)
		if f:
			if count and not complement:
				for i in range(filter_ranges[f.group(1)]):
					print line
			elif not complement:
				print line
	#if True:
		#this is probably slower than the current approach, but reduce is c compiled due to it being built in
		#if reduce(lambda f1,f2: f1 or f2.search(line),filters,False) not False:
			matched=True
			break
	if not matched and complement:
		print line

def help():
	print """
		This is a substitute for grep/egrep -f for filtering with a list of patterns from a given file, since it tends to take a long time.
		In normal (no -w) mode, the list of patterns will be searched for the target file using regexps for each pattern against each line
		in the file.  Matches spanning lines is not supported at this point.  Also, the speed may not exceed the normal speed of regexps in 
		normal mode.  If you know ahead of time that you're are doing "word" matching, please use -w.
		The following options are supported:
			-f the filter file which contains patterns to match with (required)
			-t the target file which contains the lines we want (required)
			-w matches the whole line in the pattern file against one or more lines which contain the pattern *delimited* by either the start of the line				 or a non-word character (something outside of words, digits, and the underscore "_") and/or the end of the line.
			-i case insensitive matching
			-p submit your own pattern
			-n complement of the search
			-c column of target file to use for the key search (use in conjunction with -w and possibly -p); 0-based
			-r redundant count: the same pattern exists multiple times in the filter file, report a matching line for every instance
			-h prints this help
"""

def main():
	try:
		opts,args=getopt.getopt(sys.argv[1:],"f:inp:wt:c:r",["filter_file","case_insensitive","complement_of_match","filter_pattern","match_word_mode","target_file","target_column","redundant_count"])
	except getopt.GetoptError,err:
		print str(err)
		help()
		sys.exit(-1)
	filter_file=None
	target_file=None
	case=False
	complement=False
	pattern='\t'
	word=None
	target_column=None
	redundant_count=False
	for o, a in opts:
		if o =='-f':
			filter_file=a
		elif o =='-t':
			target_file=a
		elif o == '-i':
			case=True
		elif o == '-n':
			complement=True
		elif o == '-p':
			pattern=purge(a)
		elif o == '-w':
			word=purge(a)
		elif o == '-c':
			target_column=purge(a)
		elif o == '-r':
			redundant_count=True
		else:
			assert False, "unhandled option"
			help()
			sys.exit(-1)

	if filter_file == None or target_file == None:
		print "must submit a filter file AND a target file: -f filter_file -t target_file"
		help()
		sys.exit(-1)

        #if split pattern is more complex than one character, use re.split (slower)
        if len(pattern) > 1:
            sys.stderr.write("using re.split\n")
            p = re.compile(pattern)
            global split_func
            split_func = lambda pattern,line: p.split(line)

	filters=read_filters(filter_file,case,target_column)
	fin = open(target_file,"r") 
	for line in fin:
		line=line.rstrip()
		if case:
			line = line.upper()
		if word != None or pattern != '\t':
			filter_on_words(line,case,word,pattern,complement,filters[0],filters[1],redundant_count,target_column)
		else:
			filter_on_whole_line(line,case,complement,filters[0],filters[1],redundant_count)

if __name__ == '__main__':
	main()

				
		
