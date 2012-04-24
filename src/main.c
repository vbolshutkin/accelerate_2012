#include <sys/mman.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/time.h>
#include <limits.h>

#include "bm.h"
#include "kmp.h"

#define MALLOC_IT(ptr, sz) (ptr) = malloc((sz) * sizeof(*(ptr)))
#define REALLOC_IT(ptr, sz) (ptr) = realloc(ptr, (sz) * sizeof(*(ptr)))

#define FALSE 0
#define TRUE 1

#define DECLARE_TIME_MEASUREMENT struct timeval start, finish
#define START_TIME_MEASURE gettimeofday(&start, 0)
#define FINISH_TIME_MEASURE gettimeofday(&finish, 0)
#define TIME_ELAPSED ((finish.tv_sec-start.tv_sec)*1000 + (finish.tv_usec-start.tv_usec)/1000)


typedef struct {
	char *name;
	unsigned int len;
	unsigned char *ptr;
} seq;

typedef struct {
	uint refoffset;
	uint reflimit;
	uint othoffset;
	uint othlimit;
} ans;

uint minMatchLength;
uint m;

seq refseq;

uint n_otherseqs = 0;
seq* otherseqs = NULL;

char encode_char(char c) {
	switch (c) {
	case 'G':
		return 0;
	case 'T':
		return 1;
	case 'A':
		return 2;
	case 'C':
		return 3;
	default:
		exit(-1);
	}
}

char decode_char(char c) {
	switch (c) {
	case 0:
		return 'G';
	case 1:
		return 'T';
	case 2:
		return 'A';
	case 3:
		return 'C';
	default:
		exit(-1);
	}
}

void print_seq(seq* s) {
	printf("seq %s (len=%d):", s->name, s->len);
	int i;
	for (i = 0; i < s->len; ++i) {

		char c = s->ptr[i] % 4;

		c = decode_char(c);

		printf("%c", c);
	}
	printf("\n");
}

void put_char(char c, unsigned char* ptr, uint index) {

    c = encode_char(c);

	uint step = minMatchLength / 4;

	char ix;
	for (ix = 1; (ix < 4) && (index >= step * ix); ++ix) {
		ptr[index - step * ix] += ((c == 0) ? 0 : c << (2 * ix));
	}

	ptr[index] = c;
}

void read_ref_file(char* fname) {
	FILE* fd = fopen(fname, "r");
	fseek(fd, 0, SEEK_END);
	int size = ftell(fd);
	fseek(fd, 0, SEEK_SET);

	fgetc(fd);
	char line[1024];
	char *seqname = fgets(line, 1024, fd);
	if (NULL == seqname) {
		printf("FGETS NULL");
		exit(-1);
	}
	int namelen = strlen(seqname);
	seqname[namelen - 1] = 0;

	MALLOC_IT(refseq.ptr, size - namelen - 2);
	uint seqlen = 0;

	short read;
	unsigned char buffer[1024];
	while (0 < (read = fread(buffer, 1, 1024, fd))) {
		short i;
		for (i = 0; i < read; ++i) {
			char c = buffer[i];
			if (c == 'G' || c == 'T' || c == 'A' || c == 'C') {
				put_char(c, refseq.ptr, seqlen);
				++seqlen;
			}
		}

	}

	refseq.ptr[seqlen] = 0;
	refseq.len = seqlen;

	MALLOC_IT(refseq.name, namelen - 1);
	strcpy(refseq.name, seqname);

	REALLOC_IT(refseq.ptr, seqlen + 1);
}

inline void setup_seq(unsigned char *curptr, uint seqlen, char* seqname) {

	++n_otherseqs;
	REALLOC_IT(otherseqs, n_otherseqs);
	seq *s = otherseqs + n_otherseqs - 1;

	s->ptr = curptr;
	s->ptr[seqlen] = 0;
	s->len = seqlen;
	MALLOC_IT(s->name, strlen(seqname) + 1);
	strcpy(s->name, seqname);
}

void read_input_file(char* fname) {
	FILE* fd = fopen(fname, "r");
	fseek(fd, 0, SEEK_END);
	int size = ftell(fd);
	fseek(fd, 0, SEEK_SET);

	unsigned char *globptr, *curptr;

	MALLOC_IT(globptr, size);
	curptr = globptr;

	uint globlen = 0, seqlen = 0;

	char state_reading_name = FALSE;

	short read;
	unsigned short namelen = 0, namestart;
	char *seqname = NULL;
	char buffer[1024];
	while (0 < (read = fread(buffer, 1, 1023, fd))) {
		buffer[1023] = 0;

		short i;
		for (i = 0; i < read; ++i) {
			char c = buffer[i];

			if (c == '>') {
				state_reading_name = TRUE;
				if (seqlen > 0) {
					globlen += seqlen;
					setup_seq(curptr, seqlen, seqname);
				}

				if (NULL != seqname) {
					seqname[0] = 0;
					namelen = 0;
				}
			}

			if (state_reading_name == TRUE) {
				namestart = i + 1;
				for (++i; i < read; ++i) {
					char nc = buffer[i];
					if (nc == '\n') {
						buffer[i] = 0;
						state_reading_name = FALSE;
						curptr += seqlen;
						seqlen = 0;
						break;
					}
				}

				REALLOC_IT(seqname, namelen + i - namestart + 1);
				strcat(seqname, buffer + namestart);
				namelen = strlen(seqname);
			} else {
				if (c == 'G' || c == 'T' || c == 'A' || c == 'C') {
					put_char(c, curptr, seqlen);
					++seqlen;
				}
			}
		}
	}

	if (seqlen > 0) {
		globlen += seqlen;
		setup_seq(curptr, seqlen, seqname);
	}

}

// sa - needle
// sb - haystack
char find_lcs(seq *sa, seq *sb, uint aoffset, uint alimit, uint boffset, uint blimit, ans *out) {
	unsigned char *needle = &sa->ptr[aoffset];
	unsigned char *haystack = &sb->ptr[boffset];

	const unsigned char *bmres = substring_kmp(haystack, blimit - boffset, needle, m);
	if (NULL != bmres) {
		uint hs_offset = bmres - haystack;

//		printf("bmres pos: %d\n", hs_offset);
		uint commonlen = m * 4;

		uint i;
		for (i = m; i < blimit - hs_offset; ++i) {
			if (haystack[hs_offset + i] != needle[i]) {
				break;
			}
			++commonlen;
		}

		out->refoffset = aoffset + 1;
		out->reflimit = out->refoffset + commonlen - 1;
		out->othoffset = boffset + hs_offset + 1;
		out->othlimit = out->othoffset + commonlen - 1;

//		printf("commonlen: %d\n", commonlen);
//
//		seq s;
//		s.len = commonlen;
//		s.ptr = haystack + hs_offset;
//		s.name = "substring";
//		print_seq(&s);

		return TRUE;
	}
	return FALSE;
}


int main(int argc, char* argv[]) {

	DECLARE_TIME_MEASUREMENT;
	START_TIME_MEASURE;



//	// first command line argument : number of worker threads to use
//	// not used in this program at this time
//
//	// first command line argument : window size
	minMatchLength = atol(argv[2]);
	m = minMatchLength / 4;

	read_ref_file(argv[3]);
//	print_seq(&refseq);

	read_input_file(argv[4]);
//	int i;
//	for (i = 0; i < n_otherseqs; ++i) {
//		print_seq(&otherseqs[i]);
//	}

	uint ref_iter;
	for (ref_iter = 0; ref_iter < refseq.len - minMatchLength; ++ref_iter) {
		//printf("ref_iter: %d\n", ref_iter);
		ans a;
		char f;
		uint other_shift = 0;
		do {
			f = find_lcs(&refseq, &otherseqs[0], ref_iter, refseq.len, other_shift, otherseqs[0].len, &a);
			if (f) {
				other_shift = a.othlimit;
				printf("%d %d %d %d\n", a.refoffset, a.reflimit, a.othoffset, a.othlimit);
			}
		} while(f);
	}

	FINISH_TIME_MEASURE;
	printf("OK, took %ld ms\n", TIME_ELAPSED);


	return 0;
}

//pair<string,char*> getSequence(ifstream& s){
//	string sequence;
//	string name;
//
//	if(s.get()=='>')
//		getline(s,name);
//
//	while (!s.eof()){
//		int c = s.get(); //get char from reference Sequence file
//		if(c == '>'){
//			s.unget();
//			break;
//		}
//		else if (c == 'G' || c=='T' || c=='A' || c=='C'){
//			sequence.push_back(c);
//		}
//		// ignore every other chars (newline, etc)
//	}
//
//	return make_pair(name,sequence.c_str());
//}

//	// second command line argument : reference sequence file
//	ifstream refSeqFile(argv[3],ios::in);
//	pair<string,char*> ref = getSequence(refSeqFile);
//	string refSeqName = ref.first;
//	char* refSeq = ref.second;
//
//	// following command line arguments : other files containing sequences
//	// result is stored in an associative array ordered by title
//	map<string,string> otherSequences;
//	for(int i=4; i<argc; i++){ // iterate over command arguments
//		ifstream seqFile(argv[i],ios::in);
//		while(!seqFile.eof()){
//			pair<string,string> other = getSequence(seqFile);
//			otherSequences[other.first] = other.second;
//		}
//	}
//
//	// compare other sequences to reference sequence
//	// iterate over other sequences
//	for(map<string,string>::iterator sequencesIter = otherSequences.begin(); sequencesIter!=otherSequences.end(); sequencesIter++){
//		// output sequence name
//		cout << sequencesIter->first << "\n";
//		string otherSeq = sequencesIter->second;
//
//		// L[i][j] will contain length of the longest substring
//		// ending by positions i in refSeq and j in otherSeq
//		size_t **L = new size_t*[refSeq.length()];
//		for(size_t i=0; i<refSeq.length();++i)
//			L[i] = new size_t[otherSeq.length()];
//
//		// iteration over the characters of the reference sequence
//		for(size_t i=0; i<refSeq.length();i++){
//			// iteration over the characters of the sequence to compare
//			for(size_t j=0; j<otherSeq.length();j++){
//				// if the characters are the same,
//				// increase the consecutive matching score from the previous cell
//				if(refSeq[i]==otherSeq[j]){
//					if(i==0 || j==0)
//						L[i][j]=1;
//					else
//						L[i][j] = L[i-1][j-1] + 1;
//				}
//				// or reset the matching score to 0
//				else
//					L[i][j]=0;
//			}
//		}
//
//		// output the matches for this sequence
//		// length must be at least minMatchLength
//		// and the longest possible.
//		for(size_t i=0; i<refSeq.length();i++){
//			for(size_t j=0; j<otherSeq.length();j++){
//
//				if(L[i][j]>=minMatchLength) {
//					//this match can be shifted on i and j
//					if(i+1<refSeq.length() && j+1<otherSeq.length() && L[i][j]<=L[i+1][j+1])
//						continue;
//					//this match can be shifted on j
//					if(i<refSeq.length() && j+1<otherSeq.length() && L[i][j]<=L[i][j+1])
//						continue;
//					//this match can be shifted on i
//					if(i+1<refSeq.length() && j<otherSeq.length() && L[i][j]<=L[i+1][j])
//						continue;
//					cout << i-L[i][j]+2 << " " << i+1 << " " << j-L[i][j]+2 << " " << j+1 << "\n";
//
//					// output the matching sequences for debugging :
//					//cout << refSeq.substr(i-L[i][j]+1,L[i][j]) << "\n";
//					//cout << otherSeq.substr(j-L[i][j]+1,L[i][j]) << "\n";
//				}
//			}
//		}
//
//		for(size_t i=0; i<refSeq.length();++i)
//			delete[] L[i];
//		delete[] L;
//	}
