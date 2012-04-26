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

#include <stdint.h>

#include <smmintrin.h> // SSE 4.2

#include "common.h"
#include "bm.h"
#include "kmp.h"
#include "karkkainen.h"

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

uint n_threads;

seq refseq;
int *ref_suffix_array;

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
		printf(" wrong symbol!");
		exit(-1);
	}
}

void print_seq(seq* s) {
	printf("seq %s (len=%d):", s->name, s->len);
	int i;
	for (i = 0; i < s->len; ++i) {

		unsigned char c = s->ptr[i];

		c = decode_char(c % 4);

		printf("%c", c);
	}
	printf("\n");
}

void put_char(char c, unsigned char* ptr, uint index) {

    c = encode_char(c);

	uint step = m;

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
			if (c == '>') {
				goto out;
			}
			if (c == 'G' || c == 'T' || c == 'A' || c == 'C') {
				put_char(c, refseq.ptr, seqlen);
				++seqlen;
			}
		}

	}
out:

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





typedef enum { ANY, LEFT, RIGHT } Direction;

int dihotomy(int* suffix_array,
		 int n_offset,
		 int lo, int hi,
		 const unsigned char* haystack,
	     const unsigned char* needle, Direction dir) {

	unsigned char nc = needle[n_offset];

	while (hi >= lo) {

		// check hi and lo
		if (dir == LEFT || dir == ANY) {
			if (haystack[suffix_array[lo] + n_offset] == nc) {
				return lo;
			}
		}
		if (dir == RIGHT || dir == ANY) {
			if (haystack[suffix_array[hi] + n_offset] == nc) {
				return hi;
			}
		}
		if (dir == LEFT) {
			if (haystack[suffix_array[hi-1] + n_offset] != nc) {
				return hi;
			}
		}
		if (dir == RIGHT) {
			if (haystack[suffix_array[lo+1] + n_offset] != nc) {
				return lo;
			}
		}

		int mid = (hi + lo) / 2;
		int h_index = suffix_array[mid];
		unsigned char hc = haystack[h_index + n_offset];
		if (hc > nc) {
			if (dir != ANY) {
				if (hi == mid) {
					return mid;
				}
				hi = mid;
			} else {
//				if (hi == mid + 1) {
//					break;
//				}
				hi = mid - 1;
			}
		} else if (hc < nc) {
			if (dir != ANY) {
				if (lo == mid) {
					return mid;
				}
				lo = mid;
			} else {
//				if (lo == mid - 1) {
//					break;
//				}
				lo = mid + 1;
			}
		} else {
			if (dir == LEFT) {
				if (haystack[suffix_array[mid-1] + n_offset] != nc) {
					return mid;
				}
				hi = mid;
			} else if (dir == RIGHT) {
				if (haystack[suffix_array[mid+1] + n_offset] != nc) {
					return mid;
				}
				lo = mid;
			} else { // dir == ANY
				return mid;
			}

		}
	}
    return -1;
}

void substring_found(const unsigned char* substr, const unsigned char* haystack, uint hlimit,
	     const unsigned char* needle, uint noffset) {
	uint hs_offset = substr - haystack;

	if (hs_offset > 0 && noffset > 0) {
		if (haystack[hs_offset-1] == needle[-1]) {
			return;
		}
	}

	uint commonlen = m * 4;

	uint i;
	for (i = m; i < hlimit - hs_offset - m * 3; ++i) {
		if (haystack[hs_offset + i] != needle[i]) {
			break;
		}
		++commonlen;
	}

	if (commonlen < minMatchLength) {
		return;
	}

	ans a;
	ans *out = &a;

	out->othoffset = noffset + 1;
	out->othlimit = out->othoffset + commonlen - 1;
	out->refoffset = hs_offset + 1;
	out->reflimit = out->refoffset + commonlen - 1;

	// debug output
	printf("%d %d %d %d\n", a.refoffset, a.reflimit, a.othoffset, a.othlimit);

//	seq dbgref;
//	dbgref.ptr = refseq.ptr + a.refoffset - 1;
//	dbgref.name = "dbgref";
//	dbgref.len = a.reflimit - a.refoffset + 1;
//	print_seq(&dbgref);
//
//	seq dbgoth;
//	dbgoth.ptr = otherseqs[0].ptr + a.othoffset - 1;
//	dbgoth.name = "dbgoth";
//	dbgoth.len = a.othlimit - a.othoffset + 1;
//	print_seq(&dbgoth);


}

void substring_suffix_array(int* suffix_array,
		 const unsigned char* haystack, size_t hlen,
	     const unsigned char* needle,   size_t nlen,
	     uint noffset) {

	int lo = 0, hi = hlen - 1;

	int n_offset = 0;

	int lb = lo, rb = hi;

	while (n_offset < nlen) {
		// we should return lo and hi out of here to save some time
		int d = dihotomy(suffix_array, n_offset, lo, hi, haystack, needle, ANY);
		if (d == -1) {
			return;
		}
		lb = dihotomy(suffix_array, n_offset, lo, d, haystack, needle, LEFT);
		rb = dihotomy(suffix_array, n_offset, d, hi, haystack, needle, RIGHT);
		if (lb == rb && lb != -1) {
			char eq = suffix_array[lb] >= hlen - nlen ||
					memcmp(haystack + suffix_array[lb] + n_offset, needle + n_offset, nlen - n_offset);

			if (eq == 0) {
				break;
			}
		}

		lo = lb;
		hi = rb;
		++n_offset;
	}

	// here are a lot of correct answers?
	int i;
	for (i = lb; i <= rb; ++i) {
		substring_found(haystack + suffix_array[i], haystack, hlen, needle, noffset);
	}
}










int main(int argc, char* argv[]) {

	DECLARE_TIME_MEASUREMENT;
	START_TIME_MEASURE;

	n_threads = atol(argv[1]);
	minMatchLength = atol(argv[2]);
	m = minMatchLength / 4;

	read_ref_file(argv[3]);
	read_input_file(argv[4]);

	FINISH_TIME_MEASURE;
	printf("reading took %ld ms\n", TIME_ELAPSED);
	START_TIME_MEASURE;


	MALLOC_IT(ref_suffix_array, refseq.len);
	suffixArray(refseq.ptr, ref_suffix_array, refseq.len, 256);


//	const unsigned char * ptr = substring_suffix_array(ref_suffix_array, refseq.ptr, refseq.len, otherseqs[0].ptr, 4);
//    if (ptr != NULL) {
//		seq print;
//		print.ptr = ptr;
//		print.len = 20;
//		print.name = "print";
//		print_seq(&print);
//    } else {
//    	printf("test substring not found!\n");
//    }
//
//
//
//	exit(0);

	FINISH_TIME_MEASURE;
	printf("suffix array generation took %ld ms\n", TIME_ELAPSED);
	START_TIME_MEASURE;

	seq oth = otherseqs[0];
	uint oth_iter;
	for (oth_iter = 0; oth_iter < oth.len - minMatchLength; ++oth_iter) {
		substring_suffix_array(ref_suffix_array,
				refseq.ptr, refseq.len,
				oth.ptr + oth_iter, m, oth_iter);
	}

	FINISH_TIME_MEASURE;
	printf("OK, took %ld ms\n", TIME_ELAPSED);

	free(ref_suffix_array);

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
