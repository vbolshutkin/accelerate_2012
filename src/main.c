#include <sys/mman.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <limits.h>

#include <stdint.h>

#include <pthread.h>
#include <omp.h>

#include <smmintrin.h> // SSE 4.2

#include "common.h"
#include "bm.h"
#include "kmp.h"
#include "karkkainen.h"



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

typedef struct {
	seq * ref;
	uint n_otherseqs;
	seq * otherseqs;

	uint ref_offset;
} task_t;

uint minMatchLength;
uint m;
uint m4;

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

	uint step = m4;

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







void decompose_task(task_t * parent, char parts, task_t * out) {

	int plen = parent->ref->len;
	int len = plen / parts + m4;

	int i;
	for (i = 0; i < parts; ++i) {
		int offset = i * plen / parts;

		if (len >= plen - offset) {
			len = plen - offset;
		}
		seq * MALLOC_IT(tmprefp, 1);
		seq tmpref = {"", len, parent->ref->ptr + offset};
		*tmprefp = tmpref;

		task_t tmp = {tmprefp, parent->n_otherseqs, parent->otherseqs, offset};
		out[i] = tmp;
	}

}

int main(int argc, char* argv[]) {

	DECLARE_TIME_MEASUREMENT;
	START_TIME_MEASURE;

	n_threads = atol(argv[1]);
	omp_set_dynamic(0);      // запретить библиотеке openmp менять число потоков во время исполнения
	omp_set_num_threads(n_threads); // установить число потоков в 10


	minMatchLength = atol(argv[2]);
	m4 = minMatchLength / 4;
	m = m4  + minMatchLength % 4;

	read_ref_file(argv[3]);
	read_input_file(argv[4]);

	FINISH_TIME_MEASURE;
	printf("reading took %ld ms\n", TIME_ELAPSED);
	START_TIME_MEASURE;


	seq oth = otherseqs[0];

    long offset;
    for (offset = (int) (m - oth.len + 1); offset <= refseq.len - m; ++offset) {
    	uint ref_offset;
    	uint oth_offset;
    	if (offset < 0) {
    		ref_offset = 0;
    		oth_offset = - offset;
    	} else {
    		ref_offset = offset;
    		oth_offset = 0;
    	}
    	unsigned char *ref_ptr = refseq.ptr + ref_offset;
    	unsigned char *oth_ptr = oth.ptr + oth_offset;

    	uint iter = 0;

    	int cnt = 0;
    	while (iter < refseq.len - ref_offset && iter < oth.len - oth_offset) {
    		if (ref_ptr[iter] == oth_ptr[iter]) {
    			++cnt;
    		} else {
    			if (cnt >= m) {
    				printf("%d %d %d %d\n", ref_offset + iter - cnt + 1, ref_offset + iter + 3 * m4,  oth_offset + iter - cnt + 1, oth_offset + iter + 3 * m4);
    			}
    			cnt = 0;

//    			// fast forward works slowly
//    			uint back_iter = iter + m;
//    			while (ref_ptr[back_iter] == oth_ptr[back_iter]) {
//    				--back_iter;
//    				++cnt;
//    			}

    		}
    		++iter;
    	}
    	if (cnt >= m) {
			printf("%d %d %d %d\n", ref_offset + iter - cnt + 1, ref_offset + iter,  oth_offset + iter - cnt + 1, oth_offset + iter);
		}
    }





	FINISH_TIME_MEASURE;
	printf("OK, took %ld ms\n", TIME_ELAPSED);

	return 0;
}
