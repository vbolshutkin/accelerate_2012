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

void check_sa(int * ref_suffix_array, uint len, uint ref_offset) {
	int i; long s1 = 0, s2 = 0;
	for (i= 0; i < len; ++i) {
		s1 += ref_suffix_array[i];
		s2 += i;
	}
	if (s1 == s2 ) {
		printf("sa %ld of length %d is ok\n", ref_suffix_array, len);
	} else {
		printf("WTF? %ld %ld\n", s1, s2);
	}
}

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


		//printf("hi lo %d %d\n", hi, lo);
		uint mid = (hi + lo) / 2;
		uint h_index = suffix_array[mid];
		unsigned char hc = haystack[h_index + n_offset];

		//printf("middle\n");

		if (hc > nc) {
			if (dir != ANY) {
				if (hi == mid) {
					return mid; // FIXME why not break?
				}
				hi = mid;
			} else {
				hi = mid - 1;
			}
		} else if (hc < nc) {
			if (dir != ANY) {
				if (lo == mid) {
					return mid; // FIXME why not break?
				}
				lo = mid;
			} else {
				lo = mid + 1;
			}
		} else {
			//printf("eq before\n");
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
	     const unsigned char* needle, uint noffset, uint subtask_ref_offset) {
	uint hs_offset = substr - haystack;

	if (subtask_ref_offset + hs_offset > 0 && noffset > 0) {
		if (substr[-1] == needle[-1]) {
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
	out->refoffset = subtask_ref_offset + hs_offset + 1;
	out->reflimit = out->refoffset + commonlen - 1;

	// debug output
	printf("%d %d %d %d\n", a.refoffset, a.reflimit, a.othoffset, a.othlimit);
}

void substring_suffix_array(int* suffix_array,
		 const unsigned char* haystack, size_t hlen,
	     const unsigned char* needle,   size_t nlen,
	     uint noffset, uint subtask_ref_offset) {

	int lo = 0, hi = hlen - 1;

	int n_offset = 0;

	int lb = lo, rb = hi;

	while (n_offset < nlen) {
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

	// here can be a lot of matches
	int i;
	for (i = lb; i <= rb; ++i) {
		substring_found(haystack + suffix_array[i], haystack, hlen, needle, noffset, subtask_ref_offset);
	}
}

void solve_task(task_t * st) {

	int * MALLOC_IT(ref_suffix_array, st->ref->len);
	suffixArray(st->ref->ptr, ref_suffix_array, st->ref->len, 256);

	int i;
	for (i = 0; i < st->n_otherseqs; ++i) {
		seq oth = st->otherseqs[i];
		uint oth_iter;
		for (oth_iter = 0; oth_iter < oth.len - minMatchLength;
				++oth_iter) {

			substring_suffix_array(ref_suffix_array, st->ref->ptr,
					st->ref->len, oth.ptr + oth_iter, m,
					oth_iter, st->ref_offset);
		}
	}

	free(ref_suffix_array);
}

void decompose_task(task_t * parent, char parts, task_t * out) {

	int plen = parent->ref->len;

	int i;
	for (i = 0; i < parts; ++i) {
		int offset = i * plen / parts;
		int len = plen / parts + m / 4;
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

void *thread_func_solve(void *vptr_args)
{
	task_t *task = (task_t*)vptr_args;
	solve_task(task);
    return NULL;
}

int main(int argc, char* argv[]) {

	DECLARE_TIME_MEASUREMENT;
	START_TIME_MEASURE;

	n_threads = atol(argv[1]);
	omp_set_dynamic(0);      // запретить библиотеке openmp менять число потоков во время исполнения
	omp_set_num_threads(n_threads); // установить число потоков в 10


	minMatchLength = atol(argv[2]);
	m = minMatchLength / 4;

	read_ref_file(argv[3]);
	read_input_file(argv[4]);

	FINISH_TIME_MEASURE;
	printf("reading took %ld ms\n", TIME_ELAPSED);
	START_TIME_MEASURE;

    task_t task = {&refseq, n_otherseqs, otherseqs, 0};

    char parts = n_threads;


    task_t * MALLOC_IT(subtasks, parts);
    decompose_task(&task, parts, subtasks);
    int i;

	#pragma omp parallel for private(i)
	for (i = 0; i < parts; ++i) {
		solve_task(&subtasks[i]);
	}

    free(subtasks);


	FINISH_TIME_MEASURE;
	printf("OK, took %ld ms\n", TIME_ELAPSED);

	return 0;
}
