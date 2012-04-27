#ifndef COMMON_H
#define COMMON_H

#include <sys/time.h>

#define MALLOC_IT(ptr, sz) (ptr) = malloc((sz) * sizeof(*(ptr)))
#define REALLOC_IT(ptr, sz) (ptr) = realloc(ptr, (sz) * sizeof(*(ptr)))

#define FALSE 0
#define TRUE 1

#define FREE_ARRAY(ind, sz, arr, free_func) \
	int ind; \
	for (ind = 0; ind < sz; ++ind) { \
		free_func(&arr[ind]); \
	}


#define DECLARE_TIME_MEASUREMENT struct timeval start, finish
#define START_TIME_MEASURE gettimeofday(&start, 0)
#define FINISH_TIME_MEASURE gettimeofday(&finish, 0)
#define TIME_ELAPSED ((finish.tv_sec-start.tv_sec)*1000 + (finish.tv_usec-start.tv_usec)/1000)


#endif
