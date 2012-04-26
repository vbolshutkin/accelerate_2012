#ifndef COMMON_H
#define COMMON_H

#define MALLOC_IT(ptr, sz) (ptr) = malloc((sz) * sizeof(*(ptr)))
#define REALLOC_IT(ptr, sz) (ptr) = realloc(ptr, (sz) * sizeof(*(ptr)))

#define FALSE 0
#define TRUE 1

#define FREE_ARRAY(ind, sz, arr, free_func) \
	int ind; \
	for (ind = 0; ind < sz; ++ind) { \
		free_func(&arr[ind]); \
	}

#endif
