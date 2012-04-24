#ifndef BM_H
#define BM_H

#include <string.h>


const unsigned char* substring_bm
    (const unsigned char* haystack, size_t hlen,
     const unsigned char* needle,   size_t nlen);

#endif
