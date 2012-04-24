#ifndef KMP_H
#define KMP_H

#include <string.h>

const unsigned char* substring_kmp
    (const unsigned char* haystack, size_t hlen,
     const unsigned char* needle,   size_t nlen);

#endif
