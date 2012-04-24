#include "kmp.h"

#include <stdlib.h>

const unsigned char* substring_kmp (
		 const unsigned char* haystack, size_t hlen,
	     const unsigned char* needle,   size_t nlen)
{
        int i, j;
        int *d =(int*)malloc(nlen*sizeof(int)); /* динамический массив длины М*/
        /* Вычисление префикс-функции */
        d[0]=0;
        for(i=1,j=0;i<nlen;i++)
        {
                while(j>0 && needle[j]!=needle[i]) j = d[j-1];
                if(needle[j]==needle[i])
                        j++;
                        d[i]=j;
        }
        /* поиск */
        for(i=0,j=0;i<hlen; i++)
        {
                while(j>0 && needle[j]!=haystack[i])
                        j=d[j-1];
                if(needle[j]==haystack[i])
                        j++;
                if (j==nlen)
                {
                        free (d); /* освобождение памяти массива d */
                        return &haystack[i-j+1];
                }
        }
        free (d); /* освобождение памяти массива d */
        return NULL;
}
