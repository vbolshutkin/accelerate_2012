#include <limits.h>
#include <string.h>
#include "bm.h"

/* Вход: needle+nlen - что искать
 *   offset - позиция конца суффикса; suffixlen - его длина
 * Выход: 1, если есть совпадение суффикса (и нет совпадения более длинного суффикса)
 */
static int suffix_match(const unsigned char *needle, size_t nlen, size_t offset, size_t suffixlen)
{
    if (offset > suffixlen)
        return needle[offset - suffixlen - 1] != needle[nlen - suffixlen - 1] &&
            memcmp(needle + nlen - suffixlen, needle + offset - suffixlen, suffixlen) == 0;
    else
        return memcmp(needle + nlen - offset, needle, offset) == 0;
}

static size_t max(size_t a, size_t b)
{
    return a > b ? a : b;
}

/* Вход: haystack - где искать, needle - что искать.
 *   hlen - длина haystack, nlen - длина needle
 * Выход: указатель на первое вхождение needle в haystack, либо NULL
 */
const unsigned char* substring_bm
    (const unsigned char* haystack, size_t hlen,
     const unsigned char* needle,   size_t nlen)
{
    size_t skip[nlen];          /* Таблица суффиксов */
    size_t occ[UCHAR_MAX + 1]; /* Таблица стоп-символов */

    if(nlen > hlen || nlen <= 0 || !haystack || !needle)
        return NULL;

    /* ПОСТРОЕНИЕ ТАБЛИЦЫ СТОП-СИМВОЛОВ */
    size_t a;
    /* Заполняем -1 (в Си нумерация символов с 0!!) */
    for(a=0; a<UCHAR_MAX+1; ++a)
        occ[a] = -1;

    /* В таблицу occ записывается последнее вхождение в needle  */
    /* (исключая последний символ) */
    for(a = 0; a < nlen - 1; ++a)
        occ[needle[a]] = a;

    /* ПОСТРОЕНИЕ ТАБЛИЦЫ СУФФИКСОВ */
    /* Упрощённый вариант. Можно реализовать быстрее. */
    for(a = 0; a < nlen; ++a)
    {
    	size_t offs = nlen;
        while(offs && !suffix_match(needle, nlen, offs, a))
            --offs;
        skip[nlen - a - 1] = nlen - offs;
    }

    /* ПОИСК */
    size_t hpos;
    for(hpos = 0; hpos <= hlen - nlen; )
    {
        size_t npos = nlen - 1;
        /* Сверка needle и haystack с конца */
        while(needle[npos] == haystack[npos + hpos])
        {
            if(npos == 0)
                return haystack + hpos;

            --npos;
        }
        /* Не совпало! */
        signed int second = npos - occ[haystack[npos + hpos]];
        if (second < 0) second = 0;
        hpos += max(skip[npos], second);
        /*          ^^^         ^^^^                               */
        /*        суффикс     стоп-символ                          */
    }
    return NULL;
}
