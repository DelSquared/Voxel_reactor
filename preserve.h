#include <stdio.h>
#include <stdlib.h>

/*
 * Writes a struct to a binary file.
 * Returns 0 on success, -1 on failure.
 */

#ifndef PRESERVEC
#define PRESERVEC

int preserve(const char *filename, const void *data, size_t size)
{
    FILE *fp = fopen(filename, "wb");
    if (!fp)
        return -1;

    size_t written = fwrite(data, 1, size, fp);
    fclose(fp);

    return (written == size) ? 0 : -1;
}

/*
 * Reads a struct from a binary file.
 * Returns 0 on success, -1 on failure.
 */
int unpreserve(const char *filename, void *data, size_t size)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp)
        return -1;

    size_t read = fread(data, 1, size, fp);
    fclose(fp);

    return (read == size) ? 0 : -1;
}


#endif
