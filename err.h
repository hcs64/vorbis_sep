#ifndef ERR_H_INCLUDED
#define ERR_H_INCLUDED

#include <stdlib.h>

#define expect(cond) do { if (!(cond)) { fprintf(stderr, __FILE__ " %d: " #cond " failed\n", __LINE__); exit(EXIT_FAILURE); } } while (0)

#endif // ERR_H_INCLUDED
