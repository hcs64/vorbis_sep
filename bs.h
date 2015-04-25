#ifndef BS_H_INCLUDED
#define BS_H_INCLUDED

#include <stdint.h>

#include "ogg.h"

struct bitstream_reader;

struct bitstream_reader * create_bitstream_reader(struct ogg_reader *);
uint32_t read_bits(struct bitstream_reader *, int n);
void destroy_bitstream_reader(struct bitstream_reader *);

#endif // BS_H_INCLUDED
