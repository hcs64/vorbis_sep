#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

#include "bs.h"
#include "ogg.h"
#include "err.h"

struct bitstream_reader {
  struct ogg_reader * ogg;

  unsigned int buffer;
};

struct bitstream_reader * create_bitstream_reader(struct ogg_reader * or) {
  struct bitstream_reader * br = malloc(sizeof(struct bitstream_reader));
  expect(br);

  br->ogg = or;
  br->buffer = 1;

  return br;
}

void destroy_bitstream_reader(struct bitstream_reader * br) {
  expect( get_packet_byte(br->ogg) == EOF );
  free(br);
}

uint32_t read_bits(struct bitstream_reader * br, int n) {
  expect(n <= 32);

  unsigned int t = 0;
  for (unsigned int i = 0; i < n; i++)
  {
    if (br->buffer == 1)
    {
      int c = get_packet_byte(br->ogg);
      expect(c != EOF);

      br->buffer = c | 0x100;
    }

    t |= (uint32_t)(br->buffer & 1) << i;
    br->buffer >>= 1;
  }

  return t;
}
