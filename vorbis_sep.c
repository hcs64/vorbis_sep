#include <stdlib.h>
#include <stdio.h>

#include "ogg.h"
#include "vorbis.h"
#include "bs.h"
#include "err.h"

int main(int argc, char ** argv) {
  expect(argc == 2);
  FILE * infile = fopen(argv[1], "rb");
  expect(infile);

  struct ogg_reader * or = create_ogg_reader(infile);

  struct vorbis_header v;

  {
    struct bitstream_reader * br = create_bitstream_reader(or);
    read_identification_header(br, &v);
    destroy_bitstream_reader(br);
  }
  {
    struct bitstream_reader * br = create_bitstream_reader(or);
    read_comment_header(br, &v);
    destroy_bitstream_reader(br);
  }
  {
    struct bitstream_reader * br = create_bitstream_reader(or);
    read_setup_header(br, &v);
    destroy_bitstream_reader(br);
  }

  printf("vorbis header loaded ok\n");

  // check that we don't use coupling
  bool coupling = false;
  for (unsigned int i = 0; i < v.mapping_count; i += 1)
  {
    if (v.mappings[i].coupling_steps > 0)
    {
      printf("mapping %d: %d coupling steps\n", i, v.mappings[i].coupling_steps);
      coupling = true;
    }
  }

  if (coupling)
  {
      printf("ERROR: coupling not supported\n");
      exit(EXIT_FAILURE);
  }

  struct code_reader * crs = create_code_readers(&v);

  while (!ogg_eos(or))
  {
    struct bitstream_reader * br = create_bitstream_reader(or);

    read_vorbis_audio_packet(br, &v, crs);

    destroy_bitstream_reader(br);
  }

  destroy_ogg_reader(or);

  return 0;
}

