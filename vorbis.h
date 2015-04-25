#ifndef VORBIS_H_INCLUDED
#define VORBIS_H_INCLUDED

#include <stdbool.h>
#include <stdint.h>

#include "bs.h"

struct vorbis_header {
  // identification header
  uint8_t audio_channels;
  uint32_t audio_sample_rate;
  int32_t bitrate_maximum;
  int32_t bitrate_nominal;
  int32_t bitrate_minimum;
  uint32_t blocksize_0;
  uint32_t blocksize_1;

  // comment header
  char * vendor;
  uint32_t user_comment_list_length;
  struct vorbis_comment {
    uint32_t length;
    char * comment;
  } * user_comments;

  // setup header
  unsigned int codebook_count;
  struct codebook {
    uint16_t dimensions;
    uint32_t entry_count;
    bool ordered;
    bool sparse;

    uint8_t * entry_lengths;

    uint8_t lookup_type;
    uint32_t minimum_value;
    uint32_t delta_value;
    uint8_t value_bits;
    bool sequence_p;

    uint32_t lookup_count;
    uint16_t * lookups;

  } * codebooks;

  unsigned int floor_count;
  struct floor {
    unsigned int type;

    union {
      struct floor0 {
        uint8_t order;
        uint16_t rate;
        uint16_t bark_map_size;
        uint8_t amplitude_bits;
        uint8_t amplitude_offset;
        uint8_t number_of_books;
        uint8_t * book_list;
      } type0;
      struct floor1 {
        uint8_t partitions;
        uint8_t * partition_class_list;
        int maximum_class;

        struct floor1_class {
          uint8_t dimensions;
          uint8_t subclasses;
          uint8_t masterbook;
          uint8_t subclass_books_count;
          uint8_t * subclass_books;
        } * classes;

        uint8_t multiplier;
        uint8_t rangebits;
        
        uint8_t values;
        uint16_t x_list[65];
      } type1;
    };
  } * floors;

  unsigned int residue_count;
  struct residue {
    uint16_t type;

    uint32_t begin;
    uint32_t end;
    uint32_t partition_size;
    uint8_t classifications;
    uint8_t classbook;

    uint8_t * cascade;

    uint8_t * books;
  } * residues;

  unsigned int mapping_count;
  struct mapping {
    uint16_t type;
    uint8_t submap_count;

    unsigned int coupling_steps;
    struct coupling {
      uint8_t magnitude;
      uint8_t angle;
    } * couplings;

    uint8_t * mux;
    struct submap {
      uint8_t floor;
      uint8_t residue;
    } * submaps;

  } * mappings;

  uint8_t mode_count;
  struct mode {
    bool blockflag;
    uint16_t windowtype;
    uint16_t transformtype;
    uint8_t mapping;
  } * modes;
};

struct code_tree_node {
  union {
    struct code_tree_node * children[2];

    uint32_t leaf;
  };
  bool is_leaf;
};

struct code_reader {
  unsigned int count;
  uint32_t * codewords;
  uint8_t * lengths;
  bool * used;

  struct code_tree_node * decoder;
};

unsigned int book_maptype1_quantvals(unsigned int entries, unsigned int dimensions);
unsigned int ilog(int32_t x);

void read_identification_header(struct bitstream_reader *, struct vorbis_header *);
void read_comment_header(struct bitstream_reader *, struct vorbis_header *);
void read_setup_header(struct bitstream_reader *, struct vorbis_header *);

void read_vorbis_audio_packet(struct bitstream_reader *, struct vorbis_header *, struct code_reader * cr);

struct code_reader * create_code_readers(struct vorbis_header *);
uint32_t read_code(struct bitstream_reader *, struct code_reader *);

#endif // VORBIS_H_INCLUDED
