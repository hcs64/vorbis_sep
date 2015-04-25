#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "bs.h"
#include "err.h"

#include "vorbis.h"

unsigned int ilog(int x) {
  unsigned int r = 0;

  while (x > 0)
  {
    r += 1;
    x >>= 1;
  }

  return r;
}

unsigned int book_maptype1_quantvals(unsigned int entries, unsigned int dimensions){
  /* get us a starting hint, we'll polish it below */
  int bits=ilog(entries);
  int vals=entries>>((bits-1)*(dimensions-1)/dimensions);
  
  while(1){
    unsigned long acc=1;
    unsigned long acc1=1;
    unsigned int i;
    for(i=0;i<dimensions;i++){
      acc*=vals;
      acc1*=vals+1;
    }
    if(acc<=entries && acc1>entries){
      return(vals); 
    }else{
      if(acc>entries){
        vals--;
      }else{
        vals++;
      } 
    }
  }
}

void expect_vorbis(struct bitstream_reader * br) {
  expect( read_bits(br, 8) == 'v' );
  expect( read_bits(br, 8) == 'o' );
  expect( read_bits(br, 8) == 'r' );
  expect( read_bits(br, 8) == 'b' );
  expect( read_bits(br, 8) == 'i' );
  expect( read_bits(br, 8) == 's' );
}

void read_identification_header(struct bitstream_reader * br, struct vorbis_header * v) {
  // identification header type
  expect( read_bits(br, 8) == 1 );
  expect_vorbis(br);

  // vorbis_version
  expect( read_bits(br, 32) == 0 );
  v->audio_channels = read_bits(br, 8);
  v->audio_sample_rate = read_bits(br, 32);
  v->bitrate_maximum = (int32_t)read_bits(br, 32);
  v->bitrate_nominal = (int32_t)read_bits(br, 32);
  v->bitrate_minimum = (int32_t)read_bits(br, 32);
  v->blocksize_0 = UINT32_C(1) << read_bits(br, 4);
  v->blocksize_1 = UINT32_C(1) << read_bits(br, 4);
  expect( v->blocksize_0 <= v->blocksize_1 );
  // framing flag
  expect( read_bits(br, 1) );
}

void read_comment_header(struct bitstream_reader * br, struct vorbis_header * v) {
  // comment header type
  expect( read_bits(br, 8) == 3 );
  expect_vorbis(br);

  // 
  uint32_t vendor_length = read_bits(br, 32);
  v->vendor = malloc(vendor_length);
  expect(v->vendor);

  for (unsigned int i = 0; i < vendor_length; i+= 1)
  {
    v->vendor[i] = read_bits(br, 8);
  }

  v->user_comment_list_length = read_bits(br, 32);
  v->user_comments = malloc(
    sizeof(struct vorbis_comment) * v->user_comment_list_length);
  expect(v->user_comments);

  for (unsigned int i = 0; i < v->user_comment_list_length; i += 1)
  {
    struct vorbis_comment * vc = &v->user_comments[i];

    vc->length = read_bits(br, 32);
    vc->comment = malloc(vc->length);
    expect( vc->comment );

    for (unsigned int j = 0; j < vc->length; j += 1)
    {
      vc->comment[j] = read_bits(br, 8);
    }

  }

  // framing bit
  expect( read_bits(br, 1) );
}

void read_setup_header(struct bitstream_reader * br, struct vorbis_header * v) {
  // setup header type
  expect( read_bits(br, 8) == 5 );
  expect_vorbis(br);

  v->codebook_count = read_bits(br, 8) + 1;
  v->codebooks = malloc(sizeof(struct codebook) * v->codebook_count);
  expect( v->codebooks );

  for (unsigned int cbi = 0; cbi < v->codebook_count; cbi += 1)
  {
    // VCB
    expect( read_bits(br, 24) == 0x564342 );

    struct codebook * cb = &v->codebooks[cbi];

    cb->dimensions = read_bits(br, 16);
    cb->entry_count = read_bits(br, 24);

    cb->entry_lengths = malloc(sizeof(uint8_t) * cb->entry_count);
    expect(cb->entry_lengths);


    // lengths
    cb->ordered = read_bits(br, 1);

    if (!cb->ordered)
    {
      
      cb->sparse = read_bits(br, 1);

      for (unsigned int cbe = 0; cbe < cb->entry_count; cbe += 1)
      {
        if (!cb->sparse || read_bits(br, 1))
        {
          // present
          cb->entry_lengths[cbe] = read_bits(br, 5) + 1;
        }
        else
        {
          // absent
          cb->entry_lengths[cbe] = 0;
        }
      }
    }
    else // ordered
    {
      unsigned int current_entry = 0;
      unsigned int current_length = read_bits(br, 5) + 1;

      do
      {
        unsigned int number = read_bits(br, ilog(cb->entry_count - current_entry));
        for (; number > 0; number -= 1, current_entry += 1)
        {
          expect( current_entry < cb->entry_count );
          cb->entry_lengths[current_entry] = current_length;
        }

        current_length += 1;
      }
      while (current_entry < cb->entry_count);
    }

    // lookup
    cb->lookup_type = read_bits(br, 4);

    if (cb->lookup_type == 0)
    {
      // no lookup
    }
    else if (cb->lookup_type == 1 || cb->lookup_type == 2)
    {
      cb->minimum_value = read_bits(br, 32);
      cb->delta_value = read_bits(br, 32);
      cb->value_bits = read_bits(br, 4) + 1;
      cb->sequence_p = read_bits(br, 1);

      if (cb->lookup_type == 1)
      {
        // implicitly populated value mapping (lattice VQ)

        cb->lookup_count = book_maptype1_quantvals(cb->entry_count, cb->dimensions);
      }
      else
      {
        // explicitly populated value mapping (tessellated or 'foam' VQ)

        cb->lookup_count = cb->entry_count * cb->dimensions;
      }

      cb->lookups = malloc(sizeof(uint16_t) * cb->lookup_count);
      expect(cb->lookups);

      for (uint32_t i = 0; i < cb->lookup_count; i += 1)
      {
        cb->lookups[i] = read_bits(br, cb->value_bits);
      }
    }
    else // unknown lookup type
    {
      expect(false);
    }
  } // end codebooks

  // time domain transform placeholders
  expect( read_bits(br, 6) == 0 );
  expect( read_bits(br, 16) == 0 );

  // floors
  v->floor_count = read_bits(br, 6) + 1;
  v->floors = malloc(sizeof(struct floor) * v->floor_count);
  expect(v->floors);

  for (unsigned int fi = 0; fi < v->floor_count; fi += 1)
  {
    struct floor * f = &v->floors[fi];

    f->type = read_bits(br, 16);

    if (f->type == 0)
    {
      struct floor0 * f0 = &f->type0;
      f0->order = read_bits(br, 8);
      f0->rate = read_bits(br, 16);
      f0->bark_map_size = read_bits(br, 16);
      f0->amplitude_bits = read_bits(br, 6);
      f0->amplitude_offset = read_bits(br, 8);
      f0->number_of_books = read_bits(br, 4) + 1;
      f0->book_list = malloc(sizeof(uint8_t) * f0->number_of_books);
      expect(f0->book_list);

      for (unsigned int i = 0; i < f0->number_of_books; i += 1)
      {
        f0->book_list[i] = read_bits(br, 8);
        expect( f0->book_list[i] < v->codebook_count );
      }
    }
    else if (f->type == 1)
    {
      struct floor1 * f1 = &f->type1;
      f1->partitions = read_bits(br, 5);
      f1->partition_class_list = malloc(sizeof(uint8_t) * f1->partitions);
      expect(f1->partition_class_list);

      f1->maximum_class = -1;

      for (unsigned int i = 0; i < f1->partitions; i += 1)
      {
        f1->partition_class_list[i] = read_bits(br, 4);
        if (f1->partition_class_list[i] > f1->maximum_class)
        {
          f1->maximum_class = f1->partition_class_list[i];
        }
      }

      if (f1->maximum_class > -1)
      {
        f1->classes = malloc(sizeof(struct floor1_class) * (f1->maximum_class+1));
        expect(f1->classes);
      }
      else
      {
        f1->classes = NULL;
      }

      for (unsigned int i = 0; i <= f1->maximum_class; i += 1)
      {
        struct floor1_class * c = &f1->classes[i];
        c->dimensions = read_bits(br, 3) + 1;
        c->subclasses = read_bits(br, 2);

        if (c->subclasses)
        {
          c->masterbook = read_bits(br, 8);
          expect( c->masterbook < v->codebook_count );
        }

        c->subclass_books_count = 1 << c->subclasses;
        c->subclass_books = malloc(sizeof(uint8_t) * c->subclass_books_count);
        expect(c->subclass_books);
        for (unsigned int j = 0; j < c->subclass_books_count; j += 1)
        {
          c->subclass_books[j] = read_bits(br, 8) - 1;
          expect( c->subclass_books[j] == 0xff ||
                  c->subclass_books[j] < v->codebook_count );
        }
      }

      f1->multiplier = read_bits(br, 2) + 1;
      f1->rangebits = read_bits(br, 4);

      f1->x_list[0] = 0;
      f1->x_list[1] = 1 << f1->rangebits;
      f1->values = 2;

      for (unsigned int i = 0; i < f1->partitions; i += 1)
      {
        uint8_t class_number = f1->partition_class_list[i];
        uint8_t class_dimensions = f1->classes[class_number].dimensions;

        for (unsigned int j = 0; j < class_dimensions; j += 1)
        {
          expect( f1->values < 65 );
          f1->x_list[f1->values] = read_bits(br, f1->rangebits);
          f1->values+= 1;
        }
      }
    }
    else  // not floor 0 or 1
    {
      expect(false);
    }
  } // end floors

  // residues
  v->residue_count = read_bits(br, 6) + 1;
  v->residues = malloc(sizeof(struct residue) * v->residue_count);
  expect(v->residues);

  for (unsigned int ri = 0; ri < v->residue_count; ri += 1)
  {
    struct residue * r = &v->residues[ri];
    
    r->type = read_bits(br, 16);
    expect (r->type == 0 || r->type == 1 || r->type == 2);

    r->begin = read_bits(br, 24);
    r->end = read_bits(br, 24);
    r->partition_size = read_bits(br, 24) + 1;
    r->classifications = read_bits(br, 6) + 1;
    r->classbook = read_bits(br, 8);
    expect( r->classbook < v->codebook_count );

    r->cascade = malloc(sizeof(uint8_t) * r->classifications);
    expect(r->cascade);

    for (unsigned int i = 0; i < r->classifications; i += 1)
    {
      uint8_t high_bits = 0;
      uint8_t low_bits = read_bits(br, 3);

      if (read_bits(br, 1))
      {
        high_bits = read_bits(br, 5);
      }

      r->cascade[i] = high_bits * 8 + low_bits;
    }

    r->books = malloc(sizeof(uint8_t) * r->classifications * 8);
    expect(r->books);

    for (unsigned int i = 0; i < r->classifications; i += 1)
    {
      for (unsigned int j = 0; j < 8; j += 1)
      {
        if (r->cascade[i] & (1 << j))
        {
          r->books[i*8+j] = read_bits(br, 8);
          expect( r->books[i*8+j] < v->codebook_count );
        }
      }
    }
  } // end residues

  // mappings
  v->mapping_count = read_bits(br, 6) + 1;
  v->mappings = malloc(sizeof(struct mapping) * v->mapping_count);
  expect(v->mappings);

  for (unsigned int i = 0; i < v->mapping_count; i += 1)
  {
    struct mapping * m = &v->mappings[i];

    m->type = read_bits(br, 16);
    expect(m->type == 0);

    if (read_bits(br, 1))
    {
      m->submap_count = read_bits(br, 4) + 1;
    }
    else
    {
      m->submap_count = 1;
    }

    if (read_bits(br, 1))
    {
      // square polar channel mapping
      m->coupling_steps = read_bits(br, 8) + 1;

      m->couplings = malloc(sizeof(struct coupling) * m->coupling_steps);
      expect(m->couplings);
      for (unsigned int j = 0; j < m->coupling_steps; j += 1)
      {
        struct coupling * c = &m->couplings[j];

        c->magnitude = read_bits(br, ilog(v->audio_channels-1));
        c->angle = read_bits(br, ilog(v->audio_channels-1));

        expect(c->magnitude < v->audio_channels);
        expect(c->angle < v->audio_channels);
        expect(c->magnitude != c->angle);
      }
    }
    else
    {
      m->coupling_steps = 0;
    }

    expect(read_bits(br, 2) == 0);

    if (m->submap_count > 1)
    {
      m->mux = malloc(sizeof(uint8_t) * v->audio_channels);
      expect(m->mux);
      
      for (unsigned int j = 0; j < v->audio_channels; j += 1)
      {
        m->mux[j] = read_bits(br, 4);
        expect(m->mux[j] < m->submap_count);
      }
    }

    m->submaps = malloc(sizeof(struct submap) * m->submap_count);
    expect(m->submaps);

    for (unsigned int j = 0; j < m->submap_count; j += 1)
    {
      struct submap * sm = &m->submaps[j];
      expect( read_bits(br, 8) == 0 );

      sm->floor = read_bits(br, 8);
      expect(sm->floor < v->floor_count);
      sm->residue = read_bits(br, 8);
      expect(sm->residue < v->residue_count);
    }
  } // end mappings

  // modes
  v->mode_count = read_bits(br, 6) + 1;
  v->modes = malloc(sizeof(struct mode) * v->mode_count);
  expect(v->modes);

  for (unsigned int i = 0; i < v->mode_count; i += 1)
  {
    struct mode * m = &v->modes[i];
    m->blockflag = read_bits(br, 1);
    m->windowtype = read_bits(br, 16);
    m->transformtype = read_bits(br, 16);
    m->mapping = read_bits(br, 8);

    expect(m->windowtype == 0);
    expect(m->transformtype == 0);
    expect(m->mapping < v->mapping_count);
  } // end mappings

  // framing flag
  expect( read_bits(br, 1) == 1 );
}

