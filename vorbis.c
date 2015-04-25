#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

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
          // I'm not trying to do a lookup with this later
          expect(v->codebooks[c->masterbook].lookup_type == 0);
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
    // I'm not trying to do a lookup with this later
    expect( v->codebooks[r->classbook].lookup_type == 0 );

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

    m->mux = malloc(sizeof(uint8_t) * v->audio_channels);
    expect(m->mux);
      
    if (m->submap_count > 1)
    {
      for (unsigned int j = 0; j < v->audio_channels; j += 1)
      {
        m->mux[j] = read_bits(br, 4);
        expect(m->mux[j] < m->submap_count);
      }
    }
    else
    {
      // NOTE: vorbis_mapping_mux is unspecified in the spec with only one
      // submapping, but it is clear that it needs to be all 0
      for (unsigned int j = 0; j < v->audio_channels; j += 1)
      {
        m->mux[j] = 0;
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

struct code_reader * create_code_readers(struct vorbis_header * v)
{
  struct code_reader * crs = malloc(sizeof(struct code_reader) * v->codebook_count);
  expect(crs);
  
  uint32_t masks[33];
  for (unsigned int i = 1; i <= 32; i += 1)
  {
    masks[i] = (UINT32_C(1) << (32 - i)) - 1;
    masks[i] = ~masks[i];
  }

  for (unsigned int i = 0; i < v->codebook_count; i += 1)
  {
    struct codebook * cb = &v->codebooks[i];
    struct code_reader * cr = &crs[i];

    cr->count = cb->entry_count;
    cr->codewords = malloc(sizeof(uint32_t) * cr->count);
    expect(cr->codewords);
    cr->lengths = malloc(sizeof(uint8_t) * cr->count);
    expect(cr->lengths);
    cr->used = malloc(sizeof(bool) * cr->count);
    expect(cr->used);

    #if 0
    printf("codebook %d\n", i);
    printf("type: %s\n", cb->ordered ? "ordered" : (cb->sparse ? "sparse" : "nonsparse"));
    #endif

    // determine codewords
    bool all1 = false;
    bool exhausted[33] = {0};
    uint32_t next_codewords[33] = {0};

    for (unsigned int j = 0; j < cb->entry_count; j += 1)
    {
      unsigned int length = cb->entry_lengths[j];
      expect( length <= 32 );

      if (length > 0)
      {
        cr->used[j] = true;
        cr->lengths[j] = length;

        expect(!exhausted[length]);
        uint32_t codeword = next_codewords[length];
        bool collided;
        uint32_t next_step = 0;
        do
        {
          // check for collisions
          collided = false;
          for (unsigned int k = 0; k < j; k += 1)
          {
            if (cr->used[k])
            {
              uint32_t mask = masks[length];
              if (cr ->lengths[k] < length)
              {
                mask = masks[cr->lengths[k]];
              }

              if ((cr->codewords[k] & mask) == (codeword & mask))
              {
                collided = true;
                // next attempt must at least differ within the mask
                next_step = (codeword + (~mask) + 1) & mask;
                expect(next_step > codeword);
                codeword = next_step;
                break;
              }
            }
          }
        }
        while (collided);

        next_codewords[length] = codeword + (UINT32_C(1) << (32 - length));
        if (next_codewords[length] < codeword)
        {
          exhausted[length] = true;
        }

        cr->codewords[j] = codeword;

        if (cr->codewords[j]  == (-1 & masks[length]))
        {
          expect(exhausted[length]);
          all1 = true;
        }

        #if 0
        printf("entry %3u: length %3d codeword ", j, (int)length);
        for (int k = 31; k >= 32 - length; k -= 1)
        {
          printf("%c", ( codeword & (UINT32_C(1) << k) ) ? '1' : '0');
        }
        printf("\n");
        #endif
      }
      else
      {
        cr->used[j] = false;
      }
    } // end entry loop

    // make sure we hit the end of the tree
    expect(all1);

    // build tree
    cr->decoder = NULL;
    /*
    cr->decoder = malloc(sizeof(struct code_tree_node));
    expect(cr->decoder);
    cr->decoder.is_leaf = false;
    cr->decoder.left = cr->decoder.right = NULL;
    */

    for (unsigned int j = 0; j < cr->count; j += 1)
    {
      if (!cr->used[j])
      {
        continue;
      }

      struct code_tree_node ** pcur = &cr->decoder;
      struct code_tree_node * cur = cr->decoder;

      for (int k = 31; k >= 32 - cr->lengths[j]; k -= 1)
      {
        if (!cur)
        {
          cur = malloc(sizeof(struct code_tree_node));
          expect(cur);
          *pcur = cur;
          cur->is_leaf = false;
          cur->children[0] = cur->children[1] = NULL;
        }

        unsigned int next = (cr->codewords[j] >> k) & 1;
        pcur = &cur->children[next];
        cur = cur->children[next];
      }

      expect(!cur);
      cur = malloc(sizeof(struct code_tree_node));
      expect(cur);
      *pcur = cur;
      cur->is_leaf = true;
      cur->leaf = j;
    }
  }

  return crs;
}

uint32_t read_code(struct bitstream_reader * br, struct code_reader * cr)
{
  struct code_tree_node * cur = cr->decoder;

  expect(!cur->is_leaf);
  do
  {
    cur = cur->children[read_bits(br, 1)];
    expect(cur);
  }
  while (!cur->is_leaf);

  return cur->leaf;
}


void read_vorbis_audio_packet(
  struct bitstream_reader * br,
  struct vorbis_header * v,
  struct code_reader * crs
) {

  printf("decoding packet\n");
  // audio packet type
  expect( read_bits(br, 1) == 0 );

  unsigned int mode_number = read_bits(br, ilog(v->mode_count - 1));
  expect( mode_number < v->mode_count );

  struct mode * mode = &v->modes[mode_number];

  bool blockflag = mode->blockflag;

  unsigned int blocksize;

  if (blockflag)
  {
    // long window
    bool previous_window_flag = read_bits(br, 1);
    bool next_window_flag = read_bits(br, 1);

    blocksize = v->blocksize_1;
  }
  else
  {
    blocksize = v->blocksize_0;
  }

  struct mapping * mapping = &v->mappings[mode->mapping];

  bool no_residue[v->audio_channels];
  for (unsigned int i = 0; i < v->audio_channels; i += 1)
  {
    no_residue[i] = false;
  }

  // decode floors
  for (unsigned int channel = 0; channel < v->audio_channels; channel += 1)
  {
    unsigned int submap_number = mapping->mux[channel];
    unsigned int floor_number = mapping->submaps[submap_number].floor;

    struct floor * floor = &v->floors[floor_number];
    if (floor->type == 1)
    {
      struct floor1 * f = &floor->type1;

      bool nonzero = read_bits(br, 1);
      if (!nonzero)
      {
        no_residue[channel] = true;
        continue;
      }

      static const unsigned int ranges[4] = {256, 128, 86, 64};
      unsigned int range = ranges[f->multiplier-1];

      // floor1_Y[0]
      read_bits(br, ilog(range-1));
      // floor1_Y[1]
      read_bits(br, ilog(range-1));

      for (unsigned int i = 0; i < f->partitions; i += 1)
      {
        uint8_t class = f->partition_class_list[i];
        uint8_t cdim = f->classes[class].dimensions;
        uint8_t cbits = f->classes[class].subclasses;
        uint32_t csub = (1 << cbits) - 1;

        unsigned int cval = 0;
        if (cbits > 0)
        {
          cval = read_code(br, &crs[f->classes[class].masterbook]);
          // I don't have lookup supported
          expect(v->codebooks[f->classes[class].masterbook].lookup_type == 0);
        }

        for (unsigned int j = 0; j < cdim; j += 1)
        {
          uint8_t book = f->classes[class].subclass_books[cval & csub];
          cval = cval >> cbits;

          if (book != 0xFF)
          {
            // floor1_Y[j + offset] 
            read_code(br, &crs[book]);
          }
          else
          {
            // floor1_Y[j + offset] = 0
          }
        }
      }
    }
    else
    {
      // only type 1 supported at the moment
      expect(false);
    }
  } // end floors

  // residues (per submap)
  for (unsigned int i = 0; i < mapping->submap_count; i += 1)
  {
    bool do_not_decode_flag[v->audio_channels];
    for (unsigned int j = 0; j < v->audio_channels; j += 1)
    {
      do_not_decode_flag[j] = false;
    }

    unsigned int ch = 0;
    for (unsigned int j = 0; j < v->audio_channels; j += 1)
    {

      if (mapping->mux[j] == i)
      {
        //printf("submap %u has channel %u\n", i, j);
        if (no_residue[i])
        {
          do_not_decode_flag[ch] = true;
        }

        ch += 1;
      }
    }

    unsigned int residue_number = mapping->submaps[i].residue;
    struct residue * r = &v->residues[residue_number];
    unsigned int residue_type = r->type;

    unsigned int actual_size = blocksize/2;
    if (residue_type == 2)
    {
      actual_size *= ch;
    }

    unsigned int limit_residue_begin = r->begin;
    if (actual_size > limit_residue_begin)
    {
      limit_residue_begin = actual_size;
    }
    unsigned int limit_residue_end = r->end;
    if (actual_size > limit_residue_end)
    {
      limit_residue_end = actual_size;
    }

    unsigned int classwords_per_codeword = v->codebooks[r->classbook].dimensions;
    unsigned int n_to_read = limit_residue_end - limit_residue_begin;
    unsigned int partitions_to_read = n_to_read / r->partition_size;

    if (n_to_read == 0)
    {
      continue;
    }

    for (unsigned int pass = 0; pass < 8; pass += 1)
    {
      unsigned int classifications[ch][classwords_per_codeword+partitions_to_read];
      unsigned int partition_count = 0;

      while (partition_count < partitions_to_read)
      {
        if (pass == 0)
        {
          for (unsigned int j = 0; j < ch; j += 1)
          {

            if (do_not_decode_flag[j])
            {
              continue;
            }

            // I don't have lookup supported
            expect(v->codebooks[r->classbook].lookup_type == 0);
            unsigned int temp = read_code(br, &crs[r->classbook]);

            for (int k = classwords_per_codeword-1; k >= 0; k -= 1)
            {
              classifications[j][k+partition_count] = temp % r->classifications;
              temp /= r->classifications;
            }
          }
        } // end if (pass == 0)

        for (unsigned int k = 0; k < classwords_per_codeword && partition_count < partitions_to_read; k += 1)
        {
          for (unsigned int j = 0; j < ch; j += 1)
          {
            if (do_not_decode_flag[j])
            {
              continue;
            }

            unsigned int vqclass = classifications[j][partition_count];
            unsigned int vqbook = r->books[vqclass*8 + pass];
            if (r->cascade[vqclass] & (1 << pass))
            {
              if (residue_type == 1)
              {
                unsigned int ii = 0;
                do
                {
                  for (unsigned int jj = 0; jj < v->codebooks[vqbook].dimensions; jj += 1)
                  {
                    read_code(br, &crs[vqbook]);
                    ii += 1;
                  }
                }
                while (ii < r->partition_size);
              }
              else
              {
                expect(false);
              }
            }
          }

          partition_count += 1;
        }
      } // end while (partition_count < partitions_to_read)
    }  // end of pass loop

  } // end of residue (submap loop)

  //abort();
}
