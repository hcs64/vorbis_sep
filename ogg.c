#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "err.h"

#include "ogg.h"

struct ogg_page_header {
  bool continued_flag;
  bool first_flag;
  bool last_flag;

  uint64_t granulepos;

  uint32_t stream_serial;

  uint32_t seqno;

  uint32_t checksum;

  unsigned int segments;
  uint8_t segment_table[256];
};

struct ogg_reader {
  FILE * file;

  bool page_ok;
  struct ogg_page_header page;
  unsigned int segment;
  unsigned int segment_bytes_left;

  uint32_t serial;
  bool bos;
  bool eos;
};

struct ogg_writer {
  FILE * file;
};

uint64_t get64(FILE * file) {
  uint64_t t = 0;
  for (unsigned int i = 0; i < 8; i += 1)
  {
    int c = getc(file);
    expect(c != EOF);

    t *= 0x100;
    t += c;
  }

  return t;
}

uint32_t get32(FILE * file) {
  uint32_t t = 0;
  for (unsigned int i = 0; i < 4; i += 1)
  {
    int c = getc(file);
    expect(c != EOF);

    t *= 0x100;
    t += c;
  }

  return t;
}

void get_ogg_page(FILE * file, struct ogg_page_header * h) {
  if (
    getc(file) != 'O' ||
    getc(file) != 'g' ||
    getc(file) != 'g' ||
    getc(file) != 'S' ||
    getc(file) != 0
  )
  {
    expect(false);
  }

  int c;

  c = getc(file);
  expect(c != EOF);

  h->continued_flag = c & 0x01;
  h->first_flag = c & 0x02;
  h->last_flag = c & 0x04;

  h->granulepos = get64(file);
  h->stream_serial = get32(file);
  h->seqno = get32(file);
  h->checksum = get32(file);

  c = getc(file);
  expect(c != EOF);
  h->segments = c;

  for (unsigned int i = 0; i < h->segments; i += 1) {
    c = getc(file);
    expect(c != EOF);

    h->segment_table[i] = c;
  }
}

struct ogg_reader * create_ogg_reader(FILE * file) {
  struct ogg_reader * const or = malloc(sizeof(struct ogg_reader));
  expect(or);

  or->file = file;
  or->page_ok = false;
  or->bos = true;
  or->eos = false;

  return or;
}

int get_packet_byte(struct ogg_reader * or) {

check_page:
  if (!or->page_ok)
  {
    expect(!or->eos);

    // no page is loaded, get the next one in the file
    get_ogg_page(or->file, &or->page);

    if (or->bos)
    {
      expect(or->page.first_flag);
      or->bos = false;
      or->serial = or->page.stream_serial;
    }
    else
    {
      expect (or->serial == or->page.stream_serial);
    }

    if (or->page.last_flag)
    {
      or->eos = true;
    }

    // start handling segments on this page
    or->segment = 0;

    // start handling bytes in this segment
    or->segment_bytes_left = or->page.segment_table[0];

    or->page_ok = true;
  }

  // here we have a valid page loaded

check_segment:
  if (or->segment == or->page.segments)
  {
    // we are out of segments on this page
    or->page_ok = false;
    goto check_page;
  }

  // here we have a segment available

  if (or->segment_bytes_left == 0)
  {
    or->segment += 1;
    or->segment_bytes_left = or->page.segment_table[or->segment];
    if (or->page.segment_table[or->segment-1] < 255)
    {
      return EOF;
    }

    goto check_segment;
  }

  // here we have a byte available in the segment

  int c = getc(or->file);
  expect(c != EOF);

  or->segment_bytes_left -= 1;

  return c;
}

bool ogg_eos(struct ogg_reader * or) {
  return (or->eos && or->segment == or->page.segments);
}

void destroy_ogg_reader(struct ogg_reader * or) {
  expect( or->eos );
  free(or);
}
