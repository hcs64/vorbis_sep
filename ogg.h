#ifndef INCLUDED_OGG_H
#define INCLUDED_OGG_H

#include <stdbool.h>
#include <stdio.h>

struct ogg_reader;
struct ogg_writer;

struct ogg_reader * create_ogg_reader(FILE * file);
int get_packet_byte(struct ogg_reader *);
bool ogg_eos(struct ogg_reader *);
void destroy_ogg_reader(struct ogg_reader *);


#endif // INCLUDED_OGG_H
