CFLAGS=-std=c99 -Wall

vorbis_sep: vorbis_sep.o vorbis.o ogg.o bs.o

vorbis_sep.o: vorbis_sep.c vorbis.h ogg.h err.h bs.h

vorbis.o: vorbis.c bs.h err.h

ogg.o: ogg.c ogg.h err.h

bs.o: bs.c bs.h ogg.h err.h
