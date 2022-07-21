CC = gcc
CFLAGS = -std=c99 -O3 -Wall
LIBS = -lm -lfftw3
INCL = -Isrc -Iio -Ilib -Imath

# Settings for FFTW
FFTW_DIR = 
ifneq ($(FFTW_DIR),)
  LIBS += -L$(FFTW_DIR)/lib
  INCL += -I$(FFTW_DIR)/include
endif

# Setting for single precision density fields and FFT
#LIBS += -DSINGLE_PREC

# Settings for OpenMP (comment the following line to disable OpenMP)
LIBS += -DOMP -fopenmp -lfftw3_omp

# Settings for CFITSIO (not implemented yet)

SRCS = $(wildcard src/*.c lib/*.c io/*.c math/*.c)
EXEC = POWSPEC

all:
	$(CC) $(CFLAGS) -o $(EXEC) $(SRCS) $(LIBS) $(INCL)

clean:
	rm $(EXEC)
