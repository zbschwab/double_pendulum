CC := clang
CFLAGS := -g -Wall -fsanitize=address $(shell pkg-config sdl2 SDL2_image --cflags)
LDFLAGS := $(shell pkg-config sdl2 SDL2_image --libs)

all : pendulum

clean :
	rm -rf pendulum pendulum.dSYM

pendulum : pendulum.c phys_math.c phys_math.h
	$(CC) $(CFLAGS) -o pendulum phys_math.c pendulum.c -lncurses -lm $(LDFLAGS)

phys_math.o : phys_math.c phys_math.h
	$(CC) $(CFLAGS) -c phys_math.c

.PHONY: all clean