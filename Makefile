m.PHONY : clean

all : libparse_time.so

libparse_time.so : parse_time.o
	gcc -shared -Wl -o libparse_time.so parse_time.o

parse_time.o : parse_time.c
	gcc -c -fPIC parse_time.c -o parse_time.o

clean :
	-rm -vf libparse_time.so parse_time.o parse_time.pyc