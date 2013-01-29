objs = abe-rsdec-test.o grs_encode_check.o

abe_grs: $(objs)
    gcc -o abe_grs $(objs) -L ~/.local/lib -Wl, -rpath ~/.local/lib -l pbc
abe-rsdec-test.o: abe-rsdec-test.c grs.h
    gcc -c abe-rsdec-test.c -L ~/.local/include/pbc/
grs_encode_check.o: grs_encode_check.c grs.h
    gcc -c grs_encode_check.c 

.PHONY: clean
clean:
    -rm -f *.o
