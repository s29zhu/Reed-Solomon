authorizer_keys: generate_authorizer_keys.c
	gcc -o authorizer_keys generate_authorizer_keys.c -I /home/s29zhu/.local/include/pbc/ -L /home/s29zhu/.local/lib/ -Wl,-rpath /home/s29zhu/.local/lib -l pbc
.PHONY: clean
clean:
	mv *.o
