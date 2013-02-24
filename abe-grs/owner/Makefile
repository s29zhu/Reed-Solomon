run: owner_keys encryption
	./encryption
owner_keys: generate_owner_keys.c
	gcc -o owner_keys generate_owner_keys.c  -I /home/s29zhu/.local/include/pbc/ -L /home/s29zhu/.local/lib/  -Wl,-rpath /home/s29zhu/.local/lib/ -l pbc
encryption: encryption.c
	gcc -o encryption encryption.c message_handle.c -I /home/s29zhu/.local/include/pbc/ -L /home/s29zhu/.local/lib/  -Wl,-rpath /home/s29zhu/.local/lib/ -l pbc
.PHONY: clean
clean:
	mv *.o
