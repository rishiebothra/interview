
output: main.o jsmn.o
#	gcc -g -std=c99 -O3 main.o jsmn.o -o a.out -lm
	gcc -std=c99 -O3 main.o jsmn.o -o a.out -lm

main.o: main.c
#	gcc -g -std=c99 main.c -c
	gcc -std=c99 main.c -c

jsmn.o: jsmn.c
#	gcc -g -std=c99 jsmn.c -c
	gcc -std=c99 jsmn.c -c

clean:
	rm *.o *.out
