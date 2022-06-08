run:	 qhull edf_hull_main
	./edf_hull_main < task_set.txt

qhull:
	git clone https://github.com/qhull/qhull.git qhull-src
	patch qhull-src/Makefile ./patches
	$(MAKE) -C qhull-src
	ln -s qhull-src/bin/qhull

ts_lib.o:	ts_lib.c ts_lib.h Makefile
	gcc -c -g -O0 -std=c89 -Wpedantic ts_lib.c -o ts_lib.o

edf_hull.o: edf_hull.c edf_hull.h qhull Makefile
	gcc -c -g -O0 -std=c89 -Wpedantic edf_hull.c -o edf_hull.o

edf_hull_main: edf_hull.o ts_lib.o main.c qh_init_struct.c Makefile
	gcc -g -O0 main.c edf_hull.o ts_lib.o qhull-src/lib/libqhullstatic_r.a qh_init_struct.c -lrt -lm -o edf_hull_main

clean:	
	rm -rf qhull* *.o edf_hull *~ data_in.txt data_out.txt
