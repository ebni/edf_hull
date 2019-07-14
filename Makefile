qhull:
	git clone https://github.com/qhull/qhull.git
	mv qhull qhull-src
	patch qhull-src/Makefile ./patch_M32
	$(MAKE) -C qhull-src
	ln -s qhull-src/bin/qhull

edf_hull:	edf_hull.c
	gcc edf_hull.c -o edf_hull -lm

run:	 qhull edf_hull
	./edf_hull < task_set.txt

clean:	
	rm -rf qhull* edf_hull.o edf_hull *~ data_in.txt data_out.txt
