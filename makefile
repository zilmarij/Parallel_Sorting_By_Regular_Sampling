sort: sort.cpp sort.h
	mpic++ -c -o hyb.o sort.cpp -l.
	ar rcs libpsort.a hyb.o	

   #hyb.o: Hybrid.cpp Hybrid.h
	

	