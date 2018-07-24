all:
	gcc -O3 driveGillespie.c -o gillespie.out -lm

run:
	gcc -O3 driveGillespie.c -o gillespie.out -lm

	chmod a+x batchGillespie.sh

	./batchGillespie.sh 7 v

test:
	gcc -O3 driveGillespie.c -o gillespie.out -lm

	./gillespie.out Kinase_Constant Phosphatase_Constant 6 10000000 0.2 OutputTest 0 100000


copyrepo:

	git -C ~/Documents/PolymerGit log -1 --pretty=format:%H > "CommitUsedHash.txt"

	cp -r ~/Documents/PolymerGit/src/Analytics/Gillespie/Reversible_Constant/ .



