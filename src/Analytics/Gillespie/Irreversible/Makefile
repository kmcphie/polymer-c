all:
	gcc -O3 driveGillespie.c -o gillespie.out -lm

	chmod a+x batchGillespie.sh

run:
	gcc -O3 driveGillespie.c -o gillespie.out -lm

	chmod a+x batchGillespie_Control.sh

	./batchGillespie_Control.sh 7 v

copyrepo:

	git -C ~/Documents/PolymerGit log -1 --pretty=format:%H > "CommitUsedHash.txt"

	cp -r ~/Documents/PolymerGit/src/Analytics/Gillespie/Irreversible/ .


