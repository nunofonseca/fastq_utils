

all:
	make -C src

install:
	make -C src install

tests:
	./run_tests.sh
