VERSION=0.9.0
# Requires zlib libraries


all:
	make -C src

install:
	mdkir -p bin && make -C src install

release:
	sed -i "s/#define VERSION.*/#define VERSION \"$(VERSION)\"/" src/*.c

tests:
	./run_tests.sh
