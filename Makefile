VERSION=0.9.2
# Requires zlib libraries


all:
	make -C src

install:
	mkdir -p bin && make -C src install && cp sh/*.sh bin

release:
	sed -i "s/#define VERSION.*/#define VERSION \"$(VERSION)\"/" src/*.c
	sed -i "s/#define VERSION.*/#define VERSION \"$(VERSION)\"/" src/*.h

tests: FORCE
	./run_tests.sh

FORCE: 
