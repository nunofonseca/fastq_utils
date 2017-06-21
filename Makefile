VERSION=0.9.4
# Requires zlib libraries


all:
	make -C src

install:
	mkdir -p bin && make -C src install && cp sh/*.sh bin

release:
	sed -i "s/#define VERSION.*/#define VERSION \"$(VERSION)\"/" src/*.c
	sed -i "s/#define VERSION.*/#define VERSION \"$(VERSION)\"/" src/*.h
	git pull && \
	git commit -m "New version $(VERSION)" . &&\
	git tag -a "$(VERSION)" -m "New version $(VERSION)" &&\
	git push --tags

tests: FORCE
	./run_tests.sh

FORCE: 
