VERSION=0.19.3
# Requires zlib and samtools 0.1.9
all:
	make -C src clean
	make -C src


install:
	mkdir -p bin && make -C src install && cp sh/*.sh sh/fastq2bam bin

release:
	sed -i "s/#define VERSION.*/#define VERSION \"$(VERSION)\"/" src/*.c
	sed -i "s/#define VERSION.*/#define VERSION \"$(VERSION)\"/" src/*.h
	sed -i "s/^VERSION=.*/VERSION=\"$(VERSION)\"/" sh/*
	git pull && \
	git commit -m "New version $(VERSION)" . &&\
	git push &&\
	git tag -a "$(VERSION)" -m "New version $(VERSION)" &&\
	git push --follow-tags
#--tags

tests: FORCE
	./run_tests.sh

FORCE: 


