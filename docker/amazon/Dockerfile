FROM amazonlinux:latest

LABEL maintainer="nuno.fonseca at gmail.com"
RUN yum update -y && yum install -y bzip2-devel ncurses-devel bzip2 zlib-devel git gcc wget make xz-devel tar && yum clean all
RUN git clone https://github.com/nunofonseca/fastq_utils.git && cd fastq_utils && ./install_deps.sh && make && make install && rm -rf test* && cp bin/* /usr/bin
