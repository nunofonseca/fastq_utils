#FROM amazonlinux:latest
FROM humancellatlas/upload-validator-base-alpine


LABEL maintainer="nuno.fonseca at gmail.com"
#RUN yum update -y && yum install -y bzip2-devel ncurses-devel bzip2 zlib-devel git gcc wget make xz-devel && yum clean all
RUN apk update && apk add alpine-sdk && apk add zlib-dev && apk add ncurses-dev && apk add bzip2-dev && apk add xz-dev

RUN ./install_deps.sh && make && make install && rm -rf test* && cp bin/* /usr/bin

FROM humancellatlas/upload-validator-base-alpine

COPY --from=0 /usr/bin/fastq_info /validator
RUN chmod +x /validator
ADD wrapper /wrapper
RUN chmod +x /wrapper/validation_wrapper.py

