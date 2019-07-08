FROM humancellatlas/upload-validator-base-alpine:22

LABEL maintainer="nuno.fonseca at gmail.com"

#RUN yum update -y && yum install -y bzip2-devel ncurses-devel bzip2 zlib-devel git gcc wget make xz-devel && yum clean all
RUN apk update && apk add alpine-sdk && apk add bash && apk add zlib-dev && apk add ncurses-dev && apk add bzip2-dev && apk add xz-dev
ADD deps src tests sh ./
ADD wrapper ./wrapper
COPY install_deps.sh .
RUN chmod a+x install_deps.sh

RUN ls -l
RUN ./install_deps.sh && make && make install && rm -rf test* && cp bin/* /usr/bin

ADD wrapper /wrapper
RUN chmod +x /wrapper/validation_wrapper.py
RUN ln -s /wrapper/validation_wrapper.py validator

