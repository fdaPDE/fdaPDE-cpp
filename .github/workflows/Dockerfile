# define base OS environment
FROM alpine:latest

# install base dependencies for compiling and testing fdaPDE

# base system dependencies
RUN apk --no-cache add --update \
    g++ \
    libstdc++ \
    cmake \
    git

# c++ library dependencies
RUN apk add --update \
    eigen

# test dependencies
RUN apk add --update \
    gtest-dev

RUN git clone https://github.com/AlePalu/fdaPDE/tree/develop