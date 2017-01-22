#!/bin/make -f

.PHONY:
all: build
	@cd build && cmake ..  && make

build: 
	@mkdir -p build

.PHONY:
clean:
	@cd build && make clean

