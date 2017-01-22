#!/bin/make -f

all: build
	@cd build && cmake ..  && make

build: 
	@mkdir -p build

doc:
	doxygen Doxyfile

clean:
	@cd build && make clean

.PHONY: all clean doc
