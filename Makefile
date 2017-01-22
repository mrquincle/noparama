#!/bin/make -f

all: build
	@cd build && cmake ..  && make

build: 
	@mkdir -p build

clean:
	@cd build && make clean

doc:
	@doxygen Doxyfile

.ONESHELL:
push: doc
	git add .
	git commit 
	git push
	cd doc/html
	git add .
	git commit
	git push

.PHONY: all clean doc push
