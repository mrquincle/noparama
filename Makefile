#!/bin/make -f

all: build
	@cd build && cmake ..  && make

build: 
	@mkdir -p build

clean:
	@cd build && make clean

doc:
	doxygen Doxyfile

.ONESHELL:
push-code:
	git add -u .
	git commit 
	git push

.ONESHELL:
push: push-code doc
	cd doc/html
	git add -u .
	git commit
	git push

.PHONY: all clean doc push-code push
