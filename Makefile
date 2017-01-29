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

doc-online:
	xdg-open https://mrquincle.github.com/noparama

doc-offline:
	xdg-open doc/html/index.html

.PHONY: all clean doc push doc-online
