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

list:
	 @$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs

.PHONY: all clean doc push doc-online list
