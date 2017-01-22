#!/bin/make -f

.PHONY:
all: build
	@cd build && cmake .. 

build: 
	@mkdir -p build


