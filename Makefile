#!/bin/make -f

define cmake-prepare
	@cp conf/CMakeLists.txt .
	@cp conf/test/CMakeLists.txt test
endef

define cmake-cleanup
	@rm CMakeLists.txt
	@rm test/CMakeLists.txt
endef

all: build 
	$(call cmake-prepare)
	@cd build && cmake .. && make
	@cd ..
	$(call cmake-cleanup)

verbose: build
	$(call cmake-prepare)
	@cd build && cmake .. && VERBOSE=1 make
	@cd ..
	$(call cmake-cleanup)

build: 
	@mkdir -p build

clean: build
	@cd build && make clean

doc:
	@cp conf/Doxyfile .
	@doxygen Doxyfile
	@rm Doxyfile

doc-offline:
	xdg-open doc/html/index.html

doc-online:
	xdg-open https://mrquincle.github.com/noparama

.ONESHELL:
update-version:
	@sed -i 's/^:Version: .*$$/:Version: 0.1.$(shell git rev-list master --count)/' README.rst

.ONESHELL:
push: update-version
	git add -u .
	git commit
	git push

.ONESHELL:
push-doc: push doc
	cd doc/html
	git add .
	git commit
	git push

list:
	 @$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs

# The only non-PHONY target is build, which should run (with mkdir -p build) if it does not exist

.PHONY: all clean doc doc-offline doc-online list push verbose
