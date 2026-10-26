.PHONY: all test clean
JOBS ?= 4
all:
	python3 scripts/build.py --jobs $(JOBS)
test: all
	ctest --test-dir build --output-on-failure
clean:
	cmake -E rm -rf build
