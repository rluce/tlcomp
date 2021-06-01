MATLABBIN ?= matlab

RELEASETAG = $(shell git describe --tag)
RELEASENAME = tlcomp-$(RELEASETAG)

test:
	$(MATLABBIN) -batch "addpath('matlab'); runtests test"

release:
	rm -rf build
	mkdir -p build/$(RELEASENAME)
	cp matlab/*.m build/$(RELEASENAME)
	cp README.md build/$(RELEASENAME)
	cp LICENSE build/$(RELEASENAME)
	tar -C build -czf $(RELEASENAME).tar.gz $(RELEASENAME)

clean:
	rm -rf build
	rm -f tlcomp-*.tar.gz

.PHONY: test clean
