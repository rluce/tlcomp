MATLABBIN ?= matlab

test:
	$(MATLABBIN) -batch "addpath('matlab'); runtests test"

.PHONY: test
