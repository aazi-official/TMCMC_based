SHELL=/bin/sh

SUBDIRS=bimodal_likelihood gaussian_likelihood

default: what

what:
	@echo "Specify a target."

all clean:
	@for i in $(SUBDIRS) ;\
	 do echo "Making $@ in $$i ..."; cd $$i; $(MAKE) $@; cd .. ;\
	 done
