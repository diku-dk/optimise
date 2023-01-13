FUTHARK?=futhark

.PHONY: test
test: lib/github.com/diku-dk/linalg
	$(MAKE) -C lib/github.com/diku-dk/optimise test

.PHONY: doc
doc: lib/github.com/diku-dk/linalg
	$(MAKE) -C lib/github.com/diku-dk/optimise doc

.PHONY: clean
clean:
	rm -rf *~
	$(MAKE) -C lib/github.com/diku-dk/optimise clean

.PHONY: realclean
realclean: clean
	rm -rf lib/github.com/diku-dk/cpprandom lib/github.com/diku-dk/linalg

lib/github.com/diku-dk/linalg:
	$(FUTHARK) pkg sync
