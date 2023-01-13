FUTHARK?=futhark

TESTS=simplex-test.fut

.PHONY: test
test: $(TESTS)
	$(FUTHARK) test --pass-compiler-option=--Werror $^

.PHONY: doc
doc:
	$(FUTHARK) doc -o doc .

.PHONY: clean
clean:
	rm -rf *~ *.c simplex simplex-test *.actual *.expected doc
