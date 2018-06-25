#SUBDIRS := $(wildcard */.)
SUBDIRS := graph_converter/undirected_csr \
	graph_cleaner \
	rank_by_degree

all: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
		done
	mv graph_converter/undirected_csr/text_to_bin toolkit/
	mv rank_by_degree/tc toolkit
	mv graph_cleaner/cleaner toolkit

clean: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
		done

.PHONY: all $(SUBDIRS)

