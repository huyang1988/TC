#SUBDIRS := $(wildcard */.)
SUBDIRS := graph_converter/undirected_csr \
	graph_cleaner \
	rank_by_degree

all: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
		done
	mv graph_converter/undirected_csr/tuple_to_undirected_csr.bin toolkit/
	mv rank_by_degree/rank_by_degree.bin toolkit
	mv graph_cleaner/graph_cleaner.bin toolkit

clean: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
		done

.PHONY: all $(SUBDIRS)

