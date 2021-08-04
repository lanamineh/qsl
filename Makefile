# Have a look at README.md for instructions on how to build the
# documentation. If you have doxygen and sphinx installed, you
# can just run `make docs`. The rest of the project uses cmake.
# See the installation and building instructions in the documentation
# for how to make the examples and use the library.

# Build the documentation
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = docs
BUILDDIR      = docs
.PHONY: docs
docs:
#	doxygen docs/Doxyfile
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: clean
clean:
	-rm -rf $(SOURCEDIR)/doxy $(SOURCEDIR)/html $(SOURCEDIR)/doctrees
