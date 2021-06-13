# Copyright (C) 2020 Lana Mineh and John Scott.
#
# This file is part of QSL, the quantum computer simulator.
#
# QSL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# QSL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with QSL.  If not, see <https://www.gnu.org/licenses/>.
#

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
	doxygen docs/Doxyfile
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: clean
clean:
	-rm -rf $(SOURCEDIR)/doxy $(SOURCEDIR)/html $(SOURCEDIR)/doctrees
