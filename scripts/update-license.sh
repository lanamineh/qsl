# Copyright (C) 2020 Lana Mineh and John Scott.
#
# This file is part of qsim, the quantum computer simulator.
#
# qsim is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# qsim is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with qsim.  If not, see <https://www.gnu.org/licenses/>.
#

#
# Bash script for updating the license notice in all source files
#
# Run the script to update the license notice header in all the
# source files in the src/ tree. The license notice is taking
# from templates/NOTICE. The script removes all the content before
# the first opening doxygen comment block (which contains the
# title, etc.) and replaces it with the current version of the
# license notice.
#

# Bash colors
RED="\033[1;31m"
YELLOW="\033[1;33m"
GREEN="\033[1;32m"
RESET="\033[0m"

# Project name
proj_name=QSL

# src contains the top level source file location. This
# location is recursively searched for files ending in
# .hpp, .cpp and .tpp. 
src=../src/

# tmpl contains the directory of the source file
# templates, which are used to generate new
# header, source and template files
tmpl=../other/templates/

# Loop over all files in the source tree
for file in $(find $src -name '*.cpp' -or -name '*.hpp' -or -name '*.tpp')
do

    # max contains the maximum number of lines that grep
    # should look for the doxygen /** and \file. This is
    # not for performance, it is to reduce false matches
    # later in the file. It should be about 10 longer than
    # the number of lines in the NOTICE
    max=30
    
    # Check that they have a valid doxygen comment containing \file
    # somewhere near the top
    if grep -m$max -Fq "\file" $file
    then
	# Found \file tag
	echo -e ${GREEN}"Found file tag in "$file${RESET}

	# Delete everything before that point (the old/non-existent license)
	sed -i '/\/\*\*/,$!d' $file
	
	# Insert the new license from templates/NOTICE
	temp=$tmpl/temp
	cp $tmpl/NOTICE $temp # Duplicate NOTICE to temp
	echo "" >> $temp # Add newline to temp file
	sed -i "s/%PROJNAME%/$proj_name/g" $temp
	cat $file >> $temp # Append source file to temp
	mv $temp $file # Store temp as new source file
	
    else
	# \file tag not found
	echo -e ${RED}"No file tag in "$file". Ignoring file."${RESET}
    fi
    
done

