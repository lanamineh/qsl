#
# Bash script for creating a new file from a template
#
# This script is for creating a new source file, either a
# header file (.hpp), a source file (.cpp) or a template
# file (.tpp). You can use the script in two ways. Either
# run the script with no command line arguments, in which
# case you will be prompted for a new file details; or
# you can pass the relative path of the new file (with the
# correct file extension) as the first command line argument.
#

# Bash colors
YELLOW="\033[1;33m"
GREEN="\033[1;32m"
RESET="\033[0m"

echo "This script will create a new file from one of the templates"
echo "in the templates/ directory. It will include the licence"
echo "notice and header guards if necessary."
echo
echo "You can create either a header file (with .hpp extension),"
echo "a source file (.cpp), or a template source file (.tpp)."
echo
echo "Write the relative path to the new file here, including the"
echo "filename and extension (e.g. src/somewhere/test.cpp). New"
echo "directories will be created if necessary."
echo

# tmpl contains the directory of the source file
# templates, which are used to generate new
# header, source and template files
tmpl=../templates/

# Check for command line argument
if [ -z "$1" ]
then
    read -p "Enter new file with path: " newfile
else
    newfile=$1
fi


# Separate filename and directory
name=$(basename $newfile)
path=$(dirname $newfile)

# Create new directory if necessary
mkdir -p $path

# Get file extension from filename
extension="${name##*.}"
basename="${name%.*}"

case $extension in
    hpp)
	echo -e "${YELLOW}Making header file '"$newfile"'"${RESET}
	cp $tmpl/header.hpp $newfile
	;;
    cpp)
	echo -e "${YELLOW}Making source file '"$newfile"'"${RESET}
	cp $tmpl/source.cpp $newfile
	;;
    tpp)
	echo -e "${YELLOW}Making template file '"$newfile"'"${RESET}
	cp $tmpl/template.tpp $newfile
	;;
    *)
	echo "Unrecognised file extension" $extension
	exit
	;;
esac

# Use sed to replace NAME and GUARD with the file basename
sed -i "s/NAME/$basename/g" $newfile
sed -i "s/GUARD/${basename^^}/g" $newfile
sed -i '/NOTICE/{
       r '$tmpl'/NOTICE
       d
       }' $newfile

echo

exit
