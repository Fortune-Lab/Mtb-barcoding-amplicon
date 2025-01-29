cat $1 | parallel --gnu -j 8 ./download.sh {}
