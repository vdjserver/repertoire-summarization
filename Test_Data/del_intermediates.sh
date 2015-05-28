#!/bin/bash
#GET DIR http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
#echo $DIR
rm -v $DIR/*.vdjml
rm -v $DIR/*.out
rm -v $DIR/*.txt
rm -v $DIR/*.err
rm -v $DIR/*.json
rm -v $DIR/*.tsv
rm -v $DIR/*.log
rm -v $DIR/*.sim.fna
rm -v $DIR/*.csv
