echo "destination: $1"

mkdir -p $1

mv ./data.*.* $1
mv ./pluto.0.log $1
mv ./*.out $1
cp ./init.c $1
cp ./definitions.h $1
cp ./pluto.ini $1

pvpython ~/alex_data/Diploma/dataprocessing/paraviewconverter.py $1
