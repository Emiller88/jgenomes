#!/usr/bin/env bash

# Download or regenerate manifest file

# cf https://github.com/ewels/AWS-iGenomes/pull/22

# wget https://raw.githubusercontent.com/ewels/AWS-iGenomes/refs/heads/master/ngi-igenomes_file_manifest.txt

# aws s3 --no-sign-request ls --recursive s3://ngi-igenomes/igenomes/ | cut -d "/" -f 2- > tmp
# for i in `cat tmp`; do
#     echo s3://ngi-igenomes/igenomes/$i >> manifest
# done
# mv manifest ngi-igenomes_file_manifest.txt
# rm tmp

# All source fasta
cat ngi-igenomes_file_manifest.txt | grep "\.fa" | grep -v "\.fa." | grep WholeGenomeFasta | grep -v CT_conversion | grep -v GA_conversion > all_fastas.txt

# Generate base info in species/genome/build.yml

for i in `cat all_fastas.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    mkdir -p igenomes/${species}/${genome}
    if [ -f "igenomes/${species}/${genome}/${build}.yml" ]; then
        build+="_2_"
    fi
    echo "- genome: '${build}'" > igenomes/${species}/${genome}/${build}.yml
    echo "  fasta: '${i}'" >> igenomes/${species}/${genome}/${build}.yml
    echo "  source: '${genome}'" >> igenomes/${species}/${genome}/${build}.yml
    echo "  species:  '${species}'" >> igenomes/${species}/${genome}/${build}.yml
done
