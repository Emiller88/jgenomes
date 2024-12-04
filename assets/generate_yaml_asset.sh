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

# Remove existing assets
rm -rf igenomes/

# Generate base info in species/genome/build.yml

# All source fasta.fai
# ALL fais are coming from a fasta file of the same name
# Hence I use it to generate fasta + fai (and catch with that the fasta that are not following the gemome.fa name scheme)

cat ngi-igenomes_file_manifest.txt | grep "\.fai" | grep -v "Bowtie2Index" | grep -v "fai\.gz" > all_fai.txt

for i in `cat all_fai.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    mkdir -p igenomes/${species}/${genome}
    if [ -f "igenomes/${species}/${genome}/${build}.yml" ]; then
        build+="_2_"
    fi
    echo "- genome: '${build}'" > igenomes/${species}/${genome}/${build}.yml
    echo "  fasta: '${i::-4}'" >> igenomes/${species}/${genome}/${build}.yml
    echo "  source: '${genome}'" >> igenomes/${species}/${genome}/${build}.yml
    echo "  species:  '${species}'" >> igenomes/${species}/${genome}/${build}.yml
    echo "  fasta_fai: '${i}'" >> igenomes/${species}/${genome}/${build}.yml
done

# # All source fasta.dict
# cat ngi-igenomes_file_manifest.txt | grep "\.dict" | grep -v "dict\.gz" | grep -v "dict\.old" > all_dict.txt

# # Generate base info in species/genome/build.yml

# for i in `cat all_dict.txt`;
# do
#     species=$(echo $i | cut -d "/" -f 5)
#     genome=$(echo $i | cut -d "/" -f 6)
#     build=$(echo $i | cut -d "/" -f 7)

#     echo "  fasta_dict: '${i}'" >> igenomes/${species}/${genome}/${build}.yml
# done

#  Homo_sapiens/GATK/GRCh37.yml should actually be Homo_sapiens/GATK/GRCh37decoy.yml
mv igenomes/Homo_sapiens/GATK/GRCh37.yml igenomes/Homo_sapiens/GATK/GRCh37decoy.yml
