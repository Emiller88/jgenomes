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

    echo "- genome: \"${build}\"" > igenomes/${species}/${genome}/${build}.yml
    echo "  fasta: \"${i::-4}\"" >> igenomes/${species}/${genome}/${build}.yml
    echo "  source: \"${genome}\"" >> igenomes/${species}/${genome}/${build}.yml
    echo "  species: \"${species}\"" >> igenomes/${species}/${genome}/${build}.yml
    echo "  fasta_fai: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source README
cat ngi-igenomes_file_manifest.txt | grep "README" | grep -v "Archives" | grep -v "beagle" | grep -v "plink" | grep -v "PhiX\/Illumina\/RTA\/Annotation\/README\.txt" > all_readme.txt

for i in `cat all_readme.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  readme: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source fasta.dict
cat ngi-igenomes_file_manifest.txt | grep "\.dict" | grep -v "dict\.gz" | grep -v "dict\.old" > all_dict.txt

for i in `cat all_dict.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  fasta_dict: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source BowtieIndex
cat ngi-igenomes_file_manifest.txt | grep "BowtieIndex" | grep -v "MDSBowtieIndex" | rev | cut -d "/" -f 2- | rev | sort -u > all_bowtie.txt

for i in `cat all_bowtie.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  bowtie1_index: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source Bowtie2Index
cat ngi-igenomes_file_manifest.txt | grep "Bowtie2Index" | rev | cut -d "/" -f 2- | rev | sort -u > all_bowtie2.txt

for i in `cat all_bowtie2.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  bowtie2_index: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# Just LAST source BWAmem1 (SO if we have version0.6.0, we keep it, otherwise it's version0.5.x, and if not the one with no version specified)
cat ngi-igenomes_file_manifest.txt | grep "BWAIndex" | rev | cut -d "/" -f 2- | rev | sort -u > all_bwaindex.txt

for i in `cat all_bwaindex.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    sed -i '/bwamem1_index/d' igenomes/${species}/${genome}/${build}.yml
    echo "  bwamem1_index: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source BWAmem2mem
cat ngi-igenomes_file_manifest.txt | grep "BWAmem2Index" | rev | cut -d "/" -f 2- | rev | sort -u > all_bwamem2mem.txt

for i in `cat all_bwamem2mem.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  bwamem2_index: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source Dragmap
cat ngi-igenomes_file_manifest.txt | grep "dragmap" | rev | cut -d "/" -f 2- | rev | sort -u > all_dragmap.txt

for i in `cat all_dragmap.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  dragmap_hashtable: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source star
cat ngi-igenomes_file_manifest.txt | grep "STARIndex" | rev | cut -d "/" -f 2- | rev | sort -u > all_star.txt

for i in `cat all_star.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  star_index: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

#  Homo_sapiens/GATK/GRCh37.yml should actually be Homo_sapiens/GATK/GRCh37decoy.yml
mv igenomes/Homo_sapiens/GATK/GRCh37.yml igenomes/Homo_sapiens/GATK/GRCh37decoy.yml
