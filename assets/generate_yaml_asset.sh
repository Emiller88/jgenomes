#!/usr/bin/env bash

# Download or regenerate manifest file

# Command line flags
while getopts ":adr" opt; do
    case $opt in
        a)
            NEW_MANIFEST=true
            ;;
        d)
            DOWNLOAD_MANIFEST=true
            ;;
        r)
            REGENERATE_MANIFEST=true
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            exit 1;
            ;;
    esac
done

if [[ ${DOWNLOAD_MANIFEST} ]]; then
    rm -f manifest.txt
    echo "Downloading manifest"
    wget https://raw.githubusercontent.com/ewels/AWS-iGenomes/refs/heads/master/ngi-igenomes_file_manifest.txt -O manifest.txt
fi

if [[ ${REGENERATE_MANIFEST} ]]; then
    rm -f manifest.txt
    echo "Regenerating manifest"
    # cf https://github.com/ewels/AWS-iGenomes/pull/22
    aws s3 --no-sign-request ls --recursive s3://ngi-igenomes/igenomes/ | cut -d "/" -f 2- > tmp
    for i in `cat tmp`; do
        if [[ ! $i =~ /$ ]]; then
            echo s3://ngi-igenomes/igenomes/$i >> tmp_manifest
        fi
    done
    mv tmp_manifest manifest.txt

    rm tmp
fi

if [[ ${NEW_MANIFEST} ]]; then
    rm -f manifest.txt
    echo "Regenerating NEW manifest"
    # cf https://github.com/ewels/AWS-iGenomes/pull/22

    aws s3 --profile igenomes ls --recursive s3://nf-core-references-scratch/genomes/ | grep -v "pipeline_info" | cut -d "/" -f 2- > tmp
    for i in `cat tmp`; do
        if [[ ! $i =~ /$ ]]; then
            echo s3://nf-core-references-scratch/genomes/$i >> tmp_manifest
        fi
    done
    mv tmp_manifest manifest.txt

    rm tmp
fi

total_files=$(wc -l manifest.txt | cut -d " " -f 1)

echo "Number of files in manifest: $total_files"

cp manifest.txt leftover_manifest.txt

# Remove existing assets
rm -rf igenomes/

# Generate base info in species/genome/build.yml

# All source fasta.fai
# ALL fai are coming from a fasta file of the same name
# Hence I use it to generate fasta + fai (and catch with that the fasta that are not following the gemome.fa name scheme)

cat manifest.txt | grep "\.fai" | grep -v "Bowtie2Index" | grep -v "fai\.gz" > tmp_fai.txt

echo "Populating assets for fasta and fai"

for i in `cat tmp_fai.txt`;
do
    # Create a list of fasta files
    echo "${i::-4}" >> tmp_fasta.txt

    # Grab metadata
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
cat manifest.txt | grep "README" | grep -v "Archives" | grep -v "beagle" | grep -v "plink" | grep -v "PhiX\/Illumina\/RTA\/Annotation\/README\.txt" > tmp_readme.txt

echo "Populating assets for README"

for i in `cat tmp_readme.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  readme: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source gtf (removing the onces coming from gencode)
cat manifest.txt | grep "\.gtf" | grep -v "gtf\." | grep -v "STARIndex" | grep -v "Genes\.gencode" > tmp_gtf.txt

echo "Populating assets for GTF"

for i in `cat tmp_gtf.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  gtf: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source fasta.dict
cat manifest.txt | grep "\.dict" | grep -v "dict\.gz" | grep -v "dict\.old" > tmp_dict.txt

echo "Populating assets for fasta.dict"

for i in `cat tmp_dict.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  fasta_dict: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source genes.bed
cat manifest.txt | grep "genes\.bed" > tmp_bed.txt

echo "Populating assets for genes.bed"

for i in `cat tmp_bed.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  genes_bed: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source BowtieIndex
cat manifest.txt | grep "BowtieIndex" | grep -v "MDSBowtieIndex" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_bowtie.txt

echo "Populating assets for BowtieIndex"

for i in `cat tmp_bowtie.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  bowtie1_index: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source Bowtie2Index
cat manifest.txt | grep "Bowtie2Index" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_bowtie2.txt

echo "Populating assets for Bowtie2Index"

for i in `cat tmp_bowtie2.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  bowtie2_index: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source BWAIndex (we have version0.6.0, version0.5.x, and no version specified)
cat manifest.txt | grep "BWAIndex" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_bwaindex.txt

echo "Populating assets for BWAIndex"

for i in `cat tmp_bwaindex.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    # Remove existing bwamem1_index if present
    # So that we only keep the latest version
    sed -i '\|bwamem1_index|d' igenomes/${species}/${genome}/${build}.yml
    echo "  bwamem1_index: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source BWAmem2mem
cat manifest.txt | grep "BWAmem2Index" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_bwamem2mem.txt

echo "Populating assets for BWAmem2Index"

for i in `cat tmp_bwamem2mem.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  bwamem2_index: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source Dragmap
cat manifest.txt | grep "dragmap" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_dragmap.txt

echo "Populating assets for DragmapHashtable"

for i in `cat tmp_dragmap.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  dragmap_hashtable: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source BismarkIndex
cat manifest.txt | grep "BismarkIndex\/genome\.fa" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_bismark.txt

echo "Populating assets for BismarkIndex"

for i in `cat tmp_bismark.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  bismark_index: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source star Index
cat manifest.txt | grep "STARIndex" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_star.txt

echo "Populating assets for STARIndex"

for i in `cat tmp_star.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  star_index: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source Chromosomes fasta
cat manifest.txt | grep "Chromosomes" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_chromosomes.txt

echo "Populating assets for Chromosomes fasta"

for i in `cat tmp_chromosomes.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  chromosomes_fasta: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source AbundantSequences fasta
cat manifest.txt | grep "AbundantSequences" | rev | cut -d "/" -f 2- | rev | sort -u > tmp_abundantsequences.txt

echo "Populating assets for AbundantSequences fasta"

for i in `cat tmp_abundantsequences.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  abundantsequences_fasta: \"${i}/\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source refFlat (removing the ones coming from gencode)
cat manifest.txt | grep "refFlat\.txt" | grep -v "\.gz\.bak" | grep -v "Genes\.gencode" > tmp_refflat.txt

echo "Populating assets for refFlat"

for i in `cat tmp_refflat.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  genes_refflat: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source refgene
cat manifest.txt | grep "refGene\.txt" | grep -v "Archives" | grep -v "\.gz\.bak" > tmp_refgene.txt

echo "Populating assets for refgene"

for i in `cat tmp_refgene.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  genes_refgene: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source ChromInfo.txt
cat manifest.txt | grep "ChromInfo\.txt" | grep -v "Archives" > tmp_chrominfo.txt

echo "Populating assets for ChromInfo.txt"

for i in `cat tmp_chrominfo.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  chrom_info: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source GenomeSize.xml (removing the old ones)
cat manifest.txt | grep "GenomeSize\.xml" | grep -v "\.old" > tmp_GenomeSize.txt

echo "Populating assets for GenomeSize.xml"

for i in `cat tmp_GenomeSize.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  genome_size_xml: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source hairpin.fa
cat manifest.txt | grep "SmallRNA\/hairpin\.fa" > tmp_hairpin.txt

echo "Populating assets for hairpin.fa"

for i in `cat tmp_hairpin.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  hairpin_fasta: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source mature.fa
cat manifest.txt | grep "SmallRNA\/mature\.fa" > tmp_mature.txt

echo "Populating assets for mature.fa"

for i in `cat tmp_mature.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)

    echo "  mature_fasta: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

# All source vcf
cat manifest.txt | grep "\.vcf" | grep -v "\.idx" | grep -v "\.tbi" | grep -v "\.md5" > tmp_vcf.txt

echo "Populating assets for vcf"

for i in `cat tmp_vcf.txt`;
do
    species=$(echo $i | cut -d "/" -f 5)
    genome=$(echo $i | cut -d "/" -f 6)
    build=$(echo $i | cut -d "/" -f 7)
    source_vcf=$(echo $i | cut -d "/" -f 9)
    filename=$(basename $i)
    if [ -z "${source_vcf}" ]; then
        source_vcf="unknown"
    elif [ "${source_vcf}" == "Variation" ]; then
        source_vcf="unknown"
    elif [ "${source_vcf}" == "${filename}" ]; then
        source_vcf="unknown"
    fi

    echo "- genome: \"${build}\"" >> igenomes/${species}/${genome}/${build}.yml
    echo "  source: \"${genome}\"" >> igenomes/${species}/${genome}/${build}.yml
    echo "  species: \"${species}\"" >> igenomes/${species}/${genome}/${build}.yml
    echo "  source_vcf: \"${source_vcf}\"" >> igenomes/${species}/${genome}/${build}.yml
    echo "  vcf: \"${i}\"" >> igenomes/${species}/${genome}/${build}.yml
done

echo "Fixing /Homo_sapiens/GATK/GRCh37.yml name to Homo_sapiens/GATK/GRCh37decoy.yml"

#  Homo_sapiens/GATK/GRCh37.yml should actually be Homo_sapiens/GATK/GRCh37decoy.yml
mv igenomes/Homo_sapiens/GATK/GRCh37.yml igenomes/Homo_sapiens/GATK/GRCh37decoy.yml

echo "Deleting assets from manifest"

# Delete all the lines that are used to generate the assets
for i in `cat tmp_*.txt`;
do
    sed -i "\|${i}|d" leftover_manifest.txt
done

leftover_files=$(wc -l leftover_manifest.txt | cut -d " " -f 1)
echo "Number of files in leftover manifest: $leftover_files"

let removed_files=$total_files-$leftover_files
echo "Number of files removed: $removed_files"

echo "Deleting tmp files"

rm -rf tmp_*
