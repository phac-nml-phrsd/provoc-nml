# Parsing NSGB into csvs based on lineage
# Imma try using awk. We'll see how this goes

cd ../provoc/data-raw

# GET the latest metadata
# wget --show-progress -O "data-raw/metadata_$(date +"%Y-%m-%d").tsv.gz" https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz

variants=(BA.1 BA.2 B.1.1.529 B.1.617.2 AY.4.2 AY.4 AY.25 AY.25.1 AY.24 AY.43)
for variant in ${variants[@]}
do
    echo ${variant}
    date
    zcat metadata_2022-06-20.tsv.gz | awk -F "\t" -v var="${variant}" '{ if (NR == 1 || $20 == var)  print $6 "; " $7 "; " $8 "; " $15 "; " $20 "; " $28"; "$42"; "$43"; "$44"; "$45"; "$46"; "$47"; "$48 "; "$49}' > "../../provoc-nml/output/${variant}_samples.txt" 
done

cd ../../provoc_nml
Rscript nsgb_samples.R
