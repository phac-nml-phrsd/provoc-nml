# Running tests for https://github.com/suskraem/ww_benchmark

pwd # must be gromstole

start=`date +%s`
for i in {1..100}
do
    echo "../ww_benchmark/samples/sample${i}_R1.fastq.gz"
    python scripts/minimap2.py -o ../provoc-nml/output/gromstole ../ww_benchmark/samples/sample${i}_R1.fastq.gz ../ww_benchmark/samples/sample${i}_R2.fastq.gz
done
end=`date +%s`

runtime=$((end-start))



