input="simulated_data_maf.bcf"

### list of ref panel and query inds
bcftools query -l $input | head -n200 > ref.lst
bcftools query -l $input | tail -n100 > query.lst

### samplemap for ref panel
yes "0" | head -n100 > samplemap.lst
yes "1" | head -n100 >> samplemap.lst
# yes "2" | head -n100 >> samplemap.lst
# yes "3" | head -n100 >> samplemap.lst
paste ref.lst samplemap.lst > samplemap.tsv
rm samplemap.lst

### make reference file
bcftools view $input -S ref.lst -O b -o ref.bcf.gz --threads 10
bcftools index ref.bcf.gz

### make query file
bcftools view $input -S query.lst -O b -o query.bcf.gz --threads 10
bcftools index query.bcf.gz

### make genetic map
bcftools view query.bcf.gz | grep -v ^# | awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $2 * 1.3e-6 }' > genmap.tsv

### run rfmix
rfmix -f query.bcf.gz -r ref.bcf.gz -m samplemap.tsv -g genmap.tsv -o rfmix --chromosome=1 --n-threads=64  --random-seed=1