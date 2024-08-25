K=$(cut -f2 $1/samplemap.tsv | uniq | wc -l)
mkdir $1_paramtest
cd $1_paramtest

for maf in 0 0.005 0.01 0.02 0.05
do
# gen new maf filtered bcf
bcftools view -O b -o "simulated_data_${maf}.bcf" ../$1/simulated_data.bcf -i "INFO/MAF > $maf"

for wd in 8 16 32 64 128 256
do

for lambda in 0.05 0.1 0.2
do

# hapla
hapla cluster --bcf "simulated_data_${maf}.bcf" --fixed $wd --threads 32 --out "simulated_data_${maf}_${wd}_${lambda}" -l $lambda #--min-freq 0.1
hapla admix --clusters "simulated_data_${maf}_${wd}_${lambda}.z.npy" --K $K --seed 1 --threads 32 --out "simulated_data_${maf}_${wd}_${lambda}"

for alpha in 1e-4 1e-6 1e-8 1e-10 1e-12 1e-14
do
# fwd-bwd
hapla fatash --clusters "simulated_data_${maf}_${wd}_${lambda}.z.npy" --qfile "simulated_data_${maf}_${wd}_${lambda}.K$K.s1.Q" --pfile "simulated_data_${maf}_${wd}_${lambda}.K$K.s1.P.npy" --threads 32 --out "simulated_data_${maf}_${wd}_${lambda}_${alpha}" --alpha $alpha --save-posterior

# viterbi
hapla fatash --clusters "simulated_data_${maf}_${wd}_${lambda}.z.npy" --qfile "simulated_data_${maf}_${wd}_${lambda}.K$K.s1.Q" --pfile "simulated_data_${maf}_${wd}_${lambda}.K$K.s1.P.npy" --threads 32 --out "simulated_data_${maf}_${wd}_${lambda}_${alpha}_viterbi" --alpha $alpha --viterbi

# collect match stat
Rscript ../getmatch.R "simulated_data_${maf}_${wd}_${lambda}_${alpha}.path" "simulated_data_${maf}_${wd}_${lambda}.w.info" ../$1
paste meanhaplamatch.txt <(echo -e "${maf}\t${wd}\t${lambda}\t${alpha}\tfwdbwd") >> match.lst

Rscript ../getmatch.R "simulated_data_${maf}_${wd}_${lambda}_${alpha}_viterbi.path" "simulated_data_${maf}_${wd}_${lambda}.w.info" ../$1
paste meanhaplamatch.txt <(echo -e "${maf}\t${wd}\t${lambda}\t${alpha}\tviterbi") >> match.lst

rm meanhaplamatch.txt

done
done
done
done

sort -t$'\t' -k1,1nr  match.lst > match.lst.tmp
mv match.lst.tmp match.lst
