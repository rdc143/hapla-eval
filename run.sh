mkdir $1 -p

#copy scripts to have record of settings
cp simK2.py ./$1
cp runrfmix.sh ./$1
cp runhapla.sh ./$1
cp plottracts.R ./$1

cd $1
#run sim and inference
python sim*.py
bash runrfmix.sh
bash runhapla.sh
Rscript plottracts.R