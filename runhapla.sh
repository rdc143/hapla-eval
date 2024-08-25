K=2
hapla cluster --bcf simulated_data_maf.bcf --fixed 16 --threads 32 --out hapla.simulated_data -l 0.05 #--min-freq 0.1
hapla admix --clusters hapla.simulated_data.z.npy --K $K --seed 1 --threads 32 --out hapla
hapla fatash --clusters hapla.simulated_data.z.npy --qfile hapla.K$K.s1.Q --pfile hapla.K$K.s1.P.npy --threads 32 --out hapla --alpha 1e-14 --save-posterior #--viterbi
