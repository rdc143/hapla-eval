import msprime
import tspop
import time
import subprocess
import pandas as pd
import demes
import numpy as np


t0 = time.time()

sequence_length = 1e8
seed = 1337
rho = 1.3e-8 # Recombination rate
mu = 1.25e-8  # Mutation rate

# Make the Demography object.
graph = demes.load("../aa.yaml")
demography = msprime.Demography.from_demes(graph)
demography.add_census(time=13)
demography.sort_events()


# Simulate.
ts = msprime.sim_ancestry(
    samples={"AFR": 100, "EUR": 100, "EAS": 100 ,"ADMIX" : 100},
    demography=demography,
    random_seed=seed,
    sequence_length=sequence_length,
    recombination_rate=rho
)

# Simulate mutations on the tree sequence.
mutated_ts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)

# Define the output VCF and BCF file paths.
vcf_file_path = "simulated_data.vcf"
bcf_file_path = "simulated_data_maf.bcf"

# Save the mutated tree sequence to a VCF file.
with open(vcf_file_path, "w") as vcf_file:
    mutated_ts.write_vcf(vcf_file)

# Convert VCF to BCF and filter on minor allele freq
subprocess.run(["bcftools", "+fill-tags", vcf_file_path, "-O", "b", "-o", "simulated_data.bcf", "--", "-t", "MAF"])

subprocess.run(["bcftools", "view", "-O", "b", "-o", bcf_file_path, "simulated_data.bcf", "-i", "INFO/MAF > 0.05"])

pa = tspop.get_pop_ancestry(ts, census_time=13)
# print(pa)

st = pa.squashed_table
# print(st)

# Remap numbering of populations to order of sampled pops
remap = {2: 0, 4: 1, 5: 2}
st = st.replace({"population": remap})

st_path = "simulated_true_anc.csv"
st.to_csv(st_path, sep="\t", index=False)

st["tract_len"] = st["right"] - st["left"]

st['sample'] = st['sample'] // 2 # switch to individual numbering, instead of haplotype
st_agg = st[['sample', 'population', 'tract_len']].groupby(['sample', 'population'], as_index=False).sum()
st_pivoted = st_agg.pivot(index='sample', columns='population', values='tract_len').fillna(0)
Q = st_pivoted.div(st_pivoted.sum(axis=1), axis=0).to_numpy()
np.savetxt("trueQ.csv", Q, delimiter="\t", fmt='%.4f')


t1 = time.time()
print(f"Seconds elapsed: {t1-t0}")
