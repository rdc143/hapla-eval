import msprime
import tspop
import time
import subprocess
import pandas as pd
import numpy as np

t0 = time.time()

pop_size = 5000
sequence_length = 1e8
seed = 1337
rho = 1.3e-8 # Recombination rate
mu = 1.25e-8  # Mutation rate

# Make the Demography object.
demography = msprime.Demography()
demography.add_population(name="A", initial_size=pop_size)
demography.add_population(name="B", initial_size=pop_size)
demography.add_population(name="C", initial_size=pop_size*0.5)
demography.add_population(name="D", initial_size=pop_size*2)
demography.add_population(name="ADMIX", initial_size=pop_size)
demography.add_population(name="ADMIX2", initial_size=pop_size)
demography.add_population(name="ANC_AB", initial_size=pop_size)
demography.add_population(name="ANC_CD", initial_size=pop_size)
demography.add_population(name="ANC_ABCD", initial_size=pop_size)
demography.set_migration_rate(source="B", dest="D", rate=.005)
demography.add_admixture(
    time=20, derived="ADMIX2", ancestral=["ADMIX", "C"], proportions=[0.8, 0.2])
demography.add_admixture(
    time=40, derived="ADMIX", ancestral=["A", "B"], proportions=[0.5, 0.5])
demography.add_census(time=50)
# demography.add_symmetric_migration_rate_change(time=200, populations=["D","B"], rate=100)
demography.add_population_split(
    time=1000, derived=["A", "B"], ancestral="ANC_AB")
demography.add_population_split(
    time=1400, derived=["C", "D"], ancestral="ANC_CD")
demography.add_population_split(
    time=1500, derived=["ANC_AB", "ANC_CD"], ancestral="ANC_ABCD")

# Simulate.
ts = msprime.sim_ancestry(
    samples={"A": 100, "B": 100, "C": 100, "D": 100, "ADMIX2" : 100},
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

pa = tspop.get_pop_ancestry(ts, census_time=50)
# print(pa)

st = pa.squashed_table
# print(st)

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
