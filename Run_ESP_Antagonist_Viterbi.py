import time
import argparse
from Epigenetic_Sequence_Predictor.ESP_Sim_Viterbi import run_sim

# Sample code to run this from the command line.

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--batch_name', type=str, help='Name of the batch')
parser.add_argument('--alpha_eval', type=str, help='string to evaluate as the alpha list')
parser.add_argument('--beta_eval', type=str, help='string to evaluate as the beta list')
parser.add_argument('--mu', type=float, default=0.0, help='mu value for the simulation')
parser.add_argument('--rho', type=float, default=0.5, help='rho value for the simulation')
parser.add_argument('--seed', type=int, default=42, help='Random seed')
parser.add_argument('--n_samples', type=int, default=10000, help='Number of samples')
parser.add_argument('--seq_length', type=int, default=100, help='Sequence length')
parser.add_argument('--verbose_level', type=int, default=5, help='Verbose level')

args = parser.parse_args()

full_start = time.time()

batch_name = args.batch_name
alpha_eval = args.alpha_eval
beta_eval = args.beta_eval
mu = args.mu
rho = args.rho
seed = args.seed
n_samples = args.n_samples
seq_length = args.seq_length
sim_id = batch_name
verbose_level = args.verbose_level

a_list = eval(alpha_eval)
b_list = eval(beta_eval)

with open(f'{batch_name}_BitError_Combinations.csv', 'w') as f:
    f.write('Alpha,Beta,BitError\n')

for a in a_list:
    for b in b_list:
        sim_start = time.time()
        sim_name, biterror = run_sim(sim_id=sim_id, seed=seed, alpha=round(a, 1), beta=round(b, 1), mu=mu, rho=rho,
                                     n_samples=n_samples, seq_length=seq_length, verbose_level=verbose_level)
        sim_end = time.time()
        print(f"Simulation Name = {sim_name}", end=' --- ')
        print(f"Sim Execution Time = {round(sim_end - sim_start, 4)}", end=' --- ')
        print(f"BitError = {biterror}")

        with open(f'{batch_name}_BitError_Combinations.csv', 'a') as f:
            f.write(f'{round(a, 1)},{round(b, 1)},{biterror}\n')

full_end = time.time()
print(f"Total Execution time = {full_end - full_start}")

with open(f'{batch_name}_Batch_Log.txt', 'w') as f:
    f.write(f'Log File - Batch Simulation :: {batch_name}\n')
    f.write(f"Total Execution time = {full_end - full_start}\n\n")

    f.write(f'##INPUTS##\n')
    f.write(f'sim_id = {sim_id}\n')
    f.write(f'seed = {seed}\n')
    f.write(f'alpha_eval = {alpha_eval}\n')
    f.write(f'beta_eval = {beta_eval}\n')
    f.write(f'mu = {mu}\n')
    f.write(f'rho = {rho}\n')
    f.write(f'n_samples = {n_samples}\n')
    f.write(f'seq_length = {seq_length}\n')
    f.write(f'verbose_level = {verbose_level}\n')
