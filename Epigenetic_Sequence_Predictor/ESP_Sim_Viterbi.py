import os
import numpy as np


##################################################################################################

def generate_mom(alpha, beta, mu, n=100, start_as=2):
    """
    This function generates one sequence of the mother based on alpha, beta, mu and n values.
    the alpha and beta values are used to create the A part of the mother and the mu is used to create the B.

    :param alpha: The probability that the next element node in the sequence is 1 given that the current node is 1.
    :param beta: The probability that the next element node in the sequence is 0 given that the current node is 0.
    :param mu: The probability that the element in the B sequence is 1 (only if A is 0).
    :param n: The length of the mother sequence.
    :param start_as: The first node of the mother sequence.
    :return: a numpy vector containing the sequence of the mother.
    """
    # Initializing the mother sequence
    mother_sequence_A = np.zeros([n], dtype=np.float32)
    mother_sequence_B = np.zeros([n], dtype=np.float32)
    mother_sequence = np.zeros([n], dtype=np.float32)

    # Generating the A Part.
    # if 0 or 1 is not passed as a sequence start, it is randomly assigned
    if start_as in [0, 1]:
        mother_sequence_A[0] = start_as
    else:
        mother_sequence_A[0] = np.random.randint(2)

    # Filling in the values of the sequence based on alpha and beta.
    for i in range(1, mother_sequence_A.size):
        if mother_sequence_A[i - 1] == 1:
            # The previous node is 1 so use alpha to set 1 else set 0
            if np.random.random() <= alpha:
                mother_sequence_A[i] = 1
            else:
                mother_sequence_A[i] = 0
        else:
            # The previous node is 0 so use beta to set 0 else set 1
            if np.random.random() <= beta:
                mother_sequence_A[i] = 0
            else:
                mother_sequence_A[i] = 1

    # Generating the B part.
    for i in range(mother_sequence_B.size):
        if mother_sequence_A[i] == 0:
            if np.random.random() <= mu:
                mother_sequence_B[i] = 1
            else:
                mother_sequence_B[i] = 0
        else:
            pass

            # Generating the complete mother
    for i in range(mother_sequence_A.size):
        if mother_sequence_A[i] == 0 and mother_sequence_B[i] == 0:
            mother_sequence[i] = 0
        elif mother_sequence_A[i] == 1 and mother_sequence_B[i] == 0:
            mother_sequence[i] = 1
        elif mother_sequence_A[i] == 0 and mother_sequence_B[i] == 1:
            mother_sequence[i] = 2
        # mother_sequence[i] = 3 if mother_sequence_A[i] == 1 and mother_sequence_B[i] == 1
        # Flipping the binary counting order to help the network avoid errors while having to predict just 0 and 2.
        # Now it just has to predict 0 and 1 (if mu = 0)

    return mother_sequence


def split_sequence(seq):
    seq_A = np.copy(seq)
    seq_B = np.copy(seq)

    for i in range(seq.size):
        if seq[i] == 0:
            seq_A[i] = 0
            seq_B[i] = 0
        if seq[i] == 1:
            seq_A[i] = 1
            seq_B[i] = 0
        if seq[i] == 2:
            seq_A[i] = 0
            seq_B[i] = 1

    return seq_A, seq_B


def generate_daughter(mom, rho):
    """
    This function generates a daughter sequence from a mother sequence.
    It randomly mutates a node with a 50% probability of flipping.
    :param mom: The mother sequence that is to be used.
    :param rho: The probablity that mutation occurs for a bit if a mutation is possible.
    :return: a numpy vector containing the sequence of the daughter.
    """
    daughter_sequence = np.copy(mom)

    for i in range(daughter_sequence.size):
        if daughter_sequence[i] == 0:
            pass
        elif daughter_sequence[i] == 1:
            if np.random.random() <= rho:
                daughter_sequence[i] = 0
            else:
                daughter_sequence[i] = 1
        elif daughter_sequence[i] == 2:
            if np.random.random() <= rho:
                daughter_sequence[i] = 0
            else:
                daughter_sequence[i] = 2

    return daughter_sequence


def calculate_biterror(seq1, seq2):
    """
    This function will calculate the biterror for a pair of sequences.
    """

    xor = np.logical_xor(seq1, seq2).astype(int)

    biterror = np.sum(xor) / xor.size

    return biterror


class Dataframe:
    """
    This class is a dataframe which will hold all values for one simulation.
    i.e. it will contain the mothers, corrupted daughters, corrected daughters,
    the biterror rate and all other necessary data.

    We will initialise the object with an alpha, beta, number of moms,
    the starting value of the mom sequence, number of nodes in the sequences and
    a corrupted daughter.
    """

    def __init__(self, alpha, beta, mu, rho, n_moms, mom_length, mom_start=2):
        self.alpha = alpha
        self.beta = beta
        self.mu = mu
        self.rho = rho
        self.mom_list = np.zeros([n_moms, mom_length])
        self.corrupt_daughter_list = np.zeros([n_moms, mom_length])
        self.corrected_daughter_list = np.zeros([n_moms, mom_length])

        self.biterror = np.zeros([n_moms])
        # initializing moms and making corrupt daughters
        for i in range(n_moms):
            # if we don't want all the moms to start with either 1 or 0,
            # we pass anything else to run it as a random node.
            new_mom = generate_mom(alpha=self.alpha,
                                   beta=self.beta,
                                   mu=self.mu,
                                   n=mom_length,
                                   start_as=mom_start)
            self.mom_list[i] = new_mom
            self.corrupt_daughter_list[i] = generate_daughter(new_mom, self.rho)

    def calc_errors(self):
        """
        This function will calculate the biterror for the mom_list and
        corrected_daughter_list.
        """
        for i in range(self.biterror.size):
            mom_A, mom_B = split_sequence(self.mom_list[i])
            daughter_A, daughter_B = split_sequence(self.corrected_daughter_list[i])

            self.biterror[i] = calculate_biterror(mom_A, daughter_A)

    def conclusions(self, sim_name='Untitled_Sim', level=3):
        """
        This function gives a summary of the dataframe.
        :param level: This specifies how much information is to be displayed.
        :param sim_name: Name of the simulation with which the folder is created to store the files.
        """
        # Create the log Files
        if not os.path.exists(sim_name):
            os.makedirs(sim_name)
        self.calc_errors()
        np.savetxt(sim_name + '/BitError.csv', self.biterror, delimiter=',')
        np.savetxt(sim_name + '/Mom_list.csv', self.mom_list, delimiter=',', fmt='%d')
        np.savetxt(sim_name + '/Corrupt_Daughter.csv', self.corrupt_daughter_list, delimiter=',', fmt='%d')
        np.savetxt(sim_name + '/Correct_Daughter.csv', self.corrected_daughter_list, delimiter=',', fmt='%d')
        if level > 0:
            print(f'Alpha = {self.alpha}')
            print(f'Beta = {self.beta}')
            print(f'Average BitError = {np.mean(self.biterror)}')
            if level >= 2:
                print(f'BitError = {self.biterror}')
                if level >= 3:
                    print('\nMom List')
                    print(self.mom_list)

                    print('\nCorrupt Daughter List')
                    print(self.corrupt_daughter_list)

                    print('\nCorrected Daughter List')
                    print(self.corrected_daughter_list)

        return 0


def decode_zeros(num_interim_zeros, alpha=0.2, beta=0.8, mu=0.5, rho=0.5, start_node=1, end_node=1):
    branch_metric = np.zeros(4)
    # yn = 0, xn = 0, xn-1 = 1
    branch_metric[0] = (1 - mu + (mu * rho)) * (1 - alpha)
    # yn = 0, xn = 1, xn-1 = 1
    branch_metric[1] = rho * alpha
    # yn = 0, xn = 0, xn-1 = 0
    branch_metric[2] = (1 - mu + mu * rho) * beta
    # yn = 0, xn = 1, xn-1 = 0
    branch_metric[3] = rho * (1 - beta)

    decoded_sequence = np.zeros(num_interim_zeros)

    node_path_metric_xn_0 = np.zeros((num_interim_zeros + 1, 2))
    node_path_metric_xn_1 = np.zeros((num_interim_zeros + 1, 2))

    for q in range(0, num_interim_zeros + 1):
        if q == 0:
            if start_node == 1:
                node_path_metric_xn_0[q, 0] = branch_metric[0]
                node_path_metric_xn_1[q, 0] = branch_metric[1]
                node_path_metric_xn_0[q, 1] = 1
                node_path_metric_xn_1[q, 1] = 1
            elif start_node == 2:
                node_path_metric_xn_0[q, 0] = branch_metric[2]
                node_path_metric_xn_1[q, 0] = branch_metric[3]
                node_path_metric_xn_0[q, 1] = 0
                node_path_metric_xn_1[q, 1] = 0
            else:
                print('Start Node Error')

        else:
            if end_node == 2 and q == num_interim_zeros:
                path_metric_1_0 = node_path_metric_xn_1[q - 1, 0] * (rho * mu * (1 - alpha))
                path_metric_0_0 = node_path_metric_xn_0[q - 1, 0] * (mu * (1 - rho) * beta)
            else:
                path_metric_1_0 = node_path_metric_xn_1[q - 1, 0] * branch_metric[0]
                path_metric_0_0 = node_path_metric_xn_0[q - 1, 0] * branch_metric[2]

            if path_metric_1_0 > path_metric_0_0:
            # if path_metric_1_0 >= path_metric_0_0:
                node_path_metric_xn_0[q, 0] = path_metric_1_0
                node_path_metric_xn_0[q, 1] = 1
            else:
                node_path_metric_xn_0[q, 0] = path_metric_0_0
                node_path_metric_xn_0[q, 1] = 0

            if end_node == 1 and q == num_interim_zeros:
                path_metric_0_1 = node_path_metric_xn_0[q - 1, 0] * (1 - rho) * (1 - beta)
                path_metric_1_1 = node_path_metric_xn_1[q - 1, 0] * (1 - rho) * alpha
            else:
                path_metric_0_1 = node_path_metric_xn_0[q - 1, 0] * branch_metric[3]
                path_metric_1_1 = node_path_metric_xn_1[q - 1, 0] * branch_metric[1]

            if path_metric_1_1 > path_metric_0_1:
            #if path_metric_1_1 >= path_metric_0_1:
                node_path_metric_xn_1[q, 0] = path_metric_1_1
                node_path_metric_xn_1[q, 1] = 1
            else:
                node_path_metric_xn_1[q, 0] = path_metric_0_1
                node_path_metric_xn_1[q, 1] = 0

    # print(node_path_metric_xn_0)
    # print(node_path_metric_xn_1)
    # print('The decoded sequnce is empty ::', decoded_sequence)
    if end_node == 1:
        decoded_sequence[num_interim_zeros - 1] = node_path_metric_xn_1[num_interim_zeros, 1]
    elif end_node == 2:
        decoded_sequence[num_interim_zeros - 1] = node_path_metric_xn_0[num_interim_zeros, 1]
    else:
        print("End Node Error")

    # print('The decoded sequnce\'s last bit is added based on boundary condition. ::', decoded_sequence)

    for q in range(num_interim_zeros - 1, 0, -1):
        if decoded_sequence[q] == 0:
            decoded_sequence[q - 1] = node_path_metric_xn_0[q, 1]
        elif decoded_sequence[q] == 1:
            decoded_sequence[q - 1] = node_path_metric_xn_1[q, 1]
        # print(f'adding {q}th bit ::', decoded_sequence)

    # print('The decoded sequnce is empty ::', decoded_sequence)

    return decoded_sequence


def viterbi_decode_antagonistic(daughter_chromatin, alpha, beta, mu, rho):
    # We are correcting A using B as the antagonistic part.
    # We are also assuming that the sequence always starts with a 1.

    corrected_sequence = []

    count = 1
    start_node = 1
    end_node = None

    for i in range(1, len(daughter_chromatin)):
        if daughter_chromatin[i] == daughter_chromatin[i - 1]:
            count += 1
        else:
            end_node = daughter_chromatin[i]
            if daughter_chromatin[i - 1] == 0:
                decoded_zeros = decode_zeros(num_interim_zeros=count, alpha=alpha, beta=beta, mu=mu, rho=rho,
                                             start_node=start_node, end_node=end_node)
                # print(f's({start_node})',str(daughter_chromatin[i-1])*count,f'e({end_node})',f' :: {decoded_zeros}')
                corrected_sequence.extend(decoded_zeros.astype(int))
            else:
                temp = [daughter_chromatin[i - 1] for j in range(count)]
                corrected_sequence.extend(temp)
            start_node = daughter_chromatin[i - 1]
            count = 1

    # Checking for Edge Cases
    # Assuming that if the sequence is ending with 0s the end node to be 1
    end_node = 1
    if count >= 1 and len(corrected_sequence) < len(daughter_chromatin):
        if daughter_chromatin[-1] == 0:
            decoded_zeros = decode_zeros(num_interim_zeros=count, alpha=alpha, beta=beta, mu=mu, rho=rho,
                                         start_node=start_node, end_node=end_node)
            # print(f's({start_node})',str(daughter_chromatin[i-1])*count,f'e({end_node})',f' :: {decoded_zeros}')
            # temp = [3 for j in range(count)]
            corrected_sequence.extend(decoded_zeros.astype(int))
        else:
            temp = [daughter_chromatin[-1] for j in range(count)]
            corrected_sequence.extend(temp)

    return corrected_sequence


##################################################################################################
# Viterbi Pattern Decoding Functions

import csv
import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches


def find_transition_point(alpha, beta, mu, rho, start_node, end_node):
    check_max = 100
    flag_1 = 0
    flag_0 = 0
    flag_mix = 0
    for num_interim_zeros in range(1, check_max + 1):
        result = decode_zeros(num_interim_zeros, alpha, beta, mu, rho, start_node, end_node)
        # if alpha ==0.1 and beta ==0.3:
        #     print(num_interim_zeros, result)
        #     #Check if the result is a combination of 1s and 0s
        #     print(alpha, beta, num_interim_zeros, result)

        if 0 in result and 1 in result:
            flag_mix = 1
            return None  # Stop checking for this alpha-beta pair

        # # Check if the result transitions from all 1s to all 0s
        if set(result) == {0}:
            if flag_1 == 1:
                # transition from 1s to 0s
                return num_interim_zeros
            else:
                flag_0 = 1
        if set(result) == {1}:
            if flag_0 == 1:
                # transition from 0s to 1s
                # Usually shouldnt happen hence setting it to a constant value.
                return -300
            else:
                flag_1 = 1

    # As there is no mixed or any transition, check which flag is active.
    if flag_0:
        # No Transition. All 0s
        return -200
    elif flag_1:
        # No Transition. All 1s
        return -300
    else:
        print('Some Error Happened.')
        return -500


##################################################################################################


def run_sim(sim_id, seed=42, alpha=0.8, beta=0.8, mu=0.5, rho=0.5, n_samples=10000, seq_length=100,
            verbose_level=0):
    """
    Runs the simulation with the given parameters.
    :param sim_id: A Unique ID(string) so that no 2 sims replaces the dirs created by the other unless specified.
    :param seed: Random seed for reproducibility.
    :param alpha: The probability that the next element node in the sequence is 1 given that the current node is 1.
    :param beta: The probability that the next element node in the sequence is 0 given that the current node is 0.
    :param mu: The portability that the B part of the mom is 1 if the A part of the mom is 0.
    :param rho: the probability that a mutation pushes the node to 0 if the A part of the mom is 0.
    :param n_samples: The number of samples used for training.
    :param seq_length: The length of the sequences.
    :param verbose_level: The level of verbosity during run.

    :return: The name of the simulation and The bit error rate after testing for the selected alpha and beta values.
    """
    # Set the random seed for reproducibility
    np.random.seed(seed)

    # Naming the simulation based on the sim_id, loss and alpha beta.
    # avoiding '.' in the name by using '_'
    sim_name = f'{sim_id}_A{str(alpha).replace(".", "_")}_B{str(beta).replace(".", "_")}'

    # Giving Introduction
    if verbose_level > 0:
        print("Epigenetic Sequence Predictor Viterbi_Decode Antagonistic")
        print(f"Sim Name :: {sim_name}")
        print()

    # Prepare the data
    data = Dataframe(alpha, beta, mu, rho, n_samples, seq_length)

    mother_sequences = data.mom_list
    corrupted_daughter_sequences = data.corrupt_daughter_list

    # Create the log Files
    if not os.path.exists(sim_name):
        os.makedirs(sim_name)

    for i, seq in enumerate(corrupted_daughter_sequences):
        data.corrected_daughter_list[i] = viterbi_decode_antagonistic(daughter_chromatin=seq, alpha=alpha, beta=beta,
                                                                      mu=mu, rho=rho)

    # Get conclusions
    data.conclusions(sim_name, level=verbose_level)
    # print(f'k={k}, E={np.mean(data.biterror)}')
    print(np.mean(data.biterror))

    with open(sim_name + '/log.txt', 'w') as f:
        f.write(f'Log File - Simulation :: {sim_name}\n\n')

        f.write(f'##INPUTS##\n')
        f.write(f'Random Seed for Reproducibility = {seed}\n\n')

        f.write(f'#Sample parameters#\n')
        f.write(f'alpha = {alpha}\n')
        f.write(f'beta = {beta}\n')
        f.write(f'mu = {mu}\n')
        f.write(f'rho = {rho}\n')
        f.write(f'seq_length = {seq_length}\n')
        f.write(f'Training n_samples = {n_samples}\n')

    return sim_name, np.mean(data.biterror)
