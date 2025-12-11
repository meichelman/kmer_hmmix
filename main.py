import sys
from hmm_functions import TrainModel, write_HMM_to_file, read_HMM_parameters_from_file, Calculate_Posterior_probabillities, PMAP_path, Viterbi_path, Write_posterior_probs, Emission_probs_poisson
from helper_functions import load_obs_mut
from make_mutationrate import make_mutation_rate



def decode(obs, mutrates, param, out):
    hmm_parameters = read_HMM_parameters_from_file(param)
    obs, mutrates, weights = load_obs_mut(obs, mutrates)

    print('-' * 40)
    print(hmm_parameters)  
    print(f'> number of windows: {len(obs)}. Number of kmers = {sum(obs)}')
    print('> output is', out)
    print('> Window size is', 2000, 'bp')
    print('> Decode with posterior decoding')
    print('-' * 40)
    
    emissions = Emission_probs_poisson(hmm_parameters.emissions, obs, weights, mutrates)
    posterior_probs = Calculate_Posterior_probabillities(emissions, hmm_parameters)
    pmap_path = PMAP_path(posterior_probs)
    viterbi_path = Viterbi_path(emissions, hmm_parameters)
    
    Write_posterior_probs(obs, weights, mutrates, posterior_probs, pmap_path, viterbi_path, hmm_parameters, out)



def train(obs, mutrates, param, out):
    hmm_parameters = read_HMM_parameters_from_file(param)
    obs, mutrates, weights = load_obs_mut(obs, mutrates)

    print('-' * 40)
    print(hmm_parameters)
    print(f'> number of windows: {len(obs)}. Number of kmers = {sum(obs)}')
    print('> output is', out) 
    print('> window size is', 2000, 'bp') 
    print('-' * 40)

    hmm_parameters = TrainModel(obs, mutrates, weights, hmm_parameters)
    write_HMM_to_file(hmm_parameters, out)



def main():
    args = sys.argv
    if len(args) < 2:
        sys.exit('\n\nMust input mode\n\n')

    mode = args[1]
    modes = ['train', 'decode', 'mutrate']
    if mode not in modes:
        sys.exit('\n\nERROR! Mode must be either "train", "decode", or "mutrate"\n\n')
    
    if mode == 'train':
        # python3 main.py train none HG01891_h1/vidija/HG01891_h1_counts.bed HG01891_h1/vidija/mutrate_eq1MbBin/HG01891_h1_mutrates.txt HG01891_h1/vidija/mutrate_eq1MbBin/trained.json
        if len(args) < 5:
            sys.exit('\n\nMust input:\n\tInitial guess\n\tObservations\n\tMutation rate\n\tOutput\n\n')
        init_guess, observations, mutrates, out = args[2], args[3], args[4], args[5]
        
        if init_guess == 'none':
            init_guess = None
        weights = None
        
        train(observations, mutrates, init_guess, out)
    
    elif mode == 'decode':
        # python3 main.py decode HG01891_h1/vidija/mutrate_eq1MbBin/trained.json HG01891_h1/vidija/HG01891_h1_counts.bed HG01891_h1/vidija/mutrate_eq1MbBin/HG01891_h1_mutrates.txt HG01891_h1/vidija/mutrate_eq1MbBin/probs_and_path.tsv
        if len(args) < 5:
            sys.exit('\n\nMust input:\n\Parameters\n\tObservations\n\tMutation rate\n\tOutput\n\n')
        
        params, observations, mutrates, out = args[2], args[3], args[4], args[5]
        weights = None
        
        decode(observations, mutrates, params, out)
    
    elif mode == 'mutrate':
        # python3 main.py mutrate HG01891_h1/vidija/HG01891_h1_counts.bed HG01891_h1/vidija/mutrate_eq1MbBin/HG01891_h1_mutrates.txt

        if len(args) < 3:
            sys.exit('\n\nMust input:\n\tObservations\n\tOutput\n\n')
        observations, out = args[2], args[3]
        window_size = 1_000_000
        make_mutation_rate(observations, out, window_size)



if __name__ == "__main__":
    main()
