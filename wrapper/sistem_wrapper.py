import os
import argparse
import yaml
import pickle as pkl
from typing import List

from sistem import GrowthSimulator, load_gs
from sistem.parameters import Parameters, fill_params
from sistem.selection import RandomArmLibrary, FittedArmLibrary, RandomRegionLibrary, FittedRegionLibrary, RandomHybridLibrary, FittedHybridLibrary
from sistem.anatomy import SimpleAnatomy, StaticAnatomy, GenotypeAnatomy
from sistem.data import gen_readcounts_singlecell, gen_readcounts_bulk, gen_reads

SELECTION_MODELS = {
    "RandomArmLibrary": RandomArmLibrary,
    "FittedArmLibrary": FittedArmLibrary,
    "RandomRegionLibrary": RandomRegionLibrary,
    "FittedRegionLibrary": FittedRegionLibrary,
    "RandomHybridLibrary": RandomHybridLibrary,
    "FittedHybridLibrary": FittedHybridLibrary
}

MIGRATION_MODELS = {
    "SimpleAnatomy": SimpleAnatomy,
    "StaticAnatomy": StaticAnatomy,
    "GenotypeAnatomy": GenotypeAnatomy
}

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', default=None, type=str, help='Path to configuration file in YAML format (see provided config.yaml as an example).')
    parser.add_argument('-i', '--checkpoint', default=None, type=str, help='Path to checkpoint pickle file to resume simulation.')
    args = parser.parse_args()

    if args.config is None and args.checkpoint is None:
        raise RuntimeError('Either a configuration file or checkpoint file must be passed.')
    return args

def set_parameter_types(params):
    for k in Parameters.__annotations__:
        if k not in params:
            raise ValueError('Missing parameter values in provided configuration file.')

    nonekeys = [k for k,v in params.items() if v == 'None']
    for k in nonekeys:
        params[k] = None
    
    for k,v in Parameters.__annotations__.items():
        if v == float:
            params[k] = float(params[k])
        elif v == int:
            params[k] = int(float(params[k]))

    if isinstance(params['capacities'], List):
        params['capacities'] = [int(float(x)) for x in params['capacities']]
    else:
        params['capacities'] = int(float(params['capacities']))

    return params

def process_inputs_config(config_path):
    if not os.path.isfile(config_path):
        raise FileNotFoundError("The configuration file does not exist at the specified path.")
    with open(config_path, 'r') as f:
        inputs = yaml.safe_load(f)

    # Check top-level config file format
    if 'parameters' not in inputs:
        raise ValueError("The 'parameters' top-level key does not exist in the configuration file.")
    if 'models' not in inputs:
        raise ValueError("The 'models' top-level key does not exist in the configuration file.")
    
    # Create parameters object
    raw_params = set_parameter_types(inputs['parameters'])
    params = fill_params(None, **raw_params)

    # Check if simulation has already been started at the provided simulation directory
    if os.path.isdir(params.out_dir):
        if os.path.isfile(os.path.join(params.out_dir, 'checkpoint.pkl')):
            raise RuntimeError("A simulation checkpoint already exists in the output directory specified in the provided configuration file. To resume the simulation from this checkpoint, pass the path to the checkpoint.pkl file as the --checkpoint argument. Otherwise, to run a simulation from scratch in the same directory, delete the checkpoint.pkl file and run this script again.")

    # Otherwise start from new
    options = inputs['models']
    top_level_keys = ['Selection', 'Migration', 'singlecell', 'bulk', 'readcounts', 'singlecell_reads']
    for k in top_level_keys:
        if k not in options:
            raise ValueError("The 'models' object in the provided configuration file contains missing entries.")
    if options['Selection']['model'] not in SELECTION_MODELS:
        raise ValueError("The provided Selection model does not exist. Check for typos and correct capitalization.")
    if options['Migration']['model'] not in MIGRATION_MODELS:
        raise ValueError("The provided Migration model does not exist. Check for typos and correct capitalization.")
    if 'Fitted' in options['Selection']['model']:
        if not os.path.isfile(options['Selection']['delta_file']):
            raise ValueError("The provided delta_file does not exist at the specified path.")
    if options['Migration']['model'] == 'SimpleAnatomy' and options['Migration']['create_metastatic_libraries'] == True:
        raise ValueError("Cannot create metastatic libraries with the SimpleAnatomy model.")
    
    if options['singlecell'] == options['bulk']:
        raise ValueError("One of singlecell and bulk must be True and the other False.")
    
    if options['singlecell_reads'] == True:
        if params['ref'] is None:
            raise ValueError("Cannot generate synthetic sequencing reads without providing a reference genome.")
        if options['bulk'] == True:
            raise ValueError("singlecell_reads and bulk cannot both be True.")
    
    options['step'] = 0
    print('Initializing new simulation workflow from the provided config file.')
    with open(os.path.join(params.out_dir, 'checkpoint.pkl'), 'wb') as f:
        pkl.dump((options, params), f)
    return options, params
    
def process_inputs_checkpoint(checkpoint_path):
    if not os.path.isfile(checkpoint_path):
        raise FileNotFoundError("The checkpoint file does not exist at the specified path.")
    with open(checkpoint_path, 'rb') as f:
        options, params = pkl.load(f)

    step = options['step']
    
    print(f"Resuming simulation workflow from the provided checkpoint file at stage {step}.")
    return options, params

# Step 0
def initialize_sim(options, params):
    print("Beginning stage 0: Model Initialization")
    print("Creating selection library.")
    library = SELECTION_MODELS[options['Selection']['model']](params=params)
    if 'Fitted' in options['Selection']['model']:
        library.initialize(filepath=options['Selection']['delta_file'], params=params)
    else:
        library.initialize(params=params)
    
    print("Creating anatomical landscape.")
    anatomy = MIGRATION_MODELS[options['Migration']['model']](library, params=params)
    if options['Migration']['model'] != 'SimpleAnatomy':
        anatomy.initialize_distances()
        if options['Migration']['create_metastatic_libraries'] == True:
            anatomy.create_random_metastatic_libraries(method=options['Migration']['metastatic_library_method'], params=params)
    
    print("Preparing agents.")
    gs = GrowthSimulator(anatomy)
    gs.save_checkpoint(os.path.join(params.out_dir, 'gs0.pkl'))

    options['step'] = 1
    with open(os.path.join(params.out_dir, 'checkpoint1.pkl'), 'wb') as f:
        pkl.dump((options, params), f)

# Step 1
def run_growth_sim(options, params):
    print("Beginning stage 1: Agent growth")
    if not os.path.exists(os.path.join(params.out_dir, 'gs0.pkl')):
        raise ValueError("Error: GrowthSimulator object should exist from previous step.")
    gs = load_gs(os.path.join(params.out_dir, 'gs0.pkl'))
    if gs.gen > 0:
        raise ValueError("GrowthSimulator object appears to have already run. Please restart the simulation from the previous step.")
        
    print("You can monitor population counts in the log file located in the simulation directory.")
    print("WARNING: This step may take considerable time. For quicker simulations, consider using less anatomical sites, increasing the migration rate, or reducing the minimum number of detectable cells. Additionally, decreasing the carrying capacity of all sites will lead to faster simulation cycles.")
    gs.simulate_agents(params=params)

    print("Agent growth stage completed.")
    gs.save_checkpoint(os.path.join(params.out_dir, 'gs1.pkl'))

    options['step'] = 2
    with open(os.path.join(params.out_dir, 'checkpoint2.pkl'), 'wb') as f:
        pkl.dump((options, params), f)

# Step 2
def sim_lineage(options, params):
    print("Beginning stage 2: Lineage simulation")
    if not os.path.exists(os.path.join(params.out_dir, 'gs1.pkl')):
        raise ValueError("Error: GrowthSimulator object should exist from previous step.")
    gs = load_gs(os.path.join(params.out_dir, 'gs1.pkl'))
    print("Sampling cells.")
    gs.sample_cells(params=params)

    if options['singlecell'] == True:
        print("Creating single-cell lineage")
        tree = gs.simulate_singlecell_lineage(params=params)

    if options['bulk'] == True:
        print("Creating clonal lineage")
        tree = gs.simulate_clonal_lineage(params=params)
    
    with open(os.path.join(params.out_dir, 'tree.pkl'), 'wb') as f:
        pkl.dump(tree, f)
    
    gs.save_checkpoint(os.path.join(params.out_dir, 'gs2.pkl'))
    options['step'] = 3
    with open(os.path.join(params.out_dir, 'checkpoint3.pkl'), 'wb') as f:
        pkl.dump((options, params), f)

# Step 3
def gen_data(options, params):
    if options['readcounts'] is False and options['singlecell_reads'] is False:
        print("No readcount/read generation. Skipping step 3")
    
    else:
        print("Beginning stage 3: Readcount/sequencing read generation")
        if not os.path.exists(os.path.join(params.out_dir, 'gs2.pkl')):
            raise ValueError("Error: GrowthSimulator object should exist from previous step.")
        gs = load_gs(os.path.join(params.out_dir, 'gs2.pkl'))

        if not os.path.exists(os.path.join(params.out_dir, 'tree.pkl')):
            raise ValueError("Error: Tree object should exist from previous step.")
        with open(os.path.join(params.out_dir, 'tree.pkl'), 'rb') as f:
            tree = pkl.load(f)
        
        if options['singlecell_reads']:
            print("Generating synthetic single-cell DNA sequencing reads")
            gen_reads(tree, params=params)
        
        elif options['readcounts']:
            if options['singlecell']:
                print("Generating single-cell read counts")
                gen_readcounts_singlecell(tree, params=params)
            else:
                print("Generating bulk read counts")
                gen_readcounts_bulk(tree, params=params)
    
    options['step'] = 4
    with open(os.path.join(params.out_dir, 'checkpoint4.pkl'), 'wb') as f:
        pkl.dump((options, params), f)

def main():
    args = parse_args()
    if args.checkpoint is not None:
        options, params = process_inputs_checkpoint(args.checkpoint)
    else:
        options, params = process_inputs_config(args.config)
    
    while options['step'] < 4:
        if options['step'] == 0:
            initialize_sim(options, params)
        elif options['step'] == 1:
            run_growth_sim(options, params)
        elif options['step'] == 2:
            sim_lineage(options, params)
        elif options['step'] == 3:
            gen_data(options, params)
    
    print("Simulation workflow complete.")
    return
    
if __name__ == '__main__':
    main()