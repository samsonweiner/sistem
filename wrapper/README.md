# Utilizing the SISTEM wrapper script

The `sistem_wrapper.py` script provides a wrapper for running the complete SISTEM simulation workflow end-to-end. The wrapper depends on the PyYAML package which you can install as follows:
```
# Install with pip
pip install pyyaml

# Or install with conda
conda install pyyaml
```

To use the wrapper, modify or copy the *config.yaml* configuration file with the desired model selections and parameter values, then pass this file to the wrapper as follows:

```
python sistem_wrapper.py -c config.yaml
```

The wrapper will parse the simulation options and parameters from the configuration file, and then runs SISTEM over four stages:
* Stage 0: initialize selection & migration models and prepare agents
* Stage 1: Simulate the agents until minimum detectable population sizes are reached in all anatomical sites
* Stage 2: Simulate either a single-cell or clonal lineage tree
* Stage 3: Generate readcount data or synthetic single-cell DNA sequencing reads (optional)

The wrapper will log the status of each step and, after successful completion of any step, will save the results in a file named *checkpointN.pkl* within the specified simulation output directory, where 'N' corresponds to the checkpoint before the Nth stage. Should complications arise, or the user wishes to rerun certain steps, the checkpoint file can be passed to the wrapper to restart the simulation from that stage as follows:

```
python sistem_wrapper.py -i checkpointN.pkl
```

Note that the SISTEM wrapper does not offer the full range of customization available through the SISTEM API, but rather serves as an easy starting point for running end-to-end simulations. For complete details on models and parameters, see the SISTEM documentation.