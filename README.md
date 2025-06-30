# WTA-states
Simple network architecture combining brain states and soft winner takes all

## Reproduce paper figure
To reproduce the figure from the paper, follow these steps:

1. Clone the repository [https://jugit.fz-juelich.de/w.wybo/ca-adex](https://jugit.fz-juelich.de/w.wybo/ca-adex), and follow the instructions provided therein to install the NESTML Ca-AdEx models
2. Run `python WTA_figure.py`. To select the two-compartment neuron model, in the `tune_general_config.yaml` file, set:
    ```yaml
    use_single_compartment_environment: False
    exc_neu_params_filename: "Ca-AdEx_neural_params.yaml"
    ```
    For the point neuron model, set
    ```yaml
    use_single_compartment_environment: True 
    exc_neu_params_filename: "standard_AdEx_exc_neural_params.yaml"
    ```

# Acknowledgemt of fundings
This work has been cofunded by the European Next Generation EU grants Italian MUR CUP I53C22001400006 (FAIR PE0000013 PNRR Project) and the Helmholtz Association's project-oriented funding programme (PoF 2, Topic 3).

# TRY ON EBRAINS

[![](https://nest-simulator.org/TryItOnEBRAINS.png)](https://lab.ebrains.eu/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2FAPE-group%2FWTA-states&urlpath=lab%2Ftree%2FWTA-states%2FWTA-awake-sleep.ipynb&branch=main)

# USER GUIDE (preliminary)
see file AA-UserGuide.md

# LOG OF MAIN EVENTS (from recent to past)
2025-06-28 add and tune nestml multicompartment models 

2024-1113 added a preliminary draft of User Guide
2024-1113-plastic-on-off release
- (better provenance tracking) solved issue #57 by copying neural parameters yaml in output dirs

2024-1004 to 2024-1018 solved the follwing issues:
- merged pull request #58 to solve enhancement issue 'copy neu params yaml in output dir' #57
- solved enhancement issue #55: added a directory neu_models_catalogue hosting interesting neural configurations
- solved bug #51 furter restringing the analysis window in _crop_and_plot.yaml
- removed unnecessarily tracke overwrite_output directory
- merged pull request #50 solving issue #45 add-actions-to-switch-on-plasticity-set-learning-parrameters-and-disconnect-synapses. 
- merged pull request #48 solving BUG issue #47: rates and spectral analysis for short analysis durations #47
- merged pull request #46 solving BUG issue #32: measure of firing rate when a pop has zero spikes in the time window.
- #43 specify plasticity default parameters in config file
- #39 Add first approximation of background recurrency (can set inter_assemblies connectivity to existing (True) or non existing (False))
and attibute initial values for the weights in ...general_config.yaml files
- #40 introduce flexible sparsity in synaptic matrices, added pairwise_bernolli connectivity, with programmable p_conn_exc and p_conn_inh values (defaulting to 1.0, i.e. all_to_all)

Release 2024-1002 solving issue #37  "Several issues related to Ca-AdEx calibration and plasticity management": addition of adaptation conductance, list of receptor names, learning rate settings, naming kinds of synapses and solving issue create dictionary book keeping the present kind of intra pop synapses #34

Release 2024-0912 solving issue #26 use recurrency false

Release 2024-0822 solving issue #19 switch on plasticity and #23: Added capability of disconnecting intra assemblies connections at programmable time, and what target exc pop to be disconnected

Release 2024-0820b config and plots can be saved in directories that can be specified using the file "directories_and_list_of_yamls.yaml".
Use there the "overwrite_dir" file name if you want to send last runs always to the same directory.

Release 2024-0820 Solution of issue #16. It is now possible to specify a list of events during the simulation. The first supported action is the saving of intra assembly synaptic matrices at times that can be defined in the basic_ and tune_sim_and_recording_times.yaml files

Release 2024-0817 It Solves both issue #11 and issue #12. Now the main launches a single activity analysis function. 
All parameters for the analysis and the plot of the activity can be set in the basic_ and tune_ crop_and_plot.yaml files.

Commit 2024-0817 Solving issue #12 create crop_and_plot yaml files

Release 2024-0816 Solving issue #9: "encapsulate the first part of the main in a nest_parameters_preparation() function"

Release 2024-0808 Solving issue #5: "Add basic contextual signal reaching a single assembly with a start and stop time"it is the product of branch: "#5-add-basic-contextual-signal..."

Release 2024-0807 introduced multicompartment support solving issue#3: "Support for runs with both single and multi-compartment neurons"

Release 2024-0510 Added capability of spectral analysis with error bands, by pull request #1 that solves issue #2: "Need for spectral analysis with error band and detection of maxima"

2024-0430: novel features:
- all parameters stored in reference yaml files
- all parameters can be tuned reading additional "tune_..." yaml files 
- nest network creation and simulation in separate file

2024-0428 Initial commit 
