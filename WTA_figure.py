import nest
import yaml
import nest.raster_plot
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, butter, filtfilt
from analysis_functions import *
from asserting_functions import *
from yaml_io import *
from prepare_nest_parameters import *
from nest_reset_create_connect_simulate import *
from activity_analysis import *
from synaptic_io import *
from synaptic_analysis import *


from matplotlibsettings import *


pop_colours = [
    'firebrick',
    'goldenrod',
    'limegreen',
    'blueviolet'
]
pop_ls = [
    '-',
    '--',
    '-.',
    ':'
]
pop_markers = [
    "o",
    "D",
    "*",
    "x",
] 


# Plot of rastergram of excitatories
def exc_pops_rastegram_plot(ax_raster, cropped_events, nest_pms, crop_pms, plot_pms):
    num_exc_pop = nest_pms["network"]["num_exc_pop"]

    # Plot for excitatory neurons
    for pop in range(num_exc_pop):
        pop_events = cropped_events[pop]
        ax_raster.scatter(
            pop_events["times"],
            pop_events["senders"],
            s=5,
            c=pop_colours[pop % len(pop_colours)],
            label=f"Population {pop + 1}",
        )
    
    
    ax_raster.set_xlim(crop_pms["start_ms"]-50., crop_pms["stop_ms"]+100.)
    
    ylim = ax_raster.get_ylim()
    rect = mpatches.Rectangle((5000, ylim[0]), 
        100,
        ylim[1] - ylim[0],
        facecolor=pop_colours[0],
        edgecolor='none',
        alpha=.3,
        zorder=-100,
    )
    ax_raster.add_artist(rect)
    rect = mpatches.Rectangle((5100, ylim[0]), 
        100,
        ylim[1] - ylim[0],
        facecolor=pop_colours[1],
        edgecolor='none',
        alpha=.3,
        zorder=-101,
    )
    ax_raster.add_artist(rect)
    ax_raster.axvline(3000, c='lightgrey', lw=lwidth*1.5, zorder=-101)
    
    # ax_raster.set_title("excitatory spikes")
    ax_raster.set_xticklabels([])
    # ax_raster.set_xlabel(r"$t$ [ms]", fontsize=ticksize)
    ax_raster.set_ylabel("neuron ID", fontsize=ticksize)
    myLegend(ax_raster, loc="center left", add_frame=True, fontsize=ticksize)


def plot_traces(ax_traces, multimeter, cropped_events, crop_pms, nest_pms, xlim=(4950, 5250)):

    pop_0_events = cropped_events[0]
    event_senders = np.array(pop_0_events['senders'])
    event_times = np.array(pop_0_events['times'])
    idxs = np.where(
        np.logical_and(
            event_senders == 1,
            np.logical_and(
                event_times > xlim[0],
                event_times < xlim[1],
            )
        )
    )[0]
    # breakpoint()
    for idx in idxs:
        ax_traces.axvline(event_times[idx], c=colours[0], lw=lwidth*.7, ls='--')
    
    events = nest.GetStatus(multimeter, 'events')[0]
    ax_traces.plot(events["times"], events["v_comp0"], c=colours[0], lw=lwidth, label=r"$v_{soma}$")
    if not nest_pms["use_single_compartment_environment"]:
        ax_traces.plot(events["times"], events["v_comp1"], c=colours[1], lw=lwidth, label=r"$v_{dend}$")

    ax_traces.set_ylim((-100, 100))
    ylim = ax_traces.get_ylim()
    rect = mpatches.Rectangle((5000, ylim[0]), 
        100,
        ylim[1] - ylim[0],
        facecolor=pop_colours[0],
        edgecolor='none',
        alpha=.3,
        zorder=-100,
    )
    ax_traces.add_artist(rect)
    ax_traces.axvline(3000, c='lightgrey', lw=lwidth*1.5, zorder=-101)

    myLegend(ax_traces, loc="lower center", bbox_to_anchor=(0.5, 0.99), fontsize=ticksize, handlelength=.8, labelspacing=.2, handletextpad=.3)
    
    ax_traces.set_xlim(xlim)
    ax_traces.set_xlabel(r"$t$ [ms]", fontsize=ticksize)
    ax_traces.set_ylabel(r"$v_{m}$ [mV]", fontsize=ticksize)


def plot_weight_evolution(ax_weight, list_of_syn_matrix_file_names, crop_pms):
    """
    Analyzes and plots the temporal evolution of synaptic matrices.

    Parameters:
        list_of_syn_matrix_file_names (list): A list of file names containing synaptic matrices saved at different times.
    """
    time_points_ = []  # To store the times (in ms) when the synaptic matrices were saved
    synaptic_data = (
        {}
    )  # To store synaptic weights for each source-target pair over time

    # Iterate through each file to extract data
    for file_name in list_of_syn_matrix_file_names:
        if not os.path.exists(file_name):
            print(f"File {file_name} does not exist. Skipping.")
            continue

        # Extracting the synaptic information from the pickle file
        try:
            with open(file_name, "rb") as f:
                array_of_dicts = pickle.load(f)
        except Exception as e:
            print(f"Error reading pickle file {file_name}: {e}. Skipping.")
            continue

        # Iterate through the list of synaptic data dictionaries
        for pop_dict in array_of_dicts:
            time_ms = pop_dict.get("time (ms)", None)
            if time_ms is None:
                print(
                    f"Time information not found in file {file_name}. Skipping entry."
                )
                continue

            source_target = pop_dict.get(
                "conn_index", None
            )  # Assuming 'conn_index' identifies the source-target pair
            synaptic_weights = pop_dict["connections"][
                "weight"
            ].values  # Extract the weights as a numpy array

            # Add the time point to the list if it's not already added
            if time_ms not in time_points_:
                time_points_.append(time_ms)

            # Store synaptic weights for each source-target pair
            if source_target not in synaptic_data:
                synaptic_data[source_target] = {}
            if time_ms not in synaptic_data[source_target]:
                synaptic_data[source_target][time_ms] = []
            synaptic_data[source_target][time_ms].extend(synaptic_weights)

    # Sort time points and corresponding data for proper plotting
    time_points_ = sorted(time_points_)

    marker_index = 0
    epss = [
        0,-10,10,20
    ]

    for source_target, weight_dict in synaptic_data.items():
        mean_weights = []
        std_weights = []

        # Calculate mean and standard deviation of synaptic weights across time
        for time in time_points_:
            weights_at_time = weight_dict.get(time, [])
            if weights_at_time:
                mean_weights.append(np.mean(weights_at_time))
                std_weights.append(np.std(weights_at_time))
            else:
                mean_weights.append(np.nan)
                std_weights.append(np.nan)

        # make weight evolution plot a bit clearer
        if source_target == 0:
            time_points = [time_points_[0], 2999., 3001.] + time_points_[1:]
            # time_points[3] -= 100.
            time_points[4] -= 100.
            mean_weights = [mean_weights[0], mean_weights[0], 0.] + mean_weights[1:]
            std_weights = [std_weights[0], std_weights[0], 0.] + std_weights[1:]
        elif source_target == 1:
            time_points = time_points_
            time_points[1] += 100
            mean_weights[0] = 0.0
        elif source_target in [2,3]:
            time_points = [time_points_[0], time_points_[-1]]
            mean_weights = [mean_weights[0], mean_weights[-1]]
            std_weights = [std_weights[0], std_weights[-1]]
        else:
            time_points = time_points_
        
        # breakpoint()

        # Plot with error bars and different markers for each source-target pair
        marker = pop_markers[marker_index % len(pop_markers)]
        ls = pop_ls[marker_index % len(pop_markers)]
        colour = pop_colours[marker_index % len(pop_markers)]
        t_eps = epss[marker_index % len(pop_markers)]
        ax_weight.errorbar(
            np.array(time_points) + t_eps,
            mean_weights,
            yerr=std_weights,
            label=f"population {source_target}",
            capsize=5,
            marker=marker,
            markersize=markersize,
            markeredgewidth=lwidth*1.2,
            linestyle=ls,
            lw=lwidth,
            c=colour,
        )
        marker_index += 1

        if source_target in [0,1]:
            t0 = 5000 + source_target*100
            ylim = ax_weight.get_ylim()
            rect = mpatches.Rectangle((t0, ylim[0]), 
                100,
                ylim[1] - ylim[0],
                facecolor=pop_colours[source_target],
                edgecolor='none',
                alpha=.3,
                zorder=-100,
            )
            ax_weight.add_artist(rect)
            ax_weight.axvline(3000, c='lightgrey', lw=lwidth*1.5, zorder=-101)

    ax_weight.set_xlim(crop_pms["start_ms"]-50., crop_pms["stop_ms"]+100.)
    ax_weight.set_ylim(-.2, 11.)
    ax_weight.set_xlabel(r"$t$ [ms]", fontsize=ticksize)
    ax_weight.set_ylabel(r"$w_{exc}$ [nS]", fontsize=ticksize)
    ax_weight.legend(loc="center left", fontsize=ticksize)


def sim_and_plot_WTA_awake(axes=None):

    if axes is None:
        pl.figure("fig", figsize=(10,6))
        ax_traces = pl.subplot(131)
        ax_raster = pl.subplot(132)
        ax_weight = pl.subplot(133)
    else:
        ax_traces = axes[0]
        ax_raster = axes[1]
        ax_weight = axes[2]

    is_verbose = True
    # copy configuration yamls in specified output directory
    directories_and_list_of_yamls = read_basic_directories_and_list_of_yamls(is_verbose)
    copy_yamls_in_output_dir(directories_and_list_of_yamls, is_verbose)

    is_verbose = True
    # total sim, resolution and recording times
    times = read_sim_and_recording_times_yaml(is_verbose)

    # read general network parameters
    config = read_general_config_yaml(is_verbose)

    # copy neural params files in specified output directory
    output_dir_name = directories_and_list_of_yamls["relative_output_dir"]
    copy_neu_params_yamls_in_output_dir(output_dir_name, config, is_verbose)

    # prepare all simulation parameters
    nest_pms = {}
    nest_pms = nest_parameters_preparation(times, config, is_verbose, nest_pms)
    print("nest_pms", nest_pms)

    NEST_version = nest.__version__
    if (
        NEST_version == "3.7.0"
        and nest_pms["use_single_compartment_environment"] == False
    ):
        print(
            "ASSERTION ERROR: Ca-AdEx multi-compartment neuron not supported by this NEST version",
            NEST_version,
        )
        assert False

    is_verbose = False
    num_threads = 4
    (
        sim_completed,
        spike_recorders,
        inh_spike_recorder,
        multimeters,
        list_of_syn_matrix_file_names,
    ) = nest_reset_create_connect_simulate(nest_pms, num_threads, is_verbose)
    print("sim_completed", sim_completed)

    d_inh = nest.GetStatus(inh_spike_recorder, "events")[0]

    # before analysis, preliminary sim look
    is_verbose = True
    preliminary_sim_look(
        is_verbose,
        nest_pms,
        spike_recorders,
        inh_spike_recorder,
        nest_pms["recording_pms"],
    )

    is_verbose = False
    # here we prepare all the parameters for the following analysis and print
    crop_pms, plot_pms, sampling_pms, analysis_pms = (
        prepare_crop_plot_sampling_activityAnalysis_parameters(
            directories_and_list_of_yamls, nest_pms, is_verbose
        )
    )


    # Assuming spike_recorders is a list of spike recorder IDs previously created in your NEST simulation
    cropped_events = crop_events_from_spike_recorders(crop_pms, spike_recorders)
    cropped_inh_events = crop_inh_events(crop_pms, inh_spike_recorder)


    plot_traces(ax_traces, multimeters, cropped_events, crop_pms, nest_pms, xlim=(4700, 5500))

    is_verbose = False
    # launches all analysis
    # produces from both spike-like waveforms and tissue-like responses
    # PLEASE use the basic_ and tune_crop_and_plot.yaml to select plots and parameters
    # produce_rastegrams_rates_spectra_spectrograms(
    #     nest_pms,
    #     crop_pms,
    #     plot_pms,
    #     analysis_pms,
    #     cropped_events,
    #     cropped_inh_events,
    #     is_verbose,
    # )
    exc_pops_rastegram_plot(ax_raster, cropped_events, nest_pms, crop_pms, plot_pms)

    verbose = False
    # loading and printing of info from synaptic matrices
    # for file_name in list_of_syn_matrix_file_names:
    #     array_of_dicts = load_syn(file_name, verbose)

    plot_weight_evolution(ax_weight, list_of_syn_matrix_file_names, crop_pms)


def plot_figure():
    multi_comp = False

    pl.figure("WTA_network", figsize=(7,10))

    xcoo = getXCoords([0.5,2.0,0.2,0.8,.2])
    ycoo = getYCoords([0.3,3.0,0.2,1.6,0.2])

    plabels = getAnnotations(8)

    gs0 = GridSpec(2,1)
    gs0.update(top=ycoo[-1], bottom=ycoo[-2], left=xcoo[2], right=xcoo[3], hspace=.5)
    ax0 = myAx(pl.subplot(gs0[1,0]))

    gs1 = GridSpec(7,1)
    gs1.update(top=ycoo[-3], bottom=ycoo[-4], left=xcoo[0], right=xcoo[3], hspace=.4)
    ax1 = myAx(pl.subplot(gs1[0:4,0]))
    ax2 = myAx(pl.subplot(gs1[4:7,0]))

    if multi_comp:
        ax0.add_artist(plabels[0])
        ax1.add_artist(plabels[1])
        ax2.add_artist(plabels[2])
    else:
        ax0.add_artist(plabels[3])
        ax1.add_artist(plabels[6])
        ax2.add_artist(plabels[7])

    sim_and_plot_WTA_awake(axes=[ax0, ax1, ax2])

    if multi_comp:
        pl.savefig('NESTMLpaper_WTA_multi_comp.svg', transparent=True)
    else:
        pl.savefig('NESTMLpaper_WTA_single_comp.svg', transparent=True)


if __name__ == "__main__":
    # sim_and_plot_WTA_awake()
    plot_figure()
    pl.show()