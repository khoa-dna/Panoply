import json
from os.path import basename, join
import pandas as pd
import numpy as np
import argparse
from panel_functions import *
from panel_classes import *
import time
import matplotlib.pyplot as plt

DATA_FIELDS = ["Abundance", "Brightness", "Target", "Spectra", "Database"]
PARAM_FIELDS = ['intensity_weight', 'corr_weight', 'random_prob',
                'risk_prob', 'frontier_size', 'branch_num',
                'intensity_eval_weight', 'corr_eval_weight', 'sample_size',
                'intensity_factor', 'auto_intensity', 'noise',
                'positive_fraction', 'result_num', 'unmix_num']
AFHULYM = "AfHuLym"


def arg_parser():
    """parse cmd args"""
    parser = argparse.ArgumentParser()
    parser.add_argument('config_json', type=str,
                        help="path to config json with \
                        params contains list of values want to simulate")
    parser.add_argument('output_folder', type=str,
                        help="output folder for \
                        generated panels and diagnostics plot")
    parser.add_argument('tag_name', type=str,
                        help="name unique to this run")
    parser.add_argument('--num-panel', type=int, default=100,
                        help="num panels to generate")

    args = parser.parse_args()
    return args


def load_input(path):
    """load csv or excel file input"""
    base_file_name = basename(path)
    extension = base_file_name.split(".")[-1]
    if extension == "csv":
        df = pd.read_csv(path)
    elif extension == "xlsx":
        df = pd.read_excel(path)
    else:
        raise Exception("Bad Extension: {}".format(extension))
    return df


def generate_panel_file(panels):
    """generate panel csv file from list of panel objects"""
    df = pd.DataFrame()
    count = 0
    for panel in panels:
        count += 1
        fluors = []
        markers = []
        bond_history = panel.bond_history
        for bond in bond_history:
            fluors.append(bond.fluor_end.name)
            markers.append(bond.marker_end.name)
        df.insert(loc=len(df.columns), column="Marker" +
                  str(count), value=markers)
        df.insert(loc=len(df.columns), column="Fluor"+str(count), value=fluors)
    # insert AFHULYM to all panel
    df.loc[len(df)] = [AFHULYM for i in range(len(df.columns))]
    return df


def calculate_cis(panels, auto_fluor):
    """calculate cis of panels"""
    cis = []
    for panel in panels:
        full_fluor_list = panel.fluor_history + [auto_fluor]
        cis.append(calculate_ci(full_fluor_list))
    return cis


def generate_ci_histogram(cis):
    """Diagnostics plots of panels: showing CIs"""
    plt.figure(figsize=(20, 10))
    plt.subplot(211)
    plt.title("Histogram of CIs")
    plt.ylabel("Counts")
    plt.xlabel("CI")
    # plt.xticks(np.arange(0, max(cis), int(len(cis)/10)).astype(int))
    plt.hist(cis, bins=5 * int(np.log(len(cis)) + 1))
    plt.subplot(212)
    plt.plot(sorted(cis), marker="o")
    plt.title("sorted_ci")
    plt.xlabel("panel (not index)")
    plt.ylabel("cis")
    plt.yticks(np.arange(0, max(cis), 10))
    plt.text(1, 1, "Min CI = {}".format(min(cis)))


def main():

    args = arg_parser()
    # Load data files
    with open(args.config_json, "r") as f:
        config = json.load(f)

    file_input_paths = {key: val for key,
                        val in config.items() if key in DATA_FIELDS}
    file_inputs = {key: load_input(path)
                   for key, path in file_input_paths.items()}
    params = {key: val for key, val in config.items() if key in PARAM_FIELDS}

    abundance_df, brightness_df, target_df, spectra_df, database_df =\
        file_inputs.values()

    database_df, marker_name_list, marker_abundance_list, \
        fluor_name_list, fluor_brightness_list, fluor_channel_values_list = \
        panel_data_preprocessing(
            abundance_df, brightness_df, target_df, spectra_df, database_df)

    # define autofluor
    afhulym = Fluor(name=AFHULYM, brightness=1,
                    channel_values=spectra_df[AFHULYM].values)

    # Construct objects and run beam search
    start = time.time()
    real_marker_strand = create_marker_strand(
        marker_name_list, marker_abundance_list)
    real_fluor_strand = create_fluor_strand(
        fluor_name_list, fluor_brightness_list, fluor_channel_values_list)
    real_dna = create_dna(real_fluor_strand,
                          real_marker_strand,
                          database_df["Fluor"].values,
                          database_df["Marker"].values,
                          marker_abundance_list,
                          fluor_brightness_list)
    results = constrained_beam_search(init_panel=real_dna,
                                      result_num=args.num_panel,
                                      frontier_size=params["frontier_size"],
                                      branch_num=params["branch_num"],
                                      random_prob=params["random_prob"],
                                      risk_prob=params["risk_prob"],
                                      auto_fluor=afhulym,
                                      alpha=params["intensity_weight"],
                                      beta=params["corr_weight"])

    end = time.time()
    print('{:.4f} s'.format(end-start))

    # parse panels generate panel df
    panels_flattened = [i for b in map(
        lambda x:[x] if not isinstance(x, list) else x, results) for i in b]
    panels = panels_flattened
    panel_file = generate_panel_file(panels)
    # output panels to csv
    panel_file.to_csv(
        join(args.output_folder, args.tag_name + "_panels.csv"), index=None)
    # generate ci plots
    cis = calculate_cis(panels, afhulym)
    generate_ci_histogram(cis)
    plt.savefig(join(args.output_folder, args.tag_name + "_qc.pdf"))
    # output cis to csv
    ci_df = pd.DataFrame(cis, columns=["CI"])
    ci_df["panel"] = [f"panel_{i+1}" for i in range(len(ci_df))]
    ci_df.to_csv(
        join(args.output_folder, args.tag_name + "_cis.csv"), index=None)


if __name__ == "__main__":
    main()
