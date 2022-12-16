from panel_classes import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import random
import time

from collections import defaultdict

from copy import deepcopy

from sklearn import linear_model
from sklearn.metrics import mean_squared_error

# Create marker strand: assume dataframe has Marker and Abundance columns


def create_marker_strand(marker_name_list, abundance_list) -> Marker_strand:
    marker_strand = Marker_strand()
    for i, marker_name in enumerate(marker_name_list):
        marker_strand.add_marker(
            Marker(name=marker_name, abundance=abundance_list[i]))
    marker_strand.generate_marker_dict()
    return marker_strand

# Create fluor strand: assume dataframe has Fluor and Brightness columns


def create_fluor_strand(fluor_name_list, brightness_list, channel_values) -> Fluor_strand:
    fluor_strand = Fluor_strand(ci=0)
    for i, fluor_name in enumerate(fluor_name_list):
        fluor_strand.add_fluor(Fluor(
            name=fluor_name, brightness=brightness_list[i], channel_values=channel_values[i]))
    fluor_strand.generate_fluor_dict()
    return fluor_strand

# Create all bonds from database to make DNA


def create_dna(fluor_strand, marker_strand, bond_fluor_name_list,
               bond_marker_name_list, abundance_list, brightness_list) -> DNA:
    fluor_keys = fluor_strand.fluor_strand_dict.keys()
    marker_keys = marker_strand.marker_strand_dict.keys()
    for i, val in enumerate(bond_fluor_name_list):
        if bond_fluor_name_list[i] in fluor_keys and bond_marker_name_list[i] in marker_keys:
            fluor_end = fluor_strand.fluor_strand_dict[bond_fluor_name_list[i]]
            marker_end = marker_strand.marker_strand_dict[bond_marker_name_list[i]]
            bond = Bond(fluor_end=fluor_end,
                        marker_end=marker_end, state="available")
    dna = DNA(fluor_strand=fluor_strand, marker_strand=marker_strand)
    # standardize bond strength
    strength_list = []
    for bonds in dna.all_bonds:
        for b in bonds:
            strength_list.append(b.strength)
    std_strength = (np.array(strength_list) -
                    np.mean(strength_list)) / np.std(strength_list)
    count = 0
    for bonds in dna.all_bonds:
        for b in bonds:
            b.strength = list(std_strength)[count]
            count += 1
    # add abundance/brightness to dna
    dna.abundance_list = abundance_list
    dna.brightness_list = brightness_list
    return dna


def calculate_ci(fluor_chain):
    values = []
    for f in fluor_chain:
        values.append(f.values)
    values = np.asarray(values).T
    index = max(np.linalg.svd(values)[1])/min(np.linalg.svd(values)[1])
    return index


def pearson_corr(a, b):
    return abs(np.corrcoef(np.array(a), np.array(b))[1, 0])


def find_min_marker(dna):
    m_strand = dna.marker_strand
    min_marker = m_strand.get_min_marker()
    return min_marker


def dna_found(dna):
    dna_completed = dna.bond_formed_num == len(dna.marker_strand.all_markers)
    if dna_completed:
        return True
    else:
        return False


def not_end(frontier):
    if len(frontier.get_panel_list()) == 0:
        return False
    for panel in frontier.get_panel_list():
        if dna_found(panel):
            return False
    return True


def choose_best(n, random_prob, panel_lists):
    newlist = sorted(panel_lists, key=lambda x: x.score)
    frontier_list = []
    for i in range(min(n, len(panel_lists))):
        if random.random() > random_prob:
            frontier_list.append(newlist.pop(0))
        else:
            frontier_list.append(newlist.pop(random.randrange(len(newlist))))
    return frontier_list, newlist

# form bond and commite nucleotides, and make other bonds unavilable


def form_bond(dna, strongest_bond):

    strongest_bond.state = "formed"
    strongest_bond.fluor_end.committed = True
    strongest_bond.marker_end.committed = True
    unavailable_bond_list = []
    for i in strongest_bond.fluor_end.all_bonds:
        if i.state == "available":
            i.state = "unavailable"
            unavailable_bond_list.append(i)
    for i in strongest_bond.marker_end.all_bonds:
        if i.state == "available":
            i.state = "unavailable"
            unavailable_bond_list.append(i)
    dna.bond_formed_num += 1
    dna.bond_history.append(strongest_bond)

    dna.marker_history.append(strongest_bond.marker_end)
    dna.fluor_history.append(strongest_bond.fluor_end)

    strongest_bond.fluor_end.final_bond = strongest_bond
    strongest_bond.marker_end.final_bond = strongest_bond

    return unavailable_bond_list


def unformed_bond(dna, bond, unavailable_bond_list):
    bond.state = "available"
    bond.fluor_end.committed = False
    bond.marker_end.committed = False
    for i in unavailable_bond_list:
        i.state = "available"
    dna.bond_formed_num -= 1
    dna.bond_history.pop()
    dna.marker_history.pop()
    dna.fluor_history.pop()


def calculate_composite_score(markers, fluors, new_bond, auto_fluor, **kwargs):
    ci = np.log(calculate_ci(fluors + [auto_fluor]))
    strength = new_bond.strength
    marker_abundances = []
    fluor_brightness = []
    for i in markers:
        marker_abundances.append(i.abundance)
    for i in fluors:
        fluor_brightness.append(i.brightness)
    corr = 0
    if len(markers) > 10:
        corr = pearson_corr(marker_abundances, fluor_brightness)

    # - 0.2*np.log(new_bond.fluor_end.brightness)
    return ci + kwargs["alpha"]*abs(strength) + kwargs["beta"]*corr


def choose_next_panel(panel, min_marker, branch_num, random_prob, auto_fluor, **kwargs):
    marker_list = panel.marker_history
    fluor_list = panel.fluor_history
    bond_list = panel.bond_history

    new_bonds = min_marker.get_available_bonds()
    scores = []
    for b in new_bonds:
        new_marker_list = marker_list + [b.marker_end]
        new_fluor_list = fluor_list + [b.fluor_end]

        score = calculate_composite_score(new_marker_list, new_fluor_list, b,
                                          auto_fluor, alpha=kwargs["alpha"], beta=kwargs["beta"])
        scores.append(score)

    table = tuple(zip(new_bonds, scores))
    sorted_table = list(sorted(table, key=lambda tup: tup[1]))

    best_list = []
    for i in range(min(branch_num, len(sorted_table))):
        if random.random() > random_prob:
            best_list.append(sorted_table.pop(0))
        else:
            best_list.append(sorted_table.pop(
                random.randrange(len(sorted_table))))

    children_panel_list = []
    for best in best_list:
        un_bond_list = form_bond(panel, best[0])
        children_panel = deepcopy(panel)
        children_panel.score = best[1]
        children_panel_list.append(children_panel)
        unformed_bond(panel, best[0], un_bond_list)
    return children_panel_list


def grow_panel(panel, branch_num, random_prob, auto_fluor, **kwargs):
    min_marker = find_min_marker(panel)
    if len(min_marker.get_available_bonds()) == 0:
        return None
    else:
        panels = choose_next_panel(panel, min_marker, branch_num, random_prob,
                                   auto_fluor,  alpha=kwargs["alpha"], beta=kwargs["beta"])
    return panels


def grow_frontier(frontier, frontier_size, branch_num, random_prob, risk_prob, auto_fluor, **kwargs):
    for panel in frontier.get_panel_list():
        new_panel_list = []
        new_panels = grow_panel(panel, branch_num, random_prob,
                                auto_fluor, alpha=kwargs["alpha"], beta=kwargs["beta"])
        if new_panels != None:
            for p in new_panels:
                if p != None:
                    new_panel_list.append(p)
        while (len(new_panel_list) < frontier_size) and len(frontier.reserve) > 0:
            p = frontier.reserve.pop()
            new_ps = grow_panel(p, branch_num, random_prob, auto_fluor,
                                alpha=kwargs["alpha"], beta=kwargs["beta"])
            for p in new_ps:
                new_panel_list.append(new_p)

    selected_panels, unselected_panels = choose_best(
        frontier_size, risk_prob, new_panel_list)
    new_frontier = Frontier(selected_panels)
    new_frontier.reserve = unselected_panels
    return new_frontier


def panel_data_preprocessing(abundance_df, brightness_df, target_df, spectra_df, database_df):
    abundance_df_copy, brightness_df_copy, target_df_copy, spectra_df_copy, database_df_copy = \
        deepcopy(abundance_df), deepcopy(brightness_df), deepcopy(
            target_df), deepcopy(spectra_df), deepcopy(database_df)
    database_df_copy = filter_database(
        database_df_copy, spectra_df_copy, target_df_copy)
    marker_name_list, marker_abundance_list = generate_marker_data(
        abundance_df_copy, target_df_copy)
    fluor_name_list, fluor_brightness_list, fluor_channel_values_list = generate_fluor_data(
        brightness_df_copy, database_df_copy, spectra_df_copy)
    return database_df_copy, \
        marker_name_list, marker_abundance_list, \
        fluor_name_list, fluor_brightness_list, fluor_channel_values_list


def filter_database(database_df, spectra, target_markers_df):
    # Filter database to fit in target panel and available spectra and brightness
    fluor_index_list = []
    not_included_fluors = []
    for i, val in enumerate(database_df["Fluor"].values):
        # and val in fluor_brightness_df["Fluor"].values
        if val in spectra.columns:
            fluor_index_list.append(i)
        else:
            not_included_fluors.append(val)
    print(set(not_included_fluors))
    marker_index_list = []
    for i, val in enumerate(database_df["Marker"].values):
        if val in target_markers_df["Marker"].values:
            marker_index_list.append(i)

    intersection_list = list(
        set(fluor_index_list).intersection(marker_index_list))
    database_df = database_df.iloc[intersection_list]
    database_df = database_df[["Fluor", "Marker"]]

    return database_df


def generate_marker_data(marker_abundance_df, target_markers_df):
    # marker_abundance_df["Abundance"] = np.log10(marker_abundance_df["Abundance"])
    # Impute missing value
    marker_abundance_df["Abundance"] = marker_abundance_df["Abundance"].\
        fillna(np.mean(marker_abundance_df["Abundance"]))
    marker_name_list = list(target_markers_df["Marker"].values)
    marker_abundance_list = []
    for i in marker_name_list:
        abundance = None
        if i not in marker_abundance_df["Marker"].values:
            abundance = np.mean(marker_abundance_df["Abundance"])
        else:
            abundance = marker_abundance_df[marker_abundance_df["Marker"]
                                            == i]["Abundance"].values[0]
        marker_abundance_list.append(abundance)
    return marker_name_list, marker_abundance_list


def generate_fluor_data(fluor_brightness_df, database_df, spectra):
    fluor_name_list = list(set(database_df["Fluor"].values))
    fluor_brightness_list = []
    fluor_channel_values_list = []
    # fluor_brightness_df["Brightness"] = np.log(fluor_brightness_df["Brightness"])
    fluor_brightness_df["Brightness"] = fluor_brightness_df["Brightness"].\
        fillna(np.mean(fluor_brightness_df["Brightness"]))
    for i in fluor_name_list:
        brightness = None
        if i not in fluor_brightness_df["Fluor"].values:
            brightness = np.mean(fluor_brightness_df["Brightness"])
        else:
            brightness = fluor_brightness_df[fluor_brightness_df["Fluor"]
                                             == i]["Brightness"].values[0]
        fluor_brightness_list.append(brightness)
        fluor_channel_values_list.append(spectra[i].values)
    return fluor_name_list, fluor_brightness_list, fluor_channel_values_list


def constrained_beam_search(init_panel=None, result_num=100,
                            frontier_size=3, branch_num=3, random_prob=0.2,
                            risk_prob=0.2, random_state=0, auto_fluor=None, **kwargs):
    result_list = []
    while len(result_list) < result_num:
        frontier = Frontier([init_panel])
        while not_end(frontier):
            frontier = grow_frontier(frontier, frontier_size, branch_num,
                                     random_prob, risk_prob, auto_fluor,
                                     alpha=kwargs["alpha"], beta=kwargs["beta"])
#             print(frontier)
#             clear_output(wait = True)
#             sleep(1)
        panel, l = choose_best(1, 0., frontier.get_panel_list())
        result_list.append(panel)
    return result_list


def calculate_panel_approx_score(markers, fluors, auto_fluor,
                                 intensity_eval_weight, corr_eval_weight):
    ci = np.log(calculate_ci(fluors + [auto_fluor]))
    marker_abundances = []
    fluor_brightness = []
    for i in markers:
        marker_abundances.append(i.abundance)
    for i in fluors:
        fluor_brightness.append(i.brightness)
    min_strength = min(np.log(marker_abundances)) * \
        min(np.log(fluor_brightness))
    corr = pearson_corr(marker_abundances, fluor_brightness)
    return ci + intensity_eval_weight*min_strength + corr_eval_weight*corr
