#!/usr/bin/env python3
"""
The module performs enrichment analysis of gene functions per phase.
The analysis is performed with a hypergeometric test,
the multiple testing correction is done with the benjamini hochberg method.
"""

__version__ = "1.0"
__author__ = "Tin Siroki"

from scipy.stats import hypergeom
import statsmodels.stats.multitest as multi
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from pandas import DataFrame
import re
import configparser
import os
import argparse


def hypergeometric_over(quant, sample, hit, total):
    """
     Performs hypergeometric test for overexpression.

     :param quant: number of successes in sample
     :param sample: sample size
     :param hit: number of successes in population
     :param total: population size
     :return: over p-value,
     """
    over_rep = hypergeom.sf(quant - 1, total, hit, sample)
    return over_rep


def import_profiles(path):
    """
    Import profiles(expression values) from a given file.

    :param path: path to expression profiles
    :return: gene two profile dictionary
    """
    gene_2_profile = {}
    with open(path, "r") as o:
        first = True
        for line in o:
            if first:
                first = False
                continue
            line_splitted = line.split("\t")
            if len(line_splitted) > 2:
                id = line_splitted[0]
                points = [float(x.strip()) for x in line_splitted[1:]]
                gene_2_profile[id] = points
    return gene_2_profile


def import_mappings(path, col_names, alter_id = True):
    """
    Inputs mappings of gene id to alternative gene id. Input file should be a tab delimited text file contaning
    columns orig_id and alter_id.

    :param path: path to mappings
    :param col_names: tuple containing column names of: orig_id, orig_name, alter_id, alter_name (must be in given order)
    :param alter_id: if alter id is set to False the function returns mappings orig_2_orig
    :return: orig_2_alter dictionary containing mappings
    """
    orig_2_alter = {}
    if alter_id:
        orig_index = 0
        alter_index = 0
        with open(path, "r") as input:
            first = True
            for line in input:
                if first:
                    line_splitted = line.split("\t")
                    orig_index = line_splitted.index(col_names[0])
                    alter_index = line_splitted.index(col_names[2])
                    first = False
                    continue
                line_splitted = line.split("\t")
                if line_splitted[orig_index] != "" and line_splitted[alter_index] != "":
                    orig_2_alter[line_splitted[orig_index]] = line_splitted[alter_index]
    else:
        orig_index = 0
        with open(path, "r") as input:
            first = True
            for line in input:
                if first:
                    line_splitted = line.split("\t")
                    orig_index = line_splitted.index(col_names[0])
                    first = False
                    continue
                line_splitted = line.split("\t")
                if line_splitted[orig_index] != "":
                    orig_2_alter[line_splitted[orig_index]] = line_splitted[orig_index]
    return orig_2_alter


def import_gene_2_categories(path, depth):
    """
    Imports gene categories to dictionary.
    :param path: path to file containing categories. Format:
    category id;gene;locus tag;category1;category2;category3;category4;category5
    SW 1.1.1.1;BSU01770;glmM;Cellular processes;Cell envelope and cell division;Cell wall synthesis;Biosynthesis of peptidoglycan;
    SW 1.1.1.1;BSU14040;ldt;Cellular processes;Cell envelope and cell division;Cell wall synthesis;Biosynthesis of peptidoglycan;

    :param depth: the category
    :return: gene_id to category dictionary
    """
    gene_2_cat = {}
    with open(path, "r") as o:
        first = True
        for line in o:
            if first:
                first = False
                continue
            line_splitted = line.split(",")
            id = line_splitted[1]

            categories = re.findall(r'"(.*?)"', line)
            cat_num = len(categories)
            cat = ""
            if depth < cat_num:
                cat_num = depth

            for i in range(cat_num):
                if categories[i] != "":
                    cat += categories[i].strip("\"") + "#"

            cat = cat.strip("#")
            cat = cat.strip()
            cat = cat.strip("#")
            cat = cat.strip()
            cat = cat.strip("\"")

            if id not in gene_2_cat:
                gene_2_cat[id] = set()
            gene_2_cat[id].add(cat)
    return gene_2_cat


def is_expressed(profile, indexes, cutoff):
    """
    Checks if gene(profile) is expressed in a given phase. The phase can contain multiple points in the profile, given in
    the inedexes list. The gene is expressed if its value is grater than the cutoff value in any of the points from the
    phase.

    :param profile: expression profile
    :param indexes: indexes of observed phase
    :param cutoff: if value is greater than cutoff value the gene is expressed in the observed phase
    :return: true if gene is expressed
    """
    is_exp = False
    points = [profile[i] for i in indexes]
    for p in points:
        if p > cutoff:
            is_exp = True
            break

    return is_exp


def get_total_function_counts(gene_2_profile, orig_2_alter, gene_2_cat, folder_path):
    """
    Calculates total number of annotations (phase expression combination) per functional category.
    Stores total annotations to annotations_all.txt file.

    :param gene_2_profile: dictionary with gene profiles
    :param orig_2_alter: dictionary with mappings of original to alter id-s
    :param gene_2_cat: dictionary with categories for genes
    :param folder_path: out folder path
    :return: dictionary with counts per category (function)
    """
    cat_2_tot_count = {}
    total_count = 0
    no_annot_gene_count = 0

    total_annotation_out_path = os.path.join(folder_path, "annotations_all.txt")
    with open(total_annotation_out_path, "w") as out:
        for gene_id in gene_2_profile:
            if gene_id in orig_2_alter:
                gene_bsu_id = orig_2_alter[gene_id]
            else:
                print("No orig -> alter mapping for " + gene_id)
                continue

            if gene_bsu_id in gene_2_cat:
                categories = gene_2_cat[gene_bsu_id]

                for cat in categories:
                    if cat not in cat_2_tot_count:
                        cat_2_tot_count[cat] = 0

                    out.write(gene_id + "###" + cat + "\n")
                    cat_2_tot_count[cat] += 1
                    total_count += 1
            else:
                no_annot_gene_count += 1
                cat = "no_annotation"
                if cat not in cat_2_tot_count:
                    cat_2_tot_count[cat] = 0

                out.write(gene_id + "###" + cat + "\n")
                cat_2_tot_count[cat] += 1
                total_count += 1

    print("Number of genes with no annotations: " + str(no_annot_gene_count))
    check_sum = 0
    for cat in cat_2_tot_count:
        check_sum += cat_2_tot_count[cat]

    if check_sum != total_count:
        raise Exception("Total annotations count not valid.")

    return cat_2_tot_count, total_count


def store_table_format(data, functions, phases, path, p_adj=True, sig_p=0.05):
    """
    Stores significant p-values to a excel table.

    :param data:
    :param functions:
    :param phases:
    :param path:
    :param p_adj:
    :param sig_p:
    :return:
    """
    frame = {}
    frame["functions"] = functions
    i = 0
    for phase in phases:
        if p_adj:
            frame[phase] = [x if x < sig_p else "" for x in data[:, i]]
        else:
            frame[phase] = [x if x > 0 else "" for x in data[:, i]]
        i += 1

    df = DataFrame(frame)
    df.to_excel(path, sheet_name='functions', index=False)


def store_genes_with_significant_functions(phase_annotation_path, sig_table, info_path, col_names, out_path, alter_id):
    """
    Stores all genes with significant functions.

    :param phase_annotation_path:
    :param sig_table:
    :param info_path:
    :param col_names:
    :param out_path:
    :param alter_id:
    :return:
    """
    #read gene_names
    df_info = pd.read_csv(info_path, sep='\t', lineterminator='\n')
    df = pd.read_excel(sig_table)

    sig_fun_stage = set()

    for fun in df["functions"]:
        for stage in df.columns[1:]:
            if not pd.isnull(df.loc[df["functions"] == fun][stage].item()):
                sig_fun_stage.add(fun + "#" + stage)


    phase_2_function_2_genes = {}
    current_phase = ""
    with open(phase_annotation_path) as input:
        for line in input:
            line = line.strip()
            if line == "":
                continue
            if "###" not in line:
                current_phase = line
                phase_2_function_2_genes[line] = {}
                continue
            gene_id, function = line.split("###")

            if function + "#" + current_phase not in sig_fun_stage:
                continue

            if function not in phase_2_function_2_genes[current_phase]:
                phase_2_function_2_genes[current_phase][function] = []

            phase_2_function_2_genes[current_phase][function].append(gene_id)

    frame = {}
    frame["Stage"] = []
    frame["Function"] = []
    frame["GeneID"] = []
    orig_col_id = col_names[0]
    orig_col_name = col_names[1]
    frame[orig_col_name] = []

    orig_id_2_name = dict(zip(df_info[orig_col_id], df_info[orig_col_name]))
    if alter_id:
        alter_col_name = col_names[3]
        frame[alter_col_name] = []
        orig_id_2_alter_name = dict(zip(df_info[orig_col_id], df_info[alter_col_name]))

    for phase in phase_2_function_2_genes:
        for function in phase_2_function_2_genes[phase]:
            for gene in phase_2_function_2_genes[phase][function]:
                frame["Stage"].append(phase)
                frame["Function"].append(function)
                frame["GeneID"].append(gene)
                frame[orig_col_name].append(orig_id_2_name[gene])
                if alter_id:
                    frame[alter_col_name].append(orig_id_2_alter_name[gene])

    df_out = DataFrame(frame)
    df_out.to_excel(out_path, sheet_name="functions", index=False)


def plot_store_heatmap(path_in, path_out, p_adj=True, only_save=False, sig_p=0.05):
    """
    Stores significant results to excel table, plots heatmap.

    :param path_in: input path
    :param path_out: out path
    :param p_adj: set to True if p_adj is used
    :param only_save: if false plots heatmap.
    :return:
    """
    phases = []
    functions = set()
    current_phase = ""
    phase_2_function_2_res = {}

    with open(path_in) as input:
        for line in input:
            line = line.strip()
            if line == "":
                continue
            line_splitted = line.split("\t")
            if "quant\tsample\thit\ttotal" in line:
                phases.append(line_splitted[0])
                current_phase = line_splitted[0]
                phase_2_function_2_res[line_splitted[0]] = {}
                continue
            phase_2_function_2_res[current_phase][line_splitted[0]] = [float(line_splitted[5]), float(line_splitted[6]),
                                                                       float(line_splitted[7])]
            if float(line_splitted[6]) < sig_p:
                functions.add(line_splitted[0])

    data = np.zeros((len(functions), len(phases)))
    functions = list(functions)
    functions.sort()
    i = 0
    j = 0
    for function in functions:
        for phase in phases:
            if function in phase_2_function_2_res[phase]:
                if phase_2_function_2_res[phase][function][1] < sig_p:
                    if p_adj:
                        data[i][j] = phase_2_function_2_res[phase][function][1]
                    else:
                        data[i][j] = phase_2_function_2_res[phase][function][2]
                else:
                    if p_adj:
                        data[i][j] = 1
                    else:
                        data[i][j] = 0
            else:
                if p_adj:
                    data[i][j] = 1
                else:
                    data[i][j] = 0
            j += 1
        j = 0
        i += 1

    if not only_save:
        fig, ax = plt.subplots()
        im = ax.imshow(data, cmap='hot')
        ax.set_xticks(np.arange(len(phases)))
        ax.set_yticks(np.arange(len(functions)))
        ax.set_xticklabels(phases)
        ax.set_yticklabels(functions)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")
        ax.set_title("Function enrichment per phase")
        plt.colorbar(im)
        plt.show()

    store_table_format(data, functions, phases, path_out, p_adj=p_adj, sig_p=sig_p)


def run_experiment(depth, cutoff, out_folder, expression_path, categories_path, id_names_path, col_names, phase_2_index,
                   alter_id=True, only_save=True, total_count_all = True, sig_p=0.05):
    """
    Performs enrichment analysis og gene functions.
    The analysis is performed with a hypergeometric test, the multiple testing coreccbenjamini hochberg
    correction.

    :param depth: The depth of the functions(catgeroies) used
    :param cutoff: Cutoff value for expression
    :param out_folder: Path to output folder
    :param expression_path: Path to expression values
    :param categories_path: Path to mappings of genes to functions
    :param id_names_path: Path to names and id-s
    :param col_names: Names of columns used in the id_names_path file
    :param phase_2_index: Indexes of phases in the expression value file
    :param alter_id: Alternative id
    :param only_save: If true values are stored to out directory, heatmaps are not plotted
    :param total_count_all: If true the number of successes per function in population is calculated from all genes,
    else only from expressed
    :param sig_p: The p-value considered as significant
    :return:
    """

    gene_2_profile = import_profiles(expression_path)
    gene_2_cat = import_gene_2_categories(categories_path, depth)
    orig_2_alter = import_mappings(id_names_path, col_names, alter_id)

    """
    Get expressed genes by phases.
    """
    phase_2_genes = {}
    for phase in phase_2_index:
        phase_2_genes[phase] = []

    for gene in gene_2_profile:
        profile = gene_2_profile[gene]

        for phase in phase_2_index:
            if is_expressed(profile, phase_2_index[phase], cutoff):
                phase_2_genes[phase].append(gene)

    cluster_sample_count = []
    cat_2_tot_count = {}
    total_count = 0
    clus_2_res = {}
    all_cats = set()
    no_annot_gene_count = 0

    annotation_out_path = os.path.join(out_folder, "annotations_phase_cutoff_" + str(cutoff).replace(".", "_")
                                       + "_depth_" + str(depth) + ".txt")

    with open(annotation_out_path, "w") as ann_out:
        for phase in phase_2_genes:
            print("Analyzing phase " + phase)

            ann_out.write("\n")
            ann_out.write(phase + "\n")

            clus_2_count_sample = {}
            sample_count = 0

            for gene_id in phase_2_genes[phase]:
                if gene_id in orig_2_alter:
                    alter_gene_id = orig_2_alter[gene_id]
                else:
                    print("No orig -> alter mapping for " + gene_id)
                    continue

                if alter_gene_id in gene_2_cat:
                    categories = gene_2_cat[alter_gene_id]

                    for cat in categories:
                        if cat not in clus_2_count_sample:
                            clus_2_count_sample[cat] = 0
                        if cat not in cat_2_tot_count:
                            cat_2_tot_count[cat] = 0

                        clus_2_count_sample[cat] += 1
                        sample_count += 1

                        ann_out.write(gene_id + "###" + cat + "\n")
                        all_cats.add(cat)

                        cat_2_tot_count[cat] += 1
                        total_count += 1
                else:
                    no_annot_gene_count += 1
                    cat = "no_annotation"
                    if cat not in clus_2_count_sample:
                        clus_2_count_sample[cat] = 0
                    if cat not in cat_2_tot_count:
                        cat_2_tot_count[cat] = 0

                    clus_2_count_sample[cat] += 1
                    sample_count += 1

                    ann_out.write(gene_id + "###" + cat + "\n")
                    all_cats.add(cat)

                    cat_2_tot_count[cat] += 1
                    total_count += 1

            clus_2_res[phase] = clus_2_count_sample.copy()
            cluster_sample_count.append(sample_count)

    if total_count_all:
        cat_2_tot_count, total_count = get_total_function_counts(gene_2_profile, orig_2_alter, gene_2_cat, out_folder)

    print("Total number of different categories: " + str(len(all_cats)))
    p_values = []
    total_sample_count = sum(cluster_sample_count)

    if total_sample_count != total_count and not total_count_all:
        raise Exception("Total count and total sample count must be equal!!!")
    else:
        print("OK")

    cluster_counter = 0

    for c in clus_2_res:
        print("Calculating p-values sample: " + str(cluster_counter))

        for cat in clus_2_res[c]:
            pval = hypergeometric_over(clus_2_res[c][cat], cluster_sample_count[cluster_counter],
                                       cat_2_tot_count[cat], total_count)
            p_values.append(pval)

        cluster_counter += 1

    # benjamini hochberg correction
    p_adjusted = multi.multipletests(p_values, method="fdr_bh")[1]

    p_counter = 0
    cluster_counter = 0
    res_out_path = os.path.join(out_folder, "hyper_cutoff_" + str(cutoff).replace(".", "_") \
                   + "_depth_" + str(depth) + ".txt")
    res_out_path_filtered = os.path.join(out_folder, "hyper_filtered_cutoff_" + str(cutoff).replace(".", "_") \
                            + "_depth_" + str(depth) + ".txt")


    with open(res_out_path, "w") as out:
        with open(res_out_path_filtered, "w") as out_filter:

            for c in clus_2_res:
                print("Storing sample " + str(cluster_counter))
                out.write(c + "\tquant\tsample\thit\ttotal\tp_value\tp_adj\tlog_odds" "\n")
                out_filter.write(c + "\tquant\tsample\thit\ttotal\tp_value\tp_adj\tlog_odds" "\n")

                sample_count = cluster_sample_count[cluster_counter]
                temp_out = []

                for cat in clus_2_res[c]:
                    quant = clus_2_res[c][cat]
                    sample = sample_count
                    hit = cat_2_tot_count[cat]
                    total = total_count

                    out_line = cat + "\t" + str(quant) + "\t" + str(sample) + "\t" \
                               + str(hit) + "\t" + str(total) + "\t" + \
                               str(p_values[p_counter]) + "\t" + str(p_adjusted[p_counter])

                    if quant != 0 and (sample - quant) != 0 and (total - hit - sample + quant) != 0 and (hit - quant) != 0:
                        odds_sample = quant / (sample - quant)
                        odds_rest = (hit - quant) / (total - hit - sample + quant)
                        real_log_odds = math.log2(odds_sample / odds_rest)
                        out_line += "\t" + "%.2f" % real_log_odds
                        temp_out.append((p_adjusted[p_counter], out_line))
                    else:
                        out_line += "\tnan"
                        temp_out.append((p_adjusted[p_counter], out_line))

                    p_counter += 1
                temp_out.sort()
                cluster_counter += 1

                for t in temp_out:
                    out.write(t[1] + "\n")
                    if t[0] < sig_p:
                        out_filter.write(t[1] + "\n")

                out.write("\n")
                out_filter.write("\n")

    # set paths
    table_path_log = os.path.join(out_folder, "log_odds_cutoff_" + str(cutoff).replace(".", "_") + "_depth_" \
                     + str(depth) + ".xlsx" )
    table_p_values_path = os.path.join(out_folder, "p_adj_values_cutoff_" + str(cutoff).replace(".", "_") + "_depth_" \
                     + str(depth) + ".xlsx")
    gene_sig_path = os.path.join(out_folder, "genes_functions_cutoff_" + str(cutoff).replace(".", "_") + "_depth_" \
                       + str(depth) + ".xlsx")

    plot_store_heatmap(res_out_path, table_path_log, p_adj=False, only_save=only_save, sig_p=sig_p)
    plot_store_heatmap(res_out_path, table_p_values_path, p_adj=True, only_save=only_save, sig_p=sig_p)
    store_genes_with_significant_functions(annotation_out_path, table_p_values_path, id_names_path, col_names,
                                           gene_sig_path, alter_id)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="The module performs enrichment analysis of gene functions per phase. "
                                                 "The analysis is performed with a hypergeometric test, the multiple "
                                                 "testing correction is done with the benjamini hochberg method. "
                                                 "An annotation is a combination of a gene an its function. "
                                                 "The hypergeometric test for a funcion X in a phase P_1 is performed "
                                                 "on counts: "
                                                 " number of X function annotations in phase P_1 (quant), "
                                                 "number of annotations in phase P_1 (sample), "
                                                 "number of X function annotations in all phases (hit), "
                                                 "total number of annotations (total)."
                                                 "\n")

    parser.add_argument("config_path", help="The path to config file (see example file for instructions).")
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(args.config_path)

    # load paths
    out_folder = config.get("paths", "out_folder")
    expression_path = config.get("paths", "expression_path")
    id_names_path = config.get("paths", "id_names_path")
    categories_path = config.get("paths", "categories_path")

    # load flags
    alter_id = config.getboolean("flags", "alter_id")
    only_save = config.getboolean("flags", "only_save")
    total_count_all = config.getboolean("flags", "total_count_all")

    # load column names
    orig_id_col = config.get("column_names", "orig_id")
    orig_name_col = config.get("column_names", "orig_name")
    alter_id_col = None
    alter_name_col = None

    if alter_id:
        alter_id_col = config.get("column_names", "alter_id")
        alter_name_col = config.get("column_names", "alter_name")

    depths = []
    cutoffs = []

    # load significant p_value
    try:
        for c in config["cutoff_expression"]:
            sig_p = config.getfloat("p_value", "sig_p")
    except ValueError:
        print("Significant p-value must be parsable to float.")
        exit()

    # load cutoff values
    try:
        for c in config["cutoff_expression"]:
            cutoffs.append(config.getfloat("cutoff_expression", c))
    except ValueError:
        print("Cutoff expression value must be parsable to float.")
        exit()

    # load function depth values
    try:
        for c in config["category_depth"]:
            depths.append(config.getint("category_depth", c))
    except ValueError:
        print("Depth values must be parsable to int.")
        exit()

    # load phase indexes
    phase_2_index = {}
    try:
        for phase in config["phase_indexes"]:
            indexes = []
            indexes_split = config["phase_indexes"][phase].split()
            for index in indexes_split:
                indexes.append(int(index))
            phase_2_index[phase] = indexes.copy()
    except ValueError:
        print("Indexes of phases must be parsable to int.")
        exit()

    for cutoff in cutoffs:
        cutoff_folder_path = os.path.join(out_folder, "cutoff_" + str(cutoff).replace(".", "_"))
        if not os.path.exists(cutoff_folder_path):
            print("Creating results folder: " + cutoff_folder_path)
            os.makedirs(cutoff_folder_path)
        for depth in depths:
            cutoff_depth_folder_path = os.path.join(cutoff_folder_path, "depth_" + str(depth))
            if not os.path.exists(cutoff_depth_folder_path):
                print("Creating results folder: " + cutoff_depth_folder_path)
                os.makedirs(cutoff_depth_folder_path)

            print("\n\n*******-------*******\nRunning experiment cutoff: " + str(cutoff) + " depth: " + str(depth))
            run_experiment(depth, cutoff, cutoff_depth_folder_path, expression_path, categories_path, id_names_path,
                           (orig_id_col, orig_name_col, alter_id_col, alter_name_col), phase_2_index,
                           alter_id=alter_id, only_save=only_save, total_count_all=total_count_all, sig_p=sig_p)
