#!/usr/bin/env python3
"""
The module processes gene expression time series values with replicates.
Depending on the given arguments specific preprocess steps are included (default 1 and 2).
Preprocess steps:
1. Calculate frequencies per replicate.
2. Resolving replicates to one value, replicate median is calculated from non zero replicate values.
3. Zero filtering, genes with zero values for two or more points in time series are removed.
4. Interpolation of zero values with neighbours mean or with neighbour value if zero is located at the start or end of list (connected zero values are no allowed).
5. Median log in gene transformation, if multiple zeroes are present in gene values the median is calculated with only one zero value(preventing the median to be zero).
All values are divided by the gene median, and a log2 transformation is performed.
"""

__version__ = "1.0"
__author__ = "Tin Siroki"

import csv
import numpy as np
import argparse


def process_replicates(replicates):
    """
    Processes replicates, zero values are removed and replicate median is calculated.

    >>> process_replicates([100, 200, 300])
    200.0
    >>> process_replicates([100, 0, 200])
    150.0
    >>> process_replicates([0, 0, 200])
    200.0
    >>> process_replicates([1, 2, 3, 0 , 0, 402, 10])
    3.0

    :param replicates: list of replicates
    :return: non zero median
    """
    replicates = [0 if x == "NA" else x for x in replicates]
    replicates_f = [float(x) for x in replicates]
    replicates_filtered = [x for x in replicates_f if x != 0]
    res = 0

    if (len(replicates_filtered) > 0):
        res = np.median(replicates_filtered)

    return res


def interpolate_gen_expression(points):
    """
    Interpolates zero values, with neighbours mean or with neighbour value
    if zero is located at the start or end of list. Connected zero values are no allowed.

    >>> interpolate_gen_expression([0, 3, 4, 6])
    [3, 3, 4, 6]
    >>> interpolate_gen_expression([2, 4, 0, 6, 7, 8, 9])
    [2, 4, 5.0, 6, 7, 8, 9]

    :param points: list of point
    :return: interpolated list of points
    """

    points = points.copy()

    zero_count = 0

    for p in points:
        if p == 0:
            zero_count += 1

    if zero_count >= 2:
        raise Exception("In gene interpolation zero count must not be larger than 1.")

    point_count = len(points)

    if (points[point_count-1] == 0):
        points[point_count-1] = points[point_count-2]

    if (points[0] == 0):
        points[0] = points[1]

    for i in range(1, point_count-1):
        if points[i] == 0:
            points[i] = (points[i-1] + points[i+1]) / 2

    return points


def interpolate_gen_expression_one_two(points):
    """
    Checks if zerovalue in time series (points) is interpolated with one value.

    :param points: expression points
    :return: Returns true if value is interpolated with one value.
    """
    interpolated_one_value = True
    points = points.copy()
    zero_count = 0

    for p in points:
        if p == 0:
            zero_count += 1

    if zero_count >= 2:
        raise Exception("In gene interpolation zero count must not be larger than 1.")

    point_count = len(points)

    if points[point_count-1] == 0:
        points[point_count-1] = points[point_count-2]
        interpolated_one_value = True

    if points[0] == 0:
        points[0] = points[1]
        interpolated_one_value = True

    for i in range(1, point_count-1):
        if points[i] == 0:
            points[i] = (points[i-1] + points[i+1]) / 2
            interpolated_one_value = False

    return interpolated_one_value


def raw_2_mean_stats(replicates):
    """
    Calculates stats for replicates.

    :param replicates: list with replicate values
    :return: count of zero points, mean, standard deviation
    """
    zero_count = 0
    for r in replicates:
        if r == 0:
            zero_count += 1

    return zero_count, np.mean(replicates), np.std(replicates)


def has_two_or_more_zeroes(points):
    """
    Checks if two or more points are zero.

    >>> has_two_or_more_zeroes([3, 4, 5, 0])
    False
    >>> has_two_or_more_zeroes([3, 4, 5 ,7, 0, 0])
    True

    :param points: data points
    :return: True if two or more points are zero
    """
    zero_count = 0
    for point in points:
        if point == 0:
            zero_count += 1

    if zero_count < 2:
        return False
    else:
        return True


def median_log_transform(points, smoothing_factor=0.001):
    """
    Calculates pseudo median of points in given list. If multiple zero points are present in list only one zero point is
    considered in median calculation. After median transform zero points are replaced with the smoothing factor value to
    avoid performing log operation on zero values.
    In the end the points are log2 transformed.

    >>> median_log_transform([1, 2, 3])
    [-1.0, 0.0, 0.5849625007211562]

    >>> median_log_transform([0, 1, 1, 2, 3])
    [-9.965784284662087, 0.0, 0.0, 1.0, 1.584962500721156]

    >>> median_log_transform([0, 0, 0, 1, 1, 2, 3])
    [-9.965784284662087, -9.965784284662087, -9.965784284662087, 0.0, 0.0, 1.0, 1.584962500721156]

    :param points: data points
    :param smoothing_factor: if zero value is present it is replaced with the smoothing_factor
    :return: median log transformed points
    """
    median_points = []
    zero_count = 0

    for p in points:
        if p == 0:
            zero_count += 1
        else:
            median_points.append(p)

    if zero_count > 0:
        median_points.append(0)

    median = np.median(median_points)
    res = []
    for p in points:
        r = p / median
        if r == 0:
            r = smoothing_factor
        res.append(np.log2(r))

    return res


def get_point_replicate_indexes(head):
    """
    Returns indexes of replicates in each time point. Time point id should be separated from replicate id with '_' char.

    >>> get_point_replicate_indexes(["LB24H_1" , "LB24H_2", "LB24H_3", "H12_1", "H12_2", "H12_3", "D1_1", "D1_2", "D1_3"])
    {'LB24H': [0, 1, 2], 'H12': [3, 4, 5], 'D1': [6, 7, 8]}

    >>> get_point_replicate_indexes(["LB24H_34" , "LB24H_2", "LB24H_3", "H12_1", "D1_3", "LB24H_34", "D1_1", "LB24H_34", "D4_2"])
    {'LB24H': [0, 1, 2, 5, 7], 'H12': [3], 'D1': [4, 6], 'D4': [8]}

    :param head: list of header strings
    :return: id of replicates in time point
    """
    time_point_2_replicate_index = {}

    i = 0
    for h in head:
        point_replicate = h.split("_")
        tp = point_replicate[0]

        if tp in time_point_2_replicate_index:
            time_point_2_replicate_index[tp].append(i)
        else:
            time_point_2_replicate_index[tp] = []
            time_point_2_replicate_index[tp].append(i)
        i += 1

    return time_point_2_replicate_index


def calculate_frequencies(data_matrix):
    """
    Divides each column with the sum of its elements.

    >>> calculate_frequencies([[0, 1, 2], [1, 1, 2], [1, 1, 1]])
    array([[0.        , 0.33333333, 0.4       ],
           [0.5       , 0.33333333, 0.4       ],
           [0.5       , 0.33333333, 0.2       ]])

    :param data_matrix: data matrix
    :return: frequencies by column
    """

    data_matrix = np.array(data_matrix)

    for el in data_matrix.sum(axis=0):
        if el == 0:
            raise Exception("Column sum must not be zero, division by zero in frequencies calculation!!!")

    column_sums = data_matrix.sum(axis=0)
    return data_matrix / column_sums


def preprocess_time_series_values(input_path, output_path, transform=True, double_zero_remove=True, interpolate=True):
    """
    Preprocesses time series values. Input should be a path to the file containing counts for each gene,
    the counts should be normalized in respect to the gene length.
    Example input:
    locus_tag	LB24H_1	LB24H_2	LB24H_3	X6H_1	X6H_2	X6H_3	X12H_1	X12H_2	X12H_3	X1D_1
	B4U62_test0	5	7	9	0	5	4	0	0	3	8
	B4U62_test1	0	0	0	0	0	0	5	4	3	9
	B4U62_test2	1	0	1	0	0	0	0	4	7	2
	B4U62_test3	0	1	0	5	4	4	8	7	9	5

    All values should be tab delimited, the head contains time point names and replicate names delimited with '_'. For
    example, X6H_1 is time point X6H and replicate 1.

    Preprocess steps:
    1. Calculate frequencies per replicate.
    2. Resolving replicates to one value, replicate median is calculated from non zero replicate values.
    3. Zero filtering, genes with zero values for two or more points are removed.
    4. Interpolation. Interpolates zero values, with neighbours mean or with neighbour value
    if zero is located at the start or end of list. Connected zero values are no allowed.
    5. Median log transfrom. All points in gene expression profile are divided by the gene median value, after which
    the log2 transformation is performed\n
    6. Calculated values are stored to given output file.

    Example output:
    locus_tag	LB24H	X6H	X12H	X1D
	B4U62_test0	1.023083613	0.293731203	-1.658740427	-0.36923381
	B4U62_test2	-0.376618736	0.298375965	0.756406682	-1.054690642
	B4U62_test3	-1.432959407	0.567040593	0.467504919	-0.695993813

    :param input_path: path of input file
    :param output_path: path of output file
    :param transform: flag for performing log and pseudo median transformation
    :param double_zero_remove: flag for removing genes with double zero values
    :param interpolate: flag for interpolating single zero values with median of neighbour values
    :return: stores preprocessed values to given file
    """

    if not double_zero_remove and interpolate:
        raise Exception("Invalid combination of arguments (false double zero and true interpolation!!!")

    time_point_2_replicate_index = {}
    id_2_data_index = {}
    data_matrix = []
    id_2_data = {}

    head = ""
    ids_head = ""
    sample_count = 0

    first = True
    i = 0

    with open(input_path) as tsv:
        for row in csv.reader(tsv, dialect="excel-tab"):
            if first:
                first = False
                head = row[1:]
                ids_head = row[0]
                sample_count = len(head)
                continue
            if row[0] not in id_2_data_index:
                if len(row[1:]) != sample_count:
                    raise Exception("Invalid length of data points, each gene should have same number of data points!!!")

                id_2_data_index[row[0]] = i
                data_matrix.append([float(x) if x != "NA" else float(0) for x in row[1:]])
            else:
                print(row[0])
                raise Exception("Multiple gene id's are not allowed.")
            i += 1
    #1. frequencies calculation
    data_matrix = calculate_frequencies(data_matrix)
    #Gets index of replicates in each time point
    time_point_2_replicate_index = get_point_replicate_indexes(head)

    removed_count = 0
    interpolated_count = 0
    zero_gene_input_count = 0
    zero_after_replicate_resolving_count = 0
    interpolated_one_value = 0
    interpolate_two_values = 0

    #2. Replicate resolving, 3. Zero filttering, 4. Interpolation
    for id in id_2_data_index:
        points = []
        has_zero = False
        for tp in time_point_2_replicate_index:
            indexes = time_point_2_replicate_index[tp]

            if 0 in data_matrix[id_2_data_index[id]][indexes]:
                has_zero = True

            points.append(process_replicates(data_matrix[id_2_data_index[id]][indexes]))

        if 0 in points:
            zero_after_replicate_resolving_count += 1

        if has_zero:
            zero_gene_input_count += 1

        if double_zero_remove:
            if not has_two_or_more_zeroes(points):
                id_2_data[id] = {}
                id_2_data[id]["id"] = id

                if interpolate:
                    points_before = points.copy()
                    id_2_data[id]["point_freq"] = interpolate_gen_expression(points)

                    if points_before != id_2_data[id]["point_freq"]:
                        interpolated_count += 1

                        if interpolate_gen_expression_one_two(points):
                            interpolated_one_value += 1
                        else:
                            interpolate_two_values += 1
                else:
                    id_2_data[id]["point_freq"] = points
            else:
                removed_count += 1

        else:
            id_2_data[id] = {}
            id_2_data[id]["id"] = id
            id_2_data[id]["id"] = id

            if interpolate:
                id_2_data[id]["point_freq"] = interpolate_gen_expression(points)
            else:
                id_2_data[id]["point_freq"] = points

    print("Number of genes with one or more zero value replicates: " + str(zero_gene_input_count))
    print("Number of genes removed due to two or more zero points(after replicate resolving step): " + str(removed_count))
    print("Number of genes with one or more zero after replicate resolving: " + str(zero_after_replicate_resolving_count))
    print("Number of genes with interpolated zeroes: " + str(interpolated_count))
    print("Number of genes interpolated with one value: " + str(interpolated_one_value))
    print("Number of genes interpolated with two values: " + str(interpolate_two_values))

    #5. Median log transform
    if transform:
        for id in id_2_data:
            id_2_data[id]["log_med_point"] = median_log_transform(id_2_data[id]["point_freq"])
    else:
        for id in id_2_data:
            id_2_data[id]["log_med_point"] = id_2_data[id]["point_freq"]

    #6. Store values
    with open(output_path, "w") as o:
        head = ids_head + "\t"

        for tp in time_point_2_replicate_index:
            head += tp + "\t"

        head = head.strip()
        o.write(head + "\n")

        for id in id_2_data:
            line_out = id + "\t" + "\t".join([str(x) for x in id_2_data[id]["log_med_point"]]) + "\n"
            o.write(line_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocesses gene expression time series values with replicates.\n"
                                                 "Depending on the given arguments specific preprocess steps are included (default 1 and 2).\n"
                                                 "All preprocess steps:\n"
                                                 "1. Calculate frequencies per replicate.\n"
                                                 "2. Resolving replicates to one value, "
                                                 "replicate median is calculated from non zero replicate values.\n"
                                                 "3. Zero filtering, genes with zero values for two or more"
                                                 "points are removed.\n "
                                                 "4. Interpolation. Interpolates zero values, with neighbours mean or"
                                                 " with neighbour value if zero is located at the start or end of list."
                                                 "Connected zero values are no allowed.\n"
                                                 "5. Median log transfrom. All points in gene expression profile are "
                                                 "divided by the gene median value, after which the log2 transformation"
                                                 " is performed\n"
                                                 "6. Calculated values are stored to given output file.\n")

    parser.add_argument("in_path", help="in_path should be a path to the file containing expression values for each "
                                        "gene, the counts should be normalized in respect to the gene length.\n"
                                        "All values should be tab delimited, the head contains time point names and "
                                        "replicate names delimited with '_', e.g.,"
                                        "X6H_1 is time point X6H and replicate 1.\nExample input:\n"
                                        "locus_tag	LB24H_1	LB24H_2	LB24H_3	X6H_1	X6H_2	X6H_3	X12H_1	X12H_2	X12H_3	X1D_1\n"
                                        "B4U62_test0	5	7	9	0	5	4	0	0	3	8\n"
                                        "B4U62_test1	0	0	0	0	0	0	5	4	3	9\n"
                                        "B4U62_test2	1	0	1	0	0	0	0	4	7	2\n"
                                        "B4U62_test3	0	1	0	5	4	4	8	7	9	5\n")

    parser.add_argument("out_path", help="out_path should be a path where the output is stored.")

    parser.add_argument("--transform", help="if flag is present all points in gene expression profile are divided by "
                                            " the gene median value, after which the log2 transformation is performed\n"
                        , action="store_true")

    parser.add_argument("--double_zero_remove", help="if flag is present gene expression series with two or more zeros"
                                                     " are removed.", action="store_true")

    parser.add_argument("--interpolate", help="if flag is present zero values in gene expression series is "
                                              "interpolated with neighbours mean or with neighbour values if zero is "
                                              "located at the start or end of list. Connected zero values are no "
                                              "allowed. This step could be used for preparing profiles for plotting.",
                        action="store_true")
    args = parser.parse_args()
    preprocess_time_series_values(args.in_path, args.out_path, transform=args.transform,
                                  double_zero_remove=args.double_zero_remove, interpolate=args.interpolate)
