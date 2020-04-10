#!/usr/bin/env python3
"""
The phyl_tree.py module implements the phylostratigraphy method. The method assigning a phylostratum to each focal
species gene by mapping it to a phylogenetic tree. The module imports a phylogenetic tree from nodes and names file and
maps the genes to the tree according to the given blast results.
"""

__version__ = "1.0"
__author__ = "Tin Siroki"

import sys
from pandas import DataFrame
import argparse

visited = set()


class Node:
    """
    Node class represent s phylogenetic tree node.
    """
    def __init__(self, node_id, name="unknown", is_root=False, is_leaf=False, is_phyl_node=False):
        self.id = node_id
        self.name = name
        self.children = []
        self.is_root = is_root
        self.is_leaf = is_leaf
        self.is_phyl_node = is_phyl_node
        self.parent = None
        self.phylostratum = None


class Tree:
    """
    Tree class represents a phylogenetic tree.
    """
    def __init__(self, root):
        self.root = root
        self.phyl_nodes = []
        if root.is_phyl_node:
            self.phyl_nodes.append(root)
        self.gene_2_hits = {}

    def __init__(self, path_nodes, path_names, focal_id):
        """
        Constructs a phylogenetic tree form nodes and names text files.
        The files should not contain headers and node ids should be parsable to int.

        Names file format:
        node_1_id name_1
        node_2_id name_2
        ...

        The root-id must be the parent in the first line of the input file.
        Nodes file format:
        child_id_1 parent_id_1
        ...

        :param path_nodes:
        :param path_names:
        :param focal_id:
        """

        print("Generating tree.")
        self.phyl_nodes = []
        self.focal_id = focal_id
        self.gene_2_hits = {}
        self.ph_node_id_2_PS = {}
        self.ph_node_id_2_species = {}
        self.gene_2_min_PS = {}
        self.node_id_2_PS = {}
        self.rep_2_members = {}

        node_2_name = {}
        with open(path_names) as in_file:
            for line in in_file:
                lines_splitted = line.split()
                if line.strip() == "":
                    continue
                if len(lines_splitted) < 2:
                    raise Exception("Invalid names file format. Line: " + line)
                try:
                    node_2_name[int(lines_splitted[0])] = lines_splitted[1].strip()
                except ValueError:
                    print("Node id-s must be parsable to integer. Line: " + line)

        parent_2_children = {}
        root = True
        root_id = ""

        with open(path_nodes) as input:
            for line in input:
                if line.strip() == "":
                    continue
                line_splitted = line.split()
                if len(lines_splitted) < 2:
                    raise Exception("Invalid nodes file format. Line: " + line)
                try:
                    child_id = int(line_splitted[0])
                    parent_id = int(line_splitted[1])
                except ValueError:
                    print("Node id-s must be parsable to integer. Line: " + line)
                if root:
                    root_id = parent_id
                    root = False

                if parent_id not in parent_2_children:
                    parent_2_children[parent_id] = []

                parent_2_children[parent_id].append(child_id)

        name = parent_id
        if parent_id in node_2_name:
            name = node_2_name[parent_id]
        else:
            print("ID: " + int(parent_id) + " is not in names file. Names set to id value.")

        node = Node(root_id, is_root=True, name=name)
        self.root = node
        self.__add_nodes_from_child_2_parent_dict(node, parent_2_children, node_2_name)

        # set phylostratum nodes to all nodes on the path from focal species to root
        current_node = self.get_node(focal_id)
        while not current_node.is_root:
            self.phyl_nodes.append(current_node)
            current_node.is_phyl_node = True
            current_node = current_node.parent

        # set phylostratum numbers
        PS = 1
        current_node = self.root
        while True:
            if current_node.is_phyl_node:
                self.ph_node_id_2_PS[current_node.id] = PS
                current_node.phylostratum = PS
                PS += 1
            phyl_node = None
            for c in current_node.children:
                if c.is_phyl_node:
                    if phyl_node:
                        raise ValueError("Invalid format of phylogeny. Each phyl node must contain one and only one"
                                         "phyl child.")
                    phyl_node = c
            # breaks loop when last phyl node is reached
            if not phyl_node:
                break
            current_node = phyl_node

        # set species from side branch to each phyl node
        self.ph_node_id_2_species = self.get_side_phyl_branch_leafs()

        # set species 2 PS
        for ph_node in self.ph_node_id_2_species:
            for spec in self.ph_node_id_2_species[ph_node]:
                PS = self.ph_node_id_2_PS[ph_node]
                spec.phylostratum = PS
                self.node_id_2_PS[spec.id] = PS

    def export_min_PS_map(self, path, ps_merge_path=None):
        """
        Exports genes phylostratum (ps) to given path.
        If a ps_merge_path genes with ps in given file are mapped to new ps.
        File format (tab delimited text file):
        orig_ps	new_ps
        5	4
        6	4
        7	5

        :param path: output map path
        :param ps_merge_path: path to merge file (optional)
        :return:
        """
        print("Exporting map.")
        old_ps_2_new = {}
        if ps_merge_path:
            first = True
            with open(ps_merge_path) as input:
                for line in input:
                    if first:
                        first = False
                        continue
                    l_s = line.split()
                    try:
                        old_ps = int(l_s[0])
                        new_ps = int(l_s[1])
                    except ValueError:
                        print("Phylostratum number must be parsable to integer. Line: " + line)
                    old_ps_2_new[old_ps] = new_ps

        with open(path, "w") as out:
            out.write("protein_id\tPS\n")
            for gene in self.gene_2_min_PS:
                gene_ps = self.gene_2_min_PS[gene]
                if ps_merge_path:
                    if gene_ps in old_ps_2_new:
                        gene_ps = old_ps_2_new[gene_ps]
                out.write(gene + "\t" + str(gene_ps) + "\n")

    def hgt_detection_analysis(self, path, hgt_gap=2, export_hgt_correction=False,
                                       zero_filter=None, ps_merge_path=None, additional_node_counts=None):
        """
        Performs horizontal gene transfer analysis. The genes with gaps in hits on the phylogenetic tree are candidates
        for horizontal gene transfer. This genes could be removed or corrected. The function searches for genes with gap
        grater or equal to hgt_gap. The function stores counts and mappings to given output folder.
        Result files:
        ps_species_count.xlsx               -number of species per phylostratum (ps)
        gene_hit_species_count.xlsx         -number of species hits per gene per phylostratum (here gaps can be seen)
        gene_hit_species_percent.xlsx       -percent of species hits (number_species_hits / total_number_species_ps)
        gene_hit_species_stats.xlsx         -percent hits per ps, PS (original PS), hgt assigned PS,
                                            hgt assigned PS iterative gap remove, total species hit perc (from original PS to focal species)
                                            -zero gap size
                                            -additional node counts (if given)

        map_filter_zero_filter.txt          -map with removed genes with gap grater or equal to zero filter (if zero_filter not None)
        map_hgt_correction_hgt_gap.txt      -map with genes PS corrected for genes with gaps grater or equal to hgt_gap
                                            (if export_hgt_correction)
        map_hgt_correction_hgt_gap.txt      -map with genes PS iteratively corrected for genes with gaps grater or equal
                                            to hgt_gap (if export_hgt_correction)


        :param path: path to output folder
        :param hgt_gap: gap used in horizontal gene transfer analysis and correction
        :param export_hgt_correction: if true map with corrected phylostrata is exported
        :param zero_filter: if number genes with gap grater than or equal to zero_filter are filtered from map
        :param ps_merge_path:  path to merge file (optional)
        :param additional_node_counts: list of tuples of node id and node name for which hits are counted
        :return:
        """
        old_ps_2_new = {}
        if ps_merge_path:
            first = True
            with open(ps_merge_path) as input:
                for line in input:
                    if first:
                        first = False
                        continue
                    l_s = line.split()
                    try:
                        old_ps = int(l_s[0])
                        new_ps = int(l_s[1])
                    except ValueError:
                        print("Phylostratum number must be parsable to integer. Line: " + line)
                    old_ps_2_new[old_ps] = new_ps

        # Extracts phylostrata id-s
        phylostrat = set()
        for ph_node_id in self.ph_node_id_2_PS:
            node_ps = self.ph_node_id_2_PS[ph_node_id]
            if ps_merge_path and node_ps in old_ps_2_new:
                phylostrat.add(old_ps_2_new[node_ps])
            else:
                phylostrat.add(node_ps)

        frame = {}
        frame["prot_id"] = []
        for ps in phylostrat:
            frame[ps] = []

        frame_species_count = {}
        for ps in phylostrat:
            frame_species_count[ps] = []

        spec_count = [0] * len(phylostrat)
        for node_id in self.ph_node_id_2_PS:
            ps = self.ph_node_id_2_PS[node_id]
            if ps_merge_path and ps in old_ps_2_new:
                ps = old_ps_2_new[ps]
            spec_count[ps - 1] += len(self.ph_node_id_2_species[node_id])

        for i in range(len(phylostrat)):
            frame_species_count[i + 1].append(spec_count[i])

        path_counts = path + "ps_species_count.xlsx"
        df = DataFrame(frame_species_count)
        df.to_excel(path_counts, index=False)

        for gene in self.gene_2_hits:
            frame["prot_id"].append(gene)

            species = set()
            for hit in self.gene_2_hits[gene]:
                species.add(hit[0])

            ps_hit = [0] * len(phylostrat)
            for spec in species:
                ps = self.node_id_2_PS[spec]
                if ps_merge_path and ps in old_ps_2_new:
                    ps = old_ps_2_new[ps]
                ps_hit[ps - 1] += 1

            for i in range(len(phylostrat)):
                frame[i+1].append(ps_hit[i])

        path_counts = path + "gene_hit_species_count.xlsx"
        df = DataFrame(frame)
        df.to_excel(path_counts, index=False)

        # species percentage
        frame_perc = {}
        frame_perc["prot_id"] = []
        for ps in phylostrat:
            frame_perc[ps] = []

        prot_counter = 0
        for prot_id in frame["prot_id"]:
            frame_perc["prot_id"].append(prot_id)
            for i in range(len(phylostrat)):
                if frame_species_count[i + 1][0] == 0:
                    frame_perc[i + 1].append(float("nan"))
                    continue
                frame_perc[i + 1].append(frame[i + 1][prot_counter] / frame_species_count[i + 1][0])
            prot_counter += 1

        path_perc = path + "gene_hit_species_percent.xlsx"
        df = DataFrame(frame_perc)
        df.to_excel(path_perc, index=False)


        # species percentage
        frame_stats = {}
        frame_stats["prot_id"] = []
        for ps in phylostrat:
            frame_stats[ps] = []
        frame_stats["PS"] = []
        frame_stats["hgt_PS"] = []
        frame_stats["hgt_iter_PS"] = []
        frame_stats["spec_hit_perc"] = []
        frame_stats["zero_gap_size"] = []

        if additional_node_counts:
            for count_id in additional_node_counts:
                frame_stats[count_id[1]] = []

        prot_counter = 0
        print("Exporting hits, protein count:")

        # store species for hit calculation for additional nodes
        additional_counts_2_species = {}
        if additional_node_counts:
            for add_count in additional_node_counts:
                additional_counts_2_species[add_count[1]] = set([x.id for x in self.get_leafs_subtree(self.get_node(add_count[0]))])


        for prot_id in frame["prot_id"]:
            if prot_counter % 100 == 0:
                print(str(prot_counter))
            frame_stats["prot_id"].append(prot_id)
            perc = []
            spec_hit_counts = 0
            all_spec_count = 0
            gene_ps = self.gene_2_min_PS[prot_id]

            if ps_merge_path:
                if gene_ps in old_ps_2_new:
                    gene_ps = old_ps_2_new[gene_ps]

            for i in range(len(phylostrat)):
                if frame_species_count[i + 1][0] == 0:
                    frame_stats[i + 1].append(float("nan"))
                    continue
                d = frame[i + 1][prot_counter] / frame_species_count[i + 1][0]
                frame_stats[i + 1].append(d)
                perc.append(d)

                if (i + 1) >= gene_ps:
                    all_spec_count += frame_species_count[i + 1][0]
                    spec_hit_counts += frame[i + 1][prot_counter]

            frame_stats["PS"].append(gene_ps)
            frame_stats["hgt_PS"].append(self.__resolve_horizontal_transfer(perc, min_gap_size=hgt_gap))
            frame_stats["hgt_iter_PS"].append(self.__resolve_horizontal_transfer_iter(perc, min_gap_size=hgt_gap))
            frame_stats["spec_hit_perc"].append(spec_hit_counts / all_spec_count)
            frame_stats["zero_gap_size"].append(self.__get_zero_gap_size(perc))

            if additional_node_counts:
                for count_id in additional_node_counts:
                    hit_count = 0
                    species = set()
                    for hit in self.gene_2_hits[prot_id]:
                        species.add(hit[0])
                    for spec in species:
                        if spec in additional_counts_2_species[count_id[1]]:
                            hit_count += 1
                    frame_stats[count_id[1]].append(hit_count)

            prot_counter += 1

        path_perc = path + "gene_hit_species_stats.xlsx"
        df = DataFrame(frame_stats)
        df.to_excel(path_perc, index=False)

        if zero_filter:
            map_path = path + "map_filter_" + str(zero_filter) + ".txt"
            print("Exporting map.")
            old_ps_2_new = {}
            prot_counter = 0
            with open(map_path, "w") as out:
                out.write("prot_id\tPS\n")

                for prot_id in frame["prot_id"]:
                    perc = []
                    for i in range(len(phylostrat)):
                        if frame_species_count[i + 1][0] == 0:
                            perc.append(float("nan"))
                            continue
                        d = frame[i + 1][prot_counter] / frame_species_count[i + 1][0]
                        perc.append(d)

                    prot_counter += 1
                    gene_ps = self.gene_2_min_PS[prot_id]
                    if ps_merge_path:
                        if gene_ps in old_ps_2_new:
                            gene_ps = old_ps_2_new[gene_ps]

                    if self.__get_zero_gap_size(perc) >= zero_filter:
                        continue

                    out.write(prot_id + "\t" + str(gene_ps) + "\n")

        # export map with phylostrata corrected for horizontal gene transfer
        if export_hgt_correction:
            map_path = path + "map_hgt_correction_" + str(hgt_gap) + ".txt"
            frame_hgt = {}
            frame_hgt["prot_id"] = frame_stats["prot_id"]
            frame_hgt["PS"] = frame_stats["hgt_PS"]
            df_hgt = DataFrame(frame_hgt)
            df_hgt.to_csv(map_path, index=False, sep="\t")

            map_path = path + "map_hgt_correction_iter_" + str(hgt_gap) + ".txt"
            frame_hgt = {}
            frame_hgt["prot_id"] = frame_stats["prot_id"]
            frame_hgt["PS"] = frame_stats["hgt_iter_PS"]
            df_hgt = DataFrame(frame_hgt)
            df_hgt.to_csv(map_path, index=False, sep="\t")

    def get_node(self, node_id):
        """
        Returns node with given id.

        :param node_id: node id
        :return:
        """
        return self.__get_node_by_id(self.root, node_id)

    def add_node(self, parent_id, element):
        """
        Adds node to tree.

        :param parent_id: id of nodes parent.
        :param element: node to add
        :return:
        """
        node = self.__get_node_by_id(self.root, parent_id)

        if node == None:
            print("Parent not found, ID: " + str(parent_id))
            return

        element.parent = node
        node.children.append(element)

        if element.is_phyl_node:
            self.phyl_nodes.append(element)
        return

    def remove_leaf_node(self, node_id):
        """
        Removes a leaf node from the tree.

        :param node_id: id of node to be removed
        :return:
        """
        node = self.__get_node_by_id(self.root, node_id)

        if node == None:
            raise Exception("Node id not found: " + str(node_id))
        node.parent.children.remove(node)

    def print_tree(self, path=None):
        """
        Prints tree in a readable format.

        :param path: If path is given the tree is stored to given file, else the output is printed to stdout.
        :return:
        """
        global visited
        visited = set()
        left_string = " " * 60 + "|"
        if path:
            with open(path, "w") as out:
                self.__traverse_tree_print_file(self.root, 0, left_string, out)
        else:
            self.__traverse_tree_print(self.root, 0, left_string)


    def get_leafs_subtree(self, node):
        """
        Returns list of leafs from subtree from given node.

        :param node: node
        :return:
        """
        leafs = []
        self.__traverse_tree_get_leafs(node, leafs)
        return leafs

    def get_phyl_leafs(self, node):
        """
        Gets side branch leafs of a phylostratum. Give node must be a phyl node.

        :param node: node
        :return:
        """
        if not node.is_phyl_node:
            raise Exception("Node must be phylostratum node!!!")

        leafs = []
        for child in node.children:
            if not child.is_phyl_node:
                leafs.extend(self.get_leafs_subtree(child))
        if node.is_leaf:
            leafs.append(node)
        return leafs

    def get_side_phyl_branch_leafs(self):
        """
        Returns dictionary of phylostratum nodes to corespondinf leafs.
        :return: phyl_2_node dictionary
        """
        phyl_2_leafs = {}

        for ph_node in self.phyl_nodes:
            leafs = []
            for child in ph_node.children:
                if not child.is_phyl_node:
                    leafs.extend(self.get_leafs_subtree(child))
            if ph_node.is_leaf:
                leafs.append(ph_node)
            phyl_2_leafs[ph_node.id] = leafs
        return phyl_2_leafs

    def import_hits(self, path, format="blast-6"):
        """
        Imports hits from a text file, the file should not contain a header.
        Format of genes in database:
        "pgi|0000000000535026343|ti|535026|pi|0|", where fields 0000000000535026343 is gene id and 535026 species id (node)
        used in the nodes and names file used for constructing the phylogeny.
        Supported formats: "blast-6" : tabular format blast flag -outfmt 6  [default]

        :param path: path to hits file
        :param format: format, supported blast-6
        :return:
        """
        if format == "blast-6":
            print("Importing hits from file. Blast format -6.")
            line_counter = 0
            with open(path) as input:
                for line in input:
                    if line_counter % 1000000 == 0:
                        print("Line: " + str(line_counter))
                    line_counter += 1
                    line = line.strip()
                    line_splitted = line.split()
                    if len(line_splitted) != 12:
                        raise ValueError("Import hits file must be in blast -6 format. Line: " + line)
                    query = line_splitted[0].strip()
                    hit = line_splitted[1].strip()
                    node_id = self.__get_node_id(hit)
                    try:
                        sim_value = float(line_splitted[10])
                    except ValueError:
                        print("E-value is not parsable to float. Line: " + line)
                        exit()

                    if query not in self.gene_2_hits:
                        self.gene_2_hits[query] = []

                    self.gene_2_hits[query].append((node_id, sim_value))

        print("Imported " + str(line_counter) + " hits.")
        print("Calculating best hit PS.")

        for gene in self.gene_2_hits:
            if len(self.gene_2_hits[gene]) >= 0:
                min_PS = sys.maxsize
                for hit in self.gene_2_hits[gene]:
                    PS = self.node_id_2_PS[hit[0]]
                    if PS < min_PS:
                        min_PS = PS

                self.gene_2_min_PS[gene] = min_PS

    def import_clusters(self, path, path_clusters=None):
        """
        The clustering phylostratum mapping mode maps genes from focal species to the
        species with min phylostratum in the given cluster.
        File format is a two column tsv file where first column is cluster representative gene and second column is the
        cluster member gene (the cluster member column should include also the rep gene).
        Example:
            pgi|0000000000535026343|ti|535026|pi|0| pgi|0000000000535026343|ti|535026|pi|0|
            pgi|0000000000535026343|ti|535026|pi|0| pgi|0000000000535027344|ti|535027|pi|0|
            pgi|0000000000535026343|ti|535026|pi|0| pgi|0000000000535028345|ti|535028|pi|0|
            ...
        If path_clusters is given the clusters are stored to the given path in a readable format.

        :param path: path for input file
        :param path_clusters: if given clusters are stored to path
        :return:
        """
        counter = 0
        with open(path) as input:
            for line in input:
                counter += 1
                if counter % 1000 == 0:
                    print("Loading line: " + str(counter))
                l_s = line.split()
                if len(l_s) < 2:
                    continue
                rep = l_s[0]
                mem = l_s[1]
                if rep not in self.rep_2_members:
                    self.rep_2_members[rep] = []
                self.rep_2_members[rep].append(mem)

        print("Clusters loaded.\nAssigning phylostrata...")

        if path_clusters:
            cluster_counter = 0
            with open(path_clusters, "w") as out:
                for rep in self.rep_2_members:

                    contains_focal = False
                    for mem in self.rep_2_members[rep]:
                        node_id = self.__get_node_id(mem)
                        if node_id == self.focal_id:
                            contains_focal = True

                    if not contains_focal:
                        continue

                    out.write("(" + str(cluster_counter) + ")\tCluster: " + rep + ":\n")
                    mem_counter = 0
                    for mem in self.rep_2_members[rep]:
                        node_id = self.__get_node_id(mem)
                        out.write("(" + str(mem_counter) + ")\t" + mem + "\tPS=" + str(self.node_id_2_PS[node_id]) + "\n")
                        mem_counter += 1
                    cluster_counter += 1
                    out.write("\n")

        for rep in self.rep_2_members:
            focal_genes = []
            min_PS = sys.maxsize
            for mem in self.rep_2_members[rep]:
                node_id = self.__get_node_id(mem)
                if node_id == self.focal_id:
                    focal_genes.append(mem)

                PS = self.node_id_2_PS[node_id]
                if PS < min_PS:
                    min_PS = PS

            for gene in focal_genes:
                self.gene_2_min_PS[gene] = min_PS

    def __resolve_horizontal_transfer(self, data, min_gap_size=2, hit_perc_threshold=0.1):
        orig_ps_index = -1
        for i in range(len(data)):
            if data[i] != 0:
                orig_ps_index = i
                break
        new_ps_index = orig_ps_index
        gap_counter = 0
        if data[orig_ps_index] < hit_perc_threshold:
            for i in range(orig_ps_index + 1, len(data)):
                if data[i] == 0:
                    gap_counter += 1
                else:
                    if gap_counter >= min_gap_size:
                        new_ps_index = i
                    break
        return new_ps_index + 1

    def __resolve_horizontal_transfer_iter(self, data, min_gap_size=2, hit_perc_threshold=0.1):
        tmp_data = data.copy()
        ps = self.__resolve_horizontal_transfer(tmp_data, min_gap_size=min_gap_size,
                                                hit_perc_threshold=hit_perc_threshold)
        tmp_data = tmp_data[(ps - 1):]
        while len(tmp_data) > 1:
            ps_new = ps + self.__resolve_horizontal_transfer(tmp_data, min_gap_size=min_gap_size,
                                                             hit_perc_threshold=hit_perc_threshold) - 1
            if ps_new == ps:
                ps = ps_new
                break
            else:
                ps = ps_new
                tmp_data = tmp_data[(ps - 1):]

        return ps

    def __get_zero_gap_size(self, data):
        in_gap = False
        max_gap = 0
        gap_len = 0
        first_non_zero = False

        for d in data:
            if d == 0 and not first_non_zero:
                continue
            if d != 0:
                first_non_zero = True

            if d == 0:
                if not in_gap:
                    gap_len = 1
                    in_gap = True
                else:
                    gap_len += 1
            else:
                if gap_len > max_gap:
                    max_gap = gap_len
                gap_len = 0
                in_gap = False

        return max_gap

    def __get_node_id(self, gene_id):
        gene_id = gene_id.strip()
        gene_id_split = gene_id.split("|")
        if len(gene_id_split) < 6:
            raise ValueError("Subject id must be in specified format:"
                             "\"pgi|0000000000535026343|ti|535026|pi|0|\". Line: " + gene_id)
        try:
            node_id = int(gene_id_split[3])
        except ValueError:
            print("Node id-s must be parsable to integer. Line: " + gene_id)
        return node_id

    def __add_nodes_from_child_2_parent_dict(self, node, parent_2_children, node_2_names):
        # checks if node is leaf
        if node.id not in parent_2_children:
            node.is_leaf = True
            return

        for child in parent_2_children[node.id]:
            name = child
            if child in node_2_names:
                name = node_2_names[child]
            else:
                print("ID: " + str(child) + " is not in names file. Name set to id value.")
            new_node = Node(child, name)
            new_node.parent = node
            node.children.append(new_node)
            self.__add_nodes_from_child_2_parent_dict(new_node, parent_2_children, node_2_names)

    def __traverse_tree_get_leafs(self, node, leafs):
        if node.is_leaf:
            leafs.append(node)
        for child in node.children:
            self.__traverse_tree_get_leafs(child, leafs)

    def __check_leaf_children(self, node):
        leaf_children = False
        for child in node.children:
            if child.is_leaf:
                leaf_children = True
                break
        return leaf_children

    def __check_non_leaf_brothers(self, node):
        parent = node.parent
        non_leaf_children_brother = False
        for child in parent.children:
            if child.id == node.id:
                continue
            if not self.__check_leaf_children(child):
                non_leaf_children_brother = True
                break
        return non_leaf_children_brother

    def __has_ancestor(self, node, node_id):
        # Checks if given nodes has a ancestor with node_id.
        tmp_node = node
        has_ancestor = False
        while not tmp_node.is_root:
            tmp_node = tmp_node.parent
            if tmp_node.id == node_id:
                has_ancestor = True
                break
        return has_ancestor

    def __get_node_by_id(self, node, id):
        if node.id == id:
            return node
        else:
            for child in node.children:
                result = self.__get_node_by_id(child, id)
                if result != None:
                    return result
        return None

    def __check_all_siblings_visited(self, node):
        if node.is_root:
            return True

        parent = node.parent
        non_visited_siblings = False

        for child in parent.children:
            if child.id == node.id:
                continue
            if child.id not in visited:
                non_visited_siblings = True
                break
        return not non_visited_siblings

    def __check_siblings(self, node):
        if node.is_root:
            return False
        parent = node.parent
        if len(parent.children) > 2:
            return True
        else:
            return False

    def __check_first_siblings_visited(self, node):
        if node.is_root:
            return True
        parent = node.parent
        first_visited_siblings = True

        for child in parent.children:
            if child.id == node.id:
                continue
            if child.id in visited:
                first_visited_siblings = False
                break
        return first_visited_siblings

    def __traverse_tree_print(self, node, depth, left_string):
        left_space = 60
        right_space = 140
        global visited
        visited.add(node.id)

        if node.is_leaf:
            l_str = left_string
            if self.__check_all_siblings_visited(node.parent):
                l_str = l_str.strip("|")
            sep_len = right_space - len(l_str)
            print(l_str + " " * sep_len + "|" + node.name)
        else:
            if node.is_phyl_node:
                l_str = left_string[len(node.name):]
                l_str = node.name + l_str
                print(l_str + "+|")
            else:
                if self.__check_leaf_children(node):
                    l_str = left_string[(len(node.name) + 2):]
                    l_str = "(" + node.name + ")" + l_str
                    r_length = right_space - len(l_str) - 1
                    print(l_str + "-" * r_length + "+|")
                else:
                    l_str = left_string[(len(node.name)+2):]
                    l_str = "(" + node.name + ")" + l_str
                    print(l_str + "+|")

            if not self.__check_leaf_children(node):
                if self.__check_all_siblings_visited(node):
                    left_string = left_string.strip("|") + "  |"
                else:
                    left_string = left_string + " |"

        for child in node.children:
            self.__traverse_tree_print(child, depth+1, left_string)

    def __traverse_tree_print_file(self, node, depth, left_string, out):
        left_space = 60
        right_space = 140
        global visited
        visited.add(node.id)

        if node.is_leaf:
            l_str = left_string
            if self.__check_all_siblings_visited(node.parent):
                l_str = l_str.strip("|")
            sep_len = right_space - len(l_str)
            out.write(l_str + " " * sep_len + "|" + node.name + "\n")
        else:
            if node.is_phyl_node:
                l_str = left_string[len(node.name):]
                l_str = node.name + l_str
                out.write(l_str + "+|\n")
            else:
                if self.__check_leaf_children(node):
                    l_str = left_string[(len(node.name) + 2):]
                    l_str = "(" + node.name + ")" + l_str
                    r_length = right_space - len(l_str) - 1
                    out.write(l_str + "-" * r_length + "+|\n")
                else:
                    l_str = left_string[(len(node.name)+2):]
                    l_str = "(" + node.name + ")" + l_str
                    out.write(l_str + "+|\n")

            if not self.__check_leaf_children(node):
                if self.__check_all_siblings_visited(node):
                    left_string = left_string.strip("|") + "  |"
                else:
                    left_string = left_string + " |"

        self.__sort_children_ps_last(node.children)
        for child in node.children:
            self.__traverse_tree_print_file(child, depth+1, left_string, out)

    def __sort_children_ps_last(self, children):
        ps_child = None
        for c in children:
            if c.is_phyl_node:
                ps_child = c
        if ps_child:
            children.remove(ps_child)
            children.append(ps_child)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("nodes_path", help="path to tab delimited text nodes file.\nThe file should not contain headers."
                                           "Node IDs should be parsable to int.\nThe root-id must be the parent in the "
                                           "first line of the input file.\n"
                                           "Nodes file format:\n"
                                           "child_id_1 parent_id_1\n"
                                           "child_id_3 parent_id_1\n"
                                           "...\n\n")
    parser.add_argument("names_path", help="path to tab delimited text names file.\nThe file should not contain headers."
                                           "Node IDs should be parsable to int.\n"
                                           "Names file format:\n"
                                           "node_1_id name_1\n"
                                           "node_2_id name_2\n"
                                           "...\n\n")
    parser.add_argument("focal_species_id", help="ID of the focal species which must be parsable to integer.\n\n")
    parser.add_argument("in_path", help="path for a similarity text file, the file should not contain a header.\n "
                                        "In basic mode the input is a blast result file (blast -6 forma), in cluster "
                                        "mode it should be a cluster file (mmseqs .tsv clust file)."
                                        "The cluster mode is chosen with the flag --cluster. Format of genes in"
                                        "database: \"pgi|0000000000535026343|ti|535026|pi|0|\","
                                        "where fields 0000000000535026343 is gene id and 535026 species id (node) used"
                                        "in the nodes and names file used for constructing the phylogeny.")
    parser.add_argument("out_path", help="output path for map")
    parser.add_argument("--cluster", help="sets the mode to cluster. The input file should be a mmseqs cluster in tsv"
                                          "format.", action="store_true")
    parser.add_argument("--ps_merge", help="path for text file that contains mappings for PS. Old -> new.\n"
                                                "File format:\n"
                                                "orig_ps\tnew_ps\n"
                                                "4\t6\n"
                                                "4\t7\n"
                                                "...\n")
    parser.add_argument("--print_tree", help="path for storing a tree in a readable format.")
    parser.add_argument("--hgt", help="horizontal gene transfer analysis flag, if flag is present the out path must"
                                      " be a folder.", action="store_true")
    parser.add_argument("--hgt_gap", help="horizontal gene transfer correction gap size, default 2")
    parser.add_argument("--hgt_zero_filter", help="horizontal gene transfer filter gap size, default 2")

    args = parser.parse_args()

    foc_id = 0
    try:
        foc_id = int(args.focal_species_id)
    except ValueError:
        print("Focal id-s must be parsable to integer. input: " + foc_id)
        exit()

    tree = Tree(args.nodes_path, args.names_path, foc_id)
    if args.cluster:
        tree.import_clusters(args.in_path)
    else:
        tree.import_hits(args.in_path)

    if args.hgt:
        gap = 2
        gap_filter = 2
        if args.hgt_gap:
            try:
                gap = int(args.hgt_gap)
            except ValueError:
                print("Hgt gap must be parsable to integer. input: " + gap)
                exit()

        if args.hgt_gap:
            try:
                gap_filter = int(args.hgt_zero_filter)
            except ValueError:
                print("Zero filter gap must be parsable to integer. input: " + gap_filter)
                exit()

        if args.ps_merge:
            tree.hgt_detection_analysis(args.out_path, hgt_gap=gap, export_hgt_correction=True, zero_filter=gap_filter,
                                        ps_merge_path=args.ps_merge)
        else:
            tree.hgt_detection_analysis(args.out_path, hgt_gap=gap, export_hgt_correction=True, zero_filter=gap_filter)
    else:
        if args.ps_merge:
            tree.export_min_PS_map(args.out_path, ps_merge_path=args.ps_merge)
        else:
            tree.export_min_PS_map(args.out_path)

    if args.print_tree:
        tree.print_tree(args.print_tree)
