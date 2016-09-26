import logging
import itertools
import sys
from argparse import ArgumentParser
import matplotlib_venn
from matplotlib import pyplot as plt

WRONG_FAMILIES_ERROR = """families file "NAME" does not seem to have correct format.
Should be:
Family1: Species1,Species2,Species3
Family2: Species4"""

if __name__ == "__main__":
    parser = ArgumentParser(description="Parses orthomcl groups.txt and singletons, \
            plots a Venn diagram of the number of genes in shared and non-shared clusters")
    parser.add_argument("groups", help="Path to groups.txt")
    parser.add_argument("singletons", help="Path to singletons file (use orthomclSingletons, part of orthomcl")
    parser.add_argument("families", help="Path to file detailing groups of species/families - see README")
    parser.add_argument("figure", help="Output path for final figure (add .svg for SVG picture, add .png for PNG etc.)")
    parser.add_argument("table", help="Output path for final table of groups")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    logger.info("Started parsing.")

    # first get the families/groups from families

    species_dict = {} # key: group name, value: set of species in group
    species_to_group_dict = {} # key: species name, value: group the species is in

    species_counter = 0 # how many species we see
    groups_counter = 0  # how many groups we see

    for line in open(args.families):
        ll = line.rstrip("\n\r").split(" ")
        if (":" not in ll[0]) or (len(ll) != 2):
            logger.error(WRONG_FAMILIES_ERROR.replace("NAME", args.families))
            sys.exit(1)

        group = ll[0].replace(":","")
        species_dict[group] = set()
        groups_counter += 1
        for species in ll[1].split(","):
            species_dict[group].add(species)
            species_to_group_dict[species] = group
            species_counter += 1

    logger.info("Got %s species in %s families/groups."%(species_counter, groups_counter))

    # need a data structure to count the number of overlapping genes/clusters
    counter_dict = {} # key: (group1, group2, ...., groupN), value: (number of overlapping genes, number of overlapping clusters)

    overlap_dict = {} # same key, value: list of shared genes

    all_possible_combinations = list(itertools.product(species_dict.keys(),repeat=groups_counter))
    # [('Milletoids', 'Milletoids', 'Milletoids'), ('Milletoids', 'Milletoids', 'Dalbergiods'), ('Milletoids', 'Milletoids', 'Galegoids'),]
    # take that and make it a) so that each combinaion has only unique members
    # and that b) only AB, not BA is present
    # [("Milletoids"), ("Milletoids", "Dalbergiods"), ("Milletoids","Galegoids",...]

    for combo in all_possible_combinations:
        combo = tuple(sorted(set(combo)))
        counter_dict[combo] = [0, 0]
        overlap_dict[combo] = set()

    logger.info("Parsing groups.txt.")

    with open(args.groups) as f:
        for line in f:
            ll = line.rstrip("\n\r").split(" ")
            cluster = ll[0]
            genes = ll[1:]
            # which species are present?
            species_in_cluster = set()
            number_of_genes_in_cluster = 0
            genes_in_cluster = set()

            for g in genes:
                species = g.split("|")[0]
                if species not in species_to_group_dict:
                    # Do not count species we don't have in input file
                    continue
                genes_in_cluster.add(g)
                species_in_cluster.add(species)
                number_of_genes_in_cluster += 1

            # now compare these species - is it a cluster containing only species of group A, or only species of group B, or is it connecting?
            present_groups = set()
            for g in species_dict:
                intersection = species_in_cluster.intersection(species_dict[g])
                if intersection:
                    present_groups.add(g)
            present_groups = tuple(sorted(present_groups))

            if not present_groups:
                # this happens when our families file does not have all groups that are present in groups.txt
                # and when all genes in the cluster of this line are from the "missing" family
                continue

            counter_dict[present_groups][0] += number_of_genes_in_cluster
            counter_dict[present_groups][1] += 1

            for g in genes_in_cluster:
                overlap_dict[present_groups].add(g)

    # now parse the singletons and add them to the single group clusters
    with open(args.singletons) as f:
        for line in f:
            gene = line.rstrip().split("|")
            # get the species the gene belongs to
            species = gene[0]

            # get the group the species belongs to
            try:
                relevant_group = species_to_group_dict[species]
            except KeyError:
                # again, this happens when our families file doesn't have all species
                continue

            # it's just one gene
            overlap_dict[ (relevant_group, ) ].add('|'.join(gene))
            counter_dict[ (relevant_group, ) ][0] += 1

    logger.info("Writing cluster numbers to %s"%(args.table))
    logger.info("Writing cluster members to files ending in 'shared_gene_names.txt'")
    logger.info("Also printing cluster numbers here:")
    logger.info("Group\tGenes overlap\tClusters overlap")
    with open(args.table, "w") as out:
        out.write("Group\tGenes overlap\tClusters overlap\n")
        for c in counter_dict:
            name = ":".join(c)
            shared_genes = counter_dict[c][0]
            shared_clusters = counter_dict[c][1]
            logger.info("%s\t%s\t%s"%(name, shared_genes, shared_clusters))
            out.write("%s\t%s\t%s\n"%(name, shared_genes, shared_clusters))
            this_out = open('.'.join(args.table.split('.')[:-1]) + '_'.join(list(c)) + '_shared_gene_names.txt', 'w')
            for g in overlap_dict[c]:
                this_out.write('%s\n'%g)

    # now make the Venn diagram
    labels = []
    subsets = []
    
    if groups_counter == 3:
        venn_u = matplotlib_venn.venn3_unweighted
        venn_w = matplotlib_venn.venn3

        # matplotlib_venn's order is:
        # [groupA, groupB, groupAB, groupC, groupAC, groupBC, groupABC]

        # this is SUPERUGLY since it's Friday 4pm
        # pull out the single groups 
        labels = [0, 0, 0]
        subsets = [0, 0, 0, 0, 0, 0, 0]

        singles = []

        for c in counter_dict:
            if len(c) == 1:
                singles.append(c)

        a = singles[0]
        b = singles[1]
        c = singles[2]
        subsets[0] = counter_dict[a][0]
        subsets[1] = counter_dict[b][0]
        subsets[3] = counter_dict[c][0]
        labels[0] = a[0]
        labels[1] = b[0]
        labels[2] = c[0]
        # now the doubles
        ab = tuple(sorted(a + b))
        subsets[2] = counter_dict[ab][0]
        ac = tuple(sorted(a + c))
        subsets[4] = counter_dict[ac][0]
        bc = tuple(sorted(b + c))
        subsets[5] = counter_dict[bc][0]
        # now the triples (core)
        abc = tuple(sorted(a + b + c))
        subsets[6] = counter_dict[abc][0]
        
    elif groups_counter == 2:
        venn_u = matplotlib_venn.venn2_unweighted
        venn_w = matplotlib_venn.venn2
        # order is [groupA, groupB, intersection]

        for a in counter_dict:
            if len(a) == 1:
                subsets.append(counter_dict[a][0])
                labels.append(a[0])
            else:
                # no guarantee that it comes last
                last_sub = counter_dict[a][0]
                last_label = ":".join(a)

        subsets.append(last_sub)
        labels.append(last_label)
    else:
        logger.warning("You have %s groups/families, but currently this program only \
supports Venn diagrams for up to 3 groups. Exiting."%groups_counter)
        sys.exit(0)

    figure = args.figure.split(".")
    figure_type = figure[-1]
    figure_name = ".".join(figure[:-1])
    unweighted_out = "%s_unweighted.%s"%(figure_name, figure_type)
    weighted_out = "%s_weighted.%s"%(figure_name, figure_type)

    logger.info("Plotting unweighted Venn diagram to %s and weighted Venn diagram \
to %s"%(unweighted_out, weighted_out))

    venn_u(subsets=subsets, set_labels=labels)
    plt.savefig(unweighted_out)
    plt.close()
    venn_w(subsets=subsets, set_labels=labels)
    plt.savefig(weighted_out)
    plt.close()
