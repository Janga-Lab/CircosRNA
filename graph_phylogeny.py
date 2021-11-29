import functions
import numpy as np

import math
import matplotlib.pyplot as plt

from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch

def alignment_bracket_processing(serial_brackets, serial_alignment, n):
    spine = [i[0] for i in serial_brackets]
    serial_brackets = [[i, ii[1]] for i in spine for ii in serial_brackets if i == ii[0]]
    serial_alignment = [[i, ii[1]] for i in spine for ii in serial_alignment if i == ii[0]]

    for i, j in enumerate(serial_brackets):
        serial_brackets[i][1] = functions.impure_brackets(serial_brackets[i][1], serial_alignment[i][1])

    serial_bonds = [[i[0], functions.bond_finder(i[1])] for i in serial_brackets]

    bond_collection = []
    for i, j in enumerate(serial_bonds):
        for ii, jj in enumerate(serial_bonds[i][1]):
            if jj not in bond_collection:
                bond_collection.append(jj)

    all_bonds = [ii for i in serial_bonds for ii in i[1]]

    for i, j in enumerate(bond_collection):
        count = 0
        for ii in all_bonds:
            if bond_collection[i][0] == ii[0] and bond_collection[i][1] == ii[1]:
                count += 1
                bond_collection[i].append((count - n ) / (len(serial_bonds) - n))

    bond_collection = [i for i in bond_collection if i[-1] != 0.0 and not i[-1] < 0.0]

    bond_list_backbone = [[i[0], i[1]] for i in bond_collection]
    bond_list_frequency = [i[-1] for i in bond_collection]

    bond_backbone_frequency=list(zip(bond_list_frequency,bond_list_backbone))
    bond_backbone_frequency=sorted(bond_backbone_frequency, key = lambda x: x[0])

    bond_list_backbone=[i[1] for i in bond_backbone_frequency]
    bond_list_frequency=[i[0] for i in bond_backbone_frequency]

    indices = [ii for i in bond_list_backbone for ii in i]

    low_value=1/(len(serial_bonds)-n)

    return [bond_list_backbone, bond_list_frequency, indices, low_value]


def linear_nucleotide_plot(serial_alignment, colored_nucleotides, nucleotide_color_list):

    for ii,jj in enumerate(serial_alignment):
        for i in list(range(0,len(serial_alignment[ii][1]))):
            if not colored_nucleotides:
                tp = TextPath((i, 0-ii), serial_alignment[ii][1][i], size=0.4)
                nucleotide = plt.Circle(([i+0.15, 0.15-(ii)]), 0.5, edgecolor='black', facecolor='white')
            elif colored_nucleotides:
                if serial_alignment[ii][1][i] =='A':
                    nuke_color_designation=nucleotide_color_list[0]
                    edge_designation=None
                elif serial_alignment[ii][1][i] == 'T':
                    nuke_color_designation=nucleotide_color_list[1]
                    edge_designation=None
                elif serial_alignment[ii][1][i] == 'G':
                    nuke_color_designation=nucleotide_color_list[2]
                    edge_designation=None
                elif serial_alignment[ii][1][i] == 'C':
                    nuke_color_designation=nucleotide_color_list[3]
                    edge_designation=None
                else:
                    nuke_color_designation='white'
                    edge_designation='black'


                tp = TextPath((i+0.85, -0.65 -ii), serial_alignment[ii][1][i], size=0.4)
                nucleotide=plt.Rectangle((i+0.5,-1-(ii)), 1, 1, edgecolor=edge_designation,facecolor=nuke_color_designation)
            ax = plt.gca()
            ax.add_patch(nucleotide)
            ax.add_patch(PathPatch(tp, color="black"))

def circle_graph_mod(brackets, bonds, color_spectrum, serial_nucleotides, nucleotides, colored_nucleotides,  bond_ubiquity, nucleotide_color_list):



    len_index = len(list(range(0, len(brackets) - 1)))

    abcissa_0 = []
    ordinate_0 = []

    circle_list_angles = list(np.linspace(0, 2 * math.pi, 3000))

    all_coordinates = []

    for i in bonds:
        coordinates = []

        theta_0 = ((i[0] / len_index) * 2 * math.pi) + (math.pi / 2)
        theta_1 = ((i[1] / len_index) * 2 * math.pi) + (math.pi / 2)

        coordinates.append(functions.polar_to_cartesian((len_index / (2 * math.pi)), theta_0))
        coordinates.append(theta_0)
        coordinates.append(functions.polar_to_cartesian((len_index / (2 * math.pi)), theta_1))
        coordinates.append(theta_1)
        all_coordinates.append(coordinates)

    regression_plots = []

    for i in all_coordinates:

        theta_mean = functions.middle_man(i[1], i[3])[0]

        r = (((math.pi - functions.middle_man(i[1], i[3])[1]) / (math.pi))) * (len_index / (2 * math.pi))

        i.pop(1)
        i.pop(2)
        i.append(functions.polar_to_cartesian(r, theta_mean))

        abcissa = []
        ordinate = []

        for ii in i:
            abcissa.append(ii[0])
            ordinate.append(ii[1])
        regression_plots.append([abcissa, ordinate])


    for i in range(0,len(regression_plots)):
        x_m = (regression_plots[i][0][0] + regression_plots[i][0][1]) / 2
        y_m = (regression_plots[i][1][0] + regression_plots[i][1][1]) / 2
        l = math.sqrt(((x_m - regression_plots[i][0][0]) ** 2) + ((y_m - regression_plots[i][1][0]) ** 2))
        h = math.sqrt(((x_m - regression_plots[i][0][2]) ** 2) + ((y_m - regression_plots[i][1][2]) ** 2))

        abcissa = (-l, 0, l)
        ordinate = (0, h, 0)
        polynomial = list(np.polyfit(abcissa, ordinate, 2))

        x_vector = list(np.linspace(-l, l, 1000))
        y_vector = [((ii ** 2) * polynomial[0]) + (ii * polynomial[1]) + polynomial[2] for ii in x_vector]

        maxima_index = y_vector.index(max(y_vector))

        transformed_vectors_default = functions.transrotation(x_vector, y_vector,
                                                              regression_plots[i][1][0] - regression_plots[i][1][1],
                                                              regression_plots[i][0][0] - regression_plots[i][0][1],
                                                              flip_it=False, iteration=regression_plots[i])

        transformed_vectors_test = functions.transrotation(x_vector, y_vector, regression_plots[i][1][0] - regression_plots[i][1][1], regression_plots[i][0][0] - regression_plots[i][0][1],
                                                           flip_it=True, iteration=regression_plots[i])

        d_default = math.sqrt(((transformed_vectors_default[0][maxima_index]) ** 2) + (
                    (transformed_vectors_default[1][maxima_index]) ** 2))

        d_test = math.sqrt(
            ((transformed_vectors_test[0][maxima_index]) ** 2) + ((transformed_vectors_test[1][maxima_index]) ** 2))

        if d_test < d_default:

            transformed_vectors = transformed_vectors_test
        else:

            transformed_vectors = transformed_vectors_default

        plt.plot(transformed_vectors[0], transformed_vectors[1],
                 color=tuple(color_spectrum[round(((1 - bond_ubiquity[i]) * len(color_spectrum)))])
                 , zorder=0)

    start = math.pi / 2
    jump = (2 * math.pi) / 6
    r = len_index / (2 * math.pi)
    for i in list(range(0, 6)):
        plt.text((1.15 * r) * math.cos(start + i * jump), (1.10 * r) * math.sin(start + i * jump),
                str(round((i / 6) * len_index)), horizontalalignment='center')

    plt.text((1.15 * r) * math.cos(start), (1.20 * r) * math.sin(start),
            str(len_index
                ), horizontalalignment='center')

    for i in circle_list_angles:
        abcissa_0.append(functions.polar_to_cartesian(((len_index) / (2 * math.pi)), i)[0])
        ordinate_0.append(functions.polar_to_cartesian(((len_index) / (2 * math.pi)), i)[1])

    if not nucleotides:
        plt.plot(abcissa_0, ordinate_0, 'black')

    if nucleotides or colored_nucleotides:
        for ii,jj in enumerate(serial_nucleotides):
            d_b_1=(1+ii*((2*math.pi)/len_index))/2
            d_b_2 = ((2*d_b_1 + ((2*math.pi)/len_index))) / 2

            index_coordinates = [functions.polar_to_cartesian((len_index / (2 * math.pi) + ii),(i / len_index) * 2 * math.pi + (1 / 2) * math.pi)
                                 for i in range(0, len_index)]
            outer_index_coordinates=[functions.polar_to_cartesian((len_index / (2 * math.pi) + 1 + ii),(i / len_index) * 2 * math.pi + (1 / 2) * math.pi)
                                     for i in range(0, len_index)]
            for i in range(0, len_index):

                if nucleotides:
                    tp = TextPath((index_coordinates[i][0], index_coordinates[i][1]), serial_nucleotides[ii][1][i], size=0.4)


                    nucleotide = plt.Circle((index_coordinates[i][0] + 0.15, index_coordinates[i][1] + 0.15), 0.5,
                                            edgecolor='black', facecolor='white')


                if colored_nucleotides:

                    theta_1=(i / len_index) * 2 * math.pi + (1 / 2) * math.pi
                    theta_2_1 = theta_1 + (math.pi / 2) - math.pi
                    theta_2_2 = theta_1 - (math.pi / 2) + math.pi

                    bottom_middle_coordinate=(index_coordinates[i][0], index_coordinates[i][1])
                    top_middle_coordinate=(outer_index_coordinates[i][0],outer_index_coordinates[i][1])

                    point_1=(bottom_middle_coordinate[0] + d_b_1*math.cos(theta_2_1),
                             bottom_middle_coordinate[1] + d_b_1*math.sin(theta_2_1))

                    point_2 = (bottom_middle_coordinate[0] + d_b_1 * math.cos(theta_2_2),
                               bottom_middle_coordinate[1] + d_b_1 * math.sin(theta_2_2))

                    point_3 = (top_middle_coordinate[0] + d_b_2 * math.cos(theta_2_1),
                               top_middle_coordinate[1] + d_b_2 * math.sin(theta_2_1))

                    point_4 = (top_middle_coordinate[0] + d_b_2 * math.cos(theta_2_2),
                               top_middle_coordinate[1] + d_b_2 * math.sin(theta_2_2))

                    points=[point_1, point_2, point_4, point_3]

                    center=((outer_index_coordinates[i][0]+index_coordinates[i][0])/2-0.15,(outer_index_coordinates[i][1]+index_coordinates[i][1])/2-0.15)

                    tp = TextPath(center, serial_nucleotides[ii][1][i], size=0.4)

                    if serial_nucleotides[ii][1][i] == 'A':
                        nuke_color_designation = nucleotide_color_list[0]
                        edge_designation = None
                    elif serial_nucleotides[ii][1][i] == 'T':
                        nuke_color_designation = nucleotide_color_list[1]
                        edge_designation = None
                    elif serial_nucleotides[ii][1][i] == 'G':
                        nuke_color_designation = nucleotide_color_list[2]
                        edge_designation = None
                    elif serial_nucleotides[ii][1][i] == 'C':
                        nuke_color_designation = nucleotide_color_list[3]
                        edge_designation = None
                    else:
                        nuke_color_designation = 'white'
                        edge_designation = 'black'


                    nucleotide=plt.Polygon(points, facecolor=nuke_color_designation, edgecolor=edge_designation)


                ax = plt.gca()
                ax.add_patch(nucleotide)
                ax.add_patch(PathPatch(tp, color="black"))

    plt.axis('off')

    plt.axis('square')



