import numpy as np
import functions

import matplotlib

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch

from matplotlib.legend_handler import HandlerLine2D
from matplotlib.patches import Circle

def linear_graph(linear, brackets, custom_gradient, gradient_0, gradient_color_graph, mutation_imposition, fh,
                   rgb_start_global, rgb_end_global, input_color, list_nucleotides, nucleotides, alpha_0):


    index=list(range(0,len(brackets)-1))
    len_index=len(index)

    abcissa=np.linspace(0,len(brackets),len(brackets)*5)


    if not custom_gradient:
        xy_size=('1920','1080')
        rgb_start=('255','0','0')
        rgb_end=('255','255','0')
    else:
        xy_size=('1920','1080')
        rgb_start=rgb_start_global
        rgb_end=rgb_end_global

    color_spectrum=functions.gen_gradient(xy_size, rgb_start, rgb_end)


    bonds=functions.bond_finder(brackets)




    if gradient_0:
        sqrt_prob_1=functions.sqrt_finder(bonds, fh)
        gradient_bonds=[]
        for i in list(range(0,len(bonds))):
            bonds[i].append(sqrt_prob_1[i])
            gradient_bonds.append(bonds[i])


    functions.bond_grapher(abcissa, linear, bonds,input_color, gradient_0, gradient_color_graph,
                           color_spectrum=color_spectrum, mutation_imposition=mutation_imposition,
                           alpha_0=alpha_0)

    if nucleotides:

        for i in list(range(0,len_index)):


            if '/' in list_nucleotides[i]:
                tp = TextPath((i-0.15, 0), list_nucleotides[i], size=0.4)
            else:
                tp = TextPath((i, 0), list_nucleotides[i], size=0.4)
            nucleotide = plt.Circle(([i+0.15, 0.15]), 0.5, edgecolor='black', facecolor='white')
            ax = plt.gca()
            ax.add_patch(nucleotide)
            ax.add_patch(PathPatch(tp, color="black"))





