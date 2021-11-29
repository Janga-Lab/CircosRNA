#!/usr/bin/env python

import argparse
import sys
import os

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.patches import Circle
import matplotlib.image as img
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as clr

import numpy as np
import math

from PIL import Image
import webcolors as wc

import functions
import circle_graph
import linear_graph
import graph_phylogeny

def main():
    parser=argparse.ArgumentParser(description='CLI for graphing of dotbracket notation of RNA secondary structure')
    parser.add_argument('-i', '--input', help='input file for the dot-bracket structure' , type=str)
    parser.add_argument('-i2', '--secondary_input', help="secondary input file, also signals the superimposition of one linear/circular representation onto another",type=str)
    parser.add_argument('-l', "--linear", help="Produces the linear representation of RNA secondary structure", action="store_true")
    parser.add_argument('-c', "--circular", help="Produces the circular representation of RNA secondary structure", action="store_true")
    parser.add_argument('-c1', "--color_one", help='selected color for the plot, the default is lightblue', type=str)
    parser.add_argument('-c2', "--color_two", help='selected color for the superimposed plot, the default is lightgreen',type=str)
    parser.add_argument('-c3', "--color_three",help="When graphs are superimposed, the overlapping areas should be set to a seperate color, something which contrasts well is recommended",type=str)
    parser.add_argument('-c4', "--color_four", help="overlap color of unaligned regions if -a2 is chosen", type=str)
    parser.add_argument('-a', "--align", help="Align the nucleotides before checking for structural similarities (recommended)", action="store_true")
    parser.add_argument( '-sa','--second_alignment_path', help="align and permit the overlaping regions to be some fourth color of your choice", action="store_true")
    parser.add_argument('-o', '--overlay_alpha', help='transparency value between 0 and 1 for the overlying plot in superimposed graphs', type=str)
    parser.add_argument('-u', '--underlay_alpha', help='transparency value between 0 and 1 for the underlying plot in superimposed graphs', type=str)
    parser.add_argument('-m', '--match_radii', help='by default, circular representations of secondary structure will adapt to polymer length, including this argument will cause the circular graphs to adopt uniform radii', action="store_true")
    parser.add_argument('-st','--structures', help='input files for the dot-bracket structure' ,nargs='+')
    parser.add_argument('-MSA', '--MSA_file', help='input MSA alignment output from CLUSTALW' , type=str)
    parser.add_argument('-ct', '--cutoff', help='number of homologous sequences to be ignored (start from 0, the default, and work your way up (1,2,3....) if in doubt')
    parser.add_argument('-cn', '--colored_nucleotides', help='colored nucleotide alignment for MSA', action="store_true")
    parser.add_argument('-nc', '--nucleotide_colors', help='specific colors for nucleotides given the \"colored_nucleotides\" command, the colors are ordered as A,T,G, and C, default colors are lime-green, orange-red, blue-violet, and cyan', nargs='+')
    parser.add_argument('-pm', '--p_matrix_input', help='required input *dp.ps containing the ViennaRNA generated probability matrix. Triggers gradient generation' , type=str)
    parser.add_argument('-lc', '--low_prob_color', nargs='+', help='add the low rgb values for the custom gradient', type=str)
    parser.add_argument('-hc', '--high_prob_color', nargs='+', help='add the high rgb values for the custom gradient', type=str)
    parser.add_argument('-g', '--gradient_legend', help="adds a legend to gradient graphs to show which color corresponds to a low probability and which color coresponds to a high probability", action="store_true")
    parser.add_argument('-n', '--nucleotides', help='adds nucleotides to the visualization', action="store_true")
    parser.add_argument('-d', '--dpi', help='enter the dpi needed for supderimposed graphs, there is no "one size fits all", raise or lower this value as needed, start with 96 if in doubt', type=int)


    args=parser.parse_args()

    if args.nucleotide_colors:
        nucleotide_colors=args.nucleotide_colors
    else:
        nucleotide_colors=['limegreen','orangered','blueviolet','cyan']


    if args.colored_nucleotides:
        colored_nucleotides=True
    else:
        colored_nucleotides=False

    if not args.color_four:
        color_4='purple'
    else:
        color_4=args.color_four

    if not args.structures:
        brackets=functions.bracket_reader(args.input)
    else:
        brackets='brackets'

    brackets_length_one=len(brackets)

    if not args.overlay_alpha:
        overlay_alpha=1
    else:
        overlay_alpha=float(args.overlay_alpha)

    if not args.underlay_alpha:
        underlay_alpha=1
    else:
        underlay_alpha=float(args.underlay_alpha)

    if not args.color_one:
        color_1='lightblue'
    else:
        color_1=args.color_one

    if not args.color_two:
        color_2='lightgreen'
    else:
        color_2=args.color_two


    if args.p_matrix_input:

        fh=functions.prob_matrix(args.p_matrix_input)

        gradient_0=True
    else:
        fh=0
        gradient_0=False

    if args.p_matrix_input:
        sqrt_prob=functions.sqrt_finder(functions.bond_finder(brackets), fh)
    else:
        sqrt_prob=0

    if args.low_prob_color and args.high_prob_color:
        custom_gradient=True
        try:
            rgb_start_global=tuple(args.high_prob_color)
            rgb_end_global=tuple(args.low_prob_color)
            color_spectrum = functions.gen_gradient(('1920', '1080'), args.low_prob_color, args.high_prob_color)
        except:

            s_v=list(wc.name_to_rgb(args.high_prob_color[0]))
            rgb_start_global=(s_v[0], s_v[1], s_v[2])
            s_v = list(wc.name_to_rgb(args.low_prob_color[0]))
            rgb_end_global = (s_v[0], s_v[1], s_v[2])
            color_spectrum = functions.gen_gradient(('1920', '1080'), list(rgb_start_global), list(rgb_end_global))

    else:
        custom_gradient=False
        rgb_start_global=0
        rgb_end_global=0
        color_spectrum=0

    if not custom_gradient:
        rgb_start=('255','0','0')
        rgb_end=('255','255','0')
        color_spectrum=functions.gen_gradient(('1920','1080'),rgb_start,rgb_end)


    if args.low_prob_color and args.high_prob_color:
        gradient_color_graph=True
    else:
        gradient_color_graph=False

    args=parser.parse_args()

    if args.secondary_input:

        mutation_imposition=True

    else:
        mutation_imposition=False

    if args.dpi:
        my_dpi=args.dpi
    else:
        my_dpi=150

    if args.color_three:
        color_3=args.color_three
    else:
        color_3='pink'

    if not args.structures:
        list_nucleotides=functions.nucleotide_list(args.input)

    if args.cutoff:
        n=int(args.cutoff)
    else:
        n=0

    if args.secondary_input and args.linear and not args.align and not args.second_alignment_path:

        max_bond=functions.height_finder(brackets,functions.bracket_reader(args.secondary_input))

        name_1=functions.name(args.input)


        linear_graph.linear_graph(linear=True, brackets=brackets, custom_gradient=custom_gradient,
                                  gradient_0=gradient_0, gradient_color_graph=gradient_color_graph,
                                  mutation_imposition=mutation_imposition, fh=fh,
                                  rgb_start_global=rgb_start_global, rgb_end_global=rgb_end_global,
                                  input_color=color_1, list_nucleotides=list_nucleotides,
                                  nucleotides=args.nucleotides, alpha_0=1)



        bond_zero = mpatches.Patch(color=color_1, label=name_1)

        brackets=functions.bracket_reader(args.secondary_input)
        brackets_length_two=len(brackets)

        name_2=functions.name(args.secondary_input)
        bonds_one = mpatches.Patch(color=color_2, label=name_2)

        max_brackets=max(brackets_length_one, brackets_length_two)

        plt.xlim(0,max_brackets)
        plt.ylim(0, 0.6*max_bond)


        plt.yticks([])
        overlap_patch=mpatches.Patch(color=color_3, label='OVERLAP')
        plt.legend(handles=[bond_zero, bonds_one, overlap_patch], loc='upper right', fontsize='xx-small')
        plt.savefig('plot_0000000000000000_1.png', dpi=my_dpi)
        plt.clf()

        linear_graph.linear_graph(linear=True, brackets=brackets, custom_gradient=custom_gradient,
                                          gradient_0=gradient_0, gradient_color_graph=gradient_color_graph,
                                          mutation_imposition=mutation_imposition, fh=fh,
                                          rgb_start_global=rgb_start_global, rgb_end_global=rgb_end_global,
                                            input_color=color_2, list_nucleotides=list_nucleotides,
                                            nucleotides=args.nucleotides, alpha_0=1)
        plt.xlim(0,max_brackets)
        plt.ylim(0, 0.6*max_bond)
        plt.legend(handles=[bond_zero, bonds_one, overlap_patch], loc='upper right', fontsize='xx-small')
        plt.yticks([])
        plt.savefig('plot_0000000000000000_2.png', dpi=my_dpi)


        image=functions.clean_imposition(list(wc.name_to_rgb(color_1)),
                                         list(wc.name_to_rgb(color_2)),
                                         list(wc.name_to_rgb(color_3)),
                                         'plot_0000000000000000_1.png',
                                         'plot_0000000000000000_2.png')

        os.remove('plot_0000000000000000_1.png')
        os.remove('plot_0000000000000000_2.png')

        image.show()
        save_question=input('Would you like to save the image? (Y/n)')
        if save_question == 'Y':
            name_3=(str(name_1)).strip() + '&' + (str(name_2)).strip() + '_linear.png'
            image.save(name_3)


    elif args.linear and args.second_alignment_path:
        name_1 = functions.name(args.input)

        match_brackets, first_idiosyncratic_brackets, second_idiosyncratic_brackets, first_alignment, second_alignment \
            = functions.align_redefine(
            functions.nucleotide_string(args.input), functions.nucleotide_string(args.secondary_input),

            functions.bracket_reader(args.input), functions.bracket_reader(args.secondary_input))

        list_nucleotides = functions.aligned_nucleotide_list(first_alignment, second_alignment)

        max_bond = max(functions.height_finder(first_idiosyncratic_brackets, match_brackets),
                       functions.height_finder(second_idiosyncratic_brackets, match_brackets))

        linear_graph.linear_graph(linear=True, brackets=first_idiosyncratic_brackets, custom_gradient=custom_gradient,
                                  gradient_0=gradient_0, gradient_color_graph=gradient_color_graph,
                                  mutation_imposition=mutation_imposition, fh=fh,
                                  rgb_start_global=rgb_start_global, rgb_end_global=rgb_end_global,
                                  input_color=color_1, list_nucleotides=list_nucleotides,
                                  nucleotides=first_alignment, alpha_0=1)

        linear_graph.linear_graph(linear=True, brackets=match_brackets, custom_gradient=custom_gradient,
                                  gradient_0=gradient_0, gradient_color_graph=gradient_color_graph,
                                  mutation_imposition=mutation_imposition, fh=fh,
                                  rgb_start_global=rgb_start_global, rgb_end_global=rgb_end_global,
                                  input_color=color_3, list_nucleotides=list_nucleotides,
                                  nucleotides=args.nucleotides, alpha_0=1)

        bond_zero = mpatches.Patch(color=color_1, label=name_1)

        brackets = functions.bracket_reader(args.secondary_input)
        brackets_length_two = len(brackets)

        name_2 = functions.name(args.secondary_input)
        bonds_one = mpatches.Patch(color=color_2, label=name_2)

        max_brackets = max(brackets_length_one, brackets_length_two)

        plt.xlim(0, max_brackets)
        plt.ylim(0, 0.6*max_bond)

        plt.yticks([])
        alignment_patch = mpatches.Patch(color=color_3, label='ALIGNMENT')
        overlap_patch=mpatches.Patch(color=color_4, label='OVERLAP')

        plt.legend(handles=[bond_zero, bonds_one, alignment_patch, overlap_patch], loc='upper right', fontsize='xx-small')
        plt.savefig('plot_0000000000000000_1.png', dpi=my_dpi)
        plt.clf()

        linear_graph.linear_graph(linear=True, brackets=second_idiosyncratic_brackets, custom_gradient=custom_gradient,
                                  gradient_0=gradient_0, gradient_color_graph=gradient_color_graph,
                                  mutation_imposition=mutation_imposition, fh=fh,
                                  rgb_start_global=rgb_start_global, rgb_end_global=rgb_end_global,
                                  input_color=color_2, list_nucleotides=list_nucleotides,
                                  nucleotides=second_alignment, alpha_0=1)
        linear_graph.linear_graph(linear=True, brackets=match_brackets, custom_gradient=custom_gradient,
                                  gradient_0=gradient_0, gradient_color_graph=gradient_color_graph,
                                  mutation_imposition=mutation_imposition, fh=fh,
                                  rgb_start_global=rgb_start_global, rgb_end_global=rgb_end_global,
                                  input_color=color_3, list_nucleotides=list_nucleotides,
                                  nucleotides=args.nucleotides, alpha_0=1)

        plt.xlim(0, max_brackets)
        plt.ylim(0, 0.6*max_bond)
        plt.legend(handles=[bond_zero, bonds_one, alignment_patch, overlap_patch], loc='upper right', fontsize='xx-small')
        plt.yticks([])
        plt.savefig('plot_0000000000000000_2.png', dpi=my_dpi)

        image = functions.clean_imposition(list(wc.name_to_rgb(color_1)),
                                           list(wc.name_to_rgb(color_2)),
                                           list(wc.name_to_rgb(color_4)),
                                           'plot_0000000000000000_1.png',
                                           'plot_0000000000000000_2.png')

        os.remove('plot_0000000000000000_1.png')
        os.remove('plot_0000000000000000_2.png')

        image.show()
        save_question = input('Would you like to save the image? (Y/n)')
        if save_question == 'Y':
            name_3 = (str(name_1)).strip() + '&' + (str(name_2)).strip() + '_linear.png'
            image.save(name_3)

    elif args.secondary_input and args.linear and args.align:
        match_brackets, first_idiosyncratic_brackets, second_idiosyncratic_brackets, first_alignment, second_alignment\
            =functions.align_redefine(
            functions.nucleotide_string(args.input),functions.nucleotide_string(args.secondary_input),

        functions.bracket_reader(args.input),functions.bracket_reader(args.secondary_input))

        list_nucleotides=functions.aligned_nucleotide_list(first_alignment, second_alignment)

        linear_graph.linear_graph(linear=True, brackets=first_idiosyncratic_brackets, custom_gradient=custom_gradient,
                                  gradient_0=gradient_0, gradient_color_graph=gradient_color_graph,
                                  mutation_imposition=mutation_imposition, fh=fh,
                                  rgb_start_global=rgb_start_global, rgb_end_global=rgb_end_global,
                                  input_color=color_1, list_nucleotides=list_nucleotides,
                                  nucleotides=False, alpha_0=underlay_alpha)



        linear_graph.linear_graph(linear=True, brackets=second_idiosyncratic_brackets, custom_gradient=custom_gradient,
                                  gradient_0=gradient_0, gradient_color_graph=gradient_color_graph,
                                  mutation_imposition=mutation_imposition, fh=fh,
                                  rgb_start_global=rgb_start_global, rgb_end_global=rgb_end_global,
                                  input_color=color_2, list_nucleotides=list_nucleotides,
                                  nucleotides=False, alpha_0=overlay_alpha)

        linear_graph.linear_graph(linear=True, brackets=match_brackets, custom_gradient=custom_gradient,
                                  gradient_0=gradient_0, gradient_color_graph=gradient_color_graph,
                                  mutation_imposition=mutation_imposition, fh=fh,
                                  rgb_start_global=rgb_start_global, rgb_end_global=rgb_end_global,
                                  input_color=color_3, list_nucleotides=list_nucleotides,
                                  nucleotides=args.nucleotides, alpha_0=1)

        name_1 = functions.name(args.input)
        bond_zero = mpatches.Patch(color=color_1, label=name_1)
        name_2 = functions.name(args.secondary_input)
        bonds_one = mpatches.Patch(color=color_2, label=name_2)
        overlap_patch = mpatches.Patch(color=color_3, label='OVERLAP')
        plt.legend(handles=[bond_zero, bonds_one, overlap_patch], loc='upper right', fontsize='small')
        plt.yticks([])
        if not args.nucleotides:
            abcissa_x=list(range(0,len(match_brackets)))
            abcissa_y=[0 for i in range(0,len(match_brackets))]
            plt.plot(abcissa_x,abcissa_y, 'black')
        plt.show()

    else:
        if args.linear and not args.MSA_file and not args.structures:
            linear_graph.linear_graph(linear=True, brackets=brackets, custom_gradient=custom_gradient,
                                              gradient_0=gradient_0, gradient_color_graph=gradient_color_graph,
                                              mutation_imposition=mutation_imposition, fh=fh,
                                          rgb_start_global=rgb_start_global, rgb_end_global=rgb_end_global,
                                            input_color=color_1, list_nucleotides=list_nucleotides,
                                            nucleotides=args.nucleotides, alpha_0=1)

            plt.title(functions.name(args.input))
            plt.yticks([])

            if gradient_0 and args.gradient_legend:
                attempted_colormap=[color_spectrum[-1], color_spectrum[round(len(color_spectrum)/2)],color_spectrum[0]]
                my_cmap=LinearSegmentedColormap.from_list('my_cmap', attempted_colormap)
                plt.scatter([-1,-2,-3,-4,-5,-6,-7,-8,-9,-10, -11],[1,1,1,1,1,1,1,1,1,1,1], c=[0,0.1, 0.2, 0.3 , 0.4, 0.5, 0.6, 0.7 , 0.8, 0.9, 1.0], cmap=my_cmap)
                plt.colorbar(label='bond probability (0.0 - 1.0)')
                plt.xlim(0, len(brackets))
                plt.ylim(0, len(brackets) / 2)
            plt.show()
        elif args.linear and args.MSA_file and args.structures:

            serial_brackets = [[functions.name(i), functions.bracket_reader(i)] for i in args.structures]
            serial_names=[i[0] for i in serial_brackets]
            serial_names=', '.join(serial_names)
            serial_alignment = functions.MSA_alignment_extraction(args.MSA_file)
            bond_list_backbone, bond_list_frequency, indices, low_value = graph_phylogeny.alignment_bracket_processing(serial_brackets,
                                                                                            serial_alignment , n)

            abcissa = np.linspace(0, max(indices), max(indices))
            functions.bond_grapher_2(abcissa=abcissa, linear=True, bonds=bond_list_backbone, color='purple',
                                     gradient_0=False,
                                     gradient_color_graph=False,
                                     color_spectrum=color_spectrum, mutation_imposition=False, alpha_0=1, mosaic=True,
                                     bond_ubiquity=bond_list_frequency)

            plt.yticks([])

            if args.nucleotides or args.colored_nucleotides:
                graph_phylogeny.linear_nucleotide_plot(serial_alignment,colored_nucleotides,nucleotide_color_list=nucleotide_colors)

            if args.gradient_legend:

                attempted_colormap = [color_spectrum[-1], color_spectrum[round(len(color_spectrum) / 2)], color_spectrum[0]]
                my_cmap = LinearSegmentedColormap.from_list('my_cmap', attempted_colormap)


                if n ==len(serial_alignment)-1:
                    color_list = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

                else:
                    color_list = list(np.linspace(low_value, 1.0, 11))

                plt.scatter([-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                            c=color_list, cmap=my_cmap)
                plt.colorbar(label='bond conservation probability ('+ str(round(low_value,3)) + ' - 1.0)')
            plt.xlim(0,len(serial_alignment[0][1]))
            plt.plot(list(range(0,len(serial_alignment[0][1]))),[0 for i in range(0,len(serial_alignment[0][1]))],
                     color='black', zorder=0)
            plt.title(serial_names)
            plt.show()




    if args.secondary_input and args.circular:

        if args.match_radii and not args.second_alignment_path:


            brackets_1=brackets
            brackets_2 = functions.bracket_reader(args.secondary_input)
            name_1 = functions.name(args.input)

            bond_zero = mpatches.Patch(color=color_1, label=name_1)
            name_2 = functions.name(args.secondary_input)

            bond_one = mpatches.Patch(color=color_2, label=name_2)

            if len(functions.bracket_reader(args.input)) < len(functions.bracket_reader(args.secondary_input)):
                brackets_1, brackets_2 = brackets_2, brackets_1
                name_1, name_2= name_2, name_1
                bond_zero, bond_one = bond_one, bond_zero


            len_index_input=circle_graph.circle_graph(brackets_1,sqrt_prob, color_0=color_1, gradient_0=gradient_0,
                                                      gradient_color_graph=gradient_color_graph, color_spectrum=color_spectrum,
                                                      mutation_imposition=True, first=True, len_input=0, list_nucleotides=list_nucleotides,
                                                      nucleotides=args.nucleotides, alpha=1)




            plt.xlim(-1.45*(len_index_input/(2*math.pi)),1.45*(len_index_input)/(2*math.pi))
            plt.ylim(-1.45*(len_index_input / (2 * math.pi)), 1.45*(len_index_input) / (2 * math.pi))
            overlap_patch = mpatches.Patch(color=color_3, label='OVERLAP')
            plt.legend(handles=[bond_zero, bond_one, overlap_patch], loc='upper right', fontsize='xx-small')


            plt.savefig('plot_0000000000000000_1.png', dpi=my_dpi)
            plt.clf()

            circle_graph.circle_graph(brackets_2,sqrt_prob, color_0=color_2, gradient_0=gradient_0,
                                      gradient_color_graph=gradient_color_graph, color_spectrum=color_spectrum,
                                      mutation_imposition=True, first=False, len_input=len_index_input,
                                      list_nucleotides=list_nucleotides, nucleotides=args.nucleotides, alpha=1)

            plt.xlim(-1.45*(len_index_input / (2 * math.pi)), 1.45*(len_index_input) / (2 * math.pi))
            plt.ylim(-1.45*(len_index_input / (2 * math.pi)), 1.45*(len_index_input) / (2 * math.pi))
            overlap_patch = mpatches.Patch(color=color_3, label='OVERLAP')
            plt.legend(handles=[bond_zero, bond_one, overlap_patch], loc='upper right', fontsize='xx-small')
            plt.savefig('plot_0000000000000000_2.png', dpi=my_dpi)

            image = functions.clean_imposition(list(wc.name_to_rgb(color_1)),
                                               list(wc.name_to_rgb(color_2)),
                                               list(wc.name_to_rgb(color_3)),
                                               'plot_0000000000000000_1.png',
                                               'plot_0000000000000000_2.png')

            os.remove('plot_0000000000000000_1.png')
            os.remove('plot_0000000000000000_2.png')

            image.show()
            save_question = input('Would you like to save the image? (Y/n)')
            if save_question == 'Y':
                name_3 = (str(name_1)).strip() + '&' + (str(name_2)).strip() + '_circular.png'
                image.save(name_3)

        if not args.match_radii and not args.align and not args.second_alignment_path:

            circle_graph.circle_graph(brackets, sqrt_prob, color_0=color_1, gradient_0=gradient_0,
                                                        gradient_color_graph=gradient_color_graph,
                                                        color_spectrum=color_spectrum,
                                                        mutation_imposition=False, first=True, len_input=0,
                                                        list_nucleotides=list_nucleotides, nucleotides=args.nucleotides,
                                                        alpha=1)



            brackets_2 = functions.bracket_reader(args.secondary_input)

            name_1 = functions.name(args.input)
            name_2 = functions.name(args.secondary_input)

            bond_zero = mpatches.Patch(color=color_1, label=name_1)
            bond_one = mpatches.Patch(color=color_2, label=name_2)

            max_brackets = max(len(brackets), len(brackets_2))

            plt.xlim(-1.45 * (max_brackets / (2 * math.pi)), 1.45 * (max_brackets) / (2 * math.pi))
            plt.ylim(-1.45 * (max_brackets / (2 * math.pi)), 1.45 * (max_brackets) / (2 * math.pi))

            overlap_patch = mpatches.Patch(color=color_3, label='OVERLAP')

            plt.legend(handles=[bond_zero, bond_one, overlap_patch], loc='upper right', fontsize='xx-small')

            plt.savefig('plot_0000000000000000_1.png', dpi=my_dpi)
            plt.clf()



            circle_graph.circle_graph(brackets_2, sqrt_prob, color_0=color_2, gradient_0=gradient_0,
                                      gradient_color_graph=gradient_color_graph, color_spectrum=color_spectrum,
                                      mutation_imposition=False, first=True, len_input=0, list_nucleotides=list_nucleotides,
                                      nucleotides=args.nucleotides, alpha=1)

            plt.xlim(-1.45 * (max_brackets / (2 * math.pi)), 1.45 * max_brackets / (2 * math.pi))
            plt.ylim(-1.45 * (max_brackets / (2 * math.pi)), 1.45 * max_brackets / (2 * math.pi))

            plt.legend(handles=[bond_zero, bond_one, overlap_patch], loc='upper right', fontsize='xx-small')

            plt.savefig('plot_0000000000000000_2.png', dpi=my_dpi)

            plt.legend(handles=[bond_zero, bond_one])

            image = functions.clean_imposition(list(wc.name_to_rgb(color_1)),
                                               list(wc.name_to_rgb(color_2)),
                                               list(wc.name_to_rgb(color_3)),
                                               'plot_0000000000000000_1.png',
                                               'plot_0000000000000000_2.png')

            os.remove('plot_0000000000000000_1.png')
            os.remove('plot_0000000000000000_2.png')

            image.show()
            save_question = input('Would you like to save the image? (Y/n)')
            if save_question == 'Y':
                name_3 = (str(name_1)).strip() + '&' + (str(name_2)).strip() + '_concentric_circle.png'
                image.save(name_3)

        elif args.align:

            match_brackets, first_idiosyncratic_brackets, second_idiosyncratic_brackets, first_alignment, second_alignment \
                = functions.align_redefine(
                functions.nucleotide_string(args.input), functions.nucleotide_string(args.secondary_input),

                functions.bracket_reader(args.input), functions.bracket_reader(args.secondary_input))

            list_nucleotides = functions.aligned_nucleotide_list(first_alignment, second_alignment)

            circle_graph.circle_graph(first_idiosyncratic_brackets, sqrt_prob, color_0=color_1, gradient_0=gradient_0,
                                      gradient_color_graph=gradient_color_graph,
                                      color_spectrum=color_spectrum,
                                      mutation_imposition=False, first=True, len_input=0,
                                      list_nucleotides=list_nucleotides, nucleotides=args.nucleotides, alpha=underlay_alpha)
            circle_graph.circle_graph(second_idiosyncratic_brackets, sqrt_prob, color_0=color_2, gradient_0=gradient_0,
                                      gradient_color_graph=gradient_color_graph,
                                      color_spectrum=color_spectrum,
                                      mutation_imposition=False, first=True, len_input=0,
                                      list_nucleotides=list_nucleotides, nucleotides=args.nucleotides, alpha=overlay_alpha)
            circle_graph.circle_graph(match_brackets, sqrt_prob, color_0=color_3, gradient_0=gradient_0,
                                      gradient_color_graph=gradient_color_graph,
                                      color_spectrum=color_spectrum,
                                      mutation_imposition=False, first=True, len_input=0,
                                      list_nucleotides=list_nucleotides, nucleotides=args.nucleotides, alpha=1)

            name_1 = functions.name(args.input)
            bond_zero = mpatches.Patch(color=color_1, label=name_1)
            name_2 = functions.name(args.secondary_input)
            bonds_one = mpatches.Patch(color=color_2, label=name_2)
            overlap_patch = mpatches.Patch(color=color_3, label='OVERLAP')
            plt.legend(handles=[bond_zero, bonds_one, overlap_patch], loc='upper right', fontsize='xx-small')
            plt.yticks([])
            plt.show()

        if not args.match_radii and not args.align and args.second_alignment_path:

            match_brackets, first_idiosyncratic_brackets, second_idiosyncratic_brackets, first_alignment, second_alignment \
                = functions.align_redefine(
                functions.nucleotide_string(args.input), functions.nucleotide_string(args.secondary_input),

                functions.bracket_reader(args.input), functions.bracket_reader(args.secondary_input))

            list_nucleotides = functions.aligned_nucleotide_list(first_alignment, second_alignment)

            circle_graph.circle_graph(first_idiosyncratic_brackets, sqrt_prob, color_0=color_1, gradient_0=gradient_0,
                                      gradient_color_graph=gradient_color_graph,
                                      color_spectrum=color_spectrum,
                                      mutation_imposition=False, first=True, len_input=0,
                                      list_nucleotides=list_nucleotides, nucleotides=first_alignment,
                                      alpha=1)

            circle_graph.circle_graph(match_brackets, sqrt_prob, color_0=color_3, gradient_0=gradient_0,
                                      gradient_color_graph=gradient_color_graph,
                                      color_spectrum=color_spectrum,
                                      mutation_imposition=False, first=True, len_input=0,
                                      list_nucleotides=list_nucleotides, nucleotides=args.nucleotides, alpha=1)

            brackets_2 = functions.bracket_reader(args.secondary_input)

            name_1 = functions.name(args.input)
            name_2 = functions.name(args.secondary_input)

            bond_zero = mpatches.Patch(color=color_1, label=name_1)
            bond_one = mpatches.Patch(color=color_2, label=name_2)

            max_brackets = max(len(brackets), len(brackets_2))

            plt.xlim(-1.45 * (max_brackets / (2 * math.pi)), 1.45 * (max_brackets) / (2 * math.pi))
            plt.ylim(-1.45 * (max_brackets / (2 * math.pi)), 1.45 * (max_brackets) / (2 * math.pi))

            alignment_patch = mpatches.Patch(color=color_3, label='ALIGNMENT')
            overlap_patch = mpatches.Patch(color=color_4, label='OVERLAP')

            plt.legend(handles=[bond_zero, bond_one, alignment_patch, overlap_patch], loc='upper right', fontsize='xx-small')

            plt.savefig('plot_0000000000000000_1.png', dpi=my_dpi)
            plt.clf()

            circle_graph.circle_graph(second_idiosyncratic_brackets, sqrt_prob, color_0=color_2, gradient_0=gradient_0,
                                      gradient_color_graph=gradient_color_graph, color_spectrum=color_spectrum,
                                      mutation_imposition=False, first=True, len_input=0,
                                      list_nucleotides=list_nucleotides,
                                      nucleotides=args.nucleotides, alpha=1)
            circle_graph.circle_graph(match_brackets, sqrt_prob, color_0=color_3, gradient_0=gradient_0,
                                      gradient_color_graph=gradient_color_graph,
                                      color_spectrum=color_spectrum,
                                      mutation_imposition=False, first=True, len_input=0,
                                      list_nucleotides=list_nucleotides, nucleotides=args.nucleotides, alpha=1)

            plt.xlim(-1.45 * (max_brackets / (2 * math.pi)), 1.45 * max_brackets / (2 * math.pi))
            plt.ylim(-1.45 * (max_brackets / (2 * math.pi)), 1.45 * max_brackets / (2 * math.pi))

            plt.legend(handles=[bond_zero, bond_one, alignment_patch, overlap_patch], loc='upper right', fontsize='xx-small')

            plt.savefig('plot_0000000000000000_2.png', dpi=my_dpi)

            plt.legend(handles=[bond_zero, bond_one])

            image = functions.clean_imposition(list(wc.name_to_rgb(color_1)),
                                               list(wc.name_to_rgb(color_2)),
                                               list(wc.name_to_rgb(color_4)),
                                               'plot_0000000000000000_1.png',
                                               'plot_0000000000000000_2.png')

            os.remove('plot_0000000000000000_1.png')
            os.remove('plot_0000000000000000_2.png')

            image.show()
            save_question = input('Would you like to save the image? (Y/n)')
            if save_question == 'Y':
                name_3 = (str(name_1)).strip() + '&' + (str(name_2)).strip() + '_concentric_circle.png'
                image.save(name_3)

    else:
        if args.circular and not args.MSA_file and not args.structures:

            circle_graph.circle_graph(brackets,sqrt_prob, color_0=color_1, gradient_0=gradient_0,
                                      gradient_color_graph=gradient_color_graph, color_spectrum=color_spectrum,
                                      mutation_imposition=False, first=True, len_input=0,
                                      list_nucleotides=list_nucleotides, nucleotides=args.nucleotides, alpha=1)

            name = functions.name(args.input)
            plt.title(name)
            bracket_length=len(brackets)

            if gradient_0 and args.gradient_legend:
                attempted_colormap = [color_spectrum[-1], color_spectrum[round(len(color_spectrum) / 2)], color_spectrum[0]]
                my_cmap = LinearSegmentedColormap.from_list('my_cmap', attempted_colormap)

                plt.scatter([-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11], [1000000, 1000000, 1000000, 1000000,
                                                                             1000000, 1000000, 1000000, 1000000,
                                                                             1000000, 1000000, 1000000],
                            c=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], cmap=my_cmap)
                plt.colorbar(shrink=0.5, label='bond probability (0.0 - 1.0)')


            plt.xlim(-1.25 * (bracket_length / (2 * math.pi)), 1.25 * bracket_length / (2 * math.pi))
            plt.ylim(-1.25 * (bracket_length / (2 * math.pi)), 1.25 * bracket_length / (2 * math.pi))
            plt.show()

        elif args.circular and args.MSA_file and args.structures:

            serial_brackets = [[functions.name(i), functions.bracket_reader(i)] for i in args.structures]
            bracket_length=len(serial_brackets[0][1])
            serial_alignment = functions.MSA_alignment_extraction(args.MSA_file)
            serial_names = [i[0] for i in serial_brackets]
            serial_names = ', '.join(serial_names)
            bond_list_backbone, bond_list_frequency, indices, low_value = graph_phylogeny.alignment_bracket_processing(
                serial_brackets,
                serial_alignment, n)
            graph_phylogeny.circle_graph_mod(brackets=serial_brackets[0][1], bonds=bond_list_backbone,
                                             color_spectrum=color_spectrum, serial_nucleotides=serial_alignment,
                                             nucleotides=args.nucleotides, colored_nucleotides=colored_nucleotides,
                                             bond_ubiquity=bond_list_frequency, nucleotide_color_list=nucleotide_colors)

            if args.gradient_legend:
                attempted_colormap = [color_spectrum[-1], color_spectrum[round(len(color_spectrum) / 2)], color_spectrum[0]]
                my_cmap = LinearSegmentedColormap.from_list('my_cmap', attempted_colormap)

                if n ==len(serial_alignment)-1:
                    color_list = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

                else:
                    color_list = list(np.linspace(low_value, 1.0, 11))

                plt.scatter([-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11], [1000000, 1000000, 1000000, 1000000,
                                                                             1000000, 1000000, 1000000, 1000000,
                                                                             1000000, 1000000, 1000000],
                            c=color_list, cmap=my_cmap)
                plt.colorbar(shrink=0.5, label='bond conservation probability ('+ str(low_value) +'- 1.0)')
            plt.title(serial_names)

            plt.xlim(-1.25 * (bracket_length / (2 * math.pi)), 1.25 * bracket_length / (2 * math.pi))
            plt.ylim(-1.25 * (bracket_length / (2 * math.pi)), 1.25 * bracket_length / (2 * math.pi))

            plt.show()

if __name__ == '__main__':
    main()



