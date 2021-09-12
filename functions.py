#!/bin/python3

import math

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.image as img

from PIL import Image
import webcolors as wc

import numpy as np
from Bio import pairwise2


#Description: Converts the dot-bracket input into just the brackets
#Input: dot-bracket file
#Output: brackets within file
def bracket_reader(path_to_brackets):
    fh=open(path_to_brackets)
    fh=fh.readlines()

    bracket_iteration=0
    for i,j in enumerate(fh):
        for jj in list(j):
            if jj=='(' or jj==')' or jj=='.':
                bracket_iteration+=1
        if bracket_iteration/len(j) > 0.9:
            brackets=j
            break
        else:
            brackets='null'
    brackets=brackets.split(' ')[0]
    return brackets

#Description: uses *dp.ps file to generate a list containing the 'nucleotide_index_1, nucleotide_index_2,
    #probability ubox' strings
#Input: path to *dp.ps file
#Output: list containing 'nucleotide_index_1, nucleotide_index_2, probability ubox'

def prob_matrix(path_to_prob_matrix):
        fh=open(path_to_prob_matrix)
        fh=fh.readlines()
        fh=[i for i in fh if 'ubox' in i]
        return fh

#Description: Converts a string of brackets into a set of numbers and converts into a list which pairs the first
    #and second bonds
#Input: String of brackets in the form '(.(..)...)'
#Output: List of lists containing the nucleotide indices corresponding to the hydrogen bond in the format
    #[[n_1_1,n_1_2],[n_2_1,n_2_2].........]


def bond_finder(brackets):
    tuple_list=list(enumerate(brackets))
    bonds=[]
    parentheses_indices=[]
    tuple_pop=False

    while True:

        for i in tuple_list:

            if i[1] == '(':
                bond=[]
                for tuple_bracket in tuple_list[tuple_list.index(i)+1:-1]:
                    if tuple_bracket[1] == '(':
                        break
                    elif tuple_bracket[1] == ')':
                        bond.append(i[0])
                        bond.append(tuple_bracket[0])
                        bonds.append(bond)
                        parentheses_indices.append(bond[0])
                        parentheses_indices.append(bond[1])

                        tuple_pop=True
                        break
                if tuple_pop:
                    break

        if tuple_pop:
            tuple_list.remove(tuple([bond[0],'(']))
            tuple_list.remove(tuple([bond[1],')']))
            tuple_pop=False

        pseudo_brackets=''.join([i for j, i in enumerate(brackets) if j not in parentheses_indices])
        if '(' and ')' not in pseudo_brackets:
            break
    return bonds




#Description: Breaks rgb values in a specified rage and generates a spectrum going from the intial rgb value, to
    #the final rgb value entered
#Input: xy_size is a tuple containing two strings: xy_size=('x','y'), rgb_start is the low probability rgb
    #values as a tuple containing those numbers as a string ('r','g','b'), rgb_end are the high probability
    #rgb values, also as a tuple, ('r','g','b')
#Output: list of lists containing rgb values as three seperate elements in each internal list

def gen_gradient(xy_size,rgb_start,rgb_end):

        xy_size = list(map(int, xy_size))  # convert the function's arguments
        rgb_start = list(map(int, rgb_start))  # //
        rgb_end = list(map(int, rgb_end))

        r_gap = (rgb_end[0] - rgb_start[0]) / xy_size[0]  # calculate the gap of the "R" value for every column
        g_gap = (rgb_end[1] - rgb_start[1]) / xy_size[0]  # calculate the gap of the "G" value for every column
        b_gap = (rgb_end[2] - rgb_start[2]) / xy_size[0]

        color_spectrum=[[int(rgb_start[0]+r_gap*x)/255,
                            int(rgb_start[1]+g_gap*x)/255,
                            int(rgb_start[2]+b_gap*x)/255] for x in range(xy_size[0])]
        return color_spectrum




#Description: plots the linear graph given the bonds

def bond_grapher(abcissa, linear, bonds, color, gradient_0, gradient_color_graph, color_spectrum, mutation_imposition, alpha_0):

        for i in bonds:
            s=(i[0]+i[1])/2
            r=abs(i[0]-i[1])/2
            final_ordinate=[]
            for ii in abcissa:
                try:
                    final_ordinate.append(math.sqrt((r**2)-(ii-s)**2))
                except:
                    final_ordinate.append(0)
            if gradient_0==True and linear and gradient_color_graph == False and not mutation_imposition:

                plt.plot(abcissa, final_ordinate,color=str(1-i[-1]))

            elif gradient_0==True and linear and gradient_color_graph and not mutation_imposition:
                plt.plot(abcissa, final_ordinate,
                         color=tuple(color_spectrum[round(((1 - i[-1]) * len(color_spectrum)))])
                         ,zorder=0)
            else:
                if linear:
                    plt.plot(abcissa,final_ordinate, color, zorder=0, alpha=alpha_0)

                else:
                    pass



#Description: Takes the bond numbers and uses them to determine the location of the sqrt(bond_probability) and builds
    #a new list of identical length
#Input: bonds, list of lists containing the nucleotide indices corresponding to the hydrogen bond in the format
    #[[n_1_1,n_1_2],[n_2_1,n_2_2].........]
    #fh, list containing 'nucleotide_index_1, nucleotide_index_2, probability ubox'
#Output: List of squared probabilities corresponding to H-bonds, it will be in the same order as the 'bonds' list

def sqrt_finder(bonds,fh):
        try:
            sqrt_prob=[]
            for i in bonds:
                j=0
                for ii in fh:
                    j+=1
                    if str(i[0]) and str(i[1]) in ii:
                        #Let's just let whether its squared be user specified
                        sqrt_prob.append(float(ii.split(' ')[2]))
                        break
            return sqrt_prob
        except:
            return 0
#Description: Converts polar coordinates to cartesian
#Input: distance from origin, angle relative to the positive x-axis
#Output: list containing the x and y coordinates

def polar_to_cartesian(r,theta):
        x=r*math.cos(theta)
        y=r*math.sin(theta)
        return [x,y]
#Description: Finds the angle between two angles and finds the difference between those two angles
#Input: Both of the angles in no particular order (probably)
#Ouput: A list containing the angle between the two input angles and the angle between them

def middle_man(theta_one, theta_two):

            if round(abs(theta_one-theta_two), 5) <= round(math.pi, 5):


                theta_diff=abs(theta_one-theta_two)

                return [(theta_one+theta_two)/2,theta_diff]

            else:

                theta_diff=(2*math.pi-abs(theta_one-theta_two))

                if theta_two>theta_one:
                    return [theta_one-(theta_diff/2), theta_diff]
                elif theta_two<theta_one:
                    return [theta_two-(theta_diff/2), theta_diff]


#Description: Takes the parabola generated from the quadratic regression, rotates it, and moves it along the cartesian
    #plane to the location consistent with the graph orientation based on nucleotide indices
#Input: parabola abcissa vector, parabola ordinate vector, difference between bottom or top portion of the new
    #vector and the difference between the leftmost or rightmost portion of the new vector as a basis for determining
    #the angle, and a boolean to indicate whether the final vector should be rotated another ninety degrees to
    #accomodate the arctan2 function
#Output: rotated and relocated vector coordinate abcissa and ordinate vectors

def transrotation(abcissa, ordinate, vector_diff_y, vector_diff_x, flip_it, iteration):

                abcissa_t=[]
                ordinate_t=[]
                for i in range(0,len(abcissa)):

                    d=math.sqrt(((abcissa[i])**2)+((ordinate[i])**2))
                    theta=np.arctan2(ordinate[i],abcissa[i])
                    transform_theta=np.arctan2(vector_diff_y, vector_diff_x)+math.pi

                    if flip_it:
                        transform_theta+=math.pi

                    abcissa_t.append(d*math.cos(theta+transform_theta))
                    ordinate_t.append(d*math.sin(theta+transform_theta))

                x_diff = max([abcissa_t[0], abcissa_t[-1]]) - max(iteration[0][0], iteration[0][1])
                y_diff = max([ordinate_t[0], ordinate_t[-1]]) - max(iteration[1][0], iteration[1][1])

                for ii in range(0, len(abcissa_t)):
                    abcissa_t[ii] -= x_diff
                    ordinate_t[ii] -= y_diff



                return [abcissa_t, ordinate_t]

#Description: find the name of the specimen within the txt file
def name(path_to_file):
    fh=open(path_to_file)
    fh=fh.readlines()
    name=fh[0]
    name=name.split(' ')[0]
    name_list=list(name)
    name_list.pop(0)
    name=''.join(name_list)
    return name

#Description: returns a list of nucleotides given a file containing the base-pairs
def nucleotide_list(path_to_file):
    fh=open(path_to_file)
    fh=fh.readlines()
    for i in fh:
        if 'A' and 'C' and 'U' and 'G' in i:
            list_nucleotides=i
    return list(list_nucleotides)

#Description: superimposes two graphs on one another and allows the user to specify the color of the overlapping layers
#Note:
    #The first function, 'times_255' , merely serves as a means of shortening the 'clean_imposition' function and will not be used
    #elsewhere
    #example input:
    #color_1=list(wc.name_to_rgb('lightblue'))
    #color_2=list(wc.name_to_rgb('lightgreen'))
    #color_3=list(wc.name_to_rgb('pink'))
    #clean_imposition(color_1,color_2,color_3, 'big_test_1.png', 'big_test_2.png')


def times_255(list_of_lists):
    for i,j in enumerate(list_of_lists):
        for ii, jj in enumerate(list_of_lists[i]):
            for iii, jjj in enumerate(list_of_lists[i][ii]):
                list_of_lists[i][ii][iii]=round(list_of_lists[i][ii][iii]*255)

    for i,j in enumerate(list_of_lists):
        for ii, jj in enumerate(list_of_lists[i]):
            list_of_lists[i][ii].pop(3)

    return list_of_lists

def clean_imposition(color_1, color_2, color_3, image_1, image_2):

    plot_one=img.imread(image_1)

    easy_list_1=plot_one.tolist()

    easy_list_1_255=times_255(easy_list_1)

    plot_two=img.imread(image_2)

    easy_list_2=plot_two.tolist()

    easy_list_2_255=times_255(easy_list_2)



    for i in list(range(0,len(easy_list_1_255))):
        for ii in list(range(0,len(easy_list_1_255[i]))):
            one_count_one=0
            one_count_two=0
            b_one_count_one=0
            b_one_count_two=0
            for iii in [0,1,2]:

                white_one=False
                white_two=False
                black_one=False
                black_two=False

                if easy_list_1_255[i][ii][iii] == 255:
                    one_count_one+=1
                    if one_count_one==3:
                        white_one=True

                if easy_list_2_255[i][ii][iii] == 255:
                    one_count_two+=1
                    if one_count_two==3:
                        white_two=True

                if easy_list_1_255[i][ii][iii] == 0:
                    b_one_count_one+=1
                    if b_one_count_one==3:
                        black_one=True

                if easy_list_2_255[i][ii][iii] == 0:
                    b_one_count_two+=1
                    if b_one_count_two==3:
                        black_two=True


            if easy_list_1_255[i][ii]==color_1 and easy_list_2_255[i][ii]==color_2:
                easy_list_1_255[i][ii]=color_3

            elif white_one and not white_two:
                easy_list_1_255[i][ii]=easy_list_2_255[i][ii]

            elif black_one and not black_two or not black_one and black_two:
                easy_list_1_255[i][ii]=[0,0,0]

    image_data=np.array(easy_list_1_255)

    image_data=image_data.astype(np.uint8)

    image=Image.fromarray(image_data)

    return image

#Description: returns a string of nucleotides from a provided dot_bracket file in the ViennaRNA format
#Input: path to dot-bracket file
#Output: string of nucleotides
def nucleotide_string(path_to_file):
    nucleotides=nucleotide_list(path_to_file)
    nucleotides.pop(-1)
    nuke_string=''.join(map(str, nucleotides))
    return nuke_string

#Description: takes the brackets and the alignment strings and returns the brackets with '.' in the idices of the '-'s
#Input: the brackets to be modified and the alignment to serve as the basis for modification
#Output: the modified brackets

def impure_brackets(brackets, alignment):
    dash_indices=[i for i,j in enumerate(list(alignment)) if j=='-']
    list_brackets=list(brackets)
    for i in dash_indices:
        list_brackets.insert(int(i),'.')
    string_brackets=''.join(list_brackets)
    return string_brackets

#Description: Converts bonds back into their original bracket string
#Input: the list of hydrogen bond indices and the length of the original brackets
#Output: Dot-bracket notation bonds

def bonds_to_brackets(bonds, len_of_original_brackets):
    empty_brackets=['.' for i in len_of_original_brackets]
    for i in bonds:
        empty_brackets[i[0]]='('
        empty_brackets[i[1]]=')'
    placeholder_brackets=''.join(empty_brackets)
    return placeholder_brackets

#Description: converts alignments into bracket notation, which may be converted into a visualizatoin which shows
    #matching bonds
#Input: the nucleotide strings to be aligned and both sets of brackets
#Output: Matching brackets, both sets of mismatching brackets, and both alignments

def align_redefine(first_nuke_string, second_nuke_string, first_brackets, second_brackets):

    my_alignments = pairwise2.align.globalxx(first_nuke_string, second_nuke_string)

    first_alignment_list = list(my_alignments[0][0])
    second_alignment_list = list(my_alignments[0][1])

    first_impure_brackets=impure_brackets(first_brackets,first_alignment_list)
    second_impure_brackets=impure_brackets(second_brackets,second_alignment_list)


    first_impure_bonds=bond_finder(first_impure_brackets)
    second_impure_bonds=bond_finder(second_impure_brackets)


    matches=[i for i in first_impure_bonds if i in second_impure_bonds]
    first_idiosyncratic=[i for i in first_impure_bonds if i not in matches]
    second_idiosyncratic=[i for i in second_impure_bonds if i not in matches]

    bracket_length=range(0, len(first_impure_brackets))
    match_brackets=bonds_to_brackets(matches, bracket_length)
    first_idiosyncratic_brackets=bonds_to_brackets(first_idiosyncratic, bracket_length)
    second_idiosyncratic_brackets=bonds_to_brackets(second_idiosyncratic, bracket_length)

    return [match_brackets, first_idiosyncratic_brackets, second_idiosyncratic_brackets,
            my_alignments[0][0],my_alignments[0][1]]

#Description: converts alignments into a list in which elements with matching indexes do not have matching strings
    #to a three element string with a '/' seperating the base-pair letters
#Input: pairwise2.align.globalxx outputs [0][0] and [0][1]
#Output: List of strings with a '/' between non-identical elements

def aligned_nucleotide_list(first_alignment,second_alignment):
    list_nucleotides = []
    for i in range(0, len(first_alignment or second_alignment)):
        if first_alignment[i] == second_alignment[i]:
            list_nucleotides.append(first_alignment[i])
        else:
            list_nucleotides.append(first_alignment[i] + '/' + second_alignment[i])
    return list_nucleotides

def height_finder(brackets_1, brackets_2):

    bonds_1=bond_finder(brackets_1)
    bonds_2=bond_finder(brackets_2)

    bonds=bonds_1+bonds_2
    max_bond = 0
    for i in bonds:
        if abs(i[0] - i[1]) > max_bond:
            max_bond = abs(i[0] - i[1])

    return max_bond



