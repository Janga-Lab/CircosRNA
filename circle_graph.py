import functions
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch

def circle_graph(brackets, sqrt_prob, color_0, gradient_0, gradient_color_graph, color_spectrum, mutation_imposition,
                 first, len_input, list_nucleotides, nucleotides, alpha):




        bonds=functions.bond_finder(brackets)

        if not mutation_imposition and first:
            len_index=len(list(range(0,len(brackets)-1)))

        elif not first:
            len_index=len_input
        else:
            len_index = len(list(range(0, len(brackets) - 1)))

        abcissa_0=[]
        ordinate_0=[]

        circle_list_angles=list(np.linspace(0,2*math.pi,3000))


        all_coordinates=[]

        for i in bonds:
            coordinates=[]

            theta_0 = ((i[0] / len_index) * 2 * math.pi) + (math.pi/2)
            theta_1=((i[1]/len_index)*2*math.pi)+(math.pi/2)


            coordinates.append(functions.polar_to_cartesian((len_index/(2*math.pi)),theta_0))
            coordinates.append(theta_0)
            coordinates.append(functions.polar_to_cartesian((len_index/(2*math.pi)),theta_1))
            coordinates.append(theta_1)
            all_coordinates.append(coordinates)

        regression_plots=[]

        for i in all_coordinates:

            theta_mean=functions.middle_man(i[1],i[3])[0]

            r=(((math.pi-functions.middle_man(i[1],i[3])[1])/(math.pi)))*(len_index/(2*math.pi))

            i.pop(1)
            i.pop(2)
            i.append(functions.polar_to_cartesian(r,theta_mean))

            abcissa=[]
            ordinate=[]

            for ii in i:
                abcissa.append(ii[0])
                ordinate.append(ii[1])
            regression_plots.append([abcissa,ordinate])


        if gradient_0==True:
            for i,j in enumerate(all_coordinates):
                regression_plots[i].append(sqrt_prob[i])

        for i in regression_plots:
            x_m=(i[0][0] + i[0][1])/2
            y_m=(i[1][0] + i[1][1])/2
            l=math.sqrt(((x_m-i[0][0])**2)+((y_m-i[1][0])**2))
            h = math.sqrt(((x_m - i[0][2]) ** 2) + ((y_m - i[1][2]) ** 2))

            abcissa=(-l,0,l)
            ordinate=(0,h,0)
            polynomial=list(np.polyfit(abcissa, ordinate, 2))

            x_vector=list(np.linspace(-l,l,1000))
            y_vector=[((ii**2)*polynomial[0])+(ii*polynomial[1])+polynomial[2] for ii in x_vector]

            maxima_index=y_vector.index(max(y_vector))


            transformed_vectors_default = functions.transrotation(x_vector, y_vector, i[1][0] - i[1][1], i[0][0] - i[0][1], flip_it=False, iteration=i)

            transformed_vectors_test = functions.transrotation(x_vector, y_vector, i[1][0] - i[1][1], i[0][0] - i[0][1], flip_it=True, iteration=i)

            d_default=math.sqrt(((transformed_vectors_default[0][maxima_index])**2)+((transformed_vectors_default[1][maxima_index])**2))

            d_test=math.sqrt(((transformed_vectors_test[0][maxima_index])**2)+((transformed_vectors_test[1][maxima_index])**2))

            if d_test < d_default:

                transformed_vectors=transformed_vectors_test
            else:

                transformed_vectors=transformed_vectors_default



            if not gradient_0 and not gradient_color_graph:
                plt.plot(transformed_vectors[0], transformed_vectors[1], color_0, zorder=0, alpha=alpha)


            if gradient_0 and not gradient_color_graph:
                plt.plot(transformed_vectors[0], transformed_vectors[1],color=str(1-i[-1]), zorder=0)

            try:
                if gradient_0 and gradient_color_graph:

                    plt.plot(transformed_vectors[0], transformed_vectors[1],
                        color=tuple(color_spectrum[round(((1 - i[-1]) * len(color_spectrum)))]), zorder=0)
            except:
                pass


        if not mutation_imposition:
            start=math.pi/2
            jump=(2*math.pi)/6
            r=len_index/(2*math.pi)
            for i in list(range(0,6)):
                plt.text((1.15*r)*math.cos(start+i*jump),(1.10*r)*math.sin(start+i*jump), str(round((i/6)*len_index)), horizontalalignment='center')

            plt.text((1.15 * r) * math.cos(start), (1.20 * r) * math.sin(start),
                     str(len_index
                         ), horizontalalignment='center')

        for i in circle_list_angles:
            abcissa_0.append(functions.polar_to_cartesian(((len_index)/(2*math.pi)),i)[0])
            ordinate_0.append(functions.polar_to_cartesian(((len_index)/(2*math.pi)),i)[1])

        if not nucleotides:
            plt.plot(abcissa_0,ordinate_0, 'black')

        if nucleotides:
            index_coordinates = [
                functions.polar_to_cartesian((len_index / (2 * math.pi)), (i / len_index) * 2 * math.pi + (1 / 2) * math.pi)
                for i in range(0, len_index)]

            for i in range(0,len_index):

                if '/' in list_nucleotides[i]:
                    tp = TextPath((index_coordinates[i][0]-0.15, index_coordinates[i][1]), list_nucleotides[i], size=0.4)
                else:
                    tp = TextPath((index_coordinates[i][0], index_coordinates[i][1]), list_nucleotides[i], size=0.4)
                nucleotide = plt.Circle((index_coordinates[i][0] + 0.15, index_coordinates[i][1] + 0.15), 0.5,
                                            edgecolor='black', facecolor='white')
                ax = plt.gca()
                ax.add_patch(nucleotide)
                ax.add_patch(PathPatch(tp, color="black"))

        plt.axis('off')

        plt.axis('square')


        if first:
            return len_index
