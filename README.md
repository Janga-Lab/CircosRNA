Copyright 2021 Janga Lab for Genomics and Systems Biology

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Welcome to RNACircos

INSTALLATION

Remember to change the permisions as needed to use the RNACircos.py script.

No installation is needed, but the matplotlib, pillow, webcolors, and biopython library must be installed.

arch, debian, and ubuntu all seem amendable to downloading this way:

	pip install {biopython,pillow,matplotlib,webcolors}
or
	
	pip3 install {biopython,pillow,matplotlib,webcolors}

Any recent version should work, of course, if you're out of luck, here are detailed instructions for installation for each package:

matplotlib
https://matplotlib.org/stable/users/installing.html

biopython
https://biopython.org/wiki/Packages

webcolors
https://webcolors.readthedocs.io/en/1.5/install.html

pillow
https://pillow.readthedocs.io/en/stable/installation.html

To view the command line options, use:

	./RNACircos.py -h
or
	
	python3 RNACircos.py -h
or
	
	python RNACircos.py -h

Example files have been furnished in the example_file folder included with this CLI

If you are reluctant to execute all of these commands individually, A bash file is included to run them all for you entitled RNACircos_test.sh. Simply change the permission and execute ./RNA_graph_test.sh

If you want to view all colors available to you for this CLI, they are listed in this link:
https://matplotlib.org/stable/gallery/color/named_colors.html

![image](https://media.github.iu.edu/user/17625/files/41b40d00-0672-11ec-864a-83e080b0bc7a)

HELP MENU

usage: RNACircos.py [-h] [-i INPUT] [-i2 SECONDARY_INPUT] [-l] [-c] [-c1 COLOR_ONE] [-c2 COLOR_TWO] [-c3 COLOR_THREE] [-c4 COLOR_FOUR] [-a] [-s] [-o OVERLAY_ALPHA] [-u UNDERLAY_ALPHA] [-m] [-pm P_MATRIX_INPUT]
                    [-lc LOW_PROB_COLOR [LOW_PROB_COLOR ...]] [-hc HIGH_PROB_COLOR [HIGH_PROB_COLOR ...]] [-g] [-n] [-d DPI]

CLI for graphing of dotbracket notation of RNA secondary structure

optional arguments:
  
  -h, --help            show this help message and exit
  
  -i INPUT, --input INPUT
                        input file for the dot-bracket structure
  
  -i2 SECONDARY_INPUT, --secondary_input SECONDARY_INPUT
                        secondary input file, also signals the superimposition of one linear/circular representation onto another
  
  -l, --linear          Produces the linear representation of RNA secondary structure
  
  -c, --circular        Produces the circular representation of RNA secondary structure
  
  -c1 COLOR_ONE, --color_one COLOR_ONE
                        selected color for the plot, the default is lightblue
  
  -c2 COLOR_TWO, --color_two COLOR_TWO
                        selected color for the superimposed plot, the default is lightgreen
  
  -c3 COLOR_THREE, --color_three COLOR_THREE
                        When graphs are superimposed, the overlapping areas should be set to a seperate color, something which contrasts well is recommended
  
  -c4 COLOR_FOUR, --color_four COLOR_FOUR
                        overlap color of unaligned regions if -a2 is chosen
  
  -a, --align           Align the nucleotides before checking for structural similarities (recommended)
  
  -s, --second_alignment_path
                        align and permit the overlaping regions to be some fourth color of your choice
  
  -o OVERLAY_ALPHA, --overlay_alpha OVERLAY_ALPHA
                        transparency value between 0 and 1 for the overlying plot in superimposed graphs
  
  -u UNDERLAY_ALPHA, --underlay_alpha UNDERLAY_ALPHA
                        transparency value between 0 and 1 for the underlying plot in superimposed graphs
  
  -m, --match_radii     by default, circular representations of secondary structure will adapt to polymer length, including this argument will cause the circular graphs to adopt uniform radii
  
  -pm P_MATRIX_INPUT, --p_matrix_input P_MATRIX_INPUT
                        required input *dp.ps containing the ViennaRNA generated probability matrix. Triggers gradient generation
  
  -lc LOW_PROB_COLOR [LOW_PROB_COLOR ...], --low_prob_color LOW_PROB_COLOR [LOW_PROB_COLOR ...]
                        add the low rgb values for the custom gradient
  
  -hc HIGH_PROB_COLOR [HIGH_PROB_COLOR ...], --high_prob_color HIGH_PROB_COLOR [HIGH_PROB_COLOR ...]
                        add the high rgb values for the custom gradient
  
  -g, --gradient_legend
                        adds a legend to gradient graphs to show which color corresponds to a low probability and which color coresponds to a high probability
  
  -n, --nucleotides     adds nucleotides to the visualization
  
  -d DPI, --dpi DPI     enter the dpi needed for supderimposed graphs, there is no "one size fits all", raise or lower this value as needed, start with 96 if in doubt



BASIC USAGE

Refer to the help menu for flag meaning

To generate a linear or circular structure, the input should be as follows:

	./RNACircos.py -l -i example_files/input_dot_bracket.txt

![Figure_1](https://media.github.iu.edu/user/17625/files/4b03b480-11f0-11ec-8ba0-cb1ba97a1295)

	./RNACircos.py -c -i example_files/input_dot_bracket.txt
	
![Figure_2](https://media.github.iu.edu/user/17625/files/548d1c80-11f0-11ec-8da9-d629b18a000c)

	

If one wishes to add nucleotides, add the -n flag.

	./RNACircos.py -l -n -i example_files/input_dot_bracket.txt
	
![Figure_1 25](https://media.github.iu.edu/user/17625/files/5fe04800-11f0-11ec-99a7-843059700d20)

	./RNACircos.py -c -n -i example_files/input_dot_bracket.txt
	
![Figure_1 75](https://media.github.iu.edu/user/17625/files/c9605680-11f0-11ec-98ab-2bf70cc33d9f)


the input file must be of the format:

>specimen_name
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
......((((............))).....).......(......)........(..(..

Where N denotes a nucleotide from RNA, A, U, C, and G.

To generate a gradient, a *dp.ps file containing a probability matrix from ViennaRNA's RNAfold is needed. One such file has been provided in the example_file if the user is unwilling to generate their own. In the absence of SHAPE data, the installed ViennaRNA command is

	RNAfold -p input.fasta

If the user wishes to generate a *dp.ps file using SHAPE data, the command in the intalled version of ViennaRNA is:

	RNAfold -p --shape=reacitivites.txt < input.fasta

In this command, the reactivities are organized into two vertical columns, nucleotide number and normalized reactivity score.

A sophisticated cluster is not needed to implament ViennaRNA. Using one's personal device ought to suffice. In fact, this program requires a GUI, it should run well on anything with an i3 processor, equivalent, or better, but some computations will require waiting a minute or two.

GRADIENT GENERATION

To implament the probability matrix within the *dp.ps file in this tool, use the following command.

    ./RNACircos.py -l -i example_files/gradient_dot_bracket.txt -pm example_files/example_dp.ps

![bonus_Figure_1](https://media.github.iu.edu/user/17625/files/ecd8d080-11f3-11ec-9069-8754bcbf4458)
    
    ./RNACircos.py -c -i example_files/gradient_dot_bracket.txt -pm example_files/example_dp.ps

![bonus)Figure_2png](https://media.github.iu.edu/user/17625/files/f5310b80-11f3-11ec-8f6a-c15c5662cd26)

By default, the output upon the addition of a ViennaRNA generated probability matrix will be grayscale going from 255 255 255 to 0 0 0. This is perfect for testing out the tool and considering a given output, but lacks some of the features which would communicate the meaning of the plot. Adding a gradient legend is not the default, but adding another flag will allow the user to generate an appropriate legend. When the -g flag for the gradient legend is called, a -lc and -hc must be specified. The command should look like this:

	./RNACircos.py -c -i example_files/gradient_dot_bracket.txt -pm example_files/example_dp.ps -g -lc white -hc black
	
![Figure_3](https://media.github.iu.edu/user/17625/files/e006ad80-11f0-11ec-965b-2a7957f393cf)

	./RNACircos.py -l -i example_files/gradient_dot_bracket.txt -pm example_files/example_dp.ps -g -lc white -hc black

![Figure_4](https://media.github.iu.edu/user/17625/files/560b1480-11f1-11ec-8c14-0236534abfea)

This command will still run if you enter the wrong path, so take care to ensure your path is correct.


The superimposition of two graphs is an interesting, novel application of linear and circular visualizations, so by design the flags or fine tuned for specific scenarios and user control. Two solutions were used to compare homologous RNA specimens in this work, which are described in greater detail in the article corresponding to this tool. One solution allows the user to retain matplotlib's flexible environment, the other allows the user to more precisely specify graph idiosynrasies. Let's consider all conceivable scenarios briefly in this README:

UNALIGNED IMPOSITION

If one wishes to simply take two plots, and overlap them with one another, these commands will suffice:

    ./RNACircos.py -l -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA

![Figure_5](https://media.github.iu.edu/user/17625/files/28be6680-11f1-11ec-9015-8e6d1950d61f)

    
    ./RNACircos.py -c -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA
    
![Figure_6](https://media.github.iu.edu/user/17625/files/307e0b00-11f1-11ec-9102-053c68861a09)


One of the drawbacks is the differences are not immediately observable, furthermore, this image does not adjust as the user zooms in and will be pixilated should the user try. Nucleotides cannot be displayed using this method.

MATPLOTLIB BASED ALIGNMENT

It is difficult to appreciate structural similarities this way, so the user can add an -a flag for an alignment. This addition alone to the previous command allows the user to retain matplotlib's functionalities, such as resizing the window and zooming in and out to view nucleotides, which may be shown with the -n flag.

    ./RNACircos.py -l -a -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA
    
![Figure_7](https://media.github.iu.edu/user/17625/files/3b38a000-11f1-11ec-83fa-b2c6edcc589a)


    ./RNACircos.py -c -a -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA
    
![Figure_8](https://media.github.iu.edu/user/17625/files/47246200-11f1-11ec-8ada-88c0498516eb)

or

    ./RNACircos.py -l -a -n -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA
    
![Figure_9](https://media.github.iu.edu/user/17625/files/63280380-11f1-11ec-9dcf-c9de617b0392)

    ./RNACircos.py -c -a -n -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA
    
![Figure_10](https://media.github.iu.edu/user/17625/files/6ae7a800-11f1-11ec-8924-efd34e53a879)


By default, matplotlib will show the alignment, but unaligned overlapping regions will display the second input over the first input plot. This can be resolved by calling the -o flag to determine the overlay alpha. This will make the second layer transparent and the underlying layer will show more easily, a value between 0 and 1 must be input, but the best results seem to be between 0.15 and 0.30. The user must toggle these values to find something appropriate for their specific plot. The user can also define the underlying transparency; however, it is recommended. If one must change the underlying transparency, the flag is -u, the value also is between 0 and 1.

    ./RNACircos.py -l -a -n -o 0.15 -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA
    
![Figure_11](https://media.github.iu.edu/user/17625/files/7044f280-11f1-11ec-98ba-893d7fb21894)
    

    ./RNACircos.py -c -a -n -o 0.15 -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA
    
![Figure_12](https://media.github.iu.edu/user/17625/files/8d79c100-11f1-11ec-8ba3-3003027420d5)


PILLOW BASED ALIGNMENT

Perhaps the reader feels this process is a little mundane, changing the transparency may make the plot lighter than desired. If the user simply must have a specific set of colors from the CSS color palate at a specific shade, a solution is offered using a similar scheme as the overlap plot. The alignment is again performed, but this time the user can specify precisely which color will define the unaligned overlapped regions. This is recommended if one wishes to generate publication quality visualizations, rather than experiment with alignments. The primary drawback is matplotlib's functionality is not retained. This can make resizing the visualization cumbersome. Resizing must be performed through toggling the DPI. The DPI, or dots per inch. The size is contingent on the GUI screen resolution, so no single size will be suitable for all purposes. The user must toggle this value to find their preferred size. If the user wishes to use this process, the -s flag must be added. This command should look a bit like this:

    ./RNACircos.py -l -s -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA
    
![Figure_13](https://media.github.iu.edu/user/17625/files/99fe1980-11f1-11ec-9171-ee38e4f45237)

    ./RNACircos.py -c -s -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA
    

![Figure_14](https://media.github.iu.edu/user/17625/files/b4d08e00-11f1-11ec-8a1b-9a9cba023160)

The user is free to control the superimposition graph colors with the greatest autonomy here, the flags are -c1, -c2, -c3, and -c4. These flags have different meanings in the matplotlib based plot or the unaligned imposition. By default, the imposition graph colors are lightblue, lightgreen, pink, and purple. Say the user wishes to plot the first and second specimens as blue and red, while plotting any overlapping regions as purple and any aligned regions as indigo, the commands would be:

    ./RNACircos.py -l -s -c1 blue -i example_files/Dot_Bracket_India_B16171 -c2 red -i2 example_files/Spike_NCBI_MEA -c3 indigo -c4 purple
    
![Figure_15](https://media.github.iu.edu/user/17625/files/c154e680-11f1-11ec-86ad-288bb8fb2e45)

    ./RNACircos.py -c -s -c1 blue -i example_files/Dot_Bracket_India_B16171 -c2 red -i2 example_files/Spike_NCBI_MEA -c3 indigo -c4 purple.
   
![Figure_16](https://media.github.iu.edu/user/17625/files/83a48d80-11f2-11ec-9aa1-b799ea5dc3dc)

-c1, -c2, and -c3 can all be implamented in the previous superimposed graphs. In the MATPLOTLIB BASED ALIGNMENT and UNALIGNED IMPOSITION, -c3 denotes the alignment and imposition respectively. -c4 is only implamented in the PILLOW BASED ALIGNMENT.

Again, the pillow based plot will not adjust to zooming in and out, so nucleotides cannot be viewed.






