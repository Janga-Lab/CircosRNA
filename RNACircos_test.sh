#!/bin/bash

./RNACircos.py -h		

echo //////////////////////////////////////////////////////////////
echo help 
echo //////////////////////////////////////////////////////////////

./RNACircos.py -l -i example_files/input_dot_bracket.txt

echo //////////////////////////////////////////////////////////////
echo simple linear graph
echo //////////////////////////////////////////////////////////////

./RNACircos.py -c -i example_files/input_dot_bracket.txt

echo //////////////////////////////////////////////////////////////
echo simple circular graph
echo //////////////////////////////////////////////////////////////

./RNACircos.py -l -n -i example_files/input_dot_bracket.txt

echo //////////////////////////////////////////////////////////////
echo linear graph with nucleotides
echo //////////////////////////////////////////////////////////////

./RNACircos.py -c -n -i example_files/input_dot_bracket.txt

echo //////////////////////////////////////////////////////////////
echo circular graph with nucleotides
echo //////////////////////////////////////////////////////////////

./RNACircos.py -c -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA

echo //////////////////////////////////////////////////////////////
echo superimposed circular graphs
echo //////////////////////////////////////////////////////////////

./RNACircos.py -a -o 0.15 -l -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA

echo //////////////////////////////////////////////////////////////
echo superimposed, aligned linear graphs
echo //////////////////////////////////////////////////////////////

./RNACircos.py -s -c -i example_files/Dot_Bracket_India_B16171 -i2 example_files/Spike_NCBI_MEA

echo //////////////////////////////////////////////////////////////
echo superimposed, aligned circular_graphs
echo //////////////////////////////////////////////////////////////

./RNACircos.py -c -i example_files/gradient_dot_bracket.txt -pm example_files/example_dp.ps

echo //////////////////////////////////////////////////////////////
echo circular, grayscale, gradient graph
echo //////////////////////////////////////////////////////////////


./RNACircos.py -l -i example_files/gradient_dot_bracket.txt -pm example_files/example_dp.ps -pc -lc lightblue -hc red

echo //////////////////////////////////////////////////////////////
echo linear, color, custom,  gradient graph
echo //////////////////////////////////////////////////////////////

./RNACircos.py -c -i example_files/gradient_dot_bracket.txt -pm example_files/example_dp.ps -pc -lc lightblue -hc red -g

echo //////////////////////////////////////////////////////////////
echo circular, color, custom,  gradient graph
echo //////////////////////////////////////////////////////////////
