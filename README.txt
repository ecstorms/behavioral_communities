
________________________________________________________________________________________________________________________________
0. Summary

Date: October, 2019

This is the readme file for the code and data accompanying the paper "Behavioral Communities and the Atomic Structure of Networks" by Matt Jackson and Evan Storms. This and the related files can be found in the .zip file labeled "Accompanying Data and Code." 

The accompanying code is found in the file "behavioral_communities_code_to_share.py." 

The accompanying data is found in the files "Indian Villages HH" and "Zachary Karate Club." 


________________________________________________________________________________________________________________________________
1. Code Overview

The code file is divided into 11 sections.  

Section 0 imports the Python 3 modules needed to run the rest of the code, and defines some basic utility functions used throughout.  All of the modules can be installed using Pip.  The basic network interfacing uses the Graph-Tools module, for which a tutorial can be found here: https://graph-tool.skewed.de/ Alternatively, we include the anaconda environment 'behavioral_communities.yml' which can be loaded using Anaconda running Python 3.7 or later. 

Section 1 defines basic functions for testing whether any of the elements of a parition of the nodes forms a convention. 

Section 2 implements the approximation algorithm for finding the q-atoms of a network. 

Section 3 implements a brute-force algorithm for finding the q-atoms of a network.  This is mainly used to check that the approximation algorithm generally finds the correct atomic structure. 

Section 4 adops the approximation algorithm to finding the robust Q-range atomic structures. 

Section 5 defines functions for randomly generating converntions on a network. 

Section 6 defines functions for random seeding, optimal seeding by brute force, and seeding used the atom-based heuristic described in section 6 of the paper. 

Section 7 is the analogue of section 2 for the absolute threshold (t threshold) case. 

Section 8 provides the functions for estimating a relative behavioral threshold q from a network and data on behaviors. 

Section 9 is the analogue of section 8 for the absolute threshold (t threshold) case. 

Section 10 works through an example of finding the various atomic structures for an Erdos-Renyi network. 

Section 11 reproduces the Indian village figures from section 6.1 in the paper. 

Section 12 reproduces the AddHealth data figures 5-8, 17-18, 21-24, and 28-30. 

Section 13 reproduces figures 25-27 using the data on Zachary's Karate Club. 

Section 14 reproduces the seeding simulation results from section 6 of the paper. 

_______________________________________________________________________________________________________________________________
2. Indian Village Data 

The file 'Indian Villages HH' contains the data on networks in Indian villages from Banerjee, Chandrasekhar, Duflo, and Jackson (2013)). The original readme file meant to accompany their data is included.  The figures in section 6.1 are based off of the files '26KeroGo[edges].csv' and '26KeroGo[nodes].csv' along with the file 'covariates_allvillages.csv' which contains data on nodes' caste memebership. 




________________________________________________________________________________________________________________________________
3.  Zachary's Karate Club Data

The file 'Zachary_karate.gt' contains the Zachary's Karate club network from Zachary, W. (1977).  The binary vertex property 'loyalty' records whether nodes remained in the original club or joined the splinter group. 



________________________________________________________________________________________________________________________________
4.  AddHealth Data 


The AddHealth data used to construct figures 5-8, 17-18, 21-24, and 28-20 is not made publicly available and cannot be shared. We have included the code used to produce these figures in Section 12 of the code file.  Here we include some details on where the relevant files can be found in the AddHealth data repository for those researchers with access.  

The high school social network data used to construct figures 5-8 and 17-18 can be found in the'Structure' subdirectory.  The network shown corresponds to the community with community id 1. 

The newtork and behavior data used to construt the remaining AddHealth figures can be found in the 'Waves1_2' subdirectory.  Social network data are constructed from the data in the 'Waves1_2/Friends' subdirectory, and data on behaviors are constructed from data in the 'Waves1_2/Schl_Info' subdirectory.  The figures shown are for students with school code 31. 
