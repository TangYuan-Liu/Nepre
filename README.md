# Nepre-Potential
Scoring Function based on Neighbourhood Preference Statistics  
[Documentation]() |
[Resources]() |
[Installation]() |


Introduction
-----------
Nepre is a scoring function for calculating the protein's potential energy. And it will help to predict the 3D structure of a protein.
Nepre is designed to be **efficient**, **flexible** and **protable** with **two** typical algorithm included:
* Calculate the potential energy using different cutoff (Nepre-F).
* Calculate the potential energy using the radius of amino acid (Nepre-R). 

For **Nepre-F**, users can use the cutoff between 4 angstrom and 10 angstrom to calculate the potential energy.  
For **Nepre-R**, 30,000 high-resolution protein are used to get statistical result of radius and fitted by the gaussian distribution function. Then Nepre will calculate the potential energy according to the gaussian mean data of radius of each kind of amino acid.

Copyright
-----------
Nepre is created by liulab of Beijing Computational Science Research Center.
