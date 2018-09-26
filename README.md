# Nepre-Potential
Scoring Function based on Neighbourhood Preference Statistics  
[Documentation](https://www.baidu.com) |
[Resources](https://www.baidu.com) |
[Installation](https://www.baidu.com) |


Introduction
-----------
Nepre potential include two typical method:
* Calculate the potential energy using different cutoff.
* Calculate the potential energy using the radius of amino acid. 

For the first method, users can use the cutoff between 4 angstrom and 10 angstrom to calculate the potential energy.  
For the second method, we use 30,000 high-resolution protein to get statistical result of radius and use gaussian distribution function to fit. Then Nepre will calculate the potential energy according to the gaussian mean data of radius of each kind of amino acid.
