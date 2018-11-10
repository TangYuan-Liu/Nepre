# Nepre-Potential(Cutoff)
Scoring Function based on Neighbourhood Preference Statistics  

Usage
----------
The runing folder should contain:
* AminoAcid.py(Class for establish amino acid)
* cutoff.npy(Energy matrix)

We provide 7 cutoff options with cutoff between 4 angstrom and 10 angstrom.
You can see help information by typing:
<pre><code>
Nepre@liulab:~$ python Nepre.py -h
usage: Nepre.py [-h] [-s | -m] [-o] path cutoff

Nepre-F Scoring Function Created by CSRC

positional arguments:
  path          PDB file path of folder path
  cutoff        cutoff parameter for Nepre-F

optional arguments:
  -h, --help    show this help message and exit
  -s, --single  calculate single PDB
  -m, --multi   calculate a series of PDB
  -o, --output  save the results as a text file in running folder
</code></pre>

Nepre-F program provide two kinds of calculation mode:
For single protein potential energy calculate, choose a cutoff(6 angstrom e.g) and go to linux shell and type:
<pre><code>
# Not save the results in a text file
Nepre@liulab:~$ python Nepre.py -s ./example.pdb 6
# Save the results a text file(Same folder with Nepre.py with name "latest_results.txt")
Nepre@liulab:~$ python Nepre.py -s -o ./example.pdb 6
</code></pre>

<pre><code>
# Not save the results in a text file
Nepre@liulab:~$ python Nepre.py -m ./pdb_folder/ 6
# Save the results a text file(Same folder with Nepre.py with name "latest_results.txt")
Nepre@liulab:~$ python Nepre.py -m -o ./pdb_folder/ 6
</code></pre>

You can also use Nepre.py as a module if you want to use the calculation results for other purposes:
<pre><code>
import Nepre

#choose a cutoff
cutoff = 6

#select a protein
path = "./example.pdb"
f = open(path)

#load energy matrix
Matrix = Nepre.load_EnergyMatrix(cutoff)

#calculate Nepre potential energy
E = Nepre.calculate_Energy(f,Matrix,cutoff)
</code></pre>

Extensions
----------
Nepre module also provide some useful function:
* Calculate the pearson coefficient correlation.
* Extract data from standard PDB file.
<pre><code>
"""
Pearson Coefficient
"""
import Nepre
x = [1,2,3,4]
y = [1,2,3,4]
p = Nepre.Pearson(x,y)

"""
Extract Data
"""
import Nepre
f = open("./example.pdb")
res = []
for line in f.readlines():
    res.append(Nepre.extract_Data(line))
</code></pre>
