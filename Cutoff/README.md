# Nepre-Potential(Cutoff)
Scoring Function based on Neighbourhood Preference Statistics  

Usage
----------
The runing folder should contain:
* AminoAcid.py(Class for establish amino acid)
* cutoff.npy(Energy matrix)

We provide 7 cutoff options with cutoff between 4 angstrom and 10 angstrom.

For single protein structure calculate, choose a cutoff(6 angstrom e.g)go to linux shell and type:
<pre><code>
Nepre@liulab:~$ python Nepre.py 6 ./example.pdb
</code></pre>

For multi-object job, you could speed up by using Nepre as a module:

<pre></code>
import Nepre

#choose a cutoff
cutoff = 6

#select a protein
path = "./example"

#load energy matrix
Matrix = Nepre.load_EnergyMatrix(cutoff)

#calculate Nepre potential energy
E = Nepre.Calculate(Matrix,path)
</code></pre>

