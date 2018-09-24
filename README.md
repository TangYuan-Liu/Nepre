# Nepre-Potential
Scoring Function based on Neighbourhood Preference Statistics  
[Documentation](https://xgboost.readthedocs.org) |
[Resources](demo/README.md) |
[Installation](https://xgboost.readthedocs.org/en/latest/build.html) |

Usage
----------
For **Nepre.py**
<pre><code>
Paramters: (cutoff, pdb_path)
Nepre@liulab:~$ python Nepre.py 4 ./example.pdb
</code></pre>

You could use Nepre as a module by using:
<pre><code>
import Nepre
Matrix = Nepre.load_EnergyMatrix(4)
path = "./example.pdb"
E = Nepre.Calculate(Matrix,path)
</code></pre>
