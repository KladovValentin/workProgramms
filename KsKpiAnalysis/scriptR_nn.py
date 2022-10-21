
import subprocess
import os
os.environ["PATH"] += os.pathsep + 'D:/Program Files/Graphviz/bin/'


#ptrain = subprocess.run("python trainML.py NN")
ptest = subprocess.run("python predict.py NN")
