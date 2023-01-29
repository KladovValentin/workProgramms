
import subprocess
import os
os.environ["PATH"] += os.pathsep + 'D:/Program Files/Graphviz/bin/'
import PySimpleGUI as sg


layout = [[sg.Text("train or predict")],
[sg.Button("train"),sg.Button("predictMod"),sg.Button("predictExp"),sg.Button("checkDistrs")]]
window = sg.Window("Demo", layout, margins=(500, 200))

while True:
    event, values = window.read()
    if event == "train":
        print("training")
        window.Minimize()
        ptrain = subprocess.run("python trainML.py NN")
        window.Normal()
    elif event == "predictMod":
        window.Minimize()
        ptest = subprocess.run("python predict.py Mod")
        window.Normal()
    elif event == "predictExp":
        window.Minimize()
        ptest = subprocess.run("python predict.py Exp")
        window.Normal()
    elif event == "checkDistrs":
        window.Minimize()
        pcheck = subprocess.run("python predict.py Dis")
    elif event == sg.WIN_CLOSED:
        break

window.close()

#ptrain = subprocess.run("python trainML.py NN")
#ptest = subprocess.run("python predict.py NN")
