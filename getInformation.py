import io
import math
import numpy as np


def get_inputs(main_input, generator):
    values_to_find_generator = ["Q_total", "sig_x", "sig_y"]
    values_to_find_input = ["Xrms", "Yrms", "Nue(1)", "MaxE(1)", "Phi(1)", "Nue(2)", "MaxE(2)", "Phi(2)", "MaxB(1)"]
    answer = {}
    finG = io.open(generator, 'r', encoding='utf-8', newline='\n', errors='ignore')
    for i, line in enumerate(finG):
        tokens = line.strip().split('=')
        for value in values_to_find_generator:
            if tokens[0] == value:
                answer[tokens[0]] = float(tokens[1])
    finG.close()
    finI = io.open(main_input, 'r', encoding='utf-8', newline='\n', errors='ignore')
    for i, line in enumerate(finI):
        tokens = line.strip().split('=')
        for value in values_to_find_input:
            if tokens[0] == value:
                answer[tokens[0]] = float(tokens[1])
    finI.close()
    return answer


def get_inputs_main(main_input):
    values_to_find_input = ["QBUNCH"]
    answer = {}
    finI = io.open(main_input, 'r', encoding='utf-8', newline='\n', errors='ignore')
    for i, line in enumerate(finI):
        tokens = line.strip().split('=')
        for value in values_to_find_input:
            if tokens[0] == value:
                answer[tokens[0]] = float(tokens[1])
    finI.close()
    return answer


def get_outputs(output_file):
    arr_init = np.loadtxt(output_file)
    lines_to_delete = []
    for i, line in enumerate(arr_init):
        if line[9] < 0:
            lines_to_delete.append(i)
    arr = np.delete(arr_init, lines_to_delete, 0)
    ref_is_exist = arr_init[0, 2] > 5
    answer = {}
    ref_is_exist = False
    if ref_is_exist:
        arr_ref = arr_init[0, :]
        answer["X_pos"] = float(arr_ref[0]) + np.mean(arr[:, 0])
        answer["Y_pos"] = float(arr_ref[1]) + np.mean(arr[:, 1])
        answer["Z_pos"] = float(arr_ref[2]) + np.mean(arr[:, 2])
        answer["average_kin_energy"] = float(arr_ref[5]) - 511000 + np.mean(arr[:, 5])
        answer["alfa_X"] = np.sqrt(np.mean(np.square(arr[:, 3]))) / (float(arr_ref[5]) - 511000)
        answer["alfa_Y"] = np.sqrt(np.mean(np.square(arr[:, 4]))) / (float(arr_ref[5]) - 511000)
    answer["Charge"] = np.sum(arr[:, 7]) * 0.2/100
    answer["sigZ"] = np.sqrt(np.mean(np.square(arr[:, 2])))
    answer["sigX"] = np.sqrt(np.mean(np.square(arr[:, 0])))
    answer["sigY"] = np.sqrt(np.mean(np.square(arr[:, 1])))
    answer["energy_spread"] = np.sqrt(np.mean(np.square(arr[:, 5])))
    answer["emittance_Y"] = np.sqrt(
        np.mean(np.square(arr[:, 1])) * np.mean(np.square(arr[:, 4])) - math.pow(np.mean((arr[:, 1] * arr[:, 4])), 2))
    answer["emittance_X"] = np.sqrt(
        np.mean(np.square(arr[:, 0])) * np.mean(np.square(arr[:, 3])) - math.pow(np.mean((arr[:, 0] * arr[:, 3])), 2))
    return answer


def load_dataset():
    info = ""
    f = open("information.txt", "w")
    for i in range(1):
        # x_dict = get_inputs_main("run"+str(i+1)+".in")
        # y_dict = get_outputs("run"+str(i+1)+".0528.001")
        x_dict = get_inputs_main("run1.in")
        y_dict = get_outputs("run1.0528.001")

        info += "input:"
        for j in range(len(x_dict)):
            info += ("\n" + str(list(x_dict.keys())[j]) + " " + str(list(x_dict.values())[j]))
        info += "\n" + "output:"
        for j in range(len(y_dict)):
            info += ("\n" + str(list(y_dict.keys())[j]) + " " + str(list(y_dict.values())[j]))
    f.write(info)


load_dataset()
