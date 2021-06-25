import numpy as np
import openpyxl
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
from fullfilm_model import Model
import os

def main():

    mesh_num_list = [50]
    aspect_ratio_list = [0.005, 0.01, 0.02, 0.04, 0.1]

    pcav = 133.322*10**(-6) #cavitation pressure[MPa]
    h0_list = [0.05]            #nominal film thickness [µm]
    rpm_list = [10]                      #rotational speed [rpm]
    for h0 in h0_list:
        for rpm in rpm_list:
            for mesh_num in mesh_num_list:
                P_mean_list = []

                for aspect_ratio in aspect_ratio_list:

                    model = Model(mesh_num, aspect_ratio, h0, rpm, pcav)

                    #######calulation##########################################
                    t1 = time.time()
                    H_mat = model.calc_H()

                    t2 = time.time()
                    print("calc_H:",t2-t1)

                    aN, aS, aW, aE, aO, aW2, aO2, bO = model.calc_coef(H_mat)
                    t3 = time.time()
                    print("calc_coef:",t3-t2)

                    rf_opt = model.find_rf_opt(aN, aS, aW, aE, aO, aW2, aO2, bO)
                    t4 = time.time()
                    print("find_rf_opt:",t4-t3)

                    P_mat, err_list = model.relaxation(100000, rf_opt, aN, aS, aW, aE, aO, aW2, aO2, bO)
                    t5 = time.time()
                    print("relaxation:",t5-t4)

                    model.save_data(H_mat, P_mat, err_list)
                    t6 = time.time()
                    print("save_data:",t6-t5)

                    P_mean = np.mean(P_mat)
                    P_mean_list.append(P_mean)

                create_asp_meanP_graph(aspect_ratio_list, P_mean_list, mesh_num, h0, rpm, pcav)

def create_asp_meanP_graph(aspect_ratio_list, P_mean_list, mesh_num, h0, rpm, pcav):
    fig1 = plt.figure(figsize=(10,14))
    plt.rcParams['font.size'] = 34
    #plt.tick_params(labelsize=18)
    ax1 = fig1.add_subplot(2,1,1)
    ax1.scatter(aspect_ratio_list, P_mean_list, s=200)
    ax1.set_xlabel('Aspect Ratio', fontsize=38, labelpad=12)
    ax1.set_ylabel('Dimensionless Mean P', fontsize=38, labelpad=12)

    result_dir_path = f"./my_dimple/h{h0}_{rpm}rpm_pcav{pcav}"
    os.makedirs(result_dir_path, exist_ok=True)
    fig1.savefig(result_dir_path + f"/asp-meanp_n{mesh_num}.png", bbox_inches="tight")

    fig2 = plt.figure(figsize=(10,14))
    plt.rcParams['font.size'] = 34
    #plt.tick_params(labelsize=18)
    ax2 = fig2.add_subplot(2,1,1)
    depth_list = np.array(aspect_ratio_list)*50
    ax2.scatter(depth_list, P_mean_list, s=200)
    ax2.set_xlabel('Dimple Depth [µm]', fontsize=38, labelpad=12)
    ax2.set_ylabel('Dimensionless Mean P', fontsize=38, labelpad=12)
    ax2.set_xticks(np.arange(0,5.1,1.0))


    result_dir_path = f"./my_dimple/h{h0}_{rpm}rpm_pcav{pcav}"
    os.makedirs(result_dir_path, exist_ok=True)
    fig2.savefig(result_dir_path + f"/depth-meanp_n{mesh_num}.png", bbox_inches="tight")
    plt.close()
    #./my_dimple/h{self.h0}_{self.rpm}rpm_pcav{self.p_cav}/validation_depth{self.DG}_n{self.nr}
if __name__ == '__main__':
    main()
    print('END')
