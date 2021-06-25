import numpy as np
import openpyxl
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

class Model:

    ###set parameter##########################################3

    errmax = 0.0001
    pa = 0.1                    #atmospheric pressure [MPa]

    dimple_row = 40                 #dimple array
    dimple_col = 135
    
    rout = 5.0                 #outer radius of one element [mm]
    r_range = 6.5/(dimple_row-1)    #radius range [mm]
    rin = rout - r_range        #inner radius of one element [mm]

    t_range = np.pi/dimple_col   #angle range
    eta = 0.03                 #viscosity [Pa*s]

    # RG = 0.075/2                        #groove radius [mm]
    # DG = RG*2*self.aspect_ratio*10**3   #groove depth [µm]

    # #number of elements
    # nr = self.mesh_num
    # ntheta = self.mesh_num

    # #length of elements (dimentionless)
    # dr = r_range/rout/nr
    # dtheta = t_range/ntheta

    # #groove parameter (dimentionless)
    # rg = RG/rout
    # dg = DG/h0

    def __init__(self, mesh_num, aspect_ratio, h0, rpm, pcav):
        
        self.p_cav = pcav                #cavitation pressure[MPa]
        self.h0 = h0                    #nominal film thickness [µm]
        self.rpm = rpm                   #rotational speed [rpm]

        omega = rpm*2*np.pi/60      #angular velosity [rad/s]
        self.lambd = 6*self.eta*omega/(self.pa-self.p_cav)*self.rout**2/self.h0**2


        self.RG = 0.05/2                        #groove radius [mm]
        self.DG = self.RG*2*aspect_ratio*10**3   #groove depth [µm]

        #number of elements
        self.nr = mesh_num
        self.ntheta = mesh_num

        #length of elements (dimentionless)
        self.dr = self.r_range/self.rout/self.nr
        self.dtheta = self.t_range/self.ntheta

        #groove parameter (dimentionless)
        self.rg = self.RG/self.rout
        self.dg = self.DG/self.h0


    ###function############################################
    def calc_H(self):
        #coordinate of H matrix
        nx = self.nr*2+1
        ny = self.ntheta*2+1
        H_ini = np.ones((nx, ny))           #initialize H (assume h == h0)

        #position of groove center
        xc = (self.rin + self.r_range/2)/self.rout
        yc = xc*self.t_range/2

        #position of each coodinate point
        x_mat = np.tile(np.linspace(self.rin/self.rout,1.0,nx),(ny,1)).T
        theta_mat = np.tile(np.linspace(0,self.t_range,ny),(nx,1))
        y_mat = x_mat*theta_mat

        gcd = (x_mat- xc)*(x_mat- xc) + (y_mat- yc)*(y_mat- yc)     #distance from groove center

        zg = np.where(gcd < self.rg**2, self.dg*(gcd/(self.rg**2) - 1.0), 0)             #approximated by quadratic curve
        H_mat = H_ini - zg

        return H_mat

    def calc_coef(self, H_mat):
        Hn = H_mat[2:2*self.nr-3:2, 1:2*self.ntheta:2]
        Hs = H_mat[4:2*self.nr-1:2, 1:2*self.ntheta:2]
        Hw = H_mat[3:2*self.nr-2:2, 0:2*self.ntheta-1:2]
        He = H_mat[3:2*self.nr-2:2, 2:2*self.ntheta+1:2]

        R_mat = np.tile(np.linspace(self.rin/self.rout,1.0,self.nr*2+1),((self.ntheta*2+1),1)).T          #radius of each point
        
        Rn = R_mat[2:2*self.nr-3:2, 1:2*self.ntheta:2]
        Rs = R_mat[4:2*self.nr-1:2, 1:2*self.ntheta:2]
        RO = R_mat[3:2*self.nr-2:2, 1:2*self.ntheta:2]

        aN0 = Hn**3*Rn*self.dtheta/self.dr
        aS0 = Hs**3*Rs*self.dtheta/self.dr
        aW0 = Hw**3/RO*self.dr/self.dtheta
        aE0 = He**3/RO*self.dr/self.dtheta
        aO0 = aN0 + aS0 + aW0 + aE0
        aW2 = self.lambd*self.dr*RO*Hw
        aO2 = self.lambd*self.dr*RO*He
        bO = aW2 - aO2

        return aN0, aS0, aW0, aE0, aO0, aW2, aO2, bO

    #@jit('Tuple((f8[:,:],f8[:]))(i8,f8,f8[:,:],f8[:,:],f8[:,:],f8[:,:],f8[:,:],f8[:,:],f8[:,:],f8[:,:])', nopython=True) #, nopython=True
    # @jit('f8[:,:],f8[:](i8,f8,f8[:,:],f8[:,:],f8[:,:],f8[:,:],f8[:,:],f8[:,:],f8[:,:],f8[:,:])') #, nopython=True
    def relaxation(self, itr_max, rf, aN0, aS0, aW0, aE0, aO0, aW2, aO2, bO):
        P_ini = np.ones((self.nr, self.ntheta))    #initialize pressure matrix('P'='phi')

        PO = P_ini[1:self.nr-1,:]

        err = 1.0
        itr = 0
        err_list = []
        one_arr = np.ones((1, self.ntheta))
        while err > self.errmax and itr < itr_max:
            err = 0.0
            itr += 1

            PN = np.vstack((one_arr, PO[:self.nr-3,:]))
            PS = np.vstack((PO[1:], one_arr))
            PW = np.hstack((PO[:,self.ntheta-1].reshape(self.nr-2,1), PO[:,:self.ntheta-1]))
            PE = np.hstack((PO[:,1:], PO[:,0].reshape(self.nr-2,1)))

            FO = np.where(PO<0, 0.0, 1.0)
            FN = np.where(PN<0, 0.0, 1.0)
            FS = np.where(PS<0, 0.0, 1.0)
            FW = np.where(PW<0, 0.0, 1.0)
            FE = np.where(PE<0, 0.0, 1.0)
            
            aN = aN0*FN
            aS = aS0*FS
            aW = aW0*FW + aW2*(1-FW)
            aE = aE0*FE
            aO = aO0*FO + aO2*(1-FO)
            PO_i = (aN*PN + aS*PS + aW*PW + aE*PE + bO)/aO

            #err_mat = np.where(PO_i < 0, 0, PO_i) - PO
            err_mat = PO_i - PO

            PO = PO + rf*err_mat
            #PO = np.where(PO < 0, 0, PO)

            err = np.sum(np.abs(err_mat))
            err_list.append(err)

        PO = np.where(PO<0, 0, PO)                  #modify pressure finally
        phi_mat = np.vstack([one_arr,PO,one_arr])   #phi
        P_mat = phi_mat + (1 - phi_mat)*self.p_cav/self.pa    #dimentionless pressure (p/pa)
        return P_mat, err_list

    #find optimum value of relaxation factor
    def find_rf_opt(self, aN, aS, aW, aE, aO, aW2, aO2, bO):
        rf_list = np.linspace(0, 1.0, 9)
        err_min = np.inf
        rf_opt = 1.0
        for rf in rf_list:
            _, err_current = self.relaxation(1000, rf, aN, aS, aW, aE, aO, aW2, aO2, bO)
            if err_current[-1] < err_min:
                rf_opt = rf
                err_min = err_current[-1]

        return rf_opt

    def save_data(self, H_mat, P_mat, err_list):
        
        #H_mat = H_mat.astype(np.float16)
        
        #3Dvgraph
        fig1 = plt.figure(figsize=(18,20))
        plt.rcParams['font.size'] = 24
        
        #H graph
        r_axis = np.linspace(0, self.r_range, 2*self.nr+1)
        t_axis = np.linspace(0, self.t_range, 2*self.ntheta+1)
        R_axis, T_axis  = np.meshgrid(r_axis, t_axis)

        ax0 = plt.subplot(3,2,1,projection='3d')
        ax0.plot_wireframe(T_axis, R_axis, H_mat.T, linewidth=0.8)
        ax0.set_xlabel('θ', fontsize=30, labelpad=32)
        ax0.set_ylabel('R', fontsize=30, labelpad=18)
        ax0.set_zlabel('Film Thickness', fontsize=30, labelpad=18)
        ax0.set_xticks([0, self.t_range/2, self.t_range])
        #ax0.set_xticks([0, self.t_range/4, self.t_range/4*2, self.t_range/4*3, self.t_range])
        ax0.set_yticks(np.arange(0, self.r_range*1.1, 0.1))
        theta_coef = str(int(round(np.pi/self.t_range*2, -1)))
        ax0.set_xticklabels(['0', 'π/'+theta_coef, '2π/'+theta_coef])
        ax0.tick_params(axis='x', pad=16)
        
        #error graph
        ax1 = fig1.add_subplot(3,2,2)
        itr_axis = range(len(err_list))
        ax1.plot(itr_axis, err_list)
        ax1.set_xlabel('Itteration', fontsize=30)
        ax1.set_ylabel('Error', fontsize=30)
        
        #P graph
        r_axis2 = np.linspace(0, self.r_range, self.nr)
        t_axis2 = np.linspace(0, self.t_range, self.ntheta)
        R_axis2, T_axis2  = np.meshgrid(r_axis2, t_axis2)

        ax2 = plt.subplot(2,1,2,projection='3d')
        ax2.plot_wireframe(T_axis2, R_axis2, P_mat.T, linewidth=1.2)
        
        ax2.set_xlabel('θ', fontsize=42, labelpad=84)
        ax2.set_ylabel('R', fontsize=42, labelpad=40)
        ax2.set_zlabel('Dimensionless Pressure', fontsize=42, labelpad=45)
        
        #ax2.zaxis.set_major_locator(mpl.ticker.AutoLocator())
        
        # ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        #ax2.zaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
        # ax2.xaxis.set_major_locator(mpl.ticker.LinearLocator(5))
        
        ax2.set_xticks([0, self.t_range/4, self.t_range/4*2, self.t_range/4*3, self.t_range])
        ax2.set_yticks(np.arange(0, self.r_range+0.01, 0.1))
        zmin, zmax = ax2.get_zlim()
        #ax2.set_zticks(np.arange(zmin, zmax+0.01, 0.05))
        theta_coef = str(int(round(np.pi/self.t_range*4, -1)))
        ax2.set_xticklabels(['0', 'π/'+theta_coef, '2π/'+theta_coef, '3π/'+theta_coef, '4π/'+theta_coef], rotation=10)
        #ax2.xaxis.set_label_coords(-0.2, 0.5)

        ax2.tick_params(axis='x', labelsize=36, pad=20)
        ax2.tick_params(axis='y', labelsize=36, pad=3)
        ax2.tick_params(axis='z', labelsize=36, pad=20)
        fig1.tight_layout

        result_dir_path = f"./my_dimple/h{self.h0}_{self.rpm}rpm_pcav{self.p_cav}"
        os.makedirs(result_dir_path, exist_ok=True)
        fig1.savefig(result_dir_path + f"/depth{self.DG}_n{self.nr}.png")
        plt.close()

        #2D graph
        fig2 = plt.figure(figsize=(16,22))
        plt.subplots_adjust(hspace=0.3)
        plt.rcParams['font.size'] = 38

        ax3 = fig2.add_subplot(2,1,1)
        ax3.plot(r_axis2, P_mat[:, int(self.ntheta/2)], linewidth=5.0)
        ax3.set_xlabel('R', fontsize=46)
        ax3.set_ylabel('Dimensionless Pressure', fontsize=46)

        ax4 = fig2.add_subplot(2,1,2)
        ax4.plot(t_axis2, P_mat[int(self.nr/2), :], linewidth=5.0)
        ax4.set_xlabel('θ', fontsize=46)
        ax4.set_ylabel('Dimensionless Pressure', fontsize=46)
        ax4.set_xticklabels(['0', 'π/'+theta_coef, '2π/'+theta_coef, '3π/'+theta_coef, '4π/'+theta_coef])
        ax4.xaxis.set_major_locator(mpl.ticker.LinearLocator(5))
        
        fig2.savefig(result_dir_path + f"/depth{self.DG}_n{self.nr}_2.png")
        
        plt.close()
        #save P_mat to excel
        # df = pd.DataFrame(P_mat)
        # df.to_excel(f"validation_{rpm}rpm_pcav{p_cav}_n{nr}.xlsx")

