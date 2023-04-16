print("imports")
import sys,os
import numpy as np
import matplotlib.pyplot as plt

print("setting up paths")
## edit only this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","sidmahesh\OneDrive","Documents\GitHub\Circumbinary-Accretion-Disks")

## no need to edit these paths
formulae_path = os.path.join(repo_path,"Formulae")
sys.path.insert(0,formulae_path)
results_path = os.path.join(repo_path,"Results")

import ResonantTorquePicture as rtp

print("initializing arrays")
## set choices for mass ratios and radio

mass_ratios = np.arange(0.1,0.5+0.1,0.1)
inclinations = np.arange(0,1,0.01)*np.pi/2
eccentricities = np.arange(0.0,0.5,0.01)

## create figure 7

q = 0.3
ecc = [0.0,0.2,0.5]
inc = np.arange(0,1,0.01)*np.pi
alpha = 0.01
chi = 0.05

for j in range(len(ecc)):
    T_11 = np.zeros(len(inc))
    T_21 = np.zeros(len(inc))
    T_22 = np.zeros(len(inc))
    T_31 = np.zeros(len(inc))
    T_41 = np.zeros(len(inc))
    T_51 = np.zeros(len(inc))
    
    for i in range(len(inc)):
        params_11 = [rtp.LR_location(1, 1),q,inc[i],ecc[j],1,1]
        fluidparams_11 = [params_11,alpha,chi]
        T_11[i] = rtp.zeta_T(fluidparams_11)
        
        params_22 = [rtp.LR_location(2, 2),q,inc[i],ecc[j],2,2]
        fluidparams_22 = [params_22,alpha,chi]
        T_22[i] = rtp.zeta_T(fluidparams_22)
        
        params_31 = [rtp.LR_location(3, 1),q,inc[i],ecc[j],3,1]
        fluidparams_31 = [params_31,alpha,chi]
        T_31[i] = rtp.zeta_T(fluidparams_31)
        
        if ecc[j] > 0.1:
            params_21 = [rtp.LR_location(2, 1),q,inc[i],ecc[j],2,1]
            fluidparams_21 = [params_21,alpha,chi]
            T_21[i] = rtp.zeta_T(fluidparams_21)
            
            if ecc[j] > 0.3:
                params_41 = [rtp.LR_location(4, 1),q,inc[i],ecc[j],4,1]
                fluidparams_41 = [params_41,alpha,chi]
                T_41[i] = rtp.zeta_T(fluidparams_41)
                
                params_51 = [rtp.LR_location(5, 1),q,inc[i],ecc[j],5,1]
                fluidparams_51 = [params_51,alpha,chi]
                T_51[i] = rtp.zeta_T(fluidparams_51)
    
    plt.plot(inc*180/np.pi,np.log10(T_11),color = 'blue')
    plt.plot(inc*180/np.pi,np.log10(T_22), color = 'red')
    plt.plot(inc*180/np.pi,np.log10(T_31), color = 'orange')
    if ecc[j] > 0.1:
        plt.plot(inc*180/np.pi,np.log10(T_21),color = 'green')
        if ecc[j] > 0.3:
            plt.plot(inc*180/np.pi,np.log10(T_41),color = 'purple')
            plt.plot(inc*180/np.pi,np.log10(T_51),color = 'yellow')
    plt.ylim(-2,4)
    plt.xlim(0,180)
    plt.axhline(0,linestyle='dashed',color = 'black')
    plt.savefig(os.path.join(results_path,"Fig7_ecc_"+str(j)+".png"))
    plt.show()

## create figure 8

incs = np.array([0.,22.5,45.,90.,135.])
eccs = np.arange(0.01,0.8,0.01)
colorset = ['red','cyan','purple','green','orange']

for m in range(2,4):
    for i in range(len(incs)):
        T = np.zeros(len(eccs))
        for j in range(len(eccs)):
            params = [rtp.LR_location(m, 1),q,incs[i]*np.pi/180,eccs[j],m,1]
            fluidparams = [params,alpha,chi]
            T[j] = rtp.zeta_T(fluidparams)
        plt.plot(np.log10(eccs),np.log10(T),color = colorset[i])
    plt.ylim(-2,3)
    plt.xlim(-2,np.log10(0.79))
    plt.axhline(0,linestyle = 'dashed',color = 'black')
    plt.savefig(os.path.join(results_path,'Fig9m_'+str(m)+'.png'))
    plt.show()
        
        
## create figure 9

mass_ratios_fig9 = [0.1,0.3,0.5]
eccs_fig_9 = np.arange(0.,1,0.01)
incs_figure_9 = [0,22.5,45,90,135]


for i in range(len(mass_ratios_fig9)):
    for inc in range(len(incs_figure_9)):
        xgap = []
        for j in range(len(eccs_fig_9)):
            rgap = 0.
            for m in range(1,6):
                rgapguess = rtp.LR_location(m,1)
                params = [rgapguess,mass_ratios_fig9[i],incs_figure_9[inc]*np.pi/180,eccs_fig_9[j],m,1]
            fluidparams = [params,alpha,chi]
            zeta = rtp.zeta_T(fluidparams)
            if zeta > 1:
                if rgapguess > rgap:
                    rgap = rgapguess
            rgapguess = rtp.LR_location(2,2)
            params = [rgapguess,mass_ratios_fig9[i],incs_figure_9[inc]*np.pi/180,eccs_fig_9[j],2,2]
            fluidparams = [params,alpha,chi]
            zeta = rtp.zeta_T(fluidparams)
            if zeta > 1:
                if rgapguess > rgap:
                    rgap = rgapguess
                    
            xgap.append(rgap)
        plt.plot(eccs_fig_9,xgap,color = colorset[inc])
        if incs_figure_9[inc] == 22.5 or incs_figure_9[inc] == 90 or incs_figure_9[inc] == 135:
            plt.savefig( os.path.join(results_path,'Fig9plotq_'+str(i)+'inc_'+str(inc)+'.png'))
            plt.show()
                    