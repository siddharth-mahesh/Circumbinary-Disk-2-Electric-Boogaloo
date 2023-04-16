import numpy as np
import matplotlib.pyplot as plt
import sys,os

## edit this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","sidmahesh/OneDrive","Documents\GitHub\Circumbinary-Accretion-Disks")

## no need to edit these paths
data_path = os.path.join(repo_path,"Data")

# Define MP Mass Ratios and Inclinations
mass_ratios = [0.1,0.2,0.3,0.4,0.5]
mr_labels = ["5","4","3","2","1"]
inc_labels_125 = ["000","015","020","025","030","035","040","045","050","055","060","075","090","105","120","135","150","165","180"]
incs_125 = [0,15,20,25,30,35,40,45,50,55,60,75,90,105,120,135,150,165,180]
inc_labels_34 = ["000","015","030","045","060","075","090","105","120","135","150","165","180"]
incs_34 = [0,15,30,45,60,75,90,105,120,135,150,165,180]

gaps = open(os.path.join(data_path,"numericalgapsizes.txt"),"w")

for i in range(len(mass_ratios)):
    
    if mr_labels[i] in ["1","2","5"]:
        
        for j in range(len(incs_125)):
            datfile = np.loadtxt(os.path.join(data_path,"FinalDistribution",mr_labels[i]+"_"+inc_labels_125[j]+"_06.dat"))
            radii = datfile[:,3]
            cutoff = np.argmin(np.abs(radii - 3))
            densities = datfile[:,4]
            maxdense = np.max(densities)
            sig1p = 0.02*maxdense
            if mr_labels[i] =="5":
                if inc_labels_125[j] == "180":
                    check = 1
                    #print(np.abs(densities[2:cutoff] - .02*maxdense))
            sig1p_index = np.argmin(np.abs(densities[2:cutoff] - .02*maxdense)) + 2
            print(i,j,sig1p_index)
            if densities[sig1p_index] > sig1p:
                sig1p_m_index = sig1p_index - 1
            else:
                sig1p_m_index = sig1p_index
            strtoprint = "["+ str(densities[sig1p_index-1])+","+ str(densities[sig1p_index])+","+ str(densities[sig1p_index+1])+"]"
            gaps.write("0."+str(int(10*mass_ratios[i])) + "\t" + str(incs_125[j]) + "\t" + strtoprint + "\t"+ str(sig1p) + "\t"+str(sig1p_m_index)+"\t"+str(densities[cutoff])+"\n")
    
    elif mr_labels[i] in ["3","4"]:
        
        for j in range(len(incs_34)):
            datfile = np.loadtxt(os.path.join(data_path,"FinalDistribution",mr_labels[i]+"_"+inc_labels_34[j]+"_06.dat"))
            radii = datfile[:,3]
            cutoff = np.argmin(np.abs(radii - 3))
            densities = datfile[:,4]
            maxdense = np.max(densities)
            sig1p = 0.02*maxdense
            sig1p_index = np.argmin(np.abs(densities[2:cutoff] - .02*maxdense)) + 2
            print(i,j,sig1p_index)
            if densities[sig1p_index] > sig1p:
                sig1p_m_index = sig1p_index - 1
            else:
                sig1p_m_index = sig1p_index
            strtoprint = "["+ str(densities[sig1p_index-1])+","+ str(densities[sig1p_index])+","+ str(densities[sig1p_index+1])+"]"
            gaps.write("0."+str(int(10*mass_ratios[i])) + "\t" + str(incs_34[j]) + "\t" + strtoprint + "\t"+ str(sig1p) + "\t"+str(sig1p_m_index)+"\t"+str(densities[cutoff])+"\n")

gaps.close()