import os
import math
import numpy as np
import sys

ccDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
        "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
        "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
        "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}


cdDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
        "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
        "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
        "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}

def bilinear_interpolation(matrix,emp=0):
    theta = np.shape(matrix)[0]
    phi = np.shape(matrix)[1]
    res = np.zeros([theta,phi])

    for i in range(theta):
        for j in range(phi):
            if(matrix[i][j] == 0):
                val = 0
                neighbor = []
                neighbor.append(matrix[i-1][j])
                neighbor.append(matrix[(i+1)%theta][j])
                neighbor.append(matrix[i][j-1])
                neighbor.append(matrix[i][(j+1)%phi])
                while(0 in neighbor):
                    neighbor.remove(0)
                if(neighbor):
                    res[i][j] = int(float(sum(neighbor)/len(neighbor)))
                else:
                    res[i][j] = emp
            else:
                res[i][j] = matrix[i][j]

    return res

if __name__ == "__main__":
    
    args = sys.argv[1:]
    cutoff = int(args[0])
    List = ccDict.keys()
    List.sort()
    fr = open("../../cc-10-Energy.txt")
    i = 0
    for line in fr.readlines():
        line=line.strip().split()
        for j in range(20):
            ccDict[List[i]][List[j]]=float(line[j])
        i += 1
    fr.close()
    emptyarea = {}
    print("Start to read coordinate.txt")
    fr=open("../../coordinate.txt")
    amino1=''
    amino2=''
    
    for line in fr.readlines():
        line = line.strip()
        if line=='':
            continue
        if len(line)==3:
            if amino2=='VAL' or amino1=='':
                amino1=line
                amino2=''
            else:
                amino2=line
                print amino1,amino2
        else:
            line=line.split()
            if amino2 not in cdDict[amino1]:
                cdDict[amino1][amino2]=[[float(line[0]),float(line[1]),float(line[2])]]
            else:
                cdDict[amino1][amino2].append([float(line[0]),float(line[1]),float(line[2])])
    fr.close()

    for amino1 in cdDict:
        for amino2 in cdDict[amino1]:
            emp = 0
            dualArray=np.zeros((20,20))
            print amino1,amino2
            for coor in cdDict[amino1][amino2]:
                coor = np.array(coor)
                rho = sum(coor**2)**0.5
                theta = np.arccos(coor[2]/rho)
                phi = np.arctan2(coor[1],coor[0])
                theta = min(int(theta*20/np.pi),19)
                phi = min(int((phi+np.pi)/(2*np.pi/20)),19)
                #rho = min(int(rho),10)
                
                if(rho < cutoff):
                    dualArray[theta][phi] += 1.0

            print("Start to find min data")
            mindata = None
            for theta in range(20):
                for phi in range(20):
                    if(mindata == None and dualArray[theta][phi] != 0):
                        mindata = dualArray[theta][phi]
                    if(mindata > dualArray[theta][phi] and dualArray[theta][phi] != 0):
                        mindata = dualArray[theta][phi]
            print("Start to check empty area")
            dualArray = bilinear_interpolation(dualArray,mindata)
            
            dualArray = dualArray / dualArray.sum()

            Integral = np.ones((20,20))
            for j in range(20):
                Integral[j] = Integral[j]*(np.cos(j*np.pi/20)-np.cos((j+1)*np.pi/20))*((1./3)*(2*np.pi/20))
            
            # normalized the intergral
            Integral = Integral / Integral.sum()
            # normalized the p(i,j,r,theta,phi)
            dualArray = dualArray / Integral

            cdDict[amino1][amino2] = ccDict[amino1][amino2]-np.log(dualArray)
            emptyarea[amino1+'-'+amino2] = emp

    f = file(str(cutoff) + ".npy",'wb')
    for amino1 in List:
        for amino2 in List:
            np.save(f, cdDict[amino1][amino2])

