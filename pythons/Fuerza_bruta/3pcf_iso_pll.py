import sys 
import numpy as np
import time
from joblib import Parallel, delayed

#=================================================================
#==================== Funciones a usar ===========================
#=================================================================

#=================================================================

def histo_XXX(ii, p, XXX, bn, point_max):
    
    s = ii-p 
    r = np.sqrt(s[:,0]**2+s[:,1]**2+s[:,2]**2)
        
    m = 0
    for (kk) in zip(p):
            
        r12 = np.ones_like(r[m+1:])*r[m]
        r13 = r[m+1:]
    
        ss = kk-p[m+1:]
        r32 = np.sqrt(ss[:,0]**2+ss[:,1]**2+ss[:,2]**2)
        
        dis, dim = np.histogramdd(np.array([r12,r13,r32]).T,bins = [bn,bn,bn], 
                                  range=[[0,point_max],[0,point_max],[0,point_max]])
        XXX = XXX + dis
        
        m = m+1
    
    return XXX

#=================================================================

def Histos_pll_XXX(data, bn, point_max, nuc):

    #Inicializamos el arreglo del histograma
    XXX = np.zeros((bn,bn,bn))
    return np.sum(np.array(Parallel(n_jobs=nuc)(delayed(histo_XXX)(data[i], data[i+1:], XXX, bn, point_max) for i in range(len(data)))),axis=0)

#=================================================================

def histo_XXY(ii, p, p_r, XXY, bn, point_max, n):
    
    n = n+1
    s = ii-p[n:] # Diferencia entre el punto pivote y los demas puntos siguientes 
    r = np.sqrt(s[:,0]**2+s[:,1]**2+s[:,2]**2)
        
    m = 0
    for kk in p[n:]:
            
        #XXY =====================================================
        r12 = np.ones_like(p_r.T[0])*r[m]
        
        ss = ii-p_r # Diferencia entre el punto pivote y el punto en la muestra random 
        r13 = np.sqrt(ss[:,0]**2+ss[:,1]**2+ss[:,2]**2)
        
        ss = kk-p_r # Diferencia entre el punto 2 de los datos y el punto en la muestra random 
        r32 = np.sqrt(ss[:,0]**2+ss[:,1]**2+ss[:,2]**2)
        
        # Histogramas para XXY
        dis, dim = np.histogramdd(np.array([r12,r13,r32]).T,bins = [bn,bn,bn], 
                                range=[[0,point_max],[0,point_max],[0,point_max]])
        XXY = XXY + dis
            
        m = m+1
    
    return  XXY

#=================================================================

def Histos_pll_XXY(data1, data2, bn, point_max, nuc):
    
    #Inicializamos los arreglos de los histogramas
    XXY = np.zeros((bn,bn,bn))
    return np.sum(np.array(Parallel(n_jobs=nuc)(delayed(histo_XXY)(data1[i], data1,
                                                         data2, XXY, bn, point_max, i) for i in range(len(data1)))),axis=0)

#=================================================================
#===================== Tarea principal ===========================
#=================================================================

data = np.loadtxt(sys.argv[1])[:,:3]
rand = np.loadtxt(sys.argv[2])[:,:3]
name = sys.argv[3]

nuc = 2
bins = 30

start_pll = time.perf_counter()

NDDD = Histos_XXX(data,bins,180,nuc)
NRRR = Histos_XXX(rand,bins,180,nuc)
NDDR = Histos_XXY(data,rand,bins,180,nuc)
NDRR = Histos_XXY(rand,data,bins,180,nuc)

finish_pll = time.perf_counter()
print(f'Finializó en {round(finish_pll-start_pll,2)} segundos\n')

# Simetrizamos los Histogramas:

for i in range(bins - 2):
    for j in range(i+1, bins-1):
        for k in range(j+1, bins):
            s = NDDD[i][j][k] + NDDD[i][k][j] + NDDD[j][i][k] + NDDD[j][k][i] + NDDD[k][i][j] + NDDD[k][j][i]
            NDDD[i][j][k] = s
            NDDD[i][k][j] = s
            NDDD[j][i][k] = s
            NDDD[j][k][i] = s
            NDDD[k][i][j] = s
            NDDD[k][j][i] = s
            s = NRRR[i][j][k] + NRRR[i][k][j] + NRRR[j][i][k] + NRRR[j][k][i] + NRRR[k][i][j] + NRRR[k][j][i]
            NRRR[i][j][k] = s
            NRRR[i][k][j] = s
            NRRR[j][i][k] = s
            NRRR[j][k][i] = s
            NRRR[k][i][j] = s
            NRRR[k][j][i] = s
            
for i in range(bins - 2):
    for j in range(i+1, bins-1):
        for k in range(j+1,bins):
            s = NDDR[i][j][k] + NDDR[j][i][k]+ NDDR[i][k][j]+ NDDR[j][k][i]+ NDDR[k][i][j] + NDDR[k][j][i]
            NDDR[i][j][k] = s
            NDDR[i][k][j] = s
            NDDR[j][i][k] = s
            NDDR[j][k][i] = s
            NDDR[k][i][j] = s
            NDDR[k][j][i] = s
        
            s = NDRR[i][j][k] + NDRR[j][i][k]+ NDRR[i][k][j]+ NDRR[j][k][i]+ NDRR[k][i][j]+ NDRR[k][j][i]
            NDRR[i][j][k] = s
            NDRR[i][k][j] = s
            NDRR[j][i][k] = s
            NDRR[j][k][i] = s
            NDRR[k][i][j] = s
            NDRR[k][j][i] = s

#Guaradamos los histogramas:

np.savetxt('/home/echeveste/mis_trabajos/correlation_f_2/datos/histogramas/BF_Oscar_py/DDD_iso_pll_{0}_osc.dat'.format(str(name)), NDDD.reshape((bins,bins*bins)))
np.savetxt('/home/echeveste/mis_trabajos/correlation_f_2/datos/histogramas/BF_Oscar_py/RRR_iso_pll_{0}_osc.dat'.format(str(name)), NRRR.reshape((bins,bins*bins)))
np.savetxt('/home/echeveste/mis_trabajos/correlation_f_2/datos/histogramas/BF_Oscar_py/DDR_iso_pll_{0}_osc.dat'.format(str(name)), NDDR.reshape((bins,bins*bins)))
np.savetxt('/home/echeveste/mis_trabajos/correlation_f_2/datos/histogramas/BF_Oscar_py/DRR_iso_pll_{0}_osc.dat'.format(str(name)), NDRR.reshape((bins,bins*bins)))

