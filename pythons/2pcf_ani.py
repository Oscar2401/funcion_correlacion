import sys 
import numpy as np
import time

data = np.loadtxt(sys.argv[1])[:,:3]
rand = np.loadtxt(sys.argv[2])[:,:3]
name = sys.argv[3]

#FUNCIÓN DE CORRELACIÓN Peebles-Hauser
def Histos_ani(p,p_r,bn,point_max):
    
    """ 
    Función para construir los histogramas en 2D con anisotropia
    
    p = datos
    p_r = random
    bn = tamaño de bins
    point_max = punto máximo en el histograma
    
    """
    
    # Iniciamos los arreglos los histogramas 
    NDD = np.zeros((bn,bn))
    NRR = np.zeros((bn,bn))
    NDR = np.zeros((bn,bn))
    
    n = 0
    for (ii, jj) in zip(p, p_r):
        n = n+1
        
        #DD
        s = ii-p[n:] # vectores diferencia entre dos puntos
        r_ort = np.sqrt(s[:,0]**2+s[:,1]**2)
        r_pll = np.abs(ii[2]-p[n:,2])
        
        #Hacemos los histogramas para estas distancias 
        dis, x,y = np.histogram2d(r_pll, r_ort, bins=[bn,bn], range=([[0, point_max], [0, point_max]]))
        NDD = NDD + 2*dis
        
        #RR
        s = jj-p_r[n:] # vectores diferencia entre dos puntos
        r_ort = np.sqrt(s[:,0]**2+s[:,1]**2)
        r_pll = np.abs(jj[2]-p_r[n:,2])
        dis, x,y = np.histogram2d(r_pll, r_ort, bins=[bn,bn], range=([[0, point_max], [0, point_max]]))
        NRR = NRR + 2*dis
    
    
    for ii in p:
        s = ii-p_r # vectores diferencia entre dos puntos
        r_ort = np.sqrt(s[:,0]**2+s[:,1]**2)
        r_pll = np.abs(ii[2]-p_r[:,2])
        dis, x,y  = np.histogram2d(r_pll, r_ort, bins=[bn,bn], range=([[0, point_max], [0, point_max]]))
        NDR = NDR + dis
          
    
    return  NDD, NRR, NDR


#####################################################################
#####################################################################
#####################################################################


# Corremos la función 
start = time.perf_counter()

bins = 30
NDD_2d, NRR_2d, NDR_2d, xx, yy = Histos_ani(data,rand,bins,180)

finish = time.perf_counter()

print(f'Finializó en {round(finish-start,2)} segundos = {round((finish-start)/60,2)} minutos = {round((finish-start)/60,2)} horas' )


np.savetxt('/home/echeveste/Mis_trabajos/correlation_f/mis_datos/Histos_FB/DD_{0}_osc.dat'.format(str(name)), NDD)
np.savetxt('/home/echeveste/Mis_trabajos/correlation_f/mis_datos/Histos_FB/RR_{0}_osc.dat'.format(str(name)), NRR)
np.savetxt('/home/echeveste/Mis_trabajos/correlation_f/mis_datos/Histos_FB/DR_{0}_osc.dat'.format(str(name)), NDR)
