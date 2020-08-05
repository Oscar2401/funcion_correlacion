import sys 
import numpy as np
import time

data = np.loadtxt(sys.argv[1])[:,:3]
rand = np.loadtxt(sys.argv[2])[:,:3]
name = sys.argv[3]

#FUNCIÓN DE PARA HACER HISTOGRAMAS 
def Histos(p,p_r,bn,point_max):
    """ 
    Función para construir los histogramas 
    
    p = datos
    p_r = random
    bn = tamaño de bins
    point_max = punto máximo en el histograma
    
    """
    
    #Inicializamos los arreglos de los histogramas
    NDD = np.zeros(bn)
    NRR = np.zeros(bn)
    NDR = np.zeros(bn)
    
    n = 0
    
    for (ii, jj) in zip(p, p_r):
        n = n+1
        
        # Histogramas para DD
        s = ii-p[n:] # Diferencia entre el punto pivote y los demas puntos siguientes 
        dis, r = np.histogram(np.sqrt(s[:,0]**2+s[:,1]**2+s[:,2]**2), bins=bn, range=(0, point_max))
        NDD = NDD + 2*dis
        
        # Histogramas para RR
        s = jj-p_r[n:]
        dis, r = np.histogram(np.sqrt(s[:,0]**2+s[:,1]**2+s[:,2]**2), bins=bn, range=(0, point_max))
        NRR = NRR + 2*dis   
    
    for ii in p:
        # Histogramas para DR
        s = ii-p_r
        dis, r = np.histogram(np.sqrt(s[:,0]**2+s[:,1]**2+s[:,2]**2), bins=bn, range=(0, point_max))
        NDR = NDR + dis
    
    return NDD, NRR, NDR


#####################################################################
#####################################################################
#####################################################################


# Corremos la función 
start = time.perf_counter()

bins = 30
NDD, NRR, NDR = Histos(data,rand,bins,180)

finish = time.perf_counter()

print(f'Finializó en {round(finish-start,2)} segundos = {round((finish-start)/60,2)} minutos = {round((finish-start)/60,2)} horas' )


np.savetxt('/home/echeveste/Mis_trabajos/correlation_f/mis_datos/Histos_FB/DD_{0}_osc.dat'.format(str(name)), NDD)
np.savetxt('/home/echeveste/Mis_trabajos/correlation_f/mis_datos/Histos_FB/RR_{0}_osc.dat'.format(str(name)), NRR)
np.savetxt('/home/echeveste/Mis_trabajos/correlation_f/mis_datos/Histos_FB/DR_{0}_osc.dat'.format(str(name)), NDR)

