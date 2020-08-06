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
    NDDD = np.zeros((bn,bn,bn))
    NRRR = np.zeros((bn,bn,bn))
    NDDR = np.zeros((bn,bn,bn))
    NDRR = np.zeros((bn,bn,bn))
    
    # ============================================================= 
    
    n = 0
    for (ii, jj) in zip(p, p_r):
        n = n+1
        
        s = ii-p[n:] # Diferencia entre el punto pivote y los demas puntos siguientes 
        r = np.sqrt(s[:,0]**2+s[:,1]**2+s[:,2]**2)
        
        s_ran = jj-p_r[n:] # Diferencia entre el punto pivote y los demas puntos siguientes 
        r_ran = np.sqrt(s_ran[:,0]**2+s_ran[:,1]**2+s_ran[:,2]**2)
        
        m = 0
        for (kk, ll) in zip(p[n:], p_r[n:]):
            
            #DDD ==================================================
            r12 = np.ones_like(r[m+1:])*r[m]
            r13 = r[m+1:]
        
            ss = kk-p[n:][m+1:] # Diferencia entre el punto pivote y los demas puntos siguientes 
            r32 = np.sqrt(ss[:,0]**2+ss[:,1]**2+ss[:,2]**2)
        
            # Histogramas para DDD
            dis, dim = np.histogramdd(np.array([r12,r13,r32]).T,bins = [bn,bn,bn], 
                                    range=[[0,point_max],[0,point_max],[0,point_max]])
            NDDD = NDDD + dis
            
            #RRR ==================================================
            r12 = np.ones_like(r_ran[m+1:])*r_ran[m]
            r13 = r_ran[m+1:]
        
            ss = ll-p_r[n:][m+1:] # Diferencia entre el punto pivote y los demas puntos siguientes 
            r32 = np.sqrt(ss[:,0]**2+ss[:,1]**2+ss[:,2]**2)
    
            # Histogramas para RRR
            dis, dim = np.histogramdd(np.array([r12,r13,r32]).T,bins =[bn,bn,bn], 
                                    range=[[0,point_max],[0,point_max],[0,point_max]])
            NRRR = NRRR + dis  
            
            m = m+1
    
    # =========================================================================================== 
    
    n = 0
    for (ii, jj) in zip(p, p_r):
        n = n+1
        
        s = ii-p[n:] # Diferencia entre el punto pivote y los demas puntos siguientes 
        r = np.sqrt(s[:,0]**2+s[:,1]**2+s[:,2]**2)
        
        s_ran = jj-p_r[n:] # Diferencia entre el punto pivote y los demas puntos siguientes 
        r_ran = np.sqrt(s_ran[:,0]**2+s_ran[:,1]**2+s_ran[:,2]**2)
        
        m = 0
        for (kk, ll) in zip(p[n:], p_r[n:]):
            
            #DDR ==================================================
            r12 = np.ones_like(p_r.T[0])*r[m]
            
            ss = ii-p_r # Diferencia entre el punto pivote y el punto en la muestra random 
            r13 = np.sqrt(ss[:,0]**2+ss[:,1]**2+ss[:,2]**2)
            
            ss = kk-p_r # Diferencia entre el punto 2 de los datos y el punto en la muestra random 
            r32 = np.sqrt(ss[:,0]**2+ss[:,1]**2+ss[:,2]**2)
            
            # Histogramas para DDR
            dis, dim = np.histogramdd(np.array([r12,r13,r32]).T,bins = [bn,bn,bn], 
                                    range=[[0,point_max],[0,point_max],[0,point_max]])
            NDDR = NDDR + dis
            
            #DRR ==================================================
            r12 = np.ones_like(p.T[0])*r_ran[m]
            
            ss = jj-p # Diferencia entre el punto pivote random y el punto en los datos 
            r13 = np.sqrt(ss[:,0]**2+ss[:,1]**2+ss[:,2]**2)
            
            ss = ll-p # Diferencia entre el punto 2 de la muestra random y el punto en los datos 
            r32 = np.sqrt(ss[:,0]**2+ss[:,1]**2+ss[:,2]**2)
            
            # Histogramas para DRR
            dis, dim = np.histogramdd(np.array([r12,r13,r32]).T,bins = [bn,bn,bn], 
                                    range=[[0,point_max],[0,point_max],[0,point_max]])
            NDRR = NDRR + dis
            
            m = m+1
    
    return dim, NDDD, NRRR, NDDR, NDRR


#####################################################################
#####################################################################
#####################################################################


# Corremos la función 
start = time.perf_counter()

bins = 30
dim, NDDD, NRRR, NDDR, NDRR = Histos(data,rand,bins,180)

finish = time.perf_counter()

print(f'Finializó en {round(finish-start,2)} segundos = {round((finish-start)/60,2)} minutos = {round((finish-start)/3600,2)} horas' )

# Simetrizamos los Histogramas:

for i in range(bins):
    for j in range(i, bins):
        for k in range(j, bins):
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
            
for i in range(bins):
    for j in range(i, bins):
        for k in range(j,bins):
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

np.savetxt('/home/echeveste/mis_trabajos/correlation_f_2/datos/histogramas/BF_oscar_py/DDD_iso_{0}_osc.dat'.format(str(name)), NDDD.reshape((bins,bins*bins)))
np.savetxt('/home/echeveste/mis_trabajos/correlation_f_2/datos/histogramas/BF_oscar_py/RRR_iso_{0}_osc.dat'.format(str(name)), NRRR.reshape((bins,bins*bins)))
np.savetxt('/home/echeveste/mis_trabajos/correlation_f_2/datos/histogramas/BF_oscar_py/DDR_iso_{0}_osc.dat'.format(str(name)), NDDR.reshape((bins,bins*bins))
np.savetxt('/home/echeveste/mis_trabajos/correlation_f_2/datos/histogramas/BF_oscar_py/DRR_iso_{0}_osc.dat'.format(str(name)), NDRR.reshape((bins,bins*bins)))

