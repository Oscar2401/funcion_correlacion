import matplotlib.pyplot as plt
import numpy as np
import random
import time
import sys 
import math
#Libreria para programar en parallelo
from joblib import Parallel, delayed

#=================================================================
#==================== Funciones a usar ===========================
#=================================================================

#=================================================================

def make_mesh(data, size):
	"""
	Función para dividir los datos en nodos de un malla.
	
	data = Arreglo de puntos a dividir.
	size = tamaño de los nodos que componen la malla.
	"""
	# Límites de los nodos
	x = [i for i in range(math.floor(np.min(data[:,0])), math.ceil(np.max(data[:,0]))+size, size)]
	y = [i for i in range(math.floor(np.min(data[:,1])), math.ceil(np.max(data[:,1]))+size, size)]
	# Hacemos las mallas
	mesh = []
	for i in range(len(x)-1):
		for j in range(len(y)-1):
			mesh.append(data[(data[:,0]>=x[i])&(data[:,0]<x[i+1])&(data[:,1]>=y[j])&(data[:,1]<y[j+1])])
	# Puntos medios en cada nodo
	pt_m = []	
	d_pt_m = size/2   # Distancia a punto medio 
	for i in range(int(np.sqrt(len(mesh)))):
		for j in range(int(np.sqrt(len(mesh)))):
			pt_m.append([d_pt_m+i*size, d_pt_m+j*size])
	pt_m = np.array(pt_m)
	return mesh, pt_m, d_pt_m

#=================================================================   

def histo_XX_pll(data, bn, d_max):
	"""
	Función para hacer histogramas de distancias entre un mismo arreglo de puntos.
	
	data = Arreglo de puntos a dividir. 
	bn = bin de histograma.
	d_max = distancia maxima a medir.
	"""
	XX = np.zeros(bn)
	n = 0
	for ii in data:
		n = n+1
		s = ii-data[n:]
		dis, r = np.histogram(np.sqrt(s[:,0]**2+s[:,1]**2), bins=bn, range=(0, d_max))
		XX = XX + 2*dis
	return  XX
	
#=================================================================

def histo_XY(p, p_r, bn, point_max):
	NDR = np.zeros(bn)
	for ii in p:
		s = ii-p_r
		dis, r = np.histogram(np.sqrt(s[:,0]**2+s[:,1]**2), bins=bn, range=(0, point_max))
		NDR = NDR + 2*dis
	return NDR

#=================================================================

def histo_X_XY_pll(i, pt_m, d_pt_m, mesh, bn, d_max):
	"""
	Función para hacer histogramas de distancias entre un mismo arreglo de puntos usando cómputo en parallelo.
	
	i = ciclo.
	pt_m = Arreglo de puntos medios de cada nodo. 
	d_pt_m = mitad del la longitud del nodo.
	mesh = malla.
	bn = bin de histograma.
	d_max = distancia maxima a medir.
	"""	
	# Aproximación de distancia entre nodos con sus puntos medios
	s = pt_m[i]-pt_m[i+1:] 
	dx = pt_m[i][0]-pt_m[i+1:][:,0]
	dx[dx==0] = 0.000001 
	dy = pt_m[i][1]-pt_m[i+1:][:,1]
	m = np.abs(dy/dx)
	d_nodo = np.zeros_like(s[:,0])
	d_nodo[m<1] = np.sqrt(s[m<1][:,0]**2+s[m<1][:,1]**2)*(1-(2*d_pt_m/np.abs(dx[m<1])))
	d_nodo[m>1] = np.sqrt(s[m>1][:,0]**2+s[m>1][:,1]**2)*(1-(2*d_pt_m/np.abs(dy[m>1])))
	d_nodo[m==1.] = np.sqrt(s[m==1.][:,0]**2+s[m==1.][:,1]**2)-(2*d_pt_m*np.sqrt(2))    
	
	u = mesh[i] # Tomamos el i-esimo nodo de la malla  
	v = np.concatenate(np.array(mesh[i+1:],dtype=object)[d_nodo < d_max]) # Tomamos los siguientes nodos y hacemos un arreglo único
	
	XY = np.zeros(bn)
	for ii in u:
		s = ii-v
		dis, r = np.histogram(np.sqrt(s[:,0]**2+s[:,1]**2), bins=bn, range=(0, d_max))
		XY = XY + 2*dis
	return XY

#=================================================================

def histo_mesh_XX_pll(data, size, bn, d_max, nuc=2):
	"""
	Funcion principal para hecer histogramas XX de forma paralelizada
	
	data = datos.
	size = tamaño del lado de cada nodo.
	bn = bin del histograma.
	d_max = distancia maxima a medir.	
	"""
	mesh, pt_m, d_pt_m = make_mesh(data, size)
	# Distacnias entre puntos del mismo nodo
	XX = np.sum(np.array(Parallel(n_jobs=nuc) #número de nucleos
				(delayed(histo_XX_pll) #Función 
				(ii, bn, d_max) for ii in mesh)), axis = 0 ) #Argumentos
	# Diastancias entre puntso de distinto nodo 
	XX = XX + np.sum(np.array(Parallel(n_jobs=nuc) #número de nucleos
					(delayed(histo_X_XY_pll) #Función para paralelizar
					(i, pt_m, d_pt_m,
					mesh, bn, d_max) for i in range(len(pt_m)-2))), axis = 0 ) #Argumentos
	#Distancia entre el último nodo 			
	XX = XX + histo_XY(mesh[len(pt_m)-2],mesh[len(pt_m)-1],bn,d_max)
	return XX

#=================================================================
	
def histo_XY_pll(i, pt_m, pt_m_r, d_pt_m, mesh, mesh_r, bn, point_max):
	"""
	Función para hacer histogramas de distancias entre un mismo arreglo de puntos usando cómputo en parallelo.
	
	i = ciclo.
	pt_m = Arreglo de puntos medios de cada nodo de datos. 
	pt_m_r = Arreglo de puntos medios de cada nodo de los puntos aleatorios. 
	d_pt_m = mitad del la longitud del nodo.
	mesh = malla de datos.
	mesh_r = malla de puntos aleatorios.
	bn = bin de histograma.
	d_max = distancia maxima a medir.
	"""	
	s = pt_m[i]-pt_m_r 
	dx = pt_m[i][0]-pt_m_r[:,0]
	dx[dx==0] = 0.000001 
	dy = pt_m[i][1]-pt_m_r[:,1]
	m = np.abs(dy/dx)
	m[np.isnan(m)] = 0
	d_nodo = np.zeros_like(s[:,0])
	d_nodo[m<1] = np.sqrt(s[m<1][:,0]**2+s[m<1][:,1]**2)*(1-(2*d_pt_m/np.abs(dx[m<1])))
	d_nodo[m>1] = np.sqrt(s[m>1][:,0]**2+s[m>1][:,1]**2)*(1-(2*d_pt_m/np.abs(dy[m>1])))
	d_nodo[m==1.] = np.sqrt(s[m==1.][:,0]**2+s[m==1.][:,1]**2)-(2*d_pt_m*np.sqrt(2))    
	
	u = mesh[i]
	v = np.concatenate(np.array(mesh_r,dtype=object)[d_nodo < point_max])
	
	XY = np.zeros(bn)
	for ii in u:
		s = ii-v
		dis, r = np.histogram(np.sqrt(s[:,0]**2+s[:,1]**2), bins=bn, range=(0, point_max))
		XY = XY + 2*dis
	return XY
    
#=================================================================
    
def histo_mesh_XY_pll(data,rand,size,bn,point_max,nuc=2):
	"""
	Funcion principal para hecer histogramas XY de forma paralelizada
	
	data = datos.
	rand = puntos aleatorios
	size = tamaño del lado de cada nodo.
	bn = bin del histograma.
	d_max = distancia maxima a medir.	
	"""
	malla_D, malla_R =  Parallel(n_jobs=2)(delayed(make_mesh)(ii,size) for ii in [data,rand])        
	
	XY = np.sum(np.array(Parallel(n_jobs=nuc) #número de nucleos
					(delayed(histo_XY_pll) #Función 
					(i,malla_D[1],malla_R[1],malla_D[2],malla_D[0],
					malla_R[0],bn,point_max) for i in range(len(malla_D[1])))),axis=0) #Argumentos
	return XY


#=================================================================
#===================== Tarea principal ===========================
#=================================================================

data = np.loadtxt(sys.argv[1])[:,:2]
rand = np.loadtxt(sys.argv[2])[:,:2]
name = sys.argv[3]

size = 20

start = time.perf_counter()

DD_pll = histo_mesh_XX_pll (data,size,30,180)
RR_pll = histo_mesh_XX_pll (rand,size,30,180)
DR_pll = histo_mesh_XY_pll (data,rand,size,30,180)

finish = time.perf_counter()
print(f'Finializó en {round(finish-start,2)} segundos')

np.savetxt('/home/echeveste/mis_trabajos/correlation_f_2/datos/histogramas/mallas_Oscar_py/DD_pll_iso_{0}_osc.dat'.format(str(name)), DD_pll)
np.savetxt('/home/echeveste/mis_trabajos/correlation_f_2/datos/histogramas/mallas_Oscar_py/RR_pll_iso_{0}_osc.dat'.format(str(name)), RR_pll)
np.savetxt('/home/echeveste/mis_trabajos/correlation_f_2/datos/histogramas/mallas_Oscar_py/DR_pll_iso_{0}_osc.dat'.format(str(name)), DR_pll)
