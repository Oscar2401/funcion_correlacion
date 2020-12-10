import matplotlib.pyplot as plt
import numpy as np

#Función de correlación Landy-Szalay
def estim_LS(NDDD, NRRR, NDDR, NDRR):
    return (NDDD - 3*NDDR + 3*NDRR - NRRR)/NRRR

DDD = np.loadtxt('DDDiso_mesh_3D_10k.dat')
RRR = np.loadtxt('RRRiso_mesh_3D_10k.dat')
DDR = np.loadtxt('DDRiso_mesh_3D_10k.dat')
DRR = np.loadtxt('DRRiso_mesh_3D_10k.dat')

bn = 40

DDD = np.reshape(DDD, (bn,bn,bn))
RRR = np.reshape(RRR, (bn,bn,bn))
DDR = np.reshape(DDR, (bn,bn,bn))
DRR = np.reshape(DRR, (bn,bn,bn))

eps_LS = estim_LS(DDD,RRR,DDR,DRR)
eps_LS[np.where(np.isnan(eps_LS))] = 0.0
eps_LS[np.where(np.isinf(eps_LS))] = 0.0

limt = 0.01
limt_ = -limt

d_max  = 60
for r in range(bn):
	plt.figure(figsize=(6,6), dpi=100)
	plt.imshow(eps_LS[r],origin='lower',cmap='RdBu', interpolation= 'bilinear', extent=[0,d_max,0,d_max], vmin=-limt, vmax=limt)
	cax=plt.colorbar()
	plt.contour(eps_LS[r],10,cmap=plt.cm.gray,linewidths=1,extent=[0,d_max,0,d_max],vmin=-limt, vmax=limt)
	plt.ylabel('$r_1$',fontsize = 16)
	plt.xlabel('$r_2$',fontsize = 16)
	plt.title('3PCF',fontsize = 16)
	plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
	cax.set_label('$\epsilon(r_1,r_2,r_3={0})$'.format(r*(d_max)/bn),labelpad = 15,fontsize = 15)
	plt.savefig('graphics_3PCF/3PCFiso{0}.png'.format(r))
	#plt.show();

