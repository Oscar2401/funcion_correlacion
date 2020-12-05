import matplotlib.pyplot as plt
import numpy as np

DD_mesh = np.loadtxt('DDiso_mesh_3D_full.dat')
RR_mesh = np.loadtxt('RRiso_mesh_3D_full.dat')

#Función de correlación Landy-Szalay
def estim_LS(NDD, NRR):
    return (NDD  - NRR)/NRR

DD_mesh[DD_mesh==0] = 0.00000000001
RR_mesh[RR_mesh==0] = 0.00000000001
eps_LS = estim_LS(DD_mesh, RR_mesh)
r = np.linspace(0,150,100)

fig = plt.figure(figsize=(14,8))
plt.scatter(r,eps_LS, s=50, c='g',label='LS')
plt.plot(r,eps_LS,'k-')
plt.ylim(-0.02,0.02)
plt.xlabel('r',fontsize=18)
plt.ylabel('$\epsilon(r)$',fontsize=18)
plt.legend(shadow=True, fontsize='x-large')
plt.grid();
plt.savefig('2PCFiso(BPC).png')
plt.show();


fig = plt.figure(figsize=(14,8))
plt.scatter(r,r**2*eps_LS, s=50, c='g',label='LS')
plt.plot(r,r**2*eps_LS,'k-')
plt.xlabel('r',fontsize=18)
plt.ylabel('$\epsilon(r)$',fontsize=18)
plt.legend(shadow=True, fontsize='x-large')
plt.grid();
plt.savefig('2PCFiso_2(BPC).png')
plt.show();

