import matplotlib.pyplot as plt
import numpy as np
import scipy.misc
from scipy import ndimage
import matplotlib.gridspec as gridspec
import cv2

#Función de correlación Landy-Szalay
def estim_LS(NDD, NRR, NDR):
    return (NDD - 2*NDR + NRR)/NRR

def imag(x,y,t,cmap,limt,limt_):
    
    plt.figure(figsize=(6,6), dpi=100)
    plt.imshow(x,origin='lower',cmap=cmap, extent=[0,150,0,150],
               interpolation= 'bilinear', vmin=limt_, vmax=limt)
    cax=plt.colorbar()
    #plt.contour(x,10,cmap=plt.cm.gray,linewidths=1
    #            ,extent=[0,180,0,180],vmin=-limt, vmax=limt)

    plt.ylabel('$r_{\pi}$',fontsize = 16)
    plt.xlabel('$r_{p}$',fontsize = 16)
    plt.title(y,fontsize = 16)
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax.set_label(t,labelpad = 15,fontsize = 15)
    line_colour1 = ('royalblue', 'blue', 'mediumblue', 'darkblue')
    plt.plot()

DD = np.loadtxt('DDani_mesh_3D_full.dat')
RR = np.loadtxt('RRani_mesh_3D_full.dat')
DR = np.loadtxt('DRani_mesh_3D_full.dat')

eps_LS = estim_LS(DD,RR,DR)
#eps_LS[eps_LS<1] = 1
#eps_LS = np.log10(eps_LS)
eps_LS =  cv2.GaussianBlur(eps_LS,(3,3),0.5)

limt = 0.02
limt_ = -limt

imag(eps_LS,'Función de Correlación','$\epsilon(r)$',cmap='RdBu', limt=limt, limt_=limt_)
plt.savefig('2PCFani.png')

plt.figure(figsize = (7,7))
gs1 = gridspec.GridSpec(2, 2)
gs1.update(wspace=0.0, hspace=0.0) # set the spacing between axes. 

ax1 = plt.subplot(gs1[0])
rotated_img = ndimage.rotate(np.rot90(eps_LS), 90)
plt.imshow(rotated_img, cmap='RdBu',interpolation= 'bilinear', vmin=limt_, vmax=limt)
plt.contour(rotated_img,10,cmap=plt.cm.gray,linewidths=1, vmin=limt_, vmax=limt)
plt.axis('off')
ax1 = plt.subplot(gs1[1])
rotated_img = ndimage.rotate(eps_LS.T,90)
plt.imshow(rotated_img, cmap='RdBu',interpolation= 'bilinear', vmin=limt_, vmax=limt)
plt.contour(rotated_img,10,cmap=plt.cm.gray,linewidths=1,  vmin=limt_, vmax=limt)
plt.axis('off')
ax1 = plt.subplot(gs1[2])
rotated_img = ndimage.rotate(eps_LS.T,-90)
plt.imshow(rotated_img, cmap='RdBu',interpolation= 'bilinear', vmin=limt_, vmax=limt)
plt.contour(rotated_img,10,cmap=plt.cm.gray,linewidths=1,  vmin=limt_, vmax=limt)
plt.axis('off')
ax1 = plt.subplot(gs1[3])
rotated_img = ndimage.rotate(eps_LS, 0)
plt.imshow(rotated_img, cmap='RdBu',interpolation= 'bilinear', vmin=limt_, vmax=limt)
plt.contour(rotated_img,10,cmap=plt.cm.gray,linewidths=1,  vmin=limt_, vmax=limt)
plt.axis('off')
plt.savefig('2PCFani_2.png')
plt.show()


