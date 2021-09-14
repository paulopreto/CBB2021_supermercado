# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from numpy.linalg import inv
from numpy.linalg import pinv
import sys



# %% Modo Matlab
# %% Calibration DLT
def dlt_calib(cp3d, cp2d):

# =============================================================================
#                DLT 3D  
# Calcula os parametros do DLT
# para executá-la, digite os comandos abaixos
# import rec3d
# DLT = rec3d.dlt_calib(cp3d, cd2d)
# onde:
# DLT  = vetor linha com os parametros do DLT calculados
# [L1,L2,L3...L11]
# cp3d = matriz retangular com as coordenadas 3d (X, Y, Z) dos pontos (p) do calibrador
# Xp1 Yp1 Zp1
# Xp2 Yp2 Zp2
# Xp3 Yp3 Zp3
#  .   .   .
#  .   .   .
# Xpn Ypn Zpn
# 
# cp2d = matriz retangular com as coordenadas de tela (X, Y) dos pontos (p) do calibrador
# xp1 yp1
# xp2 yp2
# xp3 yp3
#  .   .
#  .   .
# xpn ypn
# =============================================================================

    cp3d = np.asarray(cp3d)
    cp2d = np.asarray(cp2d)

    m = np.size(cp3d[:, 0], 0)
    M = np.zeros([m * 2, 11])
    N = np.zeros([m * 2, 1])

    for i in range(m):
        M[i*2,:] = [cp3d[i,0], cp3d[i,1], cp3d[i,2] ,1, 0, 0, 0, 0, -cp2d[i, 0] * cp3d[i, 0], -cp2d[i, 0] * cp3d[i, 1], -cp2d[i, 0] * cp3d[i, 2]]
        M[i*2+1,:] = [0 , 0, 0, 0, cp3d[i,0], cp3d[i,1], cp3d[i,2],1, -cp2d[i,1] * cp3d[i,0],-cp2d[i,1] * cp3d[i,1], -cp2d[i,1] * cp3d[i,2]]
        N[[i*2,i*2+1],0] = cp2d[i, :]

    Mt = M.T
    M1 = inv(Mt.dot(M))
    MN = Mt.dot(N)

    DLT = (M1).dot(MN).T

    return DLT

# %% Reconstruction 3D
def r3d(DLTs, cc2ds):
    DLTs = np.asarray(DLTs)
    cc2ds = np.asarray(cc2ds)
    
    m = len(DLTs)
    M = np.zeros([2 * m, 3])
    N = np.zeros([2 * m, 1])

    for i in range(m):
        M[i*2,:] = [DLTs[i,0]-DLTs[i,8]*cc2ds[i,0], DLTs[i,1]-DLTs[i,9]*cc2ds[i,0], DLTs[i,2]-DLTs[i,10]*cc2ds[i,0]]
        M[i*2+1,:] = [DLTs[i,4]-DLTs[i,8]*cc2ds[i,1],DLTs[i,5]-DLTs[i,9]*cc2ds[i,1],DLTs[i,6]-DLTs[i,10]*cc2ds[i,1]]
        Np1 = cc2ds[i,:].T
        Np2 = [DLTs[i,3],DLTs[i,7]]
        N[[i*2,i*2+1],0] = Np1 - Np2

    cc3d = inv(M.T.dot(M)).dot((M.T.dot(N)))
    
    return cc3d

# %% Run in IDE Python
def rec3d_ide(c1=None, c2=None, ref=None):
    # Para teste na IDE
    if c1 is None:
        dfcp2d_c1 = pd.read_csv('cp2d_c1.txt', delimiter=' ',header=None)
        dfcp2d_c2 = pd.read_csv('cp2d_c2.txt', delimiter=' ',header=None)
        dfcp3d = pd.read_csv('cp3d.txt', delimiter=' ',header=None)
    else:
        dfcp2d_c1 = c1
        dfcp2d_c2 = c2
        dfcp3d = ref
            
    cp2dc1 = np.asarray(dfcp2d_c1)
    cp2dc2 = np.asarray(dfcp2d_c2)
    cp3d = np.asarray(dfcp3d)
    
    dltc1 = dlt_calib(cp3d, cp2dc1)
    dltc2 = dlt_calib(cp3d, cp2dc2)
    
    DLTs = np.append(dltc1, dltc2, axis=0)
    
    cc3d = np.zeros([len(cp2dc1), 3])
    for i in range(len(cp2dc1)):
        cc2ds = np.append([cp2dc1[i, :]], [cp2dc2[i, :]], axis=0)
        cc3d[i, :] = r3d(DLTs, cc2ds).T
    
    # print(cc3d)
    
    # np.savetxt(input('Nome do arquivo: ')+'.3d', cc3d, fmt='%.10f')
    return cc3d


   
# %% Run CMD or Terminal Shell
if __name__ == '__main__':
    resolutionx = int(720/2) ## Arrumar a translação do sist. ref do kinovea 0.9.4
    resolutiony = int(220/2)
    
    # %% Carregando os arquivos das cameras do Kinovea
    bola1 = pd.read_csv(str(sys.argv[1]), sep='\s+', skiprows=4, header=None, decimal='.')
    bola1[1] = bola1[1] - -resolutionx
    bola1[2] = -1 * (bola1[2] - resolutiony) 
    bola1 = np.asarray(bola1[[1,2]])
    
    bola2 = pd.read_csv(str(sys.argv[2]), sep='\s+', skiprows=4, header=None, decimal='.')
    bola2[1] = bola2[1] - -resolutionx
    bola2[2] = -1 * (bola2[2] - resolutiony) 
    bola2 = np.asarray(bola2[[1,2]])
    
    plt.close('all')
    plt.subplot(2,1,1)
    plt.grid(True)
    plt.plot(bola1[:,0],bola1[:,1],'o')
    plt.ylabel('CAM 1 - Coordinate Y')
    resx = 2 * resolutionx
    resy = 2 * resolutiony
    plt.title(f'Pixels coordinates (resolution = {resx} X {resy})')
    
    plt.subplot(2,1,2)
    plt.plot(bola2[:,0],bola2[:,1],'o')
    plt.xlabel('Coordinate X')
    plt.ylabel('CAM 2 - Coordinate Y')
    plt.grid(True)
    
    
    # Carregar arquivos de calibracao 
    import pdb; pdb.set_trace()
    datcal_c1 = np.asarray(pd.read_csv(str(sys.argv[3]), sep='\s+', header=None))
    datcal_c1 = datcal_c1[0] - -resolutionx

    datcal_c2 = np.asarray(pd.read_csv(str(sys.argv[4]), sep='\s+', header=None))
    datcal_c2 = -1 * (datcal_c2[1] - resolutiony) 

    
    ref = np.asarray(pd.read_csv(sys.argv[5], sep='\s+', header=None))
    ref = ref[:,1:]
    
    
    dltc1 = dlt_calib(ref, datcal_c1)
    dltc2 = dlt_calib(ref, datcal_c2)
    dlts = np.append(dltc1, dltc2, axis=0)
    
    cc3d = np.zeros([len(bola1), 3])
    
    for i in range(len(bola1)):
        cc2ds = np.append([bola1[i, :]], [bola2[i, :]], axis=0)
        cc3d[i, :] = r3d(dlts, cc2ds).T
    
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d') 
    ax.plot3D(cc3d[:,0], cc3d[:,1], cc3d[:,2], 'ro')
    ax.plot3D(ref[:,0], ref[:,1], ref[:,2], 'b.')
    ax.set_zlabel('Z [m]')
    ax.set_ylabel('Y [m]')
    ax.set_xlabel('X [m]')
    plt.title('Calib Points and 3D Ball')
    plt.show()
    
    np.savetxt(str(sys.argv[6])+'.3d', cc3d, fmt='%.10f')