
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.fft as fft
from scipy.interpolate import griddata
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp2d

# Load data
file = open('/home/chaubrich/Documents/CBC/cbc/geqdsk_gene_comp_case5_fix.eqd', 'r')

lines = file.readlines()
mr, mz, mpsi = int(lines[1][0:9]), int(lines[1][9:18]), int(lines[1][18:])
eq_min_r, eq_max_r, eq_min_z, eq_max_z = [float(num) for num in lines[2].split()]
eq_axis_r, eq_axis_z, eq_axis_b = [float(num) for num in lines[3].split()]
eq_x_psi, eq_x_r, eq_x_z = [float(num) for num in lines[4].split()]

def read_1D_grid(mpsi, lines, start_index):
    data = []
    i = 0
    while len(data) < mpsi:
        output = [float(num) for num in lines[start_index + i].split()]
        data.extend(output)
        i += 1
    return data[:mpsi], i + start_index
#psi
psi_grid, end_index_for_psi_grid = read_1D_grid(mpsi, lines, 5)
psi_grid = np.array(psi_grid)
psi_a    = np.max(psi_grid)

I_grid, end_index_for_I_grid = read_1D_grid(mpsi, lines, end_index_for_psi_grid)
I_grid = np.array(I_grid)

def read_2D_grid(mr, mz, lines, start_index):
    data = []
    i = 0
    while len(data) < mr * mz:
        output = [float(num) for num in lines[start_index + i].split()]
        data.extend(output)
        i += 1
    return data
#
psi_rz_data = read_2D_grid(mr, mz, lines, end_index_for_I_grid)
psi_rz = np.array(psi_rz_data).reshape((mr, mz))
R = np.linspace(eq_min_r,eq_max_r,mr)
Z = np.linspace(eq_min_z,eq_max_z,mz)

#################################################################################################
#Interpolate 1d grids (R,Z) 257 --> 258
prosize = np.linspace(0,1,100)
oldsize = np.linspace(0,1,257)
newsize = np.linspace(0,1,258)

R_new   = np.interp(newsize,oldsize,R)
Z_new   = np.interp(newsize,oldsize,Z)
I_new   = np.interp(newsize,oldsize,I_grid)

# Create the interpolating function
f = interp2d(R, Z, psi_rz, kind='linear')  # You can use 'cubic' or 'quadratic' as well

# Interpolate to get the new array
psi_value = f(R_new, Z_new)

psi_max    = 0.3343107087846
psi_min    = -0.1
R_min      = 1.1
Z_min      = -0.4
Z_internal = -1.2
psi_div    = 0.305
psi_a      = 0.2

# #Normalizations
# e       = 1.6e-19
# proton  = 1.67e-27
# Bu      = 2.0
# Tu      = 2.14*1000*e  #1KeV
# omegau  = e*Bu/proton
# vu      = np.sqrt(Tu/proton)
# xu      = proton*vu/(e*Bu)
# nu      = 4.66e19
# rmaj0   = 1.67/xu
# a       = rmaj0*0.36
# lref    = rmaj0
# kappat  = 6.96
# kappate = 6.96
# kappan  = 2.23
# wt      = 0.3

# Create grids
R_grid, Z_grid = np.meshgrid(R_new, Z_new)
psi = psi_value

psilevs = np.linspace(0,0.2,100)
cs = plt.contour(R_new, Z_new, psi, levels=psilevs)

plt.colorbar()
plt.grid()
#plt.show()

mask = np.zeros_like(psi)
for i in range(len(Z_new)):
    for j in range(len(R_new)):
        if psi[i][j] < psi_max and psi[i][j] > psi_min and (Z_new[i]>Z_internal or (Z_new[i]>Z_min and psi[i][j]>psi_div)) and R_new[j]>R_min:
            mask[i][j] = 1
        else:
            mask[i][j] = 0

# def interp2d(psi1d,profile,psi):
#     #psi = np.transpose(psi)
#     output = np.zeros_like(psi)
#     for x in range(len(psi)):
#         for y in range(len(psi[0])):
#             psival = psi[x][y]
#             if psival <= psi_a:
#                 miw  = int(psival/(psi1d[2]-psi1d[1]))
#                 wmx0 = ((miw+1)*(psi1d[2]-psi1d[1])-psival)/(psi1d[2]-psi1d[1])
#                 wmx1 = 1.-wmx0
#                 output[x][y] = wmx0*profile[miw] + wmx1*profile[miw+1]
#                 if (x < 70 and psival > 0.2 and psival < 0.32) or R[y] < 1:
#                     # output[x][y] = abs(psi_a-psival) + psi_a
#                     output[x][y] = np.min(profile)
#             else:
#                 # output[x][y] = np.max(profile)*0.001
#                 output[x][y] = np.min(profile)
#                 # output[x][y] = abs(psi_a-psival) + psi_a
#     return output

def bottomleft(x_pos,y_pos,R_grid,Z_grid):
    r1  = 10
    z1  = 10
    dr1 = 10
    dr2 = 10
    dz1 = 10
    dz2 = 10
    for x in range(len(R_grid)):
        if abs(R_grid[x]-x_pos) <= dr1:
            r2  = r1
            r1  = x
            dr2 = dr1
            dr1 = abs(R_grid[x]-x_pos)
        else:
            if abs(R_grid[x]-x_pos) <= dr2:
                r2  = x
                dr2 = abs(R_grid[x]-x_pos)
    for y in range(len(Z_grid)):
        if abs(Z_grid[y]-y_pos) <= dz1:
            z2  = z1
            z1  = y
            dz2 = dz1
            dz1 = abs(Z_grid[y]-y_pos)
        else:
            if abs(Z_grid[y]-y_pos) <= dz2:
                z2  = y
                dz2 = abs(Z_grid[y]-y_pos)
    if r1 < r2:
        rmin = r1
    else:
        rmin = r2
    if z1 < z2:
        zmin = z1
    else:
        zmin = z2
    return rmin, zmin

setzero = 0

idx      = []
psilev   = []
ii       = []
jj       = []
w00      = []
w01      = []
w10      = []
w11      = []
rdl      = []
deno     = []
private  = []
tester   = []

mask = np.transpose(mask)
levs = cs.levels
j=0
counter = 0
for item in cs.collections:
    for i in item.get_paths():
        v = i.vertices
        
        rxx = v[:,0]
        zyy = v[:,1]

        rx   = []
        zy   = []
        w00t = []
        w10t = []
        w01t = []
        w11t = []
        for pos in range(len(rxx)):
            #Locate nearest "bottom left" grid point then calculate weights
            ri, zj = bottomleft(rxx[pos],zyy[pos],R_new,Z_new)
            if mask[ri][zj] == 1:
                rx.append(ri)
                zy.append(zj)
                #Calculate weights for future phi interpolation
                area = (R_new[ri+1]-R_new[ri])*(Z_new[zj+1]-Z_new[zj])
                w00t.append((R_new[ri+1]-rxx[pos])*(Z_new[zj+1]-zyy[pos])/area)
                w10t.append((rxx[pos]-R_new[ri])*(Z_new[zj+1]-zyy[pos])/area)
                w01t.append((R_new[ri+1]-rxx[pos])*(zyy[pos]-Z_new[zj])/area)
                w11t.append((rxx[pos]-R_new[ri])*(zyy[pos]-Z_new[zj])/area)
            

        # print(levs[j])
        if not rx:
            # j = j-1
            print('Outside of mask, computing next step.')
            # break
        else:
            # print(levs[j])
            psilev.append(levs[j])
            private.append(0)
    
            tester.append(len(v))
            if psilev[j-1] == psilev[j-2] and len(psilev) > 1:
                print(private)
                private[j-2] = 1
            
            idx.append(int(setzero))
            setzero = setzero + 1
    
            ii.append(rx)
            jj.append(zy)
            w00.append(w00t)
            w10.append(w10t)
            w01.append(w01t)
            w11.append(w11t)            
            drdt = np.gradient(rxx)
            dzdt = np.gradient(zyy)
            pyth = np.sqrt(drdt**2 + dzdt**2)
            rdl.append(rxx*pyth)
            deno.append(np.sum(rxx*pyth))
    j = j+1

#Create two seperate tables for above and below X-point
for one in range(len(private)):
    if private[one] == 1:
        print(one)

plt.figure()
plt.plot(psilev,deno)
plt.show()

idxx  = []
psix  = []
ix    = []
jy    = []
jac   = []
denx  = []
privx = []
w00f  = []
w10f  = []
w01f  = []
w11f  = []

for rv in range(len(ii)):

    for rvv in range(len(ii[rv])):
        idxx.append(idx[rv])
        psix.append(psilev[rv])
        privx.append(private[rv])
        ix.append(ii[rv][rvv])
        jy.append(jj[rv][rvv])
        w00f.append(w00[rv][rvv])
        w10f.append(w10[rv][rvv])
        w01f.append(w01[rv][rvv])
        w11f.append(w11[rv][rvv])
        jac.append(rdl[rv][rvv])
        denx.append(deno[rv])

for fred in range(len(w00f)):
    if (w00f[fred]+w10f[fred]+w01f[fred]+w11f[fred]) < 0.5:
        print('wrong')
        print(fred)

dataset = np.column_stack((idxx,psix,ix,jy,w00f,w10f,w01f,w11f,jac,denx,privx))

jacodata_file = "~/Documents/CBC/cbc/jacodata_cbc.dat"

dfr = pd.DataFrame(dataset)
dfr.to_csv(jacodata_file, index=False, header=False, sep='\t', float_format='%.10f')
print("Data now in folder!")



