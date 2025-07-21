
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.fft as fft
from scipy.interpolate import griddata
from scipy.interpolate import CubicSpline

# Load data
R         = np.loadtxt('/home/chaubrich/Documents/lcfsdata/R.dat')
Z         = np.loadtxt('/home/chaubrich/Documents/lcfsdata/Z.dat')
psi_value = np.loadtxt('/home/chaubrich/Documents/lcfsdata/psi_p.dat')
psi_cont  = np.loadtxt('/home/chaubrich/Documents/lcfsdata/psi_contour.dat')

pfpath    = 'TokamakFiles/DIIID/162940/p162940.02944_670'


plt.figure()
plt.contour(R,Z,psi_value,100)
plt.colorbar()
plt.show()

psi_max    = 0.34
psi_min    = -0.1
R_min      = 1.0
Z_min      = -1.5
Z_internal = -1.2
psi_div    = 0.305
psi_a      = 0.311647

#Normalizations
e       = 1.6e-19
proton  = 1.67e-27
Bu      = 2.0
Tu      = 2.14*1000*e  #1KeV
omegau  = e*Bu/proton
vu      = np.sqrt(Tu/proton)
xu      = proton*vu/(e*Bu)
nu      = 4.66e19
rmaj0   = 1.67/xu
a       = rmaj0*0.36
lref    = rmaj0
kappat  = 6.96
kappate = 6.96
kappan  = 2.23
wt      = 0.3

# Create grids
R_grid, Z_grid = np.meshgrid(R, Z)
psi = psi_value

psilevs = np.linspace(0,0.34,100)
cs = plt.contour(R, Z, psi, levels=psilevs)

plt.colorbar()
plt.grid()
#plt.show()


mask = np.zeros_like(psi)
for i in range(len(Z)):
    for j in range(len(R)):
        if psi[i][j] < psi_max and psi[i][j] > psi_min and (Z[i]>Z_internal or (Z[i]>Z_min and psi[i][j]>psi_div)) and R[j]>R_min:
            mask[i][j] = 1
        else:
            mask[i][j] = 0

# print(np.unravel_index(psi.argmin(), psi.shape))
# print(psi[254][267])
# print(psi[268][254])

plt.contour(R,Z,mask)
# plt.scatter(R[254],Z[267])
# plt.scatter(R[268],Z[254])
# plt.show()

def readpfile(pfpath):
    profiles = {}

    with open(pfpath, 'r') as ofh:
        while True:
            recs = ofh.readline().split()  # Reads Heading Line
            if len(recs) == 4:
                nrec, var0, temp = int(recs[0]), str(recs[1]).lower(), str(recs[2]).lower() # Sets resolution size, radial unit, profile unit
                var1 = temp[:temp.index("(")]  # Removes Unit Values
                ary0, ary1 = [], []
                for _ in range(nrec):
                    recs = ofh.readline().split()
                    ary0.append(float(recs[0]))
                    ary1.append(float(recs[1]))
                if var0 in profiles.keys():
                    CS = CubicSpline(ary0,ary1) 
                    profiles[var1] = CS(profiles[var0]) # Makes sure profiles are smooth and match resolution size
                else:
                    profiles[var0], profiles[var1] = ary0, ary1
            elif len(recs) > 4:
                nrec, var0, var1, var2 = int(recs[0]), str(recs[1]).upper(), str(recs[2]).upper(), str(recs[3]).upper() # Sets atomic resolution, number, charge, and mass respectivly
                ary0, ary1, ary2 = [], [], []
                for _ in range(nrec):
                    recs = ofh.readline().split()
                    ary0.append(int(float(recs[0])))
                    ary1.append(float(recs[1]))
                    ary2.append(float(recs[2]))
                profiles[var0], profiles[var1], profiles[var2] = ary0, ary1, ary2
            else:
                break
    return profiles

def interp2d(psi1d,profile,psi):
    #psi = np.transpose(psi)
    output = np.zeros_like(psi)
    for x in range(len(psi)):
        for y in range(len(psi[0])):
            psival = psi[x][y]
            if psival <= psi_a:
                miw  = int(psival/(psi1d[2]-psi1d[1]))
                wmx0 = ((miw+1)*(psi1d[2]-psi1d[1])-psival)/(psi1d[2]-psi1d[1])
                wmx1 = 1.-wmx0
                output[x][y] = wmx0*profile[miw] + wmx1*profile[miw+1]
                if (x < 70 and psival > 0.2 and psival < 0.32) or R[y] < 1:
                    # output[x][y] = abs(psi_a-psival) + psi_a
                    output[x][y] = np.min(profile)
            else:
                # output[x][y] = np.max(profile)*0.001
                output[x][y] = np.min(profile)
                # output[x][y] = abs(psi_a-psival) + psi_a
    return output

#TEMP. interpolation for profile for GEM-X#
#########################################################################################################################################
profiles = readpfile(pfpath)
psi1d    = np.array(profiles['psinorm'])*psi_a
te2d     = interp2d(psi1d,np.array(profiles['te']),psi)
ti2d     = interp2d(psi1d,np.array(profiles['ti']),psi)
ne2d     = interp2d(psi1d,np.array(profiles['ne']),psi)
ni2d     = interp2d(psi1d,np.array(profiles['ni']),psi)
#dti = interp2d(psi1d,np.array(profiles['dti/dpsiN']),psi)

#Manufactured temperature profiles for diagnostics
# tific1d = np.exp(-kappat*wt*(a/lref)*np.tanh(psi1d/(wt*a)))
# tefic1d = np.exp(-kappate*wt*(a/lref)*np.tanh(psi1d/(wt*a)))
# tific1d = 0.5*(-np.tanh(8*np.array(profiles['psinorm'])-6.5)+1)*np.max(np.array(profiles['ti'])) #ORIG: 0.5  5  2.5  1
# tefic1d = 0.5*(-np.tanh(8*np.array(profiles['psinorm'])-6.5)+1)*np.max(np.array(profiles['te']))
tific1d = 0.5*(-np.tanh(5*np.array(profiles['psinorm'])-2.5)+1)*np.max(np.array(profiles['ti'])) #ORIG: 0.5  5  2.5  1
tefic1d = 0.5*(-np.tanh(5*np.array(profiles['psinorm'])-2.5)+1)*np.max(np.array(profiles['te']))

tific = interp2d(psi1d,tific1d,psi)
tefic = interp2d(psi1d,tefic1d,psi)
grati = abs(np.gradient(tific1d))
print(psi1d[np.argmax(grati)])

plt.figure()
plt.plot(psi1d,tific1d,linewidth=3.0)
plt.xlabel('$\\psi$',fontsize=20)
plt.ylabel('T [kev]',fontsize=20)
plt.grid()
plt.show()

plt.figure()
plt.plot(profiles['psinorm'],np.array(profiles['ni'])*1E20,linewidth=2.0,label='i')
plt.plot(profiles['psinorm'],np.array(profiles['ne'])*1E20,linewidth=2.0,label='e')
plt.xlabel('$\\rho_{miller}$',fontsize=20)
plt.ylabel('$n$ $(m^{-3})$',fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=16)
plt.tight_layout()
plt.grid()
plt.show()

plt.figure()
plt.plot(profiles['psinorm'],np.array(profiles['ti']),linewidth=2.0,label='$T_i$')
plt.plot(profiles['psinorm'],np.array(profiles['te']),linewidth=2.0,label='$T_e$')
plt.xlabel('$\\rho_{miller}$',fontsize=20)
plt.ylabel('$T(keV)$',fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=16)
plt.tight_layout()
plt.grid()
plt.show()

plt.figure()
plt.contourf(R,Z,ti2d,100)
plt.colorbar()
plt.show()

# plt.plot(profil)

#send 2d arrays to folder
te0 = "~/Documents/lcfsdata/te0_profile.dat"
ti0 = "~/Documents/lcfsdata/ti0_profile.dat"
ne0 = "~/Documents/lcfsdata/ne0_profile.dat"
ni0 = "~/Documents/lcfsdata/ni0_profile.dat"
tif = "~/Documents/lcfsdata/tif_profile.dat"
tef = "~/Documents/lcfsdata/tef_profile.dat"

te = pd.DataFrame(te2d)
ti = pd.DataFrame(ti2d)
ne = pd.DataFrame(ne2d)
ni = pd.DataFrame(ni2d)

tiff = pd.DataFrame(tific)
teff = pd.DataFrame(tefic)

te.to_csv(te0, index=False, header=False, sep='\t', float_format='%.10f')
ti.to_csv(ti0, index=False, header=False, sep='\t', float_format='%.10f')
ne.to_csv(ne0, index=False, header=False, sep='\t', float_format='%.10f')
ni.to_csv(ni0, index=False, header=False, sep='\t', float_format='%.10f')
tiff.to_csv(tif, index=False, header=False, sep='\t', float_format='%.10f')
teff.to_csv(tef, index=False, header=False, sep='\t', float_format='%.10f')
print("Data Arrays in folder now in folder!")


#########################################################################################################################################

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
            ri, zj = bottomleft(rxx[pos],zyy[pos],R,Z)
            if mask[ri][zj] == 1:
                rx.append(ri)
                zy.append(zj)
                #Calculate weights for future phi interpolation
                area = (R[ri+1]-R[ri])*(Z[zj+1]-Z[zj])
                w00t.append((R[ri+1]-rxx[pos])*(Z[zj+1]-zyy[pos])/area)
                w10t.append((rxx[pos]-R[ri])*(Z[zj+1]-zyy[pos])/area)
                w01t.append((R[ri+1]-rxx[pos])*(zyy[pos]-Z[zj])/area)
                w11t.append((rxx[pos]-R[ri])*(zyy[pos]-Z[zj])/area)
            

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

# print(len(levs))
# print(psilev)
# print(private)


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

# plt.figure()
# plt.scatter(ix,jy,s=1)
# plt.show()

for fred in range(len(w00f)):
    if (w00f[fred]+w10f[fred]+w01f[fred]+w11f[fred]) < 0.5:
        print('wrong')
        print(fred)

dataset = np.column_stack((idxx,psix,ix,jy,w00f,w10f,w01f,w11f,jac,denx,privx))

jacodata_file = "~/Documents/lcfsdata/jacodata.dat"

dfr = pd.DataFrame(dataset)
dfr.to_csv(jacodata_file, index=False, header=False, sep='\t', float_format='%.10f')
print("Data now in folder!")



