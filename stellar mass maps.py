# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 09:31:55 2022

@author: shivani thakur 
"""


def get(path, params=None):
     # make HTTP GET request to path
     headers = {"api-key":"d4941cc15f92705f6596c8a61e31ce47"}
     r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
     r.raise_for_status()
     
     if r.headers['content-type'] == 'application/json':
         return r.json() # parse json responses automatically

     if 'content-disposition' in r.headers:
         filename = r.headers['content-disposition'].split("filename=")[1]
         with open(filename, 'wb') as f:
             f.write(r.content)
         return filename # return the filename string

     return r

#---------------------------------------------------------------------------------------------------------
"""
Obtains the hypotenuse or the magnitude of a 3D vector passed to it.

INPUT : X coordinate, Y coordinate, Z coordinate
Returns : hyp -> hypotenuse/magnitude (Float)
"""
def hyp3D(x, y, z): 
    hyp = np.sqrt(x**2 + y**2 + z**2)
    return hyp

#--------------------------------------------------------------------------------------------------------
"""
Creates a 3D Plot of the 3D vector data given. 

INPUT : SubFind ID, X data array, Y data array, Z data array
Returns : None
"""
def plot3D(id, xdata, ydata, zdata, sfrdata): #fignum, colour,

    fig = plt.figure() #fignum
    ax = plt.axes(projection='3d')
    p = ax.scatter3D(xdata, ydata, zdata, c=sfrdata, cmap='viridis', s = 0.5, vmin = 0, vmax = 0.001)
    fig.colorbar(p, ax=ax)
    ax.set_title('SubHalo {0} Visualized in 3D'.format(id))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

#-------------------------------------------------------------------------------------------------------
"""
Returns the redshift associated with a given snapshot number from a txt file 
(have cols : snapNo - snapshot number & redshiftNo - redshift number) 

INPUT : Snapshot Number 
Returns : Redshift associated with the snapshot (Float)
"""
def redshift(snapnum):
    snapNo, redshiftNo = np.loadtxt('snapshot_redshift_num.txt', unpack=True, usecols = (0,1))
    redIndex = np.where(snapNo == snapnum) #finds the index for given snapnum
    z = redshiftNo[redIndex]               #assigns z the value in the redshift col at the found index
    return z
    
#------------------------------------------------------------------------------------------------------fd
"""
Calculates the weighted average, given : values (xdata) and data of the weights (weightdata).
Method : Finding the weight multiplied with the corresponding values. 
Obtaining the fraction of weighted values sum and weight sum. 

Returns : totalWeightedAvg (Float)
"""
def weighted_avg(xdata, weightdata):
    particles = len(xdata)
    weighted_avg = np.zeros(particles)
    for i in np.arange(particles):
        weighted_avg[i] = xdata[i] * (weightdata[i])
    
    WeightedAvg = np.sum(weighted_avg)
    sum_weight = np.sum(weightdata)
    totWeightedAvg = WeightedAvg / sum_weight
    return totWeightedAvg
#------------------------------------------------------------------------------------------------------
"""
Calculates radius for the given mass fraction of the total stellar mass.
"""
def XmassRad (MassFrac, distArr, arrmass, totMass_star):
    #Converting r array to list
    rlist = distArr.tolist()
    
    #initializing our variables before the while loop
    totMass = 0
    ReqMass = MassFrac * totMass_star
    
    while (totMass <= ReqMass):
        min_index = rlist.index(np.nanmin(rlist))     #finding the index of the min value in r
        totMass += arrmass[min_index]     #adding the masses of the min r values 
        if(totMass > ReqMass):
            Xmass_rad = rlist[min_index]
        #print('total mass', totMass)    
        rlist[min_index] = np.nan         #replacing the previous min value with NaN
        
    print(MassFrac, " total stellar mass rad with calc => ", Xmass_rad)
    
    return Xmass_rad

#------------------------------------------------------------------------------------------------------
def momentOfInertiaTensor(distarr, x, y, z, arr_mass, halfmassrad):
    """ Calculate the moment of inertia tensor (3x3 matrix) for a subhalo-scope particle set. """

    # load required particle data for this subhalo

    wGas = np.where( (distarr <= halfmassrad) ) #& (arr_starform > 0.0) 

    masses = arr_mass[wGas]
    x = x[wGas]
    y = y[wGas]
    z = z[wGas]
    
    # construct moment of inertia
    I = np.zeros( (3,3), dtype='float32' )

    I[0,0] = np.sum( masses * (y*y + z*z) )
    I[1,1] = np.sum( masses * (x*x + z*z) )
    I[2,2] = np.sum( masses * (x*x + y*y) )
    I[0,1] = -1 * np.sum( masses * (x*y) )
    I[0,2] = -1 * np.sum( masses * (x*z) )
    I[1,2] = -1 * np.sum( masses * (y*z) )
    I[1,0] = I[0,1]
    I[2,0] = I[0,2]
    I[2,1] = I[1,2]

    return I


def rotationMatricesFromInertiaTensor(I):
    """ Calculate 3x3 rotation matrix by a diagonalization of the moment of inertia tensor.
    Note the resultant rotation matrices are hard-coded for projection with axes=[0,1] e.g. along z. """

    # get eigen values and normalized right eigenvectors
    eigen_values, rotation_matrix = np.linalg.eig(I)

    # sort ascending the eigen values
    sort_inds = np.argsort(eigen_values)
    eigen_values = eigen_values[sort_inds]

    # permute the eigenvectors into this order, which is the rotation matrix which orients the
    # principal axes to the cartesian x,y,z axes, such that if axes=[0,1] we have face-on
    new_matrix = np.matrix( (rotation_matrix[:,sort_inds[0]],
                             rotation_matrix[:,sort_inds[1]],
                             rotation_matrix[:,sort_inds[2]]) )
    
    return new_matrix



#Packages to import
import requests
import h5py
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm

#Constants
ids = [110762, 69300, 56943, 141239, 463043]   #SubFind ID
#colours = ['r', 'b']     #Colours for plotting different Subhalo 3D plots 
snapshot = 69          #Snapshot interested in 
z = redshift(snapshot)   #Redshift associated with Snapshot
sqrta_velvar = np.sqrt(1 / (1 + z))  #sqrt(a) term
h = 0.6774               #h value - redshift factor
comoving_posvar = (1/h) * (1 / ( 1 + z )) #term to account for the expanding universe (h) and comoving coordinates
  
#params = {'gas':'Coordinates,Masses,StarFormationRate,Velocities'} #flitering for the required parameters 
paramsstar = {'stars':'Coordinates,Masses,Velocities'}
# #Constants for Mass Frac
# MassFracs = [0.1, 0.3, 0.5, 0.7, 0.9]

for id in ids:
    url = "http://www.tng-project.org/api/TNG100-1/snapshots/" + str(snapshot) + "/subhalos/" + str(id)
    sub = get(url) # get json response of subhalo properties
    saved_filename = get(url + "/cutout.hdf5",paramsstar) # get and save HDF5 cutout file
    print (id), print(saved_filename)
    
    print("For SubHalo id -> ", id)
    #extracting the half mass radius - UNIT CHECK : ckpc/h to kpc (with position variable)
    halfmassrad_gas = sub['halfmassrad_gas'] * comoving_posvar    
    print("half mass radius gas = ", halfmassrad_gas)
    halfmassrad_stars = sub['halfmassrad_stars'] * comoving_posvar 
    #extracting stellar mass of subhalo - UNIT CHECK : 10^10 Mo/h to Mo
    totMass_star = sub['mass_stars'] * 10**10/h
    print("Stellar Mass = ", totMass_star)
    print("Subhalo Center in ckpc/h ", sub['pos_x'], " ", sub['pos_y'], " ", sub['pos_z'])
    #www.tng-project.org/api/Illustris-1/snapshots/80/halos/523312/cutout.hdf5?dm=Coordinates,Velocities&gas=Coordinates,Masses
    with h5py.File(saved_filename, 'r') as f:
        # NOTE! If the subhalo is near the edge of the box, you must take the periodic boundary into account! (we ignore it here)
        # Unit Check : Position Units ckpc/h (corrected with position variable) to kpc
        # Velocity Units km sqrt(a)/s (corrected with velocity variable) to km/s  
        # Mass units 10^10 Mo/h corrected to Mo
        d_x = (f['PartType4']['Coordinates'][:,0] - sub['pos_x']) * comoving_posvar
        d_y = (f['PartType4']['Coordinates'][:,1] - sub['pos_y']) * comoving_posvar
        d_z = (f['PartType4']['Coordinates'][:,2] - sub['pos_z']) * comoving_posvar
        
        v_x = (f['PartType4']['Velocities'][:,0]) * sqrta_velvar
        v_y = (f['PartType4']['Velocities'][:,1]) * sqrta_velvar
        v_z = (f['PartType4']['Velocities'][:,2]) * sqrta_velvar
        
        stellar_mass = (f['PartType4']['Masses'][:]) * 10**10/h
 #       star_mass = (f['PartType4']['Masses'][:]) * 10**10/h
    
    #extracting data from HDF5 file to array
    arrx, arry, arrz = np.array(d_x), np.array(d_y), np.array(d_z)
    arrvx, arrvy, arrvz = np.array(v_x), np.array(v_y), np.array(v_z)
    arrmass = np.array(stellar_mass)
    
       #Loading to a txt file - open() function with "w" for creating a new file or writing over an existing one
    fname5 = "D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Stellar Mass Data\Subhalo_"+ str(id) +"_" + str(snapshot) +".txt"
    with open(fname5, "w") as g:
        np.savetxt(g, np.c_[arrx, arry, arrz, arrmass], fmt = '%.5f',  delimiter=',')


 #    #dist array
 #    r = hyp3D(arrx, arry, arrz)
    
 #    #creating new arrays for positions and velocities within the half mass radius
 #    dx = []
 #    dy = []
 #    dz = []
 #    vx = []
 #    vy = []
 #    vz = []
 #    m_star = []
    
 #    cutoff_vol = 4 #Cutoff for amount of particles seem on the board.
 #    for i in np.arange(len(arrx)):
 #        if ((np.abs(r[i]) <= cutoff_vol) ):  #& ( arrsfr[i] > 0.0)checking if the absolute distance from center less than/equal to halfmass radius
 #            dx.append(arrx[i])
 #            dy.append(arry[i])
 #            dz.append(arrz[i])
            
 #            vx.append(arrvx[i])
 #            vy.append(arrvy[i])
 #            vz.append(arrvz[i])
            
 #            m_star.append(arrmass[i])
 # #           m_star.append(star_mass[i])
    
 #    #Converting the above created lists into arrays
 #    arx = np.array(dx)
 #    ary = np.array(dy)
 #    arz = np.array(dz)
    
    # arvx = np.array([vx])
    # arvy = np.array([vy])
    # arvz = np.array([vz])
    
    # mass_startxt = np.array([m])
    # mass_star = np.array(m_star)
#    mass_star = np.array(m_star)
    
     # #Loading to a txt file - open() function with "w" for creating a new file or writing over an existing one
    # fname5 = "D:\Shivani\Summer 2022 Internship\TNG Illustrus API\TNG 100-1\Output\orbit txt files\Orbit_"+ str(orbit) +"_" + subtxt +".txt"API
    # with open(fname5, "w") as g:
    # np.savetxt(g, np.c_[np.asarray(snaps), np.asarray(subhaloIDs),np.asarray(MassFracss), np.asarray(lst_vx_weightedAvg), np.asarray(lst_vy_weightedAvg), np.asarray(lst_vz_weightedAvg), np.asarray(lst_totvel_weightedavg)], fmt = '%.5f',  delimiter=',')


    # #Plot 3D Image
    # plot3D(id, arx, ary, arz, sfr_starPart)
    
    # #Obtaining the rotating matrix
    # I = momentOfInertiaTensor(r, arrx, arry, arrz, arrmass, halfmassrad_stars)
    # r_faceon = rotationMatricesFromInertiaTensor(I)
    # print("Det |Rotating Matrix| = ", np.linalg.det(r_faceon))
      
    # #Transforming old values (3D) to new face-on projected values (2D)
    # x_2D = np.zeros(len(arx))
    # y_2D = np.zeros(len(arx))
    # z_2D = np.zeros(len(arx))
    # for point in np.arange(len(arx)):
    #     v = [arx[point], ary[point], arz[point]]
    #     vector_2D = np.dot(r_faceon, v)
    #     x_2D[point] = vector_2D[0,0]
    #     y_2D[point] = vector_2D[0,1]
    #     z_2D[point] = vector_2D[0,2]
    
    # #Converting particle image to spaxel data - Stellar Mass 
    # bins = 100
    # side = 2*cutoff_vol/bins
    # area = side**2
    # sums_starmass, xbins_starmass, ybins_starmass = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side), weights=mass_star)
    # counts_starmass, _, _ = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    # fig3 = plt.figure("Spaxels_stellarmass_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    #    img3 = plt.pcolormesh(xbins_starmass, ybins_starmass, np.log(sums_starmass.T/area) , cmap='viridis') #, vmin = 0.8*10**6, vmax = 2.2*10**6
    # fig3.colorbar(img3, label = "log of stellar mass")
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.xlim([-cutoff_vol, cutoff_vol])
    # plt.ylim([-cutoff_vol, cutoff_vol])
    # plt.title('Stellar Mass Surface Density (normalized) - Mo $kpc^-2$')
    
    
    #         #Converting particle image to spaxel data
    # sums_starmass, xbins_starmass, ybins_starmass = np.histogram2d(x_2D, y_2D, bins=bins, weights=mass_gas)
    # counts_starmass, _, _ = np.histogram2d(x_2D, y_2D, bins=bins)
    # fig4 = plt.figure("Spaxels_starmass_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    #    img4 = plt.pcolormesh(xbins_starmass, ybins_starmass, sums_starmass.T , cmap='viridis') #, vmin = 0.8*10**6, vmax = 2.2*10**6
    # fig4.colorbar(img4)
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.title('X Vs Y for Star Mass (normalized)')
    
    
    # #starforming & non-starforming
    # star_form = np.where(sfr_starPart > 0.0)
    # nonstar_form = np.where(sfr_starPart == 0.0)
    
    # fig2 = plt.figure("2D_" + str(id))    
    # plt.scatter(x_2D[nonstar_form], y_2D[nonstar_form], c=sfr_starPart[nonstar_form], cmap='viridis', s = 0.5, vmin = 0, vmax = 0.001)
    # pi = plt.scatter(x_2D[star_form], y_2D[star_form], c=sfr_starPart[star_form], cmap='viridis', s = 0.5, vmin = 0, vmax = 0.001)
    # fig2.colorbar(pi)
    # plt.title('SubHalo {0} Visualized in 2D'.format(id))
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')    
        
    plt.show()
    