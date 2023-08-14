# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 12:42:48 2022

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

"""
Obtains the hypotenuse or the magnitude of a 2D vector passed to it.

INPUT : X coordinate, Y coordinate
Returns : hyp -> hypotenuse/magnitude (Float)
"""
def hyp2D(x, y): 
    hyp = np.sqrt(x**2 + y**2)
    return hyp

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
# """
# Calculates radius for the given mass fraction of the total stellar mass.
# """
# def XmassRad (MassFrac, distArr, arrmass, totMass_star):
#     #Converting r array to list
#     rlist = distArr.tolist()
    
#     #initializing our variables before the while loop
#     totMass = 0
#     ReqMass = MassFrac * totMass_star
    
#     while (totMass <= ReqMass):
#         min_index = rlist.index(np.nanmin(rlist))     #finding the index of the min value in r
#         totMass += arrmass[min_index]     #adding the masses of the min r values 
#         if(totMass > ReqMass):
#             Xmass_rad = rlist[min_index]
#         #print('total mass', totMass)    
#         rlist[min_index] = np.nan         #replacing the previous min value with NaN
        
#     print(MassFrac, " total stellar mass rad with calc => ", Xmass_rad)
    
#     return Xmass_rad

#------------------------------------------------------------------------------------------------------
def momentOfInertiaTensor(distarr, x, y, z, arr_mass, sfr, pos, halfmassrad):
    """ Calculate the moment of inertia tensor (3x3 matrix) for a subhalo-scope particle set. """

    # load required particle data for this subhalo

    wGas = np.where( (distarr <= 2 * halfmassrad) & (sfr > 0.0))

    masses = arr_mass[wGas]
    x = x[wGas] - pos[0]
    y = y[wGas] - pos[1]
    z = z[wGas] - pos[2]
    
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
                             rotation_matrix[:,sort_inds[2]]))
    
    return new_matrix

#-------------------------------------------------------------------------------------------------------
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
        curr_min_index = rlist.index(np.nanmin(rlist))     #finding the index of the min value in r
        totMass += arrmass[curr_min_index]     #adding the masses of the min r values 

        if(totMass <= ReqMass):
            prev_min_value = rlist[curr_min_index]
                
        elif (totMass > ReqMass):
            Xmass_rad = prev_min_value
            calcMass  = totMass - arrmass[curr_min_index]            
        
        rlist[curr_min_index] = np.nan         #replacing the previous min value with NaN
    
    print(MassFrac, " x total stellar mass rad with calc => ", Xmass_rad, " kpc")
    print('total mass => ', calcMass)
    
    return Xmass_rad, calcMass


#================================================================================
#Packages to import
import requests
import h5py
import os
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm

#Constants
snapshot = 73           #Snapshot interested in 

if snapshot == 69: 
    ids = [356134, 356135]   #SubFind ID
elif snapshot == 71: 
    ids = [360389, 360390]   #SubFind ID
elif snapshot == 73: 
    ids = [364600, 364601]   #SubFind ID

#colours = ['r', 'b']     #Colours for plotting different Subhalo 3D plots 

z = redshift(snapshot)   #Redshift associated with Snapshot
sqrta_velvar = np.sqrt(1 / (1 + z))  #sqrt(a) term
h = 0.6774               #h value - redshift factor
comoving_posvar = (1/h) * (1 / ( 1 + z )) #term to account for the expanding universe (h) and comoving coordinates
  
params = {'gas':'Coordinates,Masses,StarFormationRate,Velocities'} #flitering for the required parameters 
parmasstar = {'star':'Masses'}

for id in ids:
    url = "http://www.tng-project.org/api/TNG100-1/snapshots/" + str(snapshot) + "/subhalos/" + str(id)
    sub = get(url) # get json response of subhalo properties
    saved_filename = get(url + "/cutout.hdf5",params) # get and save HDF5 cutout file
    print (id), print(saved_filename)
    
    print("For SubHalo id -> ", id)
    #extracting the half mass radius - UNIT CHECK : ckpc/h to kpc (with position variable)
    halfmassrad_gas = sub['halfmassrad_gas'] * comoving_posvar    
    print("half mass radius gas = ", halfmassrad_gas)
    halfmassrad_stars = sub['halfmassrad_stars'] * comoving_posvar 
    print("half mass radius stellar = ", halfmassrad_stars)
    #extracting stellar mass of subhalo - UNIT CHECK : 10^10 Mo/h to Mo
    totMass_star = sub['mass_stars'] * 10**10/h
    print("Stellar Mass = ", totMass_star)
    print("Subhalo Center in ckpc/h ", sub['pos_x'], " ", sub['pos_y'], " ", sub['pos_z'])
    pos = np.array((sub['pos_x']* comoving_posvar, sub['pos_y']* comoving_posvar, sub['pos_z']* comoving_posvar))
    #www.tng-project.org/api/Illustris-1/snapshots/80/halos/523312/cutout.hdf5?dm=Coordinates,Velocities&gas=Coordinates,Masses
    with h5py.File(saved_filename, 'r') as f:
        # NOTE! If the subhalo is near the edge of the box, you must take the periodic boundary into account! (we ignore it here)
        # Unit Check : Position Units ckpc/h (corrected with position variable) to kpc
        # Velocity Units km sqrt(a)/s (corrected with velocity variable) to km/s  
        # Mass units 10^10 Mo/h corrected to Mo
        d_x = (f['PartType0']['Coordinates'][:,0] - sub['pos_x']) * comoving_posvar
        d_y = (f['PartType0']['Coordinates'][:,1] - sub['pos_y']) * comoving_posvar
        d_z = (f['PartType0']['Coordinates'][:,2] - sub['pos_z']) * comoving_posvar
        
        v_x = (f['PartType0']['Velocities'][:,0]) * sqrta_velvar
        v_y = (f['PartType0']['Velocities'][:,1]) * sqrta_velvar
        v_z = (f['PartType0']['Velocities'][:,2]) * sqrta_velvar
        
        gas_mass = (f['PartType0']['Masses'][:]) * 10**10/h
        sfr_gas = (f['PartType0']['StarFormationRate'][:])
         #star_mass = (f['PartType4']['Masses'][:]) * 10**10/h
    
    #extracting data from HDF5 file to array
    arrx, arry, arrz = np.array(d_x), np.array(d_y), np.array(d_z)
    arrvx, arrvy, arrvz = np.array(v_x), np.array(v_y), np.array(v_z)
    arrmass = np.array(gas_mass)
    arrsfr = np.array(sfr_gas)
    #arrHI = np.array(HI_frac)
    

    #loading subfind IDs for Snapshot 85
    directory = "D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Stellar Mass Data"
    filename = "Subhalo_"+ str(id) +"_" + str(snapshot) +".txt"
    orbitpath = os.path.join(directory, filename)

    starx, stary, starz, star_mass = np.loadtxt(orbitpath, delimiter = ',', unpack=True)
    
    #dist array - using gas coordinates 
    r = hyp3D(arrx, arry, arrz)
    rstar = hyp3D(starx, stary, starz)
    
    #creating new arrays for positions and velocities within the half mass radius
    dx = []
    dy = []
    dz = []
    vx = []
    vy = []
    vz = []
    m_gas = []
    sfr = []
    
    dxstar = []
    dystar = []
    dzstar = []
    m_star = []
    
    cutoff_vol = 4 #Cutoff for amount of particles seem on the board.
    # for i in np.arange(len(arrx)):
    #     if ((np.abs(r[i]) <= cutoff_vol) ):  #& ( arrsfr[i] > 0.0)checking if the absolute distance from center less than/equal to halfmass radius
    #         dx.append(arrx[i])
    #         dy.append(arry[i])
    #         dz.append(arrz[i])
            
    #         vx.append(arrvx[i])
    #         vy.append(arrvy[i])
    #         vz.append(arrvz[i])
            
    #         sfr.append(arrsfr[i])
    #         m_gas.append(arrmass[i])
            
            
    #         #HI.append(arrHI[i])
            
    # for i in np.arange(len(starx)):
    #     if ((np.abs(rstar[i]) <= cutoff_vol) ):  
    #         dxstar.append(starx[i])
    #         dystar.append(stary[i])
    #         dzstar.append(starz[i])
    #         m_star.append(star_mass[i])
    
    # #Converting the above created lists into arrays
    # arx = np.array(dx)
    # ary = np.array(dy)
    # arz = np.array(dz)
    
    # # mass_startxt = np.array([m])
    # sfr_starPart = np.array(sfr)
    # mass_gas = np.array(m_gas)
    
    # xstar = np.array(dxstar)
    # ystar = np.array(dystar)
    # zstar = np.array(dzstar)
    # mass_star = np.array(m_star)
    
    #frac_HI = np.array(HI)
    
    
    
    ##FACE ON VIEW OF GALAXY ------------------------------------------------
    #It is based on constructing the moment of inertia tensor for all star-forming 
    #gas cells within twice the stellar half mass radius, if there are enough, 
    #otherwise the stars within the stellar half mass radius. Diagonalizing this 
    #tensor then gives a 3x3 rotation matrix (above) which can bring a flattened 
    #system into an edge-on or face-on type of view.
    
    ## GAS PARTICLES
    #Obtaining the rotating matrix
    I = momentOfInertiaTensor(r, arrx, arry, arrz, arrmass, arrsfr, pos, halfmassrad_stars)
    r_faceon = rotationMatricesFromInertiaTensor(I)
    print("Det |Rotating Matrix| = ", np.linalg.det(r_faceon))
    

    #Transforming old values (3D) to new face-on projected values (2D)    
    x_2D = np.zeros(len(arrx))
    y_2D = np.zeros(len(arrx))
    z_2D = np.zeros(len(arrx))
    for point in np.arange(len(arrx)):
        v = [arrx[point], arry[point], arrz[point]]
        vector_2D = np.dot(r_faceon, v)
        x_2D[point] = vector_2D[0,0]
        y_2D[point] = vector_2D[0,1]
        z_2D[point] = vector_2D[0,2]
    
    x_2Dstar = np.zeros(len(starx))
    y_2Dstar = np.zeros(len(stary))
    z_2Dstar = np.zeros(len(starz))
    for point in np.arange(len(starx)):
        vstar = [starx[point], stary[point], starz[point]]
        vector_2Dstar = np.dot(r_faceon, vstar)
        x_2Dstar[point] = vector_2Dstar[0,0]
        y_2Dstar[point] = vector_2Dstar[0,1]
        z_2Dstar[point] = vector_2Dstar[0,2]
        
    #  ## STELLAR PARTICLES
    # #Obtaining the rotating matrix
    # Istar = momentOfInertiaTensor(rstar, starx, stary, starz, star_mass, halfmassrad_stars)
    # r_faceonstar = rotationMatricesFromInertiaTensor(Istar)
    # print("Det |Rotating Matrix star| = ", np.linalg.det(r_faceonstar))
      
    # #Transforming old values (3D) to new face-on projected values (2D)
    # x_2Dstar = np.zeros(len(xstar))
    # y_2Dstar = np.zeros(len(ystar))
    # z_2Dstar = np.zeros(len(zstar))
    # for point in np.arange(len(xstar)):
    #     vstar = [xstar[point], ystar[point], zstar[point]]
    #     vector_2Dstar = np.dot(r_faceonstar, vstar)
    #     x_2Dstar[point] = vector_2Dstar[0,0]
    #     y_2Dstar[point] = vector_2Dstar[0,1]
    #     z_2Dstar[point] = vector_2Dstar[0,2]
    
    
    
    #Converting particle image to spaxel data - SFR MAPPING
    #Gyr = 10**9   #Unit Conversion to G
    bins = 25
    side = 2*cutoff_vol/bins     #ctuoff_vol is the variable containing required size of teh galaxy we want to visualize
    area = side**2
    sums, xbins, ybins = np.histogram2d(x_2D, y_2D, bins = np.arange(-cutoff_vol, cutoff_vol + side, side), weights=arrsfr)
    counts, _, _ = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    with np.errstate(divide='ignore', invalid='ignore'):
        bas_SFR = sums
        den_SFR = sums / area
        SFR_mapping = np.log10(den_SFR)
    # fig1 = plt.figure("Spaxels_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings - dividing by number of particles in the spaxel
    #    img1 = plt.pcolormesh(xbins, ybins, np.log10(sums.T/area), cmap='viridis') #
    # fig1.colorbar(img1, label = "log of SFR")
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.xlim([-cutoff_vol, cutoff_vol])
    # plt.ylim([-cutoff_vol, cutoff_vol])
    # plt.title('SFR Surface Density - Mo $kpc^(-2)$ $yr^(-1)$')

    #Converting particle image to spaxel data - GAS MASS MAPPING
    sums_gasmass, xbins_gasmass, ybins_gasmass = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side), weights=arrmass)
    counts_gasmass, _, _ = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    with np.errstate(divide='ignore', invalid='ignore'):
        G_Mass = np.log10(sums_gasmass/area) 
    # fig2 = plt.figure("Spaxels_gasmass_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    #     img2 = plt.pcolormesh(xbins_gasmass, ybins_gasmass, np.log10(sums_gasmass.T/area) , cmap='viridis') #, vmin = 0, vmax = max(np.log(mass_gas/area)) / 2
    # fig2.colorbar(img2, label = "log of Gas Mass")
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.xlim([-cutoff_vol, cutoff_vol])
    # plt.ylim([-cutoff_vol, cutoff_vol])
    # plt.title('Gas Mass Surface Density - Mo $kpc^-2$')
    
    
    #Converting particle image to spaxel data - Stellar Mass 
    sums_starmass, xbins_starmass, ybins_starmass = np.histogram2d(x_2Dstar, y_2Dstar, bins=np.arange(-cutoff_vol, cutoff_vol + side, side), weights=star_mass)
    counts_starmass, _, _ = np.histogram2d(x_2Dstar, y_2Dstar, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    with np.errstate(divide='ignore', invalid='ignore'):
        S_Mass = np.log10(sums_starmass/area)
    # fig4 = plt.figure("Spaxels_stellarmass_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    #     img4 = plt.pcolormesh(xbins_starmass, ybins_starmass, np.log10(sums_starmass.T/area) , cmap='viridis') #, vmin = 0.8*10**6, vmax = 2.2*10**6
    # fig4.colorbar(img4, label = "log of stellar mass")
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.xlim([-cutoff_vol, cutoff_vol])
    # plt.ylim([-cutoff_vol, cutoff_vol])
    # plt.title('Stellar Mass Surface Density (normalized) - Mo $kpc^-2$')
    
        #Converting particle image to spaxel data - sSFR 
    with np.errstate(divide='ignore', invalid='ignore'):
        #sSFR_Mapping = SFR_mapping / S_Mass        
        sums_sSFR = sums/ sums_starmass
    with np.errstate(divide='ignore', invalid='ignore'):
        sSFR_Mapping = np.log10(sums_sSFR)
        
    # sums_sSFR, xbins_sSFR, ybins_sSFR = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side), weights=sSFR)
    # counts_sSFR, _, _ = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    # fig5 = plt.figure("Spaxels_sSFR_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    #     img5 = plt.pcolormesh(xbins_SFE, ybins_SFE, np.log10(sums_sSFR/area) , cmap='viridis') #, vmin = 0.8*10**6, vmax = 2.2*10**6
    # fig5.colorbar(img5, label = "log of sSFR")
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.xlim([-cutoff_vol, cutoff_vol])
    # plt.ylim([-cutoff_vol, cutoff_vol])
    # plt.title('sSFR Surface Density (normalized) - Mo $kpc^-2$')
    
     #Converting particle image to spaxel data - SFE MAPPING
    with np.errstate(divide='ignore', invalid='ignore'):
        #SFE_Mapping =SFR_mapping/G_Mass        
        sums_SFE = sums/sums_gasmass
    with np.errstate(divide='ignore', invalid='ignore'):
        SFE_Mapping = np.log10(sums_SFE)
        
    #     SFE = sfr_starPart / mass_gas
    # sums_SFE, xbins_SFE, ybins_SFE = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side), weights=SFE)
    # counts_SFE, _, _ = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    # # fig3 = plt.figure("Spaxels_SFE_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    #     img3 = plt.pcolormesh(xbins_SFE, ybins_SFE, np.log10(sums_SFE.T), cmap='viridis') #, vmin = 0, vmax = max(np.log(SFE) / 2)
    # fig3.colorbar(img3, label = "log of SFE")
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.xlim([-cutoff_vol, cutoff_vol])
    # plt.ylim([-cutoff_vol, cutoff_vol])
    # plt.title('SFE Surface Density - $yr^-1$')
    
    #Converting particle image to spaxel data - fH 
    with np.errstate(divide='ignore', invalid='ignore'):
        #fH_Mapping = G_Mass / S_Mass      
        sums_fH = sums_gasmass/ sums_starmass
    with np.errstate(divide='ignore', invalid='ignore'):
        fH_Mapping = np.log10(sums_fH)
    
    # sums_fH, xbins_fH, ybins_fH = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side), weights=fH)
    # counts_fH, _, _ = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    # fig6 = plt.figure("Spaxels_fH_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    #     img6 = plt.pcolormesh(xbins_SFE, ybins_SFE, np.log10(sums_fH/area) , cmap='viridis') #, vmin = 0.8*10**6, vmax = 2.2*10**6
    # fig6.colorbar(img6, label = "log of fH")
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.xlim([-cutoff_vol, cutoff_vol])
    # plt.ylim([-cutoff_vol, cutoff_vol])
    # plt.title('fH Surface Density (normalized) - Mo $kpc^-2$')
    
    # #Calculating Half Mass Radius
    # MassFrac = 0.5
    # stellar_distarr = hyp2D(x_2Dstar, y_2Dstar)
    # calc_halfmassrad, calc_halfmass = XmassRad(MassFrac, stellar_distarr, star_mass, totMass_star)
    
    #SUBPLOTS
    fig, axs = plt.subplots(2,3)
    size = 20
    title = "Surface Density Plots - Subhalo_" + str(id) + "_" + str(snapshot)
    fig.suptitle(title, fontsize=size)
    plt.grid()

    #SFR----------------------------------
    
    img1 = axs[0,0].pcolormesh(xbins, ybins, SFR_mapping, cmap='viridis', vmin = -2.5, vmax = 1)
    axs[0,0].set_xlim([-cutoff_vol, cutoff_vol])
    axs[0,0].set_ylim([-cutoff_vol, cutoff_vol])
    fig.colorbar(img1, ax=axs[0,0], label = "log of SFR")
    axs[0,0].title.set_text("log SFR - Mo $kpc^{-2}$ $yr^{-1}$")    
    axs[0,0].title.set_fontsize(size)
    axs[0,0].yaxis.label.set_fontsize(size)
    axs[0,0].xaxis.label.set_fontsize(size)
    
    #starmass----------------------------------
    
    img2 = axs[1,1].pcolormesh(xbins_starmass, ybins_starmass, S_Mass, cmap='viridis', vmin = 6.5, vmax = 10)  #, vmin = 6.5, vmax = 10
    fig.colorbar(img2, ax=axs[1,1], label = "log of Stellar Mass")
    axs[1,1].set_xlim([-cutoff_vol, cutoff_vol])
    axs[1,1].set_ylim([-cutoff_vol, cutoff_vol])
    axs[1,1].title.set_text("log Stellar Mass - Mo $kpc^{-2}$")
    axs[1,1].title.set_fontsize(size)
    axs[1,1].yaxis.label.set_fontsize(size)
    axs[1,1].xaxis.label.set_fontsize(size)
    
    #gasmass----------------------------------
    
    img3 = axs[0,2].pcolormesh(xbins_gasmass, ybins_gasmass, G_Mass, cmap='viridis', vmin = 7, vmax = 9.5) #, vmin = 0, vmax = max(np.log(mass_gas/area)) / 2 , vmin = 6.5, vmax = 9.5
    fig.colorbar(img3, ax=axs[0,2], label = "log of Gas Mass")
    axs[0,2].set_xlim([-cutoff_vol, cutoff_vol])
    axs[0,2].set_ylim([-cutoff_vol, cutoff_vol])
    axs[0,2].title.set_text("log $M_{gas}$ - Mo $kpc^{-2}$")
    axs[0,2].title.set_fontsize(size)
    axs[0,2].yaxis.label.set_fontsize(size)
    axs[0,2].xaxis.label.set_fontsize(size)
    
    #sSFR----------------------------------
    
    img4 = axs[1,0].pcolormesh(xbins, ybins, sSFR_Mapping, cmap='viridis', vmin = -10, vmax = -6.5)
    fig.colorbar(img4, ax=axs[1,0], label = "log of sSFR")
    axs[1,0].set_xlim([-cutoff_vol, cutoff_vol])
    axs[1,0].set_ylim([-cutoff_vol, cutoff_vol])
    axs[1,0].title.set_text("log sSFR - $kpc^{-2}$ $yr^{-1}$")
    axs[1,0].title.set_fontsize(size)
    axs[1,0].yaxis.label.set_fontsize(size)
    axs[1,0].xaxis.label.set_fontsize(size)
    #SFE----------------------------------

    img5 = axs[0,1].pcolormesh(xbins, ybins, SFE_Mapping, cmap='viridis', vmin = -10.25, vmax = -8.25)
    fig.colorbar(img5, ax=axs[0,1], label = "log of SFE")
    axs[0,1].set_xlim([-cutoff_vol, cutoff_vol])
    axs[0,1].set_ylim([-cutoff_vol, cutoff_vol])
    axs[0,1].title.set_text("log SFE - $yr^{-1}$")
    axs[0,1].title.set_fontsize(size)
    axs[0,1].yaxis.label.set_fontsize(size)
    axs[0,1].xaxis.label.set_fontsize(size)
    #fH----------------------------------

    img5 = axs[1,2].pcolormesh(xbins, ybins, fH_Mapping, cmap='viridis', vmin = -2, vmax = 1.5)  #, vmin = -1.5, vmax = 1.5
    fig.colorbar(img5, ax=axs[1,2], label = "log of $f_{gas}$")
    axs[1,2].set_xlim([-cutoff_vol, cutoff_vol])
    axs[1,2].set_ylim([-cutoff_vol, cutoff_vol])
    axs[1,2].title.set_text("log $f_{gas}$ (Fraction)")
    axs[1,2].title.set_fontsize(size)
    axs[1,2].yaxis.label.set_fontsize(size)
    axs[1,2].xaxis.label.set_fontsize(size)
    
    # axs[1].plot(x, y2, linewidth=2, linestyle='-', color='g')
    # axs[1].plot(x, y2orig, linewidth=2, linestyle='--')           #plotting actual Velocity data points
    # axs[1].set_ylim([0,1])
    
    # axs[2].plot(x, y3, linewidth=2, linestyle='-', color='r')
    # axs[2].plot(x, y4, linewidth=2, linestyle='-', color='b')
    
    #Plot Lables
    # common axis labels
    
    fig.text(0.5, 0.04, 'x (kpc)', ha='center', va='center',fontsize=size)
    fig.text(0.06, 0.5, 'y (kpc)', ha='center', va='center', rotation='vertical',fontsize=size)
    
    plt.show()
    figure = plt.gcf()  # get current figure
    figure.set_size_inches(20, 10) # set figure's size manually to your full screen (32x18)   
    plt.savefig('D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Second Progress Report\TNG Data Maps\Subhalo_'+ str(id) +'_' + str(snapshot) +'.png')
    
    #Loading to a txt file - open() function with "w" for creating a new file or writing over an existing one
    fname5 = "D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Txt files data for Dave\SFR_Subhalo_"+ str(id) +"_" + str(snapshot) +".txt"
    fmt = '%.5f' #'%.2f', '%.2f','%.5f'
    with open(fname5, "w") as g:
        np.savetxt(g, SFR_mapping, fmt = fmt,  header='Log Surface Density SFR (log Mo kpc^-2 yr^-1)', delimiter=',')
    
    #Loading to a txt file - open() function with "w" for creating a new file or writing over an existing one
    fname5 = "D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Txt files data for Dave\Bins_Subhalo_"+ str(id) +"_" + str(snapshot) +".txt"
    fmt = '%.2f' #'%.2f', '%.2f','%.5f'
    with open(fname5, "w") as g:
        np.savetxt(g, np.c_[xbins, ybins], fmt = fmt,  header='xbins(kpc), ybins(kpc)', delimiter=',')

    
    # #COUNTS MAP - using SFR counts data
    # fig3 = plt.figure("Spaxels_COUNT" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings - dividing by number of particles in the spaxel
    #     img3 = plt.pcolormesh(xbins, ybins, counts, cmap='viridis') #
    # fig3.colorbar(img3, label = "count of particles")
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.xlim([-cutoff_vol, cutoff_vol])
    # plt.ylim([-cutoff_vol, cutoff_vol])
    # plt.title('Mapping Count of Particles')
    
    # #XY MAPPING WITHOUT ORIENTATION
    # #Converting particle image to spaxel data - SFR MAPPING
    # #Gyr = 10**9   #Unit Conversion to G
    # bins_xy = 50
    # side_xy = 2*cutoff_vol/bins_xy
    # area_xy = side**2
    # sums_SFRxy, xbins_xy, ybins_xy = np.histogram2d(arx, ary, bins = np.arange(-cutoff_vol, cutoff_vol + side_xy, side_xy), weights=sfr_starPart)
    # counts_xy, _, _ = np.histogram2d(arx, ary, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    # fig6 = plt.figure("Spaxels_SFRXY_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    #     img6 = plt.pcolormesh(xbins_xy, ybins_xy, np.log10(sums_SFRxy/area_xy) , cmap='viridis') #, vmin = 0.8*10**6, vmax = 2.2*10**6
    # fig6.colorbar(img6, label = "log of SFR")
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.xlim([-cutoff_vol, cutoff_vol])
    # plt.ylim([-cutoff_vol, cutoff_vol])
    # plt.title('SFR - Mo $yr^{-1}$ $kpc^{-2}$')
    
    # #Stellar Mass
    # sums_Starxy, xbins_Starxy, ybins_Starxy = np.histogram2d(xstar, ystar, bins = np.arange(-cutoff_vol, cutoff_vol + side_xy, side_xy), weights=mass_star)
    # counts_Starxy, _, _ = np.histogram2d(xstar, ystar, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    # fig7 = plt.figure("Spaxels_StarXY_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    #     img7 = plt.pcolormesh(xbins_Starxy, ybins_Starxy, np.log10(sums_Starxy/area_xy) , cmap='viridis') #, vmin = 0.8*10**6, vmax = 2.2*10**6
    # fig7.colorbar(img7, label = "log of Stellar Mass")
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.xlim([-cutoff_vol, cutoff_vol])
    # plt.ylim([-cutoff_vol, cutoff_vol])
    # plt.title('Stellar Mass - Mo $kpc^{-2}$')
    
    # #Gas MAss
    # sums_Gasxy, xbins_Gasxy, ybins_Gasxy = np.histogram2d(arx, ary, bins = np.arange(-cutoff_vol, cutoff_vol + side_xy, side_xy), weights=mass_gas)
    # counts_xy, _, _ = np.histogram2d(arx, ary, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    # fig8 = plt.figure("Spaxels_GasXY_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    #     img8 = plt.pcolormesh(xbins_Gasxy, ybins_Gasxy, np.log10(sums_Gasxy/area_xy) , cmap='viridis') #, vmin = 0.8*10**6, vmax = 2.2*10**6
    # fig8.colorbar(img8, label = "log of Gas Mass")
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.xlim([-cutoff_vol, cutoff_vol])
    # plt.ylim([-cutoff_vol, cutoff_vol])
    # plt.title('Gas Mass - Mo $kpc^{-2}$')
   
    
    
# # Hide x labels and tick labels for top plots and y ticks for right plots.
# for ax in axs.flat:
#     ax.label_outer()
    # #Converting particle image to spaxel data - HI MASS MAPPING
    # sums_HImass, xbins_HImass, ybins_HImass = np.histogram2d(x_2D, y_2D, bins=bins, weights=mass_gas * frac_HI)
    # counts_HImass, _, _ = np.histogram2d(x_2D, y_2D, bins=bins)
    # fig4 = plt.figure("Spaxels_HImass_" + str(id)) #fignum
    # with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    #    img4 = plt.pcolormesh(xbins_HImass, ybins_HImass, sums_HImass.T , cmap='viridis') #, vmin = 0.8*10**6, vmax = 2.2*10**6
    # fig3.colorbar(img4)
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)')
    # plt.title('X Vs Y for HI Mass')

    
    #starforming & non-starforming
    # star_form = np.where(sfr_starPart > 0.0)
    # nonstar_form = np.where(sfr_starPart == 0.0)
    
    # fig2 = plt.figure("2D_" + str(id))    
    # plt.scatter(x_2D[nonstar_form], y_2D[nonstar_form], c=sfr_starPart[nonstar_form], cmap='viridis', s = 0.5, vmin = 0, vmax = 0.001)
    # pi = plt.scatter(x_2D[star_form], y_2D[star_form], c=sfr_starPart[star_form], cmap='viridis', s = 0.5, vmin = 0, vmax = 0.001)
    # fig2.colorbar(pi)
    # plt.title('SubHalo {0} Visualized in 2D'.format(id))
    # plt.xlabel('x (kpc)')
    # plt.ylabel('y (kpc)') 
    
    
    # fig3 = plt.figure("x VS y" + str(id))    
    # plt.scatter(arx[nonstar_form], ary[nonstar_form], c=sfr_starPart[nonstar_form], cmap='viridis', s = 0.5, vmin = 0, vmax = 0.001)
    # pi = plt.scatter(arx[star_form], ary[star_form], c=sfr_starPart[star_form], cmap='viridis', s = 0.5, vmin = 0, vmax = 0.001)
    # fig3.colorbar(pi)    
    # plt.title('SubHalo {0} - X vs Y'.format(id))
    # plt.xlabel('x')
    # plt.ylabel('y') 