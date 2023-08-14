# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 07:45:20 2023

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
from scipy.odr import ODR, Model, RealData

#Constants
ids = [110762, 69300, 56943, 141239, 463043]   #SubFind ID
#colours = ['r', 'b']     #Colours for plotting different Subhalo 3D plots 
snapshot = 69           #Snapshot interested in 
z = redshift(snapshot)   #Redshift associated with Snapshot
sqrta_velvar = np.sqrt(1 / (1 + z))  #sqrt(a) term
h = 0.6774               #h value - redshift factor
comoving_posvar = (1/h) * (1 / ( 1 + z )) #term to account for the expanding universe (h) and comoving coordinates
  
params = {'gas':'Coordinates,Masses,StarFormationRate,Velocities'} #flitering for the required parameters 
parmasstar = {'star':'Masses'}

#Lists for SFR & S Mass
lst_SFR = []
lst_SMass = []
lst_GMass = []
lst_GSMass =[]


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

    #loading subfind IDs for Snapshot 85
    directory = "D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Stellar Mass Data"
    filename = "Subhalo_"+ str(id) +"_" + str(snapshot) +".txt"
    orbitpath = os.path.join(directory, filename)

    starx, stary, starz, star_mass = np.loadtxt(orbitpath, delimiter = ',', unpack=True)
    
    #dist array - using gas coordinates 
    r = hyp3D(arrx, arry, arrz)
    rstar = hyp3D(starx, stary, starz) 
    cutoff_vol = 4 #Cutoff for amount of particles seem on the board.
    
    
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

    #Converting particle image to spaxel data - SFR MAPPING
    bins = 25
    side = 2*cutoff_vol/bins
    area = side**2
    sums, xbins, ybins = np.histogram2d(x_2D, y_2D, bins = np.arange(-cutoff_vol, cutoff_vol + side, side), weights=arrsfr)
    counts, _, _ = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    with np.errstate(divide='ignore', invalid='ignore'):
        SFR_mapping = np.log10(sums/area)

    #Converting particle image to spaxel data - GAS MASS MAPPING
    sums_gasmass, xbins_gasmass, ybins_gasmass = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side), weights=arrmass)
    counts_gasmass, _, _ = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    with np.errstate(divide='ignore', invalid='ignore'):
        G_Mass = np.log10(sums_gasmass/area) 
  
    #Converting particle image to spaxel data - Stellar Mass 
    sums_starmass, xbins_starmass, ybins_starmass = np.histogram2d(x_2Dstar, y_2Dstar, bins=np.arange(-cutoff_vol, cutoff_vol + side, side), weights=star_mass)
    counts_starmass, _, _ = np.histogram2d(x_2Dstar, y_2Dstar, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    with np.errstate(divide='ignore', invalid='ignore'):
        S_Mass = np.log10(sums_starmass/area)

    #     #Converting particle image to spaxel data - sSFR 
    # with np.errstate(divide='ignore', invalid='ignore'):
    #     sums_sSFR = sums/ sums_starmass
    # with np.errstate(divide='ignore', invalid='ignore'):
    #     sSFR_Mapping = np.log10(sums_sSFR/area)
  
    #  #Converting particle image to spaxel data - SFE MAPPING
    # with np.errstate(divide='ignore', invalid='ignore'):
    #     sums_SFE = sums/sums_gasmass
    # with np.errstate(divide='ignore', invalid='ignore'):
    #     SFE_Mapping = np.log10(sums_SFE)
 
    # #Converting particle image to spaxel data - fH 
    # with np.errstate(divide='ignore', invalid='ignore'):
    #     sums_fH = sums_gasmass/ sums_starmass
    # with np.errstate(divide='ignore', invalid='ignore'):
    #     fH_Mapping = np.log10(sums_fH)

    
    
   #Flatten the Array and filter for zero values - log SFR & log S Mass Relation
    Flat_SFR = SFR_mapping.flatten()
    Flat_SMass = S_Mass.flatten()
    
    Flat_SFR[np.isneginf(Flat_SFR)] = 0
    Flat_SMass[np.isneginf(Flat_SMass)] = 0
    
    Index_SFR = np.where(Flat_SFR != 0)
    Index_SMass = np.where(Flat_SMass != 0)
    Index_SFR_SMass = np.concatenate((Index_SFR[0],Index_SMass[0]))
    Index_SFR_SMass.sort(kind='mergesort')
    unq, count = np.unique(Index_SFR_SMass, return_counts=True)
    index_SFR_SMass = unq[count>1]
    
    flat_SFR =  Flat_SFR[index_SFR_SMass]
    flat_SMass =  Flat_SMass[index_SFR_SMass]
    print("SRF Size => ", np.size(Flat_SFR))
    print("Reduced SRF Size (removed zeros) => ", np.size(flat_SFR))
    
    lst_SFR.append(flat_SFR)
    lst_SMass.append(flat_SMass)   
    
    #Flatten the Array and filter for zero values - log G Mass & log S Mass Relation
    Flat_GMass = G_Mass.flatten()

    Flat_GMass[np.isneginf(Flat_GMass)] = 0
    
    Index_GMass = np.where(Flat_GMass != 0)
    
    Index_GMass_SMass = np.concatenate((Index_GMass[0],Index_SMass[0]))
    Index_GMass_SMass.sort(kind='mergesort')
    unq, count = np.unique(Index_GMass_SMass, return_counts=True)
    index_GMass_SMass = unq[count>1]
    
    flat_GMass =  Flat_GMass[index_GMass_SMass]
    flat_GSMass =  Flat_SMass[index_GMass_SMass]
    print("GMass Size => ", np.size(Flat_GMass))
    print("Reduced Gmass Size (removed zeros) => ", np.size(flat_GMass))
    
    lst_GMass.append(flat_GMass)
    lst_GSMass.append(flat_GSMass)   

#SFR & SMass Relation
Total_SFR = np.concatenate((lst_SFR[0:5]))
Total_SMass = np.concatenate((lst_SMass[0:5]))
#ODR Fitting 2q1
def ModelFit(beta, x):
    y_fit = beta[0] + beta[1]*x
    return y_fit
    
model = Model(ModelFit)
data = RealData(Total_SMass, Total_SFR)
odr = ODR(data, model, [7.5, 8])
odr.set_job(fit_type = 0)
output = odr.run()
output.pprint()
beta = output.beta
    
# #Plotting Surface density of SFR with Surface density of Stellar Mass, M*  
# f1 = plt.figure() 
line = "y = "+ str(round(beta[1],4))+ " x + "+ str(round(beta[0],4))
# plt.scatter(Total_SMass, Total_SFR)
# plt.plot(Total_SMass, ModelFit(beta, Total_SMass), color = 'red')
# plt.title('Scaling Relation b/w Log Surface Density of SFR & Stellar Mass')
# plt.xlabel('Log Stellar Mass (Mo/kpc^2)')
# plt.ylabel('Log SFR (Mo/yr*kpc^2)')
# plt.gca().legend((line,'datapoint'))

#GMass & SMass Relation
Total_GMass = np.concatenate((lst_GMass[0:5]))
Total_GSMass = np.concatenate((lst_GSMass[0:5]))
#ODR Fitting 2q1
    
Gdata = RealData(Total_GSMass, Total_GMass)
Godr = ODR(Gdata, model, [7.5, 8])
Godr.set_job(fit_type = 0)
Goutput = Godr.run()
Goutput.pprint()
Gbeta = Goutput.beta
    
# #Plotting Surface density of SFR with Surface density of Stellar Mass, M*  
# f2 = plt.figure()
Gline = "y = "+ str(round(Gbeta[1],4))+ " x + "+ str(round(Gbeta[0],4)) 
# plt.scatter(Total_GSMass, Total_GMass)
# plt.plot(Total_GSMass, ModelFit(Gbeta, Total_GSMass), color = 'red')
# plt.title('Scaling Relation b/w Log Surface Density of Gas Mass & Stellar Mass')
# plt.xlabel('Log Stellar Mass (Mo/kpc^2)')
# plt.ylabel('Log Gas Mass (Mo/kpc^2)')
# plt.gca().legend(( Gline,'datapoint'))

#SUBPLOTS
fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,10))
size = 20
lower_size = 10
title = "Scaling Relation"
fig.suptitle(title, fontsize=30)
#SFR----------------------------------
ax1.scatter(Total_SMass, Total_SFR)
ax1.plot(Total_SMass, ModelFit(beta, Total_SMass), color = 'red')
ax1.title.set_text(r"Log Surface Density of SFR vs $M_{\ast}$")    
ax1.set_ylabel(r"Log SFR $(\log$ $M_{\odot} kpc^{-2} yr^{-1})$")
ax1.legend(labels = ( line,'datapoint'), fontsize = 18)
# title = lg.get_label()
# title.set_fontsize('medium')
ax1.title.set_fontsize(size)
ax1.yaxis.label.set_fontsize(size)


    #starmass----------------------------------
ax2.scatter(Total_GSMass, Total_GMass)
ax2.plot(Total_GSMass, ModelFit(Gbeta, Total_GSMass), color = 'red')
ax2.title.set_text(r"Log Surface Density of $M_{gas}$ vs $M_{\ast}$")    
ax2.set_ylabel(r"Log $M_{gas}$ $(\log$ $M_{\odot} kpc^{-2})$")
ax2.legend(labels = ( Gline,'datapoint'), fontsize = 18)
ax2.title.set_fontsize(size)
ax2.yaxis.label.set_fontsize(size)

fig.text(0.5, 0.04, r"Log $M_{\ast}$ $(\log$ $M_{\odot} kpc^{-2})$", ha='center', va='center', fontsize=size)
plt.show()
    
    # X, Y = np.indices(SFR_mapping.shape)
    # mask = ~np.isinf(SFR_mapping)
    # x = X[mask]
    # y = Y[mask]
    # SFR_mapping = SFR_mapping[mask]
    
    # X, Y = np.indices(S_Mass.shape)
    # mask = ~np.isinf(S_Mass)
    # x = X[mask]
    # y = Y[mask]
    # S_Mass = S_Mass[mask]
    
  
    
    
    