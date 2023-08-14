# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:38:34 2023

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

#--------------------------------------------------------------------------------------------

def CustomCmap(cmap1, cmap2, NoOfColors1, NoOfColors2, FinalName):  #, lowerbound_1, upperbound_1, lowerbound_2 , upperbound_2
    # define top and bottom colormaps 
    top = cm.get_cmap(cmap1, NoOfColors1) # r means reversed version
    bottom = cm.get_cmap(cmap2, NoOfColors2)
    # combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, NoOfColors1)),
                       bottom(np.linspace(0.1, 0.8, NoOfColors2))))
    # create a new colormaps with a name of OrangeBlue
    FinalCmap = ListedColormap(newcolors, name=FinalName)
    
    return FinalCmap

#================================================================================
#Packages to import
import requests
import h5py
import os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
    #Color Mapping
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

#Constants
#Constants
snapshot = 69       #Snapshot interested in 

if snapshot == 69: 
    ids = [356135,356134]   #SubFind ID 
elif snapshot == 71: 
    ids = [360390,360389]   #SubFind ID  
elif snapshot == 73: 
    ids = [364601,364600]   #SubFind ID  364601,364600

# ALL CONST. ASSOCIATED WITH UNIT CORRECTION
z = redshift(snapshot)   #Redshift associated with Snapshot
sqrta_velvar = np.sqrt(1 / (1 + z))  #sqrt(a) term
h = 0.6774               #h value - redshift factor
comoving_posvar = (1/h) * (1 / ( 1 + z )) #term to account for the expanding universe (h) and comoving coordinates
  
params = {'gas':'Coordinates,Masses,StarFormationRate,Velocities'} #flitering for the required parameters 

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
    
    r = hyp3D(arrx, arry, arrz)
    
    #loading subfind IDs for Snapshot 85
    directory = "D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Stellar Mass Data"
    filename = "Subhalo_"+ str(id) +"_" + str(snapshot) +".txt"
    orbitpath = os.path.join(directory, filename)

    starx, stary, starz, star_mass = np.loadtxt(orbitpath, delimiter = ',', unpack=True)
    
    #dist array - using gas coordinates 
    rstar = hyp3D(starx, stary, starz)
    
    ## GAS PARTICLES
    #Obtaining the rotating matrix
    I = momentOfInertiaTensor(r, arrx, arry, arrz, arrmass, arrsfr, pos, halfmassrad_stars)
    r_faceon = rotationMatricesFromInertiaTensor(I)
    print("Det |Rotating Matrix| = ", np.linalg.det(r_faceon))
    
    #Transforming old values (3D) to new face-on projected values (2D)
    #Gas Particles    
    x_2D = np.zeros(len(arrx))
    y_2D = np.zeros(len(arrx))
    z_2D = np.zeros(len(arrx))
    for point in np.arange(len(arrx)):
        v = [arrx[point], arry[point], arrz[point]]
        vector_2D = np.dot(r_faceon, v)
        x_2D[point] = vector_2D[0,0]
        y_2D[point] = vector_2D[0,1]
        z_2D[point] = vector_2D[0,2]
        

    #Star Particles
    x_2Dstar = np.zeros(len(starx))
    y_2Dstar = np.zeros(len(stary))
    z_2Dstar = np.zeros(len(starz))
    for point in np.arange(len(starx)):
        vstar = [starx[point], stary[point], starz[point]]
        vector_2Dstar = np.dot(r_faceon, vstar)
        x_2Dstar[point] = vector_2Dstar[0,0]
        y_2Dstar[point] = vector_2Dstar[0,1]
        z_2Dstar[point] = vector_2Dstar[0,2]
    
    # #Calculating Half Mass Radius
    # MassFrac = 0.5
    # stellar_distarr = hyp2D(x_2Dstar, y_2Dstar)
    # calc_halfmassrad, calc_halfmass = XmassRad(MassFrac, stellar_distarr, star_mass, totMass_star)
    
    
    #SPAXEL MAPPING BEGINS ---------------------------------x
    #Converting particle image to spaxel data 
    bins = 25
    cutoff_vol = 4
    side = 2*cutoff_vol/bins
    area = side**2
    
    #SFR Mapping 
    sums, xbins, ybins = np.histogram2d(x_2D, y_2D, bins = np.arange(-cutoff_vol, cutoff_vol + side, side), weights=arrsfr)
    counts, _, _ = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    with np.errstate(divide='ignore', invalid='ignore'):
        den_SFR = sums/area
        SFR_mapping = np.log10(den_SFR)
        
    #Gas Mass Map
    sums_gasmass, xbins_gasmass, ybins_gasmass = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side), weights=arrmass)
    counts_gasmass, _, _ = np.histogram2d(x_2D, y_2D, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))
    with np.errstate(divide='ignore', invalid='ignore'):
        den_GMass = sums_gasmass/area
        G_Mass = np.log10(den_GMass) 
    
    #Converting particle image to spaxel data - Stellar Mass 
    sums_starmass, xbins_starmass, ybins_starmass = np.histogram2d(x_2Dstar, y_2Dstar, bins=np.arange(-cutoff_vol, cutoff_vol + side, side), weights=star_mass)
    counts_starmass, _, _ = np.histogram2d(x_2Dstar, y_2Dstar, bins=np.arange(-cutoff_vol, cutoff_vol + side, side))  
    #Calculating stellar mass
    with np.errstate(divide='ignore', invalid='ignore'):
        den_SMass = sums_starmass/area
        S_Mass = np.log10(den_SMass)
    
   #sSFR Map    
    with np.errstate(divide='ignore', invalid='ignore'):
        den_sSFR = den_SFR / den_SMass
        sSFR_Mapping = np.log10(den_sSFR)
    
    #SFE MAPPING       
    with np.errstate(divide='ignore', invalid='ignore'):
        #Setting SFE zero when Mgas is zero
        den_SFE = np.divide(den_SFR, den_GMass, out=np.zeros_like(den_SFR), where=den_GMass!=0)
        SFE_Mapping = np.log10(den_SFE)
    
    #fgas Mapping   
    with np.errstate(divide='ignore', invalid='ignore'):
        den_fH = den_GMass / den_SMass
        fH_Mapping = np.log10(den_fH)
    
    #Checking if SFR = 0 when Mgas = 0, by flattening array for easier looping
    Flat_SFR = den_SFR.flatten()       
    Flat_GMass = den_GMass.flatten()
    for i in np.arange(len(Flat_SFR)):
        if Flat_GMass[i] == 0:
            if Flat_SFR[i] != 0:
                print("SFR NOT ZERO!!")
    
    
    #FINDING BIN CENTERS FOR ALL BINS-(FOR FINDING CENTRAL HMR CONTRIBUTIONS)----------------------xxx
    center = []
    hypcenter = []
    
    #loop through ybins from the top and find the center
    for i in np.arange(len(ybins)-1, 0, -1):
        ycen = (ybins[i] + ybins[i-1]) / 2
        
        #loop through xbin valuse from left side and find center
        for j in np.arange(len(xbins)-1):
            xcen = (xbins[j] + xbins[j+1]) / 2
            center.append((xcen,ycen))
            #append the hypotenuse of the bin center 
            hyp = hyp2D(xcen, ycen)
            hypcenter.append(hyp)
             
    hypcenter = np.array(hypcenter)  #Centers flattened list
    
    #SCALING RELATIONS (Values found from another code) -------------------------------------xxx
    a_s = 0.9296    #SFR slope value
    b_s = -9.2065   #SFR y-intercept value
    SFR_den = a_s *  S_Mass + b_s    # units -> per yr, not Gyr
    
    a_g = 0.5228
    b_g = 3.5437
    GasMass_den = a_g *  S_Mass + b_g
    
    #Derived Log Plots 
    with np.errstate(all='ignore'):
        sSFR = SFR_den - S_Mass
        SFE = SFR_den - GasMass_den
        fGas = GasMass_den - S_Mass
    
    #for linear Plots
    exp_SFR = 10**(SFR_den)
    exp_GMass = 10**(GasMass_den)
    with np.errstate(divide='ignore', invalid='ignore'):
        exp_sSFR = exp_SFR / den_SMass
        exp_SFE = np.divide(exp_SFR, exp_GMass, out=np.zeros_like(exp_SFR), where=exp_GMass!=0)
        exp_fH = exp_GMass / den_SMass
    
    # #DATA FOR DAVE in txt file -----------------X
    #     #Loading to a txt file - open() function with "w" for creating a new file or writing over an existing one
    # fname5 = "D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Txt files data for Dave\SFR_ScalingRelations_Subhalo_"+ str(id) +"_" + str(snapshot) +".txt"
    # fmt = '%.5f' #'%.2f', '%.2f','%.5f'
    # with open(fname5, "w") as g:
    #     np.savetxt(g, SFR_den, fmt = fmt,  header='Log Surface Density SFR (log Mo kpc^-2 yr^-1) - Scaling Relation', delimiter=',')
    
    #Calculating Log Difference 
    with np.errstate(all='ignore'):
        Diff_SFR = SFR_mapping - SFR_den
        Diff_Gmass = G_Mass - GasMass_den
        Diff_sSFR = sSFR_Mapping - sSFR
        Diff_SFE = SFE_Mapping - SFE
        Diff_fGas = fH_Mapping - fGas
    
    #Calculating Linear difference
    diffSFR = den_SFR - exp_SFR
    diffGmass = den_GMass - exp_GMass
    diffsSFR = den_sSFR - exp_sSFR
    diffSFE = den_SFE - exp_SFE
    diffGas = den_fH - exp_fH
    
    #Calculating Contributions- Fgas & SFE for Spatial Map
    sSFRg = exp_SFE * diffGmass + 0.5 * diffGmass * diffSFE
    sSFRe = exp_GMass * diffSFE + 0.5 * diffGmass * diffSFE
    
    with np.errstate(divide='ignore', invalid='ignore'):
        log_sSFRg = Diff_fGas / Diff_sSFR
        log_sSFRe = Diff_SFE / Diff_sSFR
    
    
        #sSFRe contribution
    roundnum = 4
    Flat_sSFRe = sSFRe.flatten()  
    addsSFRe = np.sum(Flat_sSFRe)
        #sSFRe contribution
    Flat_sSFRg = sSFRg.flatten()  
    addsSFRg = np.sum(Flat_sSFRg)
        #sSFRe contribution
    Flat_diffSFR = diffSFR.flatten() 
    adddiffSFR = np.sum(Flat_diffSFR)

    #Print Statements
    print("NEGATIVE VALUES SFE => ", (Flat_sSFRe<0).sum())
    print("NEGATIVE VALUES GAS => ", (Flat_sSFRg<0).sum())
    print("ADD sSFRe ==> ",round(addsSFRe, roundnum))
    print("ADD sSFRg ==> ",round(addsSFRg, roundnum))
    print("ADD SFR ==> ",round(adddiffSFR, roundnum))
    
    if round(addsSFRe,roundnum) + round(addsSFRg,roundnum) == round(adddiffSFR,roundnum):
        print("TRUE")
    
    # #Loading to a txt file - open() function with "w" for creating a new file or writing over an existing one
    # fname5 = "D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Contributions Percent_Subhalo_"+ str(id) +"_" + str(snapshot) +".txt"
    # fmt = '%.5f' #'%.2f', '%.2f','%.5f'
    # with open(fname5, "w") as g:
    #     np.savetxt(g, SFR_den, fmt = fmt,  header='Log Surface Density SFR (log Mo kpc^-2 yr^-1) - Scaling Relation', delimiter=',')
    
    #CALCULATING CENTER (HMR) CONTRIBUTION ---------------------------------------------xx
    centralrad = halfmassrad_stars
    
    SFRe_central = []
    SFRg_central = []
    SFR_central = []
    for i in np.arange(len(hypcenter)):
        if (hypcenter[i] <= centralrad):
            SFRe_central.append(Flat_sSFRe[i])
            SFRg_central.append(Flat_sSFRg[i])
            SFR_central.append(Flat_diffSFR[i])
    
    #Central contribution
    SFRe_central = np.array(SFRe_central)
    addSFRe_central = np.sum(SFRe_central)
    
    SFRg_central = np.array(SFRg_central)
    addSFRg_central = np.sum(SFRg_central)
    
    SFR_central = np.array(SFR_central)
    addSFR_central = np.sum(SFR_central)
    
    print("SFRe Central ==> ", addSFRe_central, "  Percent from whole => ", addSFRe_central/addsSFRe * 100)
    print("SFRg Central ==> ", addSFRg_central, "  Percent from whole => ", addSFRg_central/addsSFRg * 100)
    print("SFR Central ==> ", addSFR_central, "  Percent from whole => ", addSFR_central/adddiffSFR * 100)
    
    #CALCULATING Outskirts CONTRIBUTION ---------------------------------------------xx
    SFRe_outskirt = []
    SFRg_outskirt = []
    SFR_outskirt = []
    for i in np.arange(len(hypcenter)):
        if (hypcenter[i] > centralrad):
            SFRe_outskirt.append(Flat_sSFRe[i])
            SFRg_outskirt.append(Flat_sSFRg[i])
            SFR_outskirt.append(Flat_diffSFR[i])
    
    #Central contribution
    SFRe_outskirt = np.array(SFRe_outskirt)
    addSFRe_outskirt = np.sum(SFRe_outskirt)
    
    SFRg_outskirt = np.array(SFRg_outskirt)
    addSFRg_outskirt = np.sum(SFRg_outskirt)
    
    SFR_outskirt = np.array(SFR_outskirt)
    addSFR_outskirt = np.sum(SFR_outskirt)
    
    print("SFRe outskirt ==> ", addSFRe_outskirt, "  Percent from whole => ", addSFRe_outskirt/addsSFRe * 100)
    print("SFRg outskirt ==> ", addSFRg_outskirt, "  Percent from whole => ", addSFRg_outskirt/addsSFRg * 100)
    print("SFR outskirt ==> ", addSFR_outskirt, "  Percent from whole => ", addSFR_outskirt/adddiffSFR * 100)
    
    # #Mimicing integrated value for the galaxy properties-----------------------x  
    # #creating our integrated values by adding up all spaxels in SFR, Gass mass and Stellar mass
    # Flat_SFR = den_SFR.flatten()
    # Flat_SMass = den_SMass.flatten()
    # Flat_GMass = den_GMass.flatten()
    # Flat_expSFR = exp_SFR.flatten()
    # Flat_expGMass = exp_GMass.flatten()
    
    # add_SFR = np.sum(Flat_SFR)
    # add_SMass = np.sum(Flat_SMass)
    # add_GMass = np.sum(Flat_GMass)
    # add_expSFR = np.sum(Flat_expSFR)
    # add_expGMass = np.sum(Flat_expGMass)
    
    # #Further calculations to find contributions
    # add_SFE = add_SFR/add_GMass
    # add_expSFE = add_expSFR/add_expGMass
    
    # add_diffSFR = add_SFR - add_expSFR
    # add_diffSFE = add_SFE - add_expSFE
    # add_diffGMass = add_GMass - add_expGMass
    
    # add_sSFRg = add_expSFE * add_diffGMass + 0.5 * add_diffGMass * add_diffSFE
    # add_sSFRe = add_expGMass * add_diffSFE + 0.5 * add_diffGMass * add_diffSFE
    
    #         #sSFRe contribution
    # roundnum = 4
    # print("Integrated sSFRe ==> ",round(add_sSFRe, roundnum))
    #     #sSFRg contribution
    # print("Integrated sSFRg ==> ",round(add_sSFRg, roundnum))
    #     #SFR contribution
    # print("Integrated SFR ==> ",round(add_diffSFR, roundnum))
    
    # if round(add_sSFRe,roundnum) + round(add_sSFRg,roundnum) == round(add_diffSFR,roundnum):
    #     print("TRUE")


    # #PLOTTING===========================================================================================
    # #SUBPLOTS - Scaling Relations
    # fig, axs = plt.subplots(2,3)
    size = 20
    # title = "Scaling Relations - Subhalo_" + str(id) + "_" + str(snapshot)
    # fig.suptitle(title, fontsize=size)
    # plt.grid()
    
    # #SFR----------------------------------
    # img1 = axs[0,0].pcolormesh(xbins_starmass, ybins_starmass, SFR_den, cmap='viridis', vmin = -3, vmax = 0) #, vmin = 6.5, vmax = 10
    # axs[0,0].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[0,0].set_ylim([-cutoff_vol, cutoff_vol])
    # fig.colorbar(img1, ax=axs[0,0])
    # axs[0,0].title.set_text(r"$\Sigma_{SFR}$ $(\log$ $M_{\odot} kpc^{-2} yr^{-1})$")    
    # axs[0,0].title.set_fontsize(size)
    # axs[0,0].yaxis.label.set_fontsize(size)
    # axs[0,0].xaxis.label.set_fontsize(size)
    
    # #starmass----------------------------------
    # img2 = axs[1,1].pcolormesh(xbins_starmass, ybins_starmass, S_Mass, cmap='viridis', vmin = 6.5, vmax = 9.75) #, vmin = 6.5, vmax = 10
    # fig.colorbar(img2, ax=axs[1,1])
    # axs[1,1].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[1,1].set_ylim([-cutoff_vol, cutoff_vol])
    # axs[1,1].title.set_text(r"$\Sigma_{\ast}$ $(\log$ $M_{\odot} kpc^{-2})$")
    # axs[1,1].title.set_fontsize(size)
    # axs[1,1].yaxis.label.set_fontsize(size)
    # axs[1,1].xaxis.label.set_fontsize(size)
    
    # #gasmass----------------------------------
    # img3 = axs[0,2].pcolormesh(xbins_starmass, ybins_starmass, GasMass_den, cmap='viridis', vmin = 7, vmax = 8.6) #, vmin = 6.5, vmax = 9.5
    # fig.colorbar(img3, ax=axs[0,2])
    # axs[0,2].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[0,2].set_ylim([-cutoff_vol, cutoff_vol])
    # axs[0,2].title.set_text(r"$\Sigma_{gas}$ $(\log$ $M_{\odot} kpc^{-2})$")
    # axs[0,2].title.set_fontsize(size)
    # axs[0,2].yaxis.label.set_fontsize(size)
    # axs[0,2].xaxis.label.set_fontsize(size)
    
    # #sSFR----------------------------------
    # img4 = axs[1,0].pcolormesh(xbins_starmass, ybins_starmass, sSFR, cmap='viridis', vmin = -9.9, vmax = -9.7) #, vmin = -1.5, vmax = 2.5
    # fig.colorbar(img4, ax=axs[1,0])
    # axs[1,0].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[1,0].set_ylim([-cutoff_vol, cutoff_vol])
    # axs[1,0].title.set_text(r"sSFR $(\log$ $yr^{-1})$")
    # axs[1,0].title.set_fontsize(size)
    # axs[1,0].yaxis.label.set_fontsize(size)
    # axs[1,0].xaxis.label.set_fontsize(size)
    
    # #SFE----------------------------------
    # img5 = axs[0,1].pcolormesh(xbins_starmass, ybins_starmass, SFE, cmap='viridis', vmin = -10, vmax = -8.8) #, vmin = -1.0, vmax = 0.6
    # fig.colorbar(img5, ax=axs[0,1])
    # axs[0,1].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[0,1].set_ylim([-cutoff_vol, cutoff_vol])
    # axs[0,1].title.set_text(r"SFE $(\log$ $yr^{-1})$")
    # axs[0,1].title.set_fontsize(size)
    # axs[0,1].yaxis.label.set_fontsize(size)
    # axs[0,1].xaxis.label.set_fontsize(size)
    
    # #fH----------------------------------
    # img6 = axs[1,2].pcolormesh(xbins_starmass, ybins_starmass, fGas, cmap='viridis', vmin = -1, vmax = 0.4) #, vmin = -1.5, vmax = 1.5
    # fig.colorbar(img6, ax=axs[1,2])
    # axs[1,2].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[1,2].set_ylim([-cutoff_vol, cutoff_vol])
    # axs[1,2].title.set_text(r"$f_{gas}$ (log Fraction)")
    # axs[1,2].title.set_fontsize(size)
    # axs[1,2].yaxis.label.set_fontsize(size)
    # axs[1,2].xaxis.label.set_fontsize(size)
    
    # # x & y labels 
    # fig.text(0.5, 0.04, 'x (kpc)', ha='center', va='center', fontsize=size)
    # fig.text(0.06, 0.5, 'y (kpc)', ha='center', va='center', rotation='vertical', fontsize=size)
    
    # plt.show()
    # figure = plt.gcf()  # get current figure
    # figure.set_size_inches(20, 10) # set figure's size manually to your full screen (32x18)   
    # plt.savefig('D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Second Progress Report\Scaling Relation Maps\ScalingSubhalo_'+ str(id) +'_' + str(snapshot) +'.png')
    
    
    
# #===========================================================================================================    
#     #SUBPLOTS - Difference Graphs
#     fig1, axs = plt.subplots(2,3)
#     title = "Difference Plots : log(a) - log(b) :- Subhalo_" + str(id) + "_" + str(snapshot)
#     fig1.suptitle(title, fontsize=size)
#     plt.grid()
    
#     #SFR----------------------------------
#     img11 = axs[0,0].pcolormesh(xbins_starmass, ybins_starmass, Diff_SFR, cmap='viridis', vmin = -1, vmax = 2) #, vmin = 6.5, vmax = 10
#     fig1.colorbar(img11, ax=axs[0,0])
#     axs[0,0].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[0,0].set_ylim([-cutoff_vol, cutoff_vol])
#     axs[0,0].title.set_text(r"$\Delta \Sigma_{SFR}$ $(\log$ $M_{\odot} kpc^{-2} yr^{-1})$")  
#     axs[0,0].title.set_fontsize(size)
#     axs[0,0].yaxis.label.set_fontsize(size)
#     axs[0,0].xaxis.label.set_fontsize(size)
    
#     #starmass----------------------------------
#     img12 = axs[0,1].pcolormesh(xbins_starmass, ybins_starmass, S_Mass, cmap='viridis', vmin = 6.5, vmax = 9.75) #, vmin = 6.5, vmax = 10
#     fig1.colorbar(img12, ax=axs[0,1])
#     axs[0,1].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[0,1].set_ylim([-cutoff_vol, cutoff_vol])
#     axs[0,1].title.set_text(r"$\Delta \Sigma_{\ast}$ $(\log$ $M_{\odot} kpc^{-2})$")
#     axs[0,1].title.set_fontsize(size)
#     axs[0,1].yaxis.label.set_fontsize(size)
#     axs[0,1].xaxis.label.set_fontsize(size)
    
#     #gasmass----------------------------------
#     img13 = axs[0,2].pcolormesh(xbins_starmass, ybins_starmass, Diff_Gmass, cmap='viridis', vmin = -1, vmax = 1.25) #, vmin = 6.5, vmax = 9.5
#     fig1.colorbar(img13, ax=axs[0,2])
#     axs[0,2].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[0,2].set_ylim([-cutoff_vol, cutoff_vol])
#     axs[0,2].title.set_text(r"$\Delta \Sigma_{gas}$ $(\log$ $M_{\odot} kpc^{-2})$")
#     axs[0,2].title.set_fontsize(size)
#     axs[0,2].yaxis.label.set_fontsize(size)
#     axs[0,2].xaxis.label.set_fontsize(size)
    
#     #sSFR----------------------------------
#     img14 = axs[1,0].pcolormesh(xbins_starmass, ybins_starmass, Diff_sSFR, cmap='viridis', vmin = -1, vmax = 2) #, vmin = -1.5, vmax = 2.5
#     fig1.colorbar(img14, ax=axs[1,0])
#     axs[1,0].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[1,0].set_ylim([-cutoff_vol, cutoff_vol])
#     axs[1,0].title.set_text(r"$\Delta sSFR$ $(\log$ $yr^{-1})$")
#     axs[1,0].title.set_fontsize(size)
#     axs[1,0].yaxis.label.set_fontsize(size)
#     axs[1,0].xaxis.label.set_fontsize(size)
    
#     #SFE----------------------------------
#     img15 = axs[1,1].pcolormesh(xbins_starmass, ybins_starmass, Diff_SFE, cmap='viridis', vmin = -1, vmax = 0.8) #, vmin = -1.0, vmax = 0.6
#     fig1.colorbar(img15, ax=axs[1,1])
#     axs[1,1].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[1,1].set_ylim([-cutoff_vol, cutoff_vol])
#     axs[1,1].title.set_text(r"$\Delta SFE$ $(\log$ $yr^{-1})$")
#     axs[1,1].title.set_fontsize(size)
#     axs[1,1].yaxis.label.set_fontsize(size)
#     axs[1,1].xaxis.label.set_fontsize(size)
    
#     #fH----------------------------------
#     img16 = axs[1,2].pcolormesh(xbins_starmass, ybins_starmass, Diff_fGas, cmap='viridis', vmin = -1, vmax = 1.25) #, vmin = -1.5, vmax = 1.5
#     fig1.colorbar(img16, ax=axs[1,2])
#     axs[1,2].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[1,2].set_ylim([-cutoff_vol, cutoff_vol])
#     axs[1,2].title.set_text(r"$\Delta f_{gas}$ (log Fraction)")
#     axs[1,2].title.set_fontsize(size)
#     axs[1,2].yaxis.label.set_fontsize(size)
#     axs[1,2].xaxis.label.set_fontsize(size)
    
#     # x & y labels 
#     fig1.text(0.5, 0.04, 'x (kpc)', ha='center', va='center',fontsize=size)
#     fig1.text(0.06, 0.5, 'y (kpc)', ha='center', va='center', rotation='vertical', fontsize=size)
    
#     #plt.show()
#     figure2 = plt.gcf()  # get current figure
#     figure2.set_size_inches(20, 10) # set figure's size manually to your full screen (32x18)   
#     plt.savefig('D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Second Progress Report\Difference Maps\DiffSubhalo_'+ str(id) +'_' + str(snapshot) +'.png')
    
    
#     #=======================+++++++++++++++++++++++++++++++++++========================================
    
#     #TEST - LATE LOG------------------------------------------------------------------------------
#     fig3, axs = plt.subplots(2,3)
#     title = "Linear Difference Plots (a - b) :- Subhalo_" + str(id) + "_" + str(snapshot)
#     fig3.suptitle(title, fontsize=size)
#     plt.grid()
    
#     #SFR----------------------------------
#     imglatelog1 = axs[0,0].pcolormesh(xbins_starmass, ybins_starmass, diffSFR, cmap='viridis') #, vmin = -4, vmax = 0.5
#     axs[0,0].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[0,0].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog1, ax=axs[0,0])
#     axs[0,0].title.set_text(r"$\delta \Sigma_{SFR}$ $(M_{\odot} kpc^{-2} yr^{-1})$")    
#     axs[0,0].title.set_fontsize(size)
#     axs[0,0].yaxis.label.set_fontsize(size)
#     axs[0,0].xaxis.label.set_fontsize(size)
    
#     #starmass----------------------------------
#     imglatelog2 = axs[0,1].pcolormesh(xbins_starmass, ybins_starmass, den_SMass, cmap='viridis') #, vmin = 6.5, vmax = 9.5
#     axs[0,1].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[0,1].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog2, ax=axs[0,1])
#     axs[0,1].title.set_text(r"$ \Sigma_{\ast}$ $(M_{\odot} kpc^{-2})$")    
#     axs[0,1].title.set_fontsize(size)
#     axs[0,1].yaxis.label.set_fontsize(size)
#     axs[0,1].xaxis.label.set_fontsize(size)
    
#     #Gas Mass----------------------------------
#     imglatelog3 = axs[0,2].pcolormesh(xbins_starmass, ybins_starmass, diffGmass, cmap='viridis') #, vmin = 5, vmax = 9
#     axs[0,2].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[0,2].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog3, ax=axs[0,2])
#     axs[0,2].title.set_text(r"$\delta \Sigma_{gas}$ $(M_{\odot} kpc^{-2})$")    
#     axs[0,2].title.set_fontsize(size)
#     axs[0,2].yaxis.label.set_fontsize(size)
#     axs[0,2].xaxis.label.set_fontsize(size)
    
#     #sSFR----------------------------------
#     imglatelog4 = axs[1,0].pcolormesh(xbins_starmass, ybins_starmass, diffsSFR, cmap='viridis') #, vmin = -13, vmax = -8
#     axs[1,0].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[1,0].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog4, ax=axs[1,0])
#     axs[1,0].title.set_text(r"$\delta sSFR$ $(yr^{-1})$")    
#     axs[1,0].title.set_fontsize(size)
#     axs[1,0].yaxis.label.set_fontsize(size)
#     axs[1,0].xaxis.label.set_fontsize(size)
    
#     #SFE---------------------------------
#     imglatelog5 = axs[1,1].pcolormesh(xbins_starmass, ybins_starmass, diffSFE, cmap='viridis') #, vmin = -13, vmax = -8.5
#     axs[1,1].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[1,1].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog5, ax=axs[1,1])
#     axs[1,1].title.set_text(r"$\delta SFE$ $(yr^{-1})$")    
#     axs[1,1].title.set_fontsize(size)
#     axs[1,1].yaxis.label.set_fontsize(size)
#     axs[1,1].xaxis.label.set_fontsize(size)
    
#     #fgas----------------------------------
#     imglatelog6 = axs[1,2].pcolormesh(xbins_starmass, ybins_starmass, diffGas, cmap='viridis') #, vmin = -4, vmax = 1.5
#     axs[1,2].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[1,2].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog6, ax=axs[1,2])
#     axs[1,2].title.set_text(r"$\delta f_{gas}$ (Fraction)")    
#     axs[1,2].title.set_fontsize(size)
#     axs[1,2].yaxis.label.set_fontsize(size)
#     axs[1,2].xaxis.label.set_fontsize(size)
#     #-------------------------------------------------------------------------------------------------
    
#     # x & y labels 
#     fig3.text(0.5, 0.04, 'x (kpc)', ha='center', va='center',fontsize=size)
#     fig3.text(0.06, 0.5, 'y (kpc)', ha='center', va='center', rotation='vertical', fontsize=size)
    
#     #plt.show()
#     figure3 = plt.gcf()  # get current figure
#     figure3.set_size_inches(20, 10) # set figure's size manually to your full screen (32x18)   
#     plt.savefig('D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Second Progress Report\Linear Difference Maps\LinearDiffSubhalo_'+ str(id) +'_' + str(snapshot) +'.png')
    
#     #=======================+++++++++++++++++++++++++++++++++++========================================
    
#     #ALL in ONE SFR Map ------------------------------------------------------------------------------
#     fig3, axs = plt.subplots(2,3)
#     title = "ALL in ONE SFR Plots- Subhalo_" + str(id) + "_" + str(snapshot)
#     fig3.suptitle(title, fontsize=size)
#     plt.grid()
    
#     #SFR TNG Data----------------------------------
#     imglatelog1 = axs[0,0].pcolormesh(xbins_starmass, ybins_starmass, SFR_mapping, cmap='viridis', vmin = -3, vmax = 1) #, vmin = 6.5, vmax = 10   , vmin = -2.5, vmax = 0.5
#     axs[0,0].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[0,0].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog1, ax=axs[0,0])
#     axs[0,0].title.set_text("log SFR TNG Data")    
#     axs[0,0].title.set_fontsize(size)
#     axs[0,0].yaxis.label.set_fontsize(size)
#     axs[0,0].xaxis.label.set_fontsize(size)
    
#     #SFR Scaling relation----------------------------------
#     imglatelog2 = axs[0,1].pcolormesh(xbins_starmass, ybins_starmass, SFR_den, cmap='viridis', vmin = -3, vmax = 1) #, vmin = 6.5, vmax = 10 , vmin = -2.5, vmax = 0.5
#     axs[0,1].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[0,1].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog2, ax=axs[0,1])
#     axs[0,1].title.set_text("log SFR Scaling Relation")    
#     axs[0,1].title.set_fontsize(size)
#     axs[0,1].yaxis.label.set_fontsize(size)
#     axs[0,1].xaxis.label.set_fontsize(size)
    
#     #DELTA SFR ----------------------------------
#     Diff = SFR_mapping - SFR_den
#     imglatelog3 = axs[0,2].pcolormesh(xbins_starmass, ybins_starmass, Diff, cmap='viridis') #, vmin = 6.5, vmax = 10
#     axs[0,2].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[0,2].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog3, ax=axs[0,2])
#     axs[0,2].title.set_text(r"$\Delta SFR$ (log difference)")    
#     axs[0,2].title.set_fontsize(size)
#     axs[0,2].yaxis.label.set_fontsize(size)
#     axs[0,2].xaxis.label.set_fontsize(size)
    
#     #linear SFR TNG ----------------------------------
#     imglatelog4 = axs[1,0].pcolormesh(xbins_starmass, ybins_starmass, den_SFR, cmap='viridis', vmin = 0, vmax = 10) #, vmin = 6.5, vmax = 10 , vmin = 0, vmax = 3.5
#     axs[1,0].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[1,0].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog4, ax=axs[1,0])
#     axs[1,0].title.set_text("Linear SFR TNG Data")    
#     axs[1,0].title.set_fontsize(size)
#     axs[1,0].yaxis.label.set_fontsize(size)
#     axs[1,0].xaxis.label.set_fontsize(size)
    
#     #linear SFR Scaling Relation----------------------------------
#     imglatelog5 = axs[1,1].pcolormesh(xbins_starmass, ybins_starmass, exp_SFR, cmap='viridis', vmin = 0, vmax = 10) #, vmin = 6.5, vmax = 10  , vmin = 0, vmax = 3.5
#     axs[1,1].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[1,1].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog5, ax=axs[1,1])
#     axs[1,1].title.set_text("Linear SFR Scaling Relation")    
#     axs[1,1].title.set_fontsize(size)
#     axs[1,1].yaxis.label.set_fontsize(size)
#     axs[1,1].xaxis.label.set_fontsize(size)
    
#     #linear delta SFR (Difference)----------------------------------
#     diff10 = den_SFR - exp_SFR
#     imglatelog6 = axs[1,2].pcolormesh(xbins_starmass, ybins_starmass, diff10, cmap='viridis') #, vmin = 6.5, vmax = 10
#     axs[1,2].set_xlim([-cutoff_vol, cutoff_vol])
#     axs[1,2].set_ylim([-cutoff_vol, cutoff_vol])
#     fig3.colorbar(imglatelog6, ax=axs[1,2])
#     axs[1,2].title.set_text(r"$\delta SFR$ (difference)")    
#     axs[1,2].title.set_fontsize(size)
#     axs[1,2].yaxis.label.set_fontsize(size)
#     axs[1,2].xaxis.label.set_fontsize(size)
#     #-------------------------------------------------------------------------------------------------
    
#     # x & y labels 
#     fig3.text(0.5, 0.04, 'x (kpc)', ha='center', va='center',fontsize=size)
#     fig3.text(0.06, 0.5, 'y (kpc)', ha='center', va='center', rotation='vertical', fontsize=size)
    
#     #plt.show()
#     figure3 = plt.gcf()  # get current figure
#     figure3.set_size_inches(20, 10) # set figure's size manually to your full screen (32x18)   
#     plt.savefig('D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Second Progress Report\All in One SFR Maps\AllinOneSFRSubhalo_'+ str(id) +'_' + str(snapshot) +'.png')
    
    #=======================+++++++++++++++++++++++++++++++++++========================================
    #Custom Colormap testing
    FinalCmap = CustomCmap('viridis', 'autumn_r', 50, 256, 'FinalCmap')
    
    #Spatial Maps for Contributions of SFE & fgas ------------------------------------------------------------------------------
    fig3, axs = plt.subplots(2,3)
    title = "Spatial Maps for Contributions of SFE & fgas - Subhalo_" + str(id) + "_" + str(snapshot)
    fig3.suptitle(title, fontsize=size)
    plt.grid()

    # if snapshot == 69:
    #     emin, emax = -0.2 , 0.8
    # elif snapshot == 71: 
    #     emin, emax = 0, 6
    # elif snapshot == 73:
    #     emin, emax = -0.5, 1.75
    
    # #Central Contribution (Within stellar HMR)
    # for spaxel in 
    # if(xbins_starmass[])
    circle1 = plt.Circle((0,0),halfmassrad_stars , color = 'y',fill = False)
    circle2 = plt.Circle((0,0),halfmassrad_stars , color = 'y',fill = False)
    circle3 = plt.Circle((0,0),halfmassrad_stars , color = 'y',fill = False)
    
    #Efficiency contribution----------------------------------
    #Histogram doesn't follow cartesian convention, take transpose for visualization purposes
    imglatelog1 = axs[0,1].pcolormesh(xbins_starmass, ybins_starmass, sSFRe.T, cmap=FinalCmap, vmin = -0.5, vmax = 6) #, vmin = -0.5, vmax = 6   
    axs[0,1].add_patch( circle1 )
    axs[0,1].set_xlim([-cutoff_vol, cutoff_vol])
    axs[0,1].set_ylim([-cutoff_vol, cutoff_vol])
    fig3.colorbar(imglatelog1, ax=axs[0,1])
    axs[0,1].title.set_text(r"$\delta SFR_{e}$ $(yr^{-1})$")    
    axs[0,1].title.set_fontsize(size)
    axs[0,1].yaxis.label.set_fontsize(size)
    axs[0,1].xaxis.label.set_fontsize(size)
    
    #Gas Fraction contribution----------------------------------
    imglatelog2 = axs[0,2].pcolormesh(xbins_starmass, ybins_starmass, sSFRg.T, cmap=FinalCmap, vmin = -0.5, vmax = 6) #, vmin = -0.5, vmax = 6
    axs[0,2].add_patch( circle2 )
    axs[0,2].set_xlim([-cutoff_vol, cutoff_vol])
    axs[0,2].set_ylim([-cutoff_vol, cutoff_vol])
    fig3.colorbar(imglatelog2, ax=axs[0,2])
    axs[0,2].title.set_text(r"$\delta SFR_{g}$ $(M_{\odot} kpc^{-2})$")    
    axs[0,2].title.set_fontsize(size)
    axs[0,2].yaxis.label.set_fontsize(size)
    axs[0,2].xaxis.label.set_fontsize(size)
    
    #delta sSFR----------------------------------
    imglatelog3 = axs[0,0].pcolormesh(xbins_starmass, ybins_starmass, diffSFR.T, cmap=FinalCmap, vmin = -0.737691, vmax = 9) #, vmin = -0.5, vmax = 9
    axs[0,0].add_patch( circle3 )
    axs[0,0].set_xlim([-cutoff_vol, cutoff_vol])
    axs[0,0].set_ylim([-cutoff_vol, cutoff_vol])
    fig3.colorbar(imglatelog3, ax=axs[0,0])
    axs[0,0].title.set_text(r"$\delta \Sigma_{SFR}$ $(M_{\odot} kpc^{-2} yr^{-1})$") #r"$\delta sSFR$ $(yr^{-1})$"
    axs[0,0].title.set_fontsize(size)
    axs[0,0].yaxis.label.set_fontsize(size)
    axs[0,0].xaxis.label.set_fontsize(size)
    
    #Effeciency contribution Percent ----------------------------------
    imglatelog4 = axs[1,0].pcolormesh(xbins_starmass, ybins_starmass, Diff_SFE.T, cmap='viridis', vmin = -1, vmax = 0.8) #, vmin = 6.5, vmax = 10 , vmin = 0, vmax = 3.5
    #axs[1,0].add_patch( circle1 )
    axs[1,0].set_xlim([-cutoff_vol, cutoff_vol])
    axs[1,0].set_ylim([-cutoff_vol, cutoff_vol])
    fig3.colorbar(imglatelog4, ax=axs[1,0])
    axs[1,0].title.set_text(r"$\Delta SFE$ $(yr^{-1})$")    
    axs[1,0].title.set_fontsize(size)
    axs[1,0].yaxis.label.set_fontsize(size)
    axs[1,0].xaxis.label.set_fontsize(size)
    
    #linear SFR Scaling Relation----------------------------------
    imglatelog5 = axs[1,1].pcolormesh(xbins_starmass, ybins_starmass, Diff_fGas.T, cmap='viridis', vmin = -1, vmax = 1.25) #, vmin = 6.5, vmax = 10  , vmin = 0, vmax = 3.5
    #axs[1,1].add_patch( circle1 )
    axs[1,1].set_xlim([-cutoff_vol, cutoff_vol])
    axs[1,1].set_ylim([-cutoff_vol, cutoff_vol])
    fig3.colorbar(imglatelog5, ax=axs[1,1])
    axs[1,1].title.set_text(r"$\Delta fGas$ (Fraction)")    
    axs[1,1].title.set_fontsize(size)
    axs[1,1].yaxis.label.set_fontsize(size)
    axs[1,1].xaxis.label.set_fontsize(size)
    
    #linear delta SFR (Difference)----------------------------------
    imglatelog6 = axs[1,2].pcolormesh(xbins_starmass, ybins_starmass, Diff_sSFR.T, cmap='viridis', vmin = -1, vmax = 2) #, vmin = 6.5, vmax = 10
    #axs[1,2].add_patch( circle1 )
    axs[1,2].set_xlim([-cutoff_vol, cutoff_vol])
    axs[1,2].set_ylim([-cutoff_vol, cutoff_vol])
    fig3.colorbar(imglatelog6, ax=axs[1,2])
    axs[1,2].title.set_text(r"$\Delta sSFR$ $(yr^{-1})$")    
    axs[1,2].title.set_fontsize(size)
    axs[1,2].yaxis.label.set_fontsize(size)
    axs[1,2].xaxis.label.set_fontsize(size)
    # #-------------------------------------------------------------------------------------------------
    
    # x & y labels 
    fig3.text(0.5, 0.04, 'x (kpc)', ha='center', va='center',fontsize=size)
    fig3.text(0.06, 0.5, 'y (kpc)', ha='center', va='center', rotation='vertical', fontsize=size)
    
    #plt.show()
    figure3 = plt.gcf()  # get current figure
    figure3.set_size_inches(20, 10) # set figure's size manually to your full screen (32x18)   
    plt.savefig('D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Second Progress Report\Contributions Spatial Map\ContributionsSpatialSubhalo_'+ str(id) +'_' + str(snapshot) +'.png')
    
    # #=======================+++++++++++++++++++++++++++++++++++========================================
 # #   SFE vs Gmass maps
    
 #    #SUBPLOTS - Ratio Log Graphs
 #    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,10))
 #    size = 20
 #    title = "SFE vs fgas"
 #    fig.suptitle(title, fontsize=size)
 #    cm = plt.cm.get_cmap('viridis')
    

 #    #Log difference----------------------------------
 #    sc1 = ax1.scatter(sSFRg, sSFRe, c = diffSFR, cmap = cm)
 #    fig.colorbar(sc1, ax=ax1)
 #    ax1.title.set_text("Log difference SFE vs Mgas - Subhalo_" + str(id) + "_" + str(snapshot))
 #    ax1.set_ylabel("SFRe")
 #    ax1.set_xlabel("SFRg")
 #    ax1.legend(labels = "Points", loc = 'top right', fontsize = 10)
 #    # ax1.set_ylabel(r"$\Delta SFE$ $(\log$ $yr^{-1})$")
 #    # ax1.set_xlabel(r"$\Delta M_{gas}$ (Fraction)")
 #    ax1.title.set_fontsize(size)
 #    ax1.yaxis.label.set_fontsize(size)
 #    ax1.xaxis.label.set_fontsize(size)
 #    # #now plot both limits against eachother
 #    # Diff_fGas1 = Diff_fGas
 #    # Diff_SFE1 = Diff_SFE
 #    # Diff_fGas1[np.isinf(Diff_fGas1)] = 0
 #    # Diff_SFE1[np.isinf(Diff_SFE1)] = 0
 #    # lims1 = [
 #    # np.min([np.nanmin(Diff_fGas1), np.nanmin(Diff_SFE1)]),  # min of both axes
 #    # np.max([np.nanmax(Diff_fGas1), np.nanmax(Diff_SFE1)]),  # max of both axes
 #    # ]
 #    # ax1.plot(lims1, lims1, 'k-', alpha=0.75, zorder=0)
 #    # ax1.set_aspect('equal')
    
 #    #Difference----------------------------------
 #    sc2 = ax2.scatter(diffGmass, diffSFE, c = diffSFR, cmap = cm)
 #    fig.colorbar(sc2, ax=ax2)
 #    #ax2.plot(Total_GSMass, ModelFit(Gbeta, Total_GSMass), color = 'red')
 #    ax2.title.set_text("Difference SFE vs Mgas")    
 #    ax2.set_ylabel(r"$\delta SFE$ $(\log$ $yr^{-1})$")
 #    ax2.set_xlabel(r"$\delta M_{gas}$ (Fraction)")
 #    #ax2.legend(( Gline,'datapoint'))
 #    ax2.title.set_fontsize(size)
 #    ax2.yaxis.label.set_fontsize(size)
 #    ax2.xaxis.label.set_fontsize(size)
    
 #    # #now plot both limits against eachother
 #    # diffGas1 = diffGas
 #    # diffSFE1 = diffSFE
 #    # diffGas1[np.isinf(diffGas1)] = 0
 #    # diffSFE1[np.isinf(diffSFE1)] = 0
 #    # lims2 = [
 #    # np.min([np.nanmin(diffGas1), np.nanmin(diffSFE1)]),  # min of both axes
 #    # np.max([np.nanmax(diffGas1), np.nanmax(diffSFE1)]),  # max of both axes
 #    # ]
 #    # ax2.plot(lims2, lims2, 'k-', alpha=0.75, zorder=0)
 #    # ax2.set_aspect('equal')
 #    # ax2.set_xlim(lims2)
 #    # ax2.set_ylim(lims2)
    
 #    plt.show()
 #    figure4 = plt.gcf()  # get current figure
 #    figure4.set_size_inches(20, 10) # set figure's size manually to your full screen (32x18)   
 #    plt.savefig('D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Second Progress Report\SFE vs Mgas Maps\SFEvsMgasSubhalo_'+ str(id) +'_' + str(snapshot) +'.png')
    
   
    
    
    # #SFE vs fgas maps
    
    # #SUBPLOTS - Ratio Log Graphs
    # fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,10))
    # size = 20
    # title = "SFE vs fgas"
    # fig.suptitle(title, fontsize=size)
    # cm = plt.cm.get_cmap('viridis')
    

    # #SFR----------------------------------
    # sc1 = ax1.scatter(Diff_fGas, Diff_SFE, c = diffSFR, cmap = cm)
    # fig.colorbar(sc1, ax=ax1)
    # ax1.title.set_text("Ratio Log SFE vs fgas - Subhalo_" + str(id) + "_" + str(snapshot))
    # ax1.set_ylabel(r"$\Delta SFE$ $(\log$ $yr^{-1})$")
    # ax1.set_xlabel(r"$\Delta f_{gas}$ (Fraction)")
    # ax1.title.set_fontsize(size)
    # ax1.yaxis.label.set_fontsize(size)
    # ax1.xaxis.label.set_fontsize(size)
    # #now plot both limits against eachother
    # Diff_fGas1 = Diff_fGas
    # Diff_SFE1 = Diff_SFE
    # Diff_fGas1[np.isinf(Diff_fGas1)] = 0
    # Diff_SFE1[np.isinf(Diff_SFE1)] = 0
    # lims1 = [
    # np.min([np.nanmin(Diff_fGas1), np.nanmin(Diff_SFE1)]),  # min of both axes
    # np.max([np.nanmax(Diff_fGas1), np.nanmax(Diff_SFE1)]),  # max of both axes
    # ]
    # ax1.plot(lims1, lims1, 'k-', alpha=0.75, zorder=0)
    # ax1.set_aspect('equal')
    # # ax1.set_xlim(lims1)
    # # ax1.set_ylim(lims1)
    
    # # #Percentage SFE
    # # per = Diff_SFE / Diff_fGas
    
    # # per_SFE = Diff_SFE / lims1[1] * 100
    # # sum_SFE = sum(per_SFE) / len(per_SFE)
    # # str_SFE = str(round(sum_SFE,2))
    # # ax1.text(np.nanmin(Diff_fGas1), np.nanmax(Diff_SFE1), str_SFE, ha='center', va='center',fontsize=size)
    
    # # per_fGas = Diff_fGas / lims1[1] * 100
    # # sum_fGas = sum(per_fGas) / len(per_fGas)
    # # str_fGas = str(round(sum_fGas,2))
    # # ax1.text(np.nanmin(Diff_SFE1), np.nanmax(Diff_fGas1), str_fGas, ha='center', va='center',fontsize=size)
    
    # #starmass----------------------------------
    # sc2 = ax2.scatter(diffGas, diffSFE, c = diffSFR, cmap = cm)
    # fig.colorbar(sc2, ax=ax2)
    # #ax2.plot(Total_GSMass, ModelFit(Gbeta, Total_GSMass), color = 'red')
    # ax2.title.set_text("Difference Log SFE vs fgas")    
    # ax2.set_ylabel(r"$\Delta SFE$ $(\log$ $yr^{-1})$")
    # ax2.set_xlabel(r"$\Delta f_{gas}$ (Fraction)")
    # #ax2.legend(( Gline,'datapoint'))
    # ax2.title.set_fontsize(size)
    # ax2.yaxis.label.set_fontsize(size)
    # ax2.xaxis.label.set_fontsize(size)
    
    # #now plot both limits against eachother
    # diffGas1 = diffGas
    # diffSFE1 = diffSFE
    # diffGas1[np.isinf(diffGas1)] = 0
    # diffSFE1[np.isinf(diffSFE1)] = 0
    # lims2 = [
    # np.min([np.nanmin(diffGas1), np.nanmin(diffSFE1)]),  # min of both axes
    # np.max([np.nanmax(diffGas1), np.nanmax(diffSFE1)]),  # max of both axes
    # ]
    # ax2.plot(lims2, lims2, 'k-', alpha=0.75, zorder=0)
    # ax2.set_aspect('equal')
    # ax2.set_xlim(lims2)
    # ax2.set_ylim(lims2)
    
    #plt.show()
    # figure4 = plt.gcf()  # get current figure
    # figure4.set_size_inches(20, 10) # set figure's size manually to your full screen (32x18)   
    # plt.savefig('D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Second Progress Report\SFE vs fgas Maps\SFEvsfgasSubhalo_'+ str(id) +'_' + str(snapshot) +'.png')
    
    
    
    
    
    
    #     #=======================+++++++++++++++++++++++++++++++++++========================================
    
    # #TEST - LATE LOG------------------------------------------------------------------------------
    # fig3, axs = plt.subplots(2,3)
    # title = "Late Log Difference Plots log(a - b) :- Subhalo_" + str(id) + "_" + str(snapshot)
    # fig3.suptitle(title, fontsize=size)
    # plt.grid()
    
    # #SFR----------------------------------
    # imglatelog1 = axs[0,0].pcolormesh(xbins_starmass, ybins_starmass, diffSFR, cmap='viridis', vmin = -4, vmax = 0.5) #, vmin = 6.5, vmax = 10
    # axs[0,0].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[0,0].set_ylim([-cutoff_vol, cutoff_vol])
    # fig3.colorbar(imglatelog1, ax=axs[0,0])
    # axs[0,0].title.set_text(r"$\Delta \Sigma_{SFR}$ $(\log$ $M_{\odot} kpc^{-2} yr^{-1})$")    
    # axs[0,0].title.set_fontsize(size)
    # axs[0,0].yaxis.label.set_fontsize(size)
    # axs[0,0].xaxis.label.set_fontsize(size)
    
    # #starmass----------------------------------
    # imglatelog2 = axs[0,1].pcolormesh(xbins_starmass, ybins_starmass, S_Mass, cmap='viridis', vmin = 6.5, vmax = 9.5) #, vmin = 6.5, vmax = 10
    # axs[0,1].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[0,1].set_ylim([-cutoff_vol, cutoff_vol])
    # fig3.colorbar(imglatelog2, ax=axs[0,1])
    # axs[0,1].title.set_text(r"$\Delta \Sigma_{\ast}$ $(\log$ $M_{\odot} kpc^{-2})$")    
    # axs[0,1].title.set_fontsize(size)
    # axs[0,1].yaxis.label.set_fontsize(size)
    # axs[0,1].xaxis.label.set_fontsize(size)
    
    # #Gas Mass----------------------------------
    # imglatelog3 = axs[0,2].pcolormesh(xbins_starmass, ybins_starmass, diffGmass, cmap='viridis', vmin = 5, vmax = 9) #, vmin = 6.5, vmax = 10
    # axs[0,2].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[0,2].set_ylim([-cutoff_vol, cutoff_vol])
    # fig3.colorbar(imglatelog3, ax=axs[0,2])
    # axs[0,2].title.set_text(r"$\Delta \Sigma_{gas}$ $(\log$ $M_{\odot} kpc^{-2})$")    
    # axs[0,2].title.set_fontsize(size)
    # axs[0,2].yaxis.label.set_fontsize(size)
    # axs[0,2].xaxis.label.set_fontsize(size)
    
    # #sSFR----------------------------------
    # imglatelog4 = axs[1,0].pcolormesh(xbins_starmass, ybins_starmass, diffsSFR, cmap='viridis', vmin = -13, vmax = -8) #, vmin = 6.5, vmax = 10
    # axs[1,0].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[1,0].set_ylim([-cutoff_vol, cutoff_vol])
    # fig3.colorbar(imglatelog4, ax=axs[1,0])
    # axs[1,0].title.set_text(r"$\Delta sSFR$ $(\log$ $yr^{-1})$")    
    # axs[1,0].title.set_fontsize(size)
    # axs[1,0].yaxis.label.set_fontsize(size)
    # axs[1,0].xaxis.label.set_fontsize(size)
    
    # #SFE----------------------------------
    # diffSFE = den_SFE - exp_SFE
    # imglatelog5 = axs[1,1].pcolormesh(xbins_starmass, ybins_starmass, np.log10(np.abs(diffSFE)), cmap='viridis', vmin = -13, vmax = -8.5) #, vmin = 6.5, vmax = 10
    # axs[1,1].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[1,1].set_ylim([-cutoff_vol, cutoff_vol])
    # fig3.colorbar(imglatelog5, ax=axs[1,1])
    # axs[1,1].title.set_text(r"$\Delta SFE$ $(\log$ $yr^{-1})$")    
    # axs[1,1].title.set_fontsize(size)
    # axs[1,1].yaxis.label.set_fontsize(size)
    # axs[1,1].xaxis.label.set_fontsize(size)
    
    # #fgas----------------------------------
    # diffGas = den_fH - exp_fH
    # imglatelog6 = axs[1,2].pcolormesh(xbins_starmass, ybins_starmass, np.log10(np.abs(diffGas)), cmap='viridis', vmin = -4, vmax = 1.5) #, vmin = 6.5, vmax = 10
    # axs[1,2].set_xlim([-cutoff_vol, cutoff_vol])
    # axs[1,2].set_ylim([-cutoff_vol, cutoff_vol])
    # fig3.colorbar(imglatelog6, ax=axs[1,2])
    # axs[1,2].title.set_text(r"$\Delta f_{gas}$ (log Fraction)")    
    # axs[1,2].title.set_fontsize(size)
    # axs[1,2].yaxis.label.set_fontsize(size)
    # axs[1,2].xaxis.label.set_fontsize(size)
    # #-------------------------------------------------------------------------------------------------
    
    # # x & y labels 
    # fig3.text(0.5, 0.04, 'x (kpc)', ha='center', va='center',fontsize=size)
    # fig3.text(0.06, 0.5, 'y (kpc)', ha='center', va='center', rotation='vertical', fontsize=size)
    
    # #plt.show()
    # figure3 = plt.gcf()  # get current figure
    # figure3.set_size_inches(20, 10) # set figure's size manually to your full screen (32x18)   
    # plt.savefig('D:\Shivani\Academic Trent 2022-23\Fall 2022\Project Course\Second Progress Report\Late Log Difference Maps\LateLogDiffSubhalo_'+ str(id) +'_' + str(snapshot) +'.png')
    
    # #=======================+++++++++++++++++++++++++++++++++++========================================
    
    
    
    
    
    