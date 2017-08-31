import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import linecache
import scipy.io
import glob
import operator
import sys

## WARM LHC REGIONS (can be put into an external file, but for convenience we leave it here)

lhc_warm=np.array([[  0.00000000e+00,   2.25365000e+01],
       [  5.48530000e+01,   1.52489000e+02],
       [  1.72165500e+02,   1.92400000e+02],
       [  1.99484700e+02,   2.24300000e+02],
       [  3.09545428e+03,   3.15562858e+03],
       [  3.16774008e+03,   3.18843308e+03],
       [  3.21144458e+03,   3.26386758e+03],
       [  3.30990008e+03,   3.35497408e+03],
       [  3.40100558e+03,   3.45342858e+03],
       [  3.47644008e+03,   3.49406558e+03],
       [  3.50588528e+03,   3.56831858e+03],
       [  6.40540880e+03,   6.45791380e+03],
       [  6.46877850e+03,   6.85951380e+03],
       [  6.87037850e+03,   6.92353380e+03],
       [  9.73590702e+03,   9.82473052e+03],
       [  9.83083202e+03,   9.86173052e+03],
       [  9.87873202e+03,   9.93998552e+03],
       [  9.95054802e+03,   1.00434620e+04],
       [  1.00540245e+04,   1.01152780e+04],
       [  1.01322795e+04,   1.01639705e+04],
       [  1.01700720e+04,   1.02576030e+04],
       [  1.31049892e+04,   1.31298045e+04],
       [  1.31368892e+04,   1.31571237e+04],
       [  1.31768002e+04,   1.32716472e+04],
       [  1.33067527e+04,   1.33518257e+04],
       [  1.33869312e+04,   1.34817782e+04],
       [  1.35014547e+04,   1.35227845e+04],
       [  1.35298692e+04,   1.35546845e+04],
       [  1.63946378e+04,   1.64508713e+04],
       [  1.64569728e+04,   1.64872713e+04],
       [  1.64933728e+04,   1.68308713e+04],
       [  1.68369728e+04,   1.68672713e+04],
       [  1.68733728e+04,   1.69282948e+04],
       [  1.97348504e+04,   1.97606997e+04],
       [  1.97715644e+04,   2.02179087e+04],
       [  2.02287734e+04,   2.02529744e+04],
       [  2.30899797e+04,   2.31385770e+04],
       [  2.31503967e+04,   2.31713755e+04],
       [  2.31943870e+04,   2.32468100e+04],
       [  2.32928425e+04,   2.33379155e+04],
       [  2.33839480e+04,   2.34363710e+04],
       [  2.34593825e+04,   2.34800825e+04],
       [  2.34921940e+04,   2.35531160e+04],
       [  2.64334879e+04,   2.64583032e+04],
       [  2.64653879e+04,   2.64867177e+04],
       [  2.65063942e+04,   2.66012412e+04],
       [  2.66363467e+04,   2.66588832e+04]])

hllhc_warm=np.array([[  0.00000000e+00,   2.25000000e+01],
       [  8.31530000e+01,   1.36689000e+02],
       [  1.82965500e+02,   2.01900000e+02],
       [  2.10584700e+02,   2.24300000e+02],
       [  3.09545428e+03,   3.15562858e+03],
       [  3.16774008e+03,   3.18843308e+03],
       [  3.21144458e+03,   3.26386758e+03],
       [  3.30990008e+03,   3.35497408e+03],
       [  3.40100558e+03,   3.45342858e+03],
       [  3.47644008e+03,   3.49406558e+03],
       [  3.50588528e+03,   3.56831858e+03],
       [  6.40540880e+03,   6.45791380e+03],
       [  6.46877850e+03,   6.85951380e+03],
       [  6.87037850e+03,   6.92353380e+03],
       [  9.73590702e+03,   9.82473052e+03],
       [  9.83083202e+03,   9.86173052e+03],
       [  9.87873202e+03,   9.93998552e+03],
       [  9.95054802e+03,   1.00434620e+04],
       [  1.00540245e+04,   1.01152780e+04],
       [  1.01322795e+04,   1.01639705e+04],
       [  1.01700720e+04,   1.02576030e+04],
       [  1.31036000e+04,   1.31200300e+04],
       [  1.31238892e+04,   1.31471237e+04],
       [  1.31918002e+04,   1.32476472e+04],
       [  1.33067940e+04,   1.33520892e+04],
       [  1.34110312e+04,   1.34670082e+04],
       [  1.35114547e+04,   1.35357845e+04],
       [  1.35388592e+04,   1.35552845e+04],
       [  1.63946378e+04,   1.64508713e+04],
       [  1.64569728e+04,   1.64872713e+04],
       [  1.64933728e+04,   1.68308713e+04],
       [  1.68369728e+04,   1.68672713e+04],
       [  1.68733728e+04,   1.69282948e+04],
       [  1.97348504e+04,   1.97606997e+04],
       [  1.97715644e+04,   2.02179087e+04],
       [  2.02287734e+04,   2.02529744e+04],
       [  2.30899797e+04,   2.31385770e+04],
       [  2.31503967e+04,   2.31713755e+04],
       [  2.31943870e+04,   2.32468100e+04],
       [  2.32928425e+04,   2.33379155e+04],
       [  2.33839480e+04,   2.34363710e+04],
       [  2.34593825e+04,   2.34800825e+04],
       [  2.34921940e+04,   2.35531160e+04],
       [  2.64334879e+04,   2.64483032e+04],
       [  2.64569832e+04,   2.64759232e+04],
       [  2.65221932e+04,   2.65757332e+04],
       [  2.66363832e+04,   2.66588832e+04]])


# USER INPUT

simulations = [{
    'beam' : 1,
    'batch_path' : '/Users/swretbor/Workspace/work_afs/sixtrack/simulations/oml_study/v26/',
    'pathtosim' : 'run*/',
    'pathtoLPIs' : 'LP*_BLP_out.s',
    'pathtoCollSum' : 'coll_summary.dat',
    'collPosFileName' : 'clean_input/CollPositions.b1.dat',
},{
    'beam' : 2,
    'batch_path' : '/Users/swretbor/Workspace/work_afs/sixtrack/simulations/oml_study/v27/',
    'pathtosim' : 'run*/',
    'pathtoLPIs' : 'LP*_BLP_out.s',
    'pathtoCollSum' : 'coll_summary.dat',
    'collPosFileName' : 'clean_input/CollPositions.b2.dat',
}]
SimulationName  = "Lossmap injection, ToyModel 11s + SixTrack 10k. Linear distribution."     # title of the simulation
dataFileName    = 'losses.dat'           # name of file where to save data
                                         # set it to None in case you don't want data to be saved in a file
# beam4 = beam == 2
# b4inB1RefSys = beam == 2
plotCounts = False
clhc = 26658.8832 # [m]
# in case lattices do not start at IP1, set this value to the original value of the
#    starting element as if IP1 was at 0.0
# NB: in case of B4, do not put here the value in the B1 ref sys
s0 = 0.0          # [m]


simu1co = []
totcollosses1=0                                                           # total # of lost particles at collim
ap_loss_cold = np.array([])
ap_loss_warm = np.array([])
for s in simulations:

    ##################################### INITIALIZATION #########################################################
    print(s)
    pathtosim = s['batch_path'] + s['pathtosim']
    pathtoLPIs = pathtosim + s['pathtoLPIs']
    pathtoCollSum = pathtosim + s['pathtoCollSum']
    collPosFileName = s['batch_path'] + s['collPosFileName']
    beam4 = s['beam'] == 2
    b4inB1RefSys = s['beam'] == 2

    ## Collimator List (mandatory!)

    copos=np.loadtxt(collPosFileName, skiprows=1, usecols=(0, 2))
    if ( beam4 ):
        warm_parts = clhc - lhc_warm
        warm_parts = np.flipud(np.fliplr(warm_parts))
    else:
        warm_parts = lhc_warm
    
    
    ##################################### PROCESSING OF APERTURE LOSSES ###########################################
    
    
    ######## SIMULATION 1
    
    simu1ap = []
    # aperture losses
    n = 0
    print("Reading LPI files")
    for fname in glob.glob(pathtoLPIs):                          # loop over directories
        array=np.loadtxt(fname,ndmin=2)                          # ndmin=2: minimum 2 dimensions (single line files)
        if len(array)!=0:
            simu1ap.append(array[:,2])
        print("LPI files: {}".format(n), end="\r")
        n += 1
    print("LPI files: {}".format(n))
    simu1ap=[item for sublist in simu1ap for item in sublist]    # flatten array
    
    if ( s0 != 0.0 ):
        for ii in range( len( simu1ap ) ):
            simu1ap[ii]+=s0
            if ( simu1ap[ii] > clhc ):
                simu1ap[ii]-=clhc
    
    
    ######## DIVIDE INTO COLD AND WARM LOSSES
    
    aplosses = simu1ap
    
    ap_losses_wc = np.zeros(shape=(len(aplosses),2))
    
    for i in range(len(ap_losses_wc)):
        mark = 0
        for j in range(len(warm_parts)):
            if (warm_parts[j,0] < aplosses[i] and aplosses[i] < warm_parts[j,1]):
                mark = 1
        ap_losses_wc[i] = [aplosses[i],mark]
    
    ap_losses_cold = []
    ap_losses_warm = []
    
    
    for i in range(len(ap_losses_wc)):
        if (ap_losses_wc[i,1] == 0):
            ap_losses_cold.append(ap_losses_wc[i,0])
        else:
            ap_losses_warm.append(ap_losses_wc[i,0])
    
    ap_losses_cold1=ap_losses_cold
    ap_losses_warm1=ap_losses_warm
    
    ApertureLossesColdSample1=len(ap_losses_cold)
    ApertureLossesWarmSample1=len(ap_losses_warm)
    
    
    #### COLLIMATOR LOSSES
    collimators = []
    for j in range(0,len(copos)):
        collimators.append([int(copos[j][0]),float(copos[j][1]),0,0])
        if b4inB1RefSys:
            collimators[-1][1] = clhc - collimators[-1][1]
    
    n = 0
    print("Reading coll_summary files")
    for fname in glob.glob(pathtoCollSum):
        # print fname
        array=np.loadtxt(fname, usecols=(0, 3, 6))
        for j in range(len(array)):
            id=np.where(int(array[j][0])==np.array(collimators)[:,0])[0][0]
            collimators[id][2] += int(array[j][1])
            collimators[id][3] = float(array[j][2])                               # write collimator length
            totcollosses1 += int(array[j][1])                                   # increase total lossnumber
        print("Coll summary files: {}".format(n), end="\r")
        n += 1
    print("Coll summary files: {}".format(n))

    simu1co += collimators
    
    if b4inB1RefSys:
        ap_loss_cold = np.concatenate((ap_loss_cold, clhc - np.array(ap_losses_cold1)))
        ap_loss_warm = np.concatenate((ap_loss_warm, clhc - np.array(ap_losses_warm1)))
    else:
        ap_loss_cold = np.concatenate((ap_loss_cold, np.array(ap_losses_cold1)))
        ap_loss_warm = np.concatenate((ap_loss_warm, np.array(ap_losses_warm1)))

ap_losses_cold1 = ap_loss_cold
ap_losses_warm1 = ap_loss_warm

################################## PLOT #################################

#### CREATE HISTOGRAM DATA

fig, (ax0) = plt.subplots()


nbins=int(clhc/0.10)
if ( plotCounts ):
    ymax = 1.0E+08
    ymin = 0.1
else:
    ymax = 3.0
    ymin = 1.0E-07
# bar width=80% of bin width 
width=0.80*clhc/nbins
binWidth=clhc/nbins
# total number of tracked particles
normfac = float(totcollosses1+ApertureLossesColdSample1+ApertureLossesWarmSample1)                 


## Cold Losses

#  - Generate histogram
hist_cold, bins = np.histogram(ap_losses_cold1, bins=nbins,range=(0,clhc))
#  - Remove empty entries
bins_cold = np.delete((bins[:-1] + bins[1:]) / 2,np.where(hist_cold==0)[0])
if ( plotCounts ):
    hist_cold = np.delete(hist_cold,np.where(hist_cold==0)[0])
else:
    hist_cold = np.delete(hist_cold,np.where(hist_cold==0)[0])/(binWidth*normfac)
#  - lossmap of cold elements
ax0.bar(bins_cold, hist_cold, width=width , color='blue', edgecolor='blue', label='Cold')


## Warm Losses

#  - Generate histogram
hist_warm, bins = np.histogram(ap_losses_warm1, bins=nbins,range=(0,clhc))
#  - Remove empty entries
bins_warm = np.delete( (bins[:-1] + bins[1:]) / 2, np.where(hist_warm==0)[0])
if ( plotCounts ):
    hist_warm = np.delete(hist_warm,np.where(hist_warm==0)[0])
else:
    hist_warm = np.delete(hist_warm,np.where(hist_warm==0)[0])/(binWidth*normfac)
#  - lossmap of warm elements
ax0.bar(bins_warm, hist_warm, width=width , color='red', edgecolor='red', label='Warm')


## Collimator Losses
#  - hack to get the correct label
if ( plotCounts ):
    ax0.bar(simu1co[0][1], simu1co[0][2], width=simu1co[0][3] , color='black', edgecolor='black', label='Collimator')
else:
    if simu1co[0][3] > 0.0:
        ax0.bar(simu1co[0][1], simu1co[0][2]/(normfac*simu1co[0][3]), width=simu1co[0][3] , color='black', edgecolor='black', label='Collimator')
#  - lossmap of collimators
if ( plotCounts ):
    for k in range(1,len(simu1co)):
        if ( simu1co[k][3] > 0.0 ):
            ax0.bar(simu1co[k][1], simu1co[k][2], width=simu1co[k][3] , color='black', edgecolor='black')
else:
    for k in range(1,len(simu1co)):
        if ( simu1co[k][3] > 0.0 ):
            ax0.bar(simu1co[k][1], simu1co[k][2]/(normfac*simu1co[k][3]), width=simu1co[k][3] , color='black', edgecolor='black')


## PLOT PROPERTIES

ax0.set_yscale('log', nonposy='clip')
ax0.set_xlim(0,clhc)
ax0.set_ylim(ymin,ymax)
ax0.set_xlabel(r'Longitudinal Coordinate (m)')
if ( plotCounts ):
    ax0.set_ylabel(r'counts ()')
else:
    ax0.set_ylabel(r'Local cleaning inefficiency $\eta$ (1/m)')
ax0.xaxis.grid(True)
ax0.yaxis.grid(False)
ax0.grid()
ax0.legend(fontsize = 10, loc=1, borderaxespad=0.)

ax0.set_title(SimulationName)

# save plot
plt.savefig('LM_LHC.pdf',bbox_inches='tight')

plt.show()

# in case, save data in txt file:
if ( dataFileName is not None ):
    print(' saving histograms in file %s ...' % ( dataFileName ))
    oFile = open( dataFileName, 'w' )
    oFile.write( '# normalisation: %i \n' % ( normfac ) )
    oFile.write( '# cold/warm nbins, bin width [m]: %i %.8e \n' % ( nbins, binWidth ) )
    #
    oFile.write( '\n\n' )
    oFile.write( '# cold losses - total: %i \n' % ( ApertureLossesColdSample1 ) )
    if ( plotCounts ):
        oFile.write( '# s_mean [m], counts [] \n' )
        for tmpBin,tmpVal in zip(bins_cold,hist_cold):
            oFile.write( '%13.4f %13i \n' % (tmpBin,tmpVal) )
    else:
        oFile.write( '# s_mean [m], pdf [m-1] \n' )
        for tmpBin,tmpVal in zip(bins_cold,hist_cold):
            oFile.write( '%13.4f %13.4E \n' % (tmpBin,tmpVal) )
    #
    oFile.write( '\n\n' )
    oFile.write( '# warm losses - total: %i \n' % ( ApertureLossesWarmSample1 ) )
    if ( plotCounts ):
        oFile.write( '# s_mean [m], counts [] \n' )
        for tmpBin,tmpVal in zip(bins_warm,hist_warm):
            oFile.write( '%13.4f %13i \n' % (tmpBin,tmpVal) )
    else:
        oFile.write( '# s_mean [m], pdf [m-1] \n' )
        for tmpBin,tmpVal in zip(bins_warm,hist_warm):
            oFile.write( '%13.4f %13.4E \n' % (tmpBin,tmpVal) )
    #
    oFile.write( '\n\n' )
    oFile.write( '# collimator losses - total: %i \n' % ( totcollosses1 ) )
    if ( plotCounts ):
        oFile.write( '# s_min [m], s_max[m], counts [] \n' )
        for k in range(len(simu1co)):
            if ( simu1co[k][3] > 0.0 ):
                oFile.write( '%13.4f %13.4f %13i \n' % (
                    simu1co[k][1]-0.5*simu1co[k][3], simu1co[k][1]+0.5*simu1co[k][3], simu1co[k][2] ) )
    else:
        oFile.write( '# s_min [m], s_max[m], pdf [m-1] \n' )
        for k in range(len(simu1co)):
            if ( simu1co[k][3] > 0.0 ):
                oFile.write( '%13.4f %13.4f %13.4E \n' % (
                    simu1co[k][1]-0.5*simu1co[k][3], simu1co[k][1]+0.5*simu1co[k][3], simu1co[k][2]/(normfac*simu1co[k][3] ) ) )
    #
    oFile.close()

# done
print('...done.')
