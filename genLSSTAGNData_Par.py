import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys


###############################################################################
#                        Set Up Command Line Args                             #
###############################################################################
import argparse

def coords(s):
    try:
        x, y = map(float, s.split(','))
        return x, y
    except:
        raise argparse.ArgumentTypeError("Coordinate Tuples must be X1,Y1 X2,Y2 etc.")

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('-R', '--Radius', help='Aperature Radius (in degrees; LSST = 1.75)', type=float, nargs=1, required=True)
parser.add_argument('-RA', '--SearchRA', help='RA Tuple of the Search Region', type=coords, nargs=1, required=True)
parser.add_argument('-Dec', '--SearchDec', help='Dec Tuple of the Search Region', type=coords, nargs=1, required=True)
parser.add_argument('-AM', '--Airmass', help='Airmass Range Tuple', type=coords, nargs=1, default=[(1.0,1.5)])
parser.add_argument('-N', '--NumOfCores', help='Number of CPU Cores to Use', type=int, nargs=1, default=[1])
args = vars(parser.parse_args())


###############################################################################
#                       Set Up Ad-Hoc Timing Class                            #
###############################################################################
import time
class Timer(object):
    def __init__(self, verbose=False):
        self.verbose = verbose

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000  # millisecs
        if self.verbose:
            print 'Elapsed Time: %f ms' % self.msecs


###############################################################################
#                      Set Up Observation Variables                           #
###############################################################################
global AperatureRadius

AperatureType = 'circle'
#AperatureRadius = float(raw_input("Please insert a value for the Aperature Radius: ")) #FoV in degrees
#AperatureRadius = 1.75 #LSST's Actual FoV in degrees
AperatureRadius = args['Radius'][0]
#SearchRegionRA = (20.0,30.0)
#SearchRegionDec = (-30.0,-20.0)
SearchRegionRA = args['SearchRA'][0]
SearchRegionDec = args['SearchDec'][0]
SearchAirmass = args['Airmass'][0]
DesiredFilter = None


###############################################################################
#                     Search DB for the Observations                          #
###############################################################################
import eups
import os
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catalogs.generation.db import ObservationMetaData
#help(ObservationMetaDataGenerator)

with Timer() as t:
    opsimdb = os.path.join(eups.productDir('sims_data'),'OpSimData','enigma_1189_sqlite.db')
    gen = ObservationMetaDataGenerator(driver='sqlite', database=opsimdb)
    SimObData = gen.getObservationMetaData(boundType=AperatureType, boundLength=AperatureRadius, 
                                       fieldRA=SearchRegionRA, fieldDec=SearchRegionDec,
                                       airmass=SearchAirmass, telescopeFilter=DesiredFilter)
print 'Time to search the region defined by,', (SearchRegionRA,SearchRegionDec),', for observations: %s seconds' %t.secs
#print SimsObData[0].__dict__


###############################################################################
#                  Output Basic Table of Returned Query                       #
###############################################################################
from prettytable import PrettyTable

NumOfObservations = len(SimObData)
print 'Number of Observations:', NumOfObservations
UniquePointings = list({(o.unrefractedRA,o.unrefractedDec) for o in SimObData})
NumOfPointings = len(UniquePointings)
print 'Number of Unique Pointings:', NumOfPointings,'\n'
table2 = PrettyTable(["Pointing RA","Pointing Dec"])
for x in UniquePointings:
    table2.add_row([x[0], x[1]])
print table2


###############################################################################
#              Restructure Returned Query to Correct Format                   #
###############################################################################
ObMetaData = [[] for _ in xrange(NumOfPointings)]
for i in xrange(NumOfPointings):
    for o in SimObData:
        if UniquePointings[i][0] == o.unrefractedRA and UniquePointings[i][1] == o.unrefractedDec:
            ObMetaData[i].append(o)


###############################################################################
#                   Connect to the UoW Galaxy Database                        #
###############################################################################
#from lsst.sims.catUtils.baseCatalogModels import GalaxyTileObj
#galaxyDB = GalaxyTileObj()
#galaxyDB.show_mapped_columns()


###############################################################################
#                 Define Instance Catalog Class for AGNs                      #
###############################################################################
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.photUtils import PhotometryGalaxies, VariabilityGalaxies

class FastVarAgnCat(InstanceCatalog, PhotometryGalaxies, VariabilityGalaxies):
    cannot_be_null = ['uAgn'] #Locating only the AGNs in the FoV.
    column_outputs = ['AGNID', 'redshift', 'raJ2000', 'decJ2000', 'ObsAgn', 'sigma_ObsAgn']
    transformations = {'raJ2000':np.degrees, 'decJ2000':np.degrees}
    
    #Only Append the Column of AGN Magnitudes in the Observed Bandpass
    @compound('ObsAgn', 'sigma_ObsAgn')
    def get_PhotometryCols(self):
        mag = None
        sigma = None
        if self.obs_metadata.bandpass is not None:
            if not hasattr(self.obs_metadata.bandpass, '__iter__'):
                Filter = self.obs_metadata.bandpass
                mag = self.column_by_name(Filter+'Agn')
                sigma = self.column_by_name('sigma_'+Filter+'Agn')
        return np.array([mag, sigma])
    
    #Don't Calculate Photometry for the Galaxy Bulge or Disk
    @compound('sedFilenameBulge', 'sedFilenameDisk')
    def get_QuickSED(self):
        ra = self.column_by_name('raJ2000') #Finding how many objects are in the column.
        names = []
        for r in ra:
            names.append('None') #Tricking the catalog into thinking these galaxies don't have bulges or disks.
        return np.array([names, names])


###############################################################################
#                  Define the Observation Header Writer                       #
###############################################################################
def ObsHeader(IC, file_handle):
    ObsHeaderTransformations = {'Unrefracted_RA':np.degrees, 'Unrefracted_Dec':np.degrees,
                               'Opsim_moonra':np.degrees, 'Opsim_moondec':np.degrees,
                               'Opsim_rotskypos':np.degrees, 'Opsim_rottelpos':np.degrees,
                               'Opsim_sunalt':np.degrees, 'Opsim_moonalt':np.degrees,
                               'Opsim_dist2moon':np.degrees, 'Opsim_altitude':np.degrees,
                               'Opsim_azimuth':np.degrees
                               }
    md = IC.obs_metadata.phoSimMetaData
    for k in md:
        if k in ObsHeaderTransformations.keys():
            file_handle.write(str(k)+" "+str(ObsHeaderTransformations[k](md[k][0]))+"\n")
        else:
            file_handle.write(str(k)+" "+str(md[k][0])+"\n")
            

###############################################################################
#            Define the Catalog Writer for Parallel Processes                 #
###############################################################################
from lsst.sims.catUtils.baseCatalogModels import GalaxyTileObj
def CreateCatParallel(Observation):
    galaxyDB = GalaxyTileObj()
    #galaxyDB.show_mapped_columns()
    PointingDir = 'RA[%.2f]_DEC[%.2f]' %(Observation.unrefractedRA, Observation.unrefractedDec)
    WorkingDir = PrimaryDir+'/'+CatDir+'/'+PointingDir
    if not os.path.exists(WorkingDir):
        os.makedirs(WorkingDir)
    CatFileName = 'LSSTAGN_%.5f.txt' %(Observation.mjd)
    if os.path.isfile(WorkingDir+'/'+CatFileName) == False:
        SQLRules = '(sedname_agn IS NOT NULL) AND (magnorm_agn < 24.0)'
        variableAgn = FastVarAgnCat(galaxyDB, obs_metadata=Observation, constraint=SQLRules)
        # Writing the Observational Meta-Data Header
        with open(WorkingDir+'/'+CatFileName, "w") as fh:
            ObsHeader(variableAgn, fh)
            fh.write("\n")
            fh.close()
        variableAgn.write_catalog(WorkingDir+'/'+CatFileName, write_mode='a')


###############################################################################
#              Create the Catalogs with Parallel Processing                   #
###############################################################################
from joblib import Parallel, delayed  
import multiprocessing

global PrimaryDir, CatDir, WorkingDir
#NumUsedCores = int(raw_input("Please select the number of cores to use: "))
NumUsedCores = args['NumOfCores'][0]
PrimaryDir = 'TestCats/TimeTestSurvey_R%.2f_C%i' %(AperatureRadius, NumUsedCores)
CatDir = 'Catalogs'

if not os.path.exists(PrimaryDir+'/'+CatDir):
        os.makedirs(PrimaryDir+'/'+CatDir)
with Timer() as t:
    for Pointing in ObMetaData[0:1]:
        PointingDir = 'RA[%.2f]_DEC[%.2f]' %(Pointing[0].unrefractedRA, Pointing[0].unrefractedDec)
        WorkingDir = PrimaryDir+'/'+CatDir+'/'+PointingDir
        if not os.path.exists(WorkingDir):
            os.makedirs(WorkingDir)
        with Timer() as t1:
            #NumUsedCores = multiprocessing.cpu_count()
            #Parallel(n_jobs=NumUsedCores, verbose=5)(delayed(CreateCatParallel)(o) for o in Pointing)
            pool = multiprocessing.Pool(processes=NumUsedCores)
            pool.map(CreateCatParallel, Pointing)
            pool.close()
        print 'Done with one Pointing! It took %.2f seconds!' %t1.secs
print "Creating all of the catalogs took %s seconds." %t.secs
