import numpy as np
import scipy
import matplotlib.pyplot as plt

AperatureType = 'circle'
AperatureRadius = 0.5 #FoV in degrees
#AperatureRadius = 1.75 #LSST's Actual FoV in degrees
SearchRegionRA = (-30.0,-20.0)
SearchRegionDec = (-30.0,-20.0)
SearchAirmass = (1.0,1.5)
DesiredFilter = None

import eups
import os
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catalogs.generation.db import ObservationMetaData
#help(ObservationMetaDataGenerator)

opsimdb = os.path.join(eups.productDir('sims_data'),'OpSimData','enigma_1189_sqlite.db')
gen = ObservationMetaDataGenerator(driver='sqlite', database=opsimdb)
SimObData = gen.getObservationMetaData(boundType=AperatureType, boundLength=AperatureRadius, 
                                       fieldRA=SearchRegionRA, fieldDec=SearchRegionDec,
                                       airmass=SearchAirmass, telescopeFilter=DesiredFilter)
NumOfObservations = len(SimObData)
#print SimsObData[0].__dict__

from prettytable import PrettyTable
UniquePointings = list({(o.unrefractedRA,o.unrefractedDec) for o in SimObData})
NumOfPointings = len(UniquePointings)
print 'Number of Unique Pointings:', NumOfPointings,'\n'

table2 = PrettyTable(["Pointing RA","Pointing Dec"])
for x in UniquePointings:
    table2.add_row([x[0], x[1]])
print table2

ObMetaData = [[] for _ in xrange(NumOfPointings)]
for i in xrange(NumOfPointings):
    for o in SimObData:
        if UniquePointings[i][0] == o.unrefractedRA and UniquePointings[i][1] == o.unrefractedDec:
            ObMetaData[i].append(o)
            
##########

from lsst.sims.catUtils.baseCatalogModels import GalaxyTileObj
galaxyDB = GalaxyTileObj()
#galaxyDB.show_mapped_columns()

from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.photUtils import PhotometryGalaxies, VariabilityGalaxies

class FastVarAgnCat(InstanceCatalog, PhotometryGalaxies, VariabilityGalaxies):
    
    cannot_be_null = ['uAgn'] #Locating only the AGNs in the FoV.
    
    column_outputs = ['AGNID', 'redshift', 'raJ2000', 'decJ2000', 'ObsAgn', 'sigma_ObsAgn']
    
    transformations = {'raJ2000':numpy.degrees, 'decJ2000':numpy.degrees}
    
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
    
    @compound('sedFilenameBulge', 'sedFilenameDisk')
    def get_QuickSED(self):
        ra = self.column_by_name('raJ2000') #Finding how many objects are in the column.
        names = []
        for r in ra:
            names.append('None') #Tricking the catalog into thinking these galaxies don't have bulges or disks.
        return np.array([names, names])
        
def ObsHeader(IC, file_handle):
    ObsHeaderTransformations = {'Unrefracted_RA':numpy.degrees, 'Unrefracted_Dec':numpy.degrees,
                               'Opsim_moonra':numpy.degrees, 'Opsim_moondec':numpy.degrees,
                               'Opsim_rotskypos':numpy.degrees, 'Opsim_rottelpos':numpy.degrees,
                               'Opsim_sunalt':numpy.degrees, 'Opsim_moonalt':numpy.degrees,
                               'Opsim_dist2moon':numpy.degrees, 'Opsim_altitude':numpy.degrees,
                               'Opsim_azimuth':numpy.degrees
                               }
    md = IC.obs_metadata.phoSimMetaData
    for k in md:
        if k in ObsHeaderTransformations.keys():
            file_handle.write(str(k)+" "+str(ObsHeaderTransformations[k](md[k][0]))+"\n")
        else:
            file_handle.write(str(k)+" "+str(md[k][0])+"\n")

###########

import time
PrimaryDir = 'TestCats/QuickSurvey_R%.2f' %(AperatureRadius)
CatDir = 'Catalogs'
if not os.path.exists(PrimaryDir+'/'+CatDir):
        os.makedirs(PrimaryDir+'/'+CatDir)
start = time.clock()
for Pointing in ObMetaData:
    PointingDir = 'RA[%.2f]_DEC[%.2f]' %(Pointing[0].unrefractedRA, Pointing[0].unrefractedDec)
    WorkingDir = PrimaryDir+'/'+CatDir+'/'+PointingDir
    if not os.path.exists(WorkingDir):
        os.makedirs(WorkingDir)
    for Observation in Pointing:
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
end = time.clock()
print end - start

