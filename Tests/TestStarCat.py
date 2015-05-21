#this is just a method to let us see the contents of catalogs written to
#disk without leaving this notebook
def printCatalogToScreen(catName):
    data = open(catName,'r')
    lines = data.readlines()
    data.close()
    for line in lines:
        print line

import numpy
from lsst.sims.catalogs.measures.instance import InstanceCatalog

class simpleStarCatalog(InstanceCatalog):
    column_outputs = ['raJ2000', 'decJ2000', 'sedFilename']
    
    transformations = {'raJ2000':numpy.degrees, 'decJ2000':numpy.degrees}

from lsst.sims.catalogs.generation.db import ObservationMetaData

myObsMetadata = ObservationMetaData(unrefractedRA=45.0, unrefractedDec=-10.0,
                                    boundType='circle', boundLength=0.02)

from lsst.sims.catUtils.baseCatalogModels import StarObj

starTableConnection = StarObj()

myCatalog = simpleStarCatalog(starTableConnection,
                              obs_metadata = myObsMetadata)

myCatalog.write_catalog('test_catalog.txt')

readCatalog = open('test_catalog.txt', 'r').readlines()

printCatalogToScreen('test_catalog.txt')