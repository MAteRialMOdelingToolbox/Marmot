# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.14-1 replay file
# Internal Version: 2014_06_04-21.37.49 134264
# Run by c8441146 on Thu May  4 08:54:13 2017
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=194.5166015625, 
    height=287.337493896484)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
openMdb(
    pathName='/media/c8441146/MagdalenaBackUp/Dropbox (Green Energy Center)/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax.cae')
#: The model database "/media/c8441146/MagdalenaBackUp/Dropbox (Green Energy Center)/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p = mdb.models['Model-1'].parts['Part-1']
p.deleteMesh()
p = mdb.models['Model-1'].parts['Part-1']
p.deleteSeeds()
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#4 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=14, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
p.generateMesh()
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, 
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
o1 = session.openOdb(
    name='/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          11
#: Number of Steps:              2
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S33'), )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(
    plotState=CONTOURS_ON_DEF)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='initialStep', frame=0)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
o1 = session.openOdb(
    name='/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#* OdbError: Cannot open file 
#* /home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax.odb. 
#* ***ERROR: The Abaqus database file is corrupt. If this file was transferred 
#* from another machine using FTP or equivalent, ensure that the file was 
#* copied using the binary mode instead of ASCII mode.
o1 = session.openOdb(
    name='/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          11
#: Number of Steps:              2
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='initialStep', frame=0)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S11'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S22'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S33'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S11'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S22'), )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S33'), )
o1 = session.openOdb(
    name='/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax_nonlocal.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax_nonlocal.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          11
#: Number of Steps:              2
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_UNDEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S11'), )
session.viewports['Viewport: 1'].view.setValues(nearPlane=189.505, 
    farPlane=293.817, width=125.959, height=68.0743, viewOffsetX=3.86529, 
    viewOffsetY=4.51965)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='initialStep', frame=0)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S22'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S33'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U1'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].view.setValues(nearPlane=193.358, 
    farPlane=289.963, width=74.4864, height=40.2561, viewOffsetX=6.24153, 
    viewOffsetY=15.1664)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT, 
    'Mises'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S11'), )
session.viewports['Viewport: 1'].view.setValues(nearPlane=186.818, 
    farPlane=296.504, width=168.127, height=90.864, viewOffsetX=24.2313, 
    viewOffsetY=8.4238)
session.viewports['Viewport: 1'].view.setValues(nearPlane=186.601, 
    farPlane=296.721, width=167.932, height=90.7585, viewOffsetX=42.5479, 
    viewOffsetY=-6.75524)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S22'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S33'), )
session.viewports['Viewport: 1'].view.setValues(nearPlane=184.572, 
    farPlane=298.75, width=199.228, height=107.673, viewOffsetX=46.4585, 
    viewOffsetY=-5.07522)
session.viewports['Viewport: 1'].view.setValues(nearPlane=184.3, 
    farPlane=299.022, width=198.934, height=107.514, viewOffsetX=29.4313, 
    viewOffsetY=-9.94089)
session.viewports['Viewport: 1'].view.setValues(nearPlane=186.379, 
    farPlane=296.943, width=167.732, height=90.6508, viewOffsetX=23.5913, 
    viewOffsetY=-4.59949)
o7 = session.odbs['/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=o7)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_UNDEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=189.321, 
    farPlane=294, width=128.405, height=69.3965, viewOffsetX=2.64357, 
    viewOffsetY=0.293658)
session.viewports['Viewport: 1'].view.setValues(nearPlane=189.495, 
    farPlane=293.826, width=128.523, height=69.4603, viewOffsetX=25.6739, 
    viewOffsetY=0.884243)
session.viewports['Viewport: 1'].view.setValues(nearPlane=187.27, 
    farPlane=296.052, width=161.86, height=87.4772, viewOffsetX=27.5899, 
    viewOffsetY=0.679152)
session.viewports['Viewport: 1'].view.setValues(nearPlane=187.046, 
    farPlane=296.276, width=161.666, height=87.3725, viewOffsetX=47.9405, 
    viewOffsetY=-3.03437)
o1 = session.openOdb(
    name='/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S11'), )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
o7 = session.odbs['/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax_nonlocal.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=o7)
#: Warning: The output database '/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax_nonlocal.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax_nonlocal.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax_nonlocal.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          11
#: Number of Steps:              2
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_UNDEF, ))
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='initialStep', frame=0)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S22'), )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S11'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S33'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=0 )
#: Warning: The output database '/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax_nonlocal.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
session.viewports['Viewport: 1'].view.setValues(nearPlane=224.882, 
    farPlane=257.44, width=210.728, height=113.888, viewOffsetX=27.7614, 
    viewOffsetY=-4.7704)
o1 = session.openOdb(
    name='/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_UNDEF, ))
o1 = session.openOdb(
    name='/home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax_nonlocal.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/c8441146/Dropbox/PHD/ABAqus/2017-03-28_TriaxDavid/jobBiax_nonlocal.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          11
#: Number of Steps:              2
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_UNDEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='initialStep', frame=0)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT, 
    'Mises'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S33'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S22'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 
    'S11'), )
