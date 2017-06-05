#!/sw/local/anaconda2/bin/python
import os,sys,re
from sys import exit as q
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

def get_latlon_index(lon0,lat0,lon2,lat2):
	return np.argmin(np.abs(lon2[0,:]-lon0)),np.argmin(np.abs(lat2[:,0]-lat0))

def write_river_input(riveroutfile,lon2,lat2,time_counter,runoffGrid,temperatureGrid,salinityGrid,depthGrid, bathy):

	# write temporary netCDF:
	ncfile=Dataset(riveroutfile,'w')
	# create the x and y dimensions.
	ncfile.createDimension('x',len(lon2[0,:]))
	ncfile.createDimension('y',len(lat2[:,0]))
	ncfile.createDimension('time_counter',None) # None defines unlimited dimension

	x = ncfile.createVariable('x','f8',('x'))
	x.standard_name = "longitude" ;
	x.long_name = "Longitude East" ;
	x.units = "degrees_east" ;
	x.valid_max = str(np.max(lon2[0,:])) ;
	x.valid_min = str(np.min(lon2[0,:])) ;

	y = ncfile.createVariable('y','f8',('y'))
	y.standard_name = "longitude" ;
	y.long_name = "Longitude East" ;
	y.units = "degrees_north" ;
	y.valid_max = str(np.max(lat2[:,0])) ;
	y.valid_min = str(np.min(lat2[:,0])) ;	

	t = ncfile.createVariable('time_counter','f8',('time_counter'))
	t.standard_name = "time"
	
	rorunoff = ncfile.createVariable('rorunoff','f8',('time_counter','y','x'))
	rorunoff.axis = "XYT" ;
	rorunoff.coordinates = "lon lat  t" ;

	rosaline = ncfile.createVariable('rosaline','f8',('time_counter','y','x'))
	rosaline.axis = "XYT" ;
	rosaline.coordinates = "lon lat  t" ;	

	rotemper = ncfile.createVariable('rotemper','f8',('time_counter','y','x'))
	rotemper.axis = "XYT" ;
	rotemper.coordinates = "lon lat  t" ;	
	
	rodepth = ncfile.createVariable('rodepth','f8',('y','x'))
	rodepth.axis = "XYT" ;
	rodepth.coordinates = "lon lat  t" ;	

	bathymetry = ncfile.createVariable('bathymetry','f8',('y','x'))
	bathymetry.standard_name = "depth in meters" ;

	# write data to variable.
	x[:]=lon2[0,:][:]
	y[:]=lat2[:,0][:]
	t[:]=time_counter[:]
	rorunoff[:]=np.transpose(runoffGrid)[:]
	rosaline[:] =np.transpose(temperatureGrid)[:]
	rotemper[:] = np.transpose(salinityGrid)[:]
	rodepth[:] = np.transpose(depthGrid)[:]
	bathymetry[:] = bathy[:]

def test_river_nc(fname):
	ncfile=Dataset(fname,'r')
	x = ncfile.variables['x']
	y = ncfile.variables['y']
	depth =  ncfile.variables['rodepth']
	r = ncfile.variables['rorunoff']

	x2,y2 = np.meshgrid(x,y)
	print x2.shape,y2.shape,r.shape

	plt.figure(2)
	levels=[0,10]
	plt.pcolor(x2,y2,np.array(r[10,:,:]))
	plt.contour(x2,y2,np.array(depth))
	plt.colorbar()
	plt.show()

def read_grid_nc(fname):
	ncfile=Dataset(fname,'r')
	return ncfile.variables['nav_lon'],ncfile.variables['nav_lat'],ncfile.variables['Bathymetry']

def read_params_from_namelist_cfg(fname):
	f=open(fname,'r')
	lines = f.readlines()
	for line in lines:
		if  'nn_itend' in line:
			mt = re.findall(r'.*=\s+(\d+)',line)
		elif 'rn_rdt ' in  line:
			mr = re.findall(r'.*=\s+(\d+)',line)
	if mt and mr:
		return int(mt[0]),int(mr[0])
	else:
		print "params not found: mt, mr = ", mt, mr
		q()

def main():

	# input files:
	gridfile = '/home/momo/NEMOGCM/TOOLS/VENTUS_TOOLS/create_bathy/bathy_meter.nc'
	namelistfile = '/home/momo/NEMOGCM/CONFIG/MY_GYRE_ncbathy/EXP00/namelist_cfg'
	
	# output file:
	riveroutfile = '/home/momo/NEMOGCM/TOOLS/VENTUS_TOOLS/create_rivers/adriatic_rivers.nc'

	# read grid file:
	lon2,lat2,bathy = read_grid_nc(gridfile)

	# read simulation duration from namelist:
	nt,rn_rdt = read_params_from_namelist_cfg(namelistfile)
	nHourlySteps = nt*rn_rdt/3600
	tmp=np.ones(nHourlySteps)

	# ----------------------------------------set river data dictionary:
	rivers = {\
	'Po1':{'lon':12.6665756,'lat':44.9863457,'depth':1,'temperature':10.0*tmp,'salinity':2.0*tmp,'runoff':1.0e-3*tmp},\
	'Po2':{'lon':12.7197013,'lat':44.8968036,'depth':1,'temperature':10.0*tmp,'salinity':2.0*tmp,'runoff':1.0e-3*tmp},\
	'Po3':{'lon':12.5708825,'lat':44.8525362,'depth':1,'temperature':10.0*tmp,'salinity':2.0*tmp,'runoff':1.0e-3*tmp},\
	'Po4':{'lon':12.5534698,'lat':45.0257579,'depth':1,'temperature':10.0*tmp,'salinity':2.0*tmp,'runoff':1.0e-3*tmp},\
	}
	#-----------------------------------------
	# loop over river dictionary and put river data on ocean grid:
	runoffGrid = np.zeros([lon2.shape[1],lon2.shape[0],nHourlySteps])
	temperatureGrid = np.zeros([lon2.shape[1],lon2.shape[0],nHourlySteps])
	salinityGrid = np.zeros([lon2.shape[1],lon2.shape[0],nHourlySteps])
	depthGrid = -np.ones([lon2.shape[1],lon2.shape[0]])

	# set time_counter:
	time_counter = [i*rn_rdt for i in range(nHourlySteps)]

	for river in rivers:
		runoffGrid[get_latlon_index(rivers[river]['lon'],rivers[river]['lat'],lon2,lat2)] = rivers[river]['runoff']
		temperatureGrid[get_latlon_index(rivers[river]['lon'],rivers[river]['lat'],lon2,lat2)] = rivers[river]['temperature']
		salinityGrid[get_latlon_index(rivers[river]['lon'],rivers[river]['lat'],lon2,lat2)] = rivers[river]['salinity']
		depthGrid[get_latlon_index(rivers[river]['lon'],rivers[river]['lat'],lon2,lat2)] = rivers[river]['depth']
		
	write_river_input(riveroutfile,lon2,lat2, time_counter,runoffGrid,temperatureGrid,salinityGrid,depthGrid,bathy)

	# test_river_nc(riveroutfile)



if __name__=="__main__":
	main()



