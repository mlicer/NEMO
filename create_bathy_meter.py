#!/sw/local/anaconda2/bin/python
from sys import exit as q
import os, numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter as gsf
from scipy.interpolate import RectBivariateSpline as rbs

def write_bathymetry(fname,lon1,lat1,bathy):
	print "---",fname
	print lat1
	print "---"
	[lon2,lat2] = np.meshgrid(lon1,lat1)
	dlon = lon2[1,1]-lon2[1,0]
	dlat = lat2[1,1]-lat2[0,1]
	# write temporary netCDF:
	ncfile=Dataset(fname,'w')
	# create the x and y dimensions.
	ncfile.createDimension('nav_lon',len(lon1))
	ncfile.createDimension('nav_lat',len(lat1))
	nclon = ncfile.createVariable('nav_lon','f8',('nav_lat','nav_lon'))
	nclon.resolution = str(dlon)
	nclat = ncfile.createVariable('nav_lat','f8',('nav_lat','nav_lon'))
	nclat.resolution = str(dlat)

	data = ncfile.createVariable('Bathymetry','f8',('nav_lat','nav_lon'))
	data1 = ncfile.createVariable('Bathy_level','f8',('nav_lat','nav_lon'))


	# write data to variable.
	nclon[:]=lon2[:]
	nclat[:]=lat2[:]
	ncdlon = dlon
	ncdlat = dlat
	data[:] = bathy[:]
	data[:] = bathy[:]
	ncfile.close()

def generate_synthetic_bathy(lonmin,lonmax,latmin,latmax,nLon,nLat):
	# generate grid:
	lon1 = np.linspace(lonmin,lonmax,nLon)
	lat1 = np.linspace(latmin,latmax,nLat)

	dlon = lon1[1]-lon1[0]
	dlat = lat1[1]-lat1[0]
	print "nLon,nLat: ",nLon,nLat
	print "lonmin,lonmax,latmin,latmax: ", lonmin,lonmax,latmin,latmax
	print "dLon = ",dlon
	print "dLat = ",dlat
	meanlon = np.mean(lon1)
	meanlat = np.mean(lat1)
	
	[lon2,lat2] = np.meshgrid(lon1,lat1)

	bathy =-80. + 500 * (1 - np.exp(-((lon2-meanlon+3)**2+(lat2-meanlat-1)**2)/10)- np.exp(-((lon2-meanlon-3.)**2+(lat2-meanlat+1.)**4)/10))
	bathy = np.where(bathy<0.,0,bathy)
	return np.array(lon1),np.array(lat1),np.array(bathy)

def read_bathymetry_from_nc(fname):
	ncfile=Dataset(fname,'r')
	lon2 = ncfile.variables['east_c']
	lat2 = ncfile.variables['north_c']
	bathy =  ncfile.variables['h']
	bathy = np.where(bathy<0.,0.,bathy)

	return np.array(lon2[0,:]),np.array(lat2[:,0]),np.array(bathy)

def main():
	lonmin = 12.
	lonmax = 25.
	latmin = 39.
	latmax = 46.
	nLon = 130
	nLat = 90
	
	# generate artificial analytic bathymetry:
	# lon1,lat1,bathy = generate_synthetic_bathy(lonmin,lonmax,latmin,latmax,nLon,nLat)

	# generate real bathymetry:
	fname='/home/momo_si/bathymetries/sbadripom_grid.nc'
	lon1, lat1, bathy = read_bathymetry_from_nc(fname)
	lon2,lat2 = np.meshgrid(lon1,lat1)
	
	resampleToSyntheticGrid=False
	if resampleToSyntheticGrid:
		# interpolate real bathy to synthetic grid:
		## read synth grid:
		lonsy1,latsy1,bathysy = generate_synthetic_bathy(lonmin,lonmax,latmin,latmax,nLon,nLat)
		# create interpolant:
		print "interpolating..."
		print bathy.shape
		f = rbs(lon1,lat1,bathy.transpose())
		# evaluate interpolant on synthetic grid:
		bathy = f(lonsy1,latsy1).transpose()
		print "interpolation done..."
		
		bathy = np.where(bathy<0.0,-999.,bathy)
		lon1 = lonsy1
		lat1 = latsy1
		lon2,lat2 = np.meshgrid(lon1,lat1)

	# smooth bathymetry:
	smoothBathy=False
	if smoothBathy:
		smoothed_bathy = gsf(bathy,2)
		smoothed_bathy = \
		np.where(smoothed_bathy<0.0,-999.,smoothed_bathy)
		

	plt.figure(1,figsize=(12,36))
	plt.subplot(2,1,1)
	plt.pcolor(lon2,lat2,bathy)
	plt.colorbar()
	# plt.contour(lon2,lat2,bathy,[0,1])

	if smoothBathy:
		plt.subplot(2,1,2)
		plt.pcolor(lon2,lat2,smoothed_bathy)
		plt.colorbar()
		# plt.contour(lon2,lat2,smoothed_bathy,[0,0.01])
	
		write_bathymetry('bathy_meter.nc',lon1,lat1,smoothed_bathy)
		write_bathymetry('bathy_level.nc',lon1,lat1,smoothed_bathy)
	else:
		write_bathymetry('bathy_meter.nc',lon1,lat1,bathy)
		write_bathymetry('bathy_level.nc',lon1,lat1,bathy)		

	plt.show()

	# create rivers on this grid:
	os.system('/home/momo/NEMOGCM/TOOLS/VENTUS_TOOLS/create_rivers/create_rivers.py')
	os.system('cp /home/momo/NEMOGCM/TOOLS/VENTUS_TOOLS/create_rivers/adriatic_rivers.nc /home/momo/NEMOGCM/CONFIG/MY_GYRE_ADR/EXP00/')
	os.system('ncdump -h '+'bathy_meter.nc')

if __name__=="__main__":
	main()



