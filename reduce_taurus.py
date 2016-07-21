from astropy.io import fits as afits
import pyfits as pf
import numpy as np
import matplotlib as mpl
#import matplotlib.pyplot as plt
import pickle
import subprocess
import math
import scipy
from scipy.optimize import curve_fit
from scipy.ndimage.interpolation import shift
from scipy.ndimage.fourier import fourier_shift as fshift
from os.path import exists
import threading
import random
import gc
from multiprocessing.dummy import Pool as ThreadPool

#sky()
#create_reduced_images()
#centroid()

#pathname to ds9 in taurus
ds9_taurus = '/Users/gduchene/Applications/Ureka/bin/ds9 '

#change the following when dealing with different date
#or different telescope
#or changes to img files e.g. new stars added
directory = 'ShaneAO'
date = 'Jun15'
arr_date_ShaneAO = ['Jun15', 'Nov15', 'Mar16']
struct_ShaneAO_total_sets = {'Jun15':28, 'Nov15':25, 'Mar16':25}
total_sets = struct_ShaneAO_total_sets[date]

#Showing date when script is run
print 'Working with', date, 'files'

#Set parameters of region where light hits telescope
arr_region_bounds = open(directory + '/' + date + '/telescope_region.txt', 'rb').read().splitlines()
region_x_min_index = int(arr_region_bounds[0])
region_x_max_index = int(arr_region_bounds[1])
region_y_min_index = int(arr_region_bounds[2])
region_y_max_index = int(arr_region_bounds[3])

print region_x_min_index
print region_x_max_index
print region_y_min_index
print region_y_max_index

#Correlation matrix parameters:
max_numb_imgs = 5000

#------------------
#General guidelines
#------------------
#make flat.txt and dark.txt files for dark and flat file numbers
#create_txt() to create txt files with info of different stars
#Use check_img_center() and change the parameters where light hits telescope
#run_create_sky() creates the sky imgs for all sets.
#create_dark()
#print_filters to view available filters for all imgs
#sort_filters create pickle files
#create_flats() to median flats and create master flat
#cleanup_flat() to fix broken pixels
#run_create_reduced_imgs() to subtract imgs by sky img & divide imgs by flat
#run_center_imgs():
#------------------

def check_img_center():
        outputfilename = 'test_imgcenter.fits'
        img = pf.open('ShaneAO/Nov15/s0101.fits')[0].data
        img_center = img[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index]
        print 'img shape', img_center.shape
        hdu = pf.PrimaryHDU(img_center)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(outputfilename ,clobber=True)
        subprocess.call('/home/apps/ds9/ds9 ' + outputfilename, shell = True)

def create_txt(date = 'Mar16'):
	path_to_file = directory + '/' + date + '/'
	setnumb = raw_input('set number:')
	startnumb = raw_input('starting number:')
	endnumb = raw_input('ending number:')
	star_id = raw_input('Star ID:')

	thatfile = open(path_to_file + 'set_' + setnumb + '_info.txt', 'wb')
	thatfile.write(startnumb + '\n')
	thatfile.write(endnumb + '\n')
	thatfile.write(star_id + '\n')
	thatfile.close

def print_filters():
	arr_filternames = []
	for i in np.arange(1, 5000): #CHANGE ACCORDINGLY
		if i > 999:
			filename = directory + '/' + date + '/' + 's' + str(i) + '.fits'
			if exists(filename):
				fits = pf.open(filename)
			else:
				continue
		else:
			filename = directory + '/' + date + '/' + 's0' + str(i) + '.fits'
			if exists(filename):
				fits = pf.open(filename)
			else:
				continue
		hdr = fits[0].header
		try:
			filtername = hdr['filt1nam']
		except:
                        print "file", i, "does not have the tag 'FILT1NAME'"
			print hdr.keys()
			print len(hdr.keys())
		if filtername not in arr_filternames:
			print filtername
			arr_filternames.append(filtername)
        print 'filters are:', arr_filternames

def sort_filters():
        arr_ks_flat = []
        arr_brg_flat = []
        arr_j_flat = []
        arr_startend = open(directory + '/' + date + '/flat.txt', 'rb').read().splitlines()
        start = int(arr_startend[0])
        end = int(arr_startend[1])
        arr_flats = np.arange(start, end+1)
	for i in arr_flats:
		print i
		if i > 999:
			filename = directory + '/' + date + '/' + 's' + str(i) + '.fits'
			if exists(filename):
				fits = pf.open(filename)
			else:
				print filename, 'doesnt exist'
				continue
		else:
			filename = directory + '/' + date + '/' + 's0' + str(i) + '.fits'
			if exists(filename):
				fits = pf.open(filename)
			else:
				print filename, 'doesnt exist'
				continue

		img = fits[0].data
		hdr = fits[0].header
		filtername = hdr['filt1nam']
		print filtername
		if filtername[0] == 'K':
			arr_ks_flat.append(i)
		elif filtername[0] == 'B':
			arr_brg_flat.append(i)
		if filtername[0] == 'J':
			arr_j_flat.append(i)
	arr_ks_flat = np.array(arr_ks_flat)
        print arr_ks_flat
	arr_brg_flat = np.array(arr_brg_flat)
        print arr_brg_flat
	arr_j_flat = np.array(arr_j_flat)
        print arr_j_flat
	pickle.dump(arr_ks_flat, open(directory + '/' + date + '/' + 'arr_ks_flat.p', 'wb'))
	pickle.dump(arr_brg_flat, open(directory + '/' + date + '/' + 'arr_brg_flat.p', 'wb'))
	pickle.dump(arr_j_flat, open(directory + '/' + date + '/' + 'arr_j_flat.p', 'wb'))
                

def create_sky(setnumb1):
        #Get sky image by medianing all images of sky
        #setnumb1 = raw_input('Enter set number (1,2,3,4, etc.):')
	setnumb1 = str(int(setnumb1))
        with open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb') as f_temp:
		arr_startend1 = f_temp.read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
        arr_sky = np.arange(start1, end1+1)
        output = median_it(arr_sky)
	hdu = pf.PrimaryHDU(output)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto(directory + '/' + date + '/' + 'img_sky_' + str(setnumb1) +'.fits',clobber=True)

def create_dark():
        arr_startend = open(directory + '/' + date + '/dark.txt', 'rb').read().splitlines()
        start = int(arr_startend[0])
        end = int(arr_startend[1])
        arr_dark = np.arange(start, end+1)
        img_output = median_it(arr_dark)
        hdu = pf.PrimaryHDU(img_output)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(directory + '/' + date + '/' + 'img_dark.fits',clobber=True)

def median_it(arr_filenumbs):
        arr_mdn = []
        print arr_filenumbs
        for i in arr_filenumbs:
		print i
                if i > 999:
			filename = directory + '/' + date + '/' + 's' + str(i) + '.fits'
			if exists(filename):
				#with pf.open(filename, memmap = False) as fits:
				fits = pf.open(filename, memmap = False)
				img = fits[0].data
			else:
				print filename, 'doesnt exist'
		else:
			filename = directory + '/' + date + '/' + 's0' + str(i) + '.fits'
			if exists(filename):
				#with pf.open(filename, memmap = False) as fits:
				fits = pf.open(filename, memmap = False)
				img = fits[0].data
			else:
				print filename, 'doesnt exist'
                arr_mdn.append(img)
		del img
		del fits
        arr_output = np.median(np.array(arr_mdn), axis = 0)
        return arr_output

def norm_median_subdark(arr_filenumbs):

        arr_mdn = []
        print arr_filenumbs
        img_dark = pf.open(directory + '/' + date + '/' + 'img_dark.fits')[0].data
        for i in arr_filenumbs:
		print i
                if i > 999:
			filename = directory + '/' + date + '/' + 's' + str(i) + '.fits'
			if exists(filename):
				#with pf.open(filename, memmap = False) as fits:
				fits = pf.open(filename, memmap = False)
				img = fits[0].data
			else:
				print filename, 'doesnt exist'
                                continue
		else:
			filename = directory + '/' + date + '/' + 's0' + str(i) + '.fits'
			if exists(filename):
				#with pf.open(filename, memmap = False) as fits:
				fits = pf.open(filename, memmap = False)
				img = fits[0].data
			else:
				print filename, 'doesnt exist'
                                continue
                img -= img_dark
                img_small = img[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index]
                mdn_imgsmall = np.median(img_small)
                img /= mdn_imgsmall
                arr_mdn.append(img)
		del img
		del fits
        arr_output = np.median(np.array(arr_mdn), axis = 0)
        return arr_output


def create_flats():

        arr_ks_flat = pickle.load(open(directory + '/' + date + '/' + "arr_ks_flat.p", "rb" ))
	print 'arr_ks_flat', arr_ks_flat
        arr_brg_flat = pickle.load(open(directory + '/' + date + '/' + "arr_brg_flat.p", "rb" ))
	print 'arr_brg_flat', arr_brg_flat
        arr_j_flat = pickle.load(open(directory + '/' + date + '/' + "arr_j_flat.p", "rb" ))
	print 'arr_j_flat', arr_j_flat
        


        try:

                img_ks_flat = norm_median_subdark(arr_ks_flat)
                hdu = pf.PrimaryHDU(img_ks_flat)
                hdulist = pf.HDUList([hdu])
                hdulist.writeto(directory + '/' + date + '/' + 'img_flat_ks.fits',clobber=True)
                print 'img_flat_ks.fits created'
        except:
                print 'error with img_flat_ks'
        try:
                img_brg_flat = norm_median_subdark(arr_brg_flat)
                hdu = pf.PrimaryHDU(img_brg_flat)
                hdulist = pf.HDUList([hdu])
                hdulist.writeto(directory + '/' + date + '/' + 'img_flat_brg.fits',clobber=True)
                print 'img_flat_brg.fits created'                
        except:
                print 'error with img_flat_brg'                
        try:
                img_j_flat = norm_median_subdark(arr_j_flat)
                hdu = pf.PrimaryHDU(img_j_flat)
                hdulist = pf.HDUList([hdu])
                hdulist.writeto(directory + '/' + date + '/' + 'img_flat_j.fits',clobber=True)
                print 'img_flat_j.fits created'                
        except:
                print 'error with img_flat_j'                


def cleanup_flat():
        #### ___Change the next few lines as necessary!!!___
	threshold_min = 0.7
	threshold_max = 1.3
	array_flatname = ['img_flat_']#*3

        array_flatname[0] += 'ks'
	#array_flatname[1] += 'j'
        #array_flatname[2] += 'brg'
        
	array_flatname = [x+'.fits' for x in array_flatname]
        #-----------------------------
        
	for i in array_flatname:
		img = pf.open(directory + '/' + date + '/' + i, padding = True, ignore_missing_end=True)[0].data
		img_search_small = img[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index]
                
		#search for values where: threshold_min < value < threshold_max
		index_small = np.argwhere(np.logical_or(img_search_small < threshold_min, img_search_small > threshold_max))
		index_small[:,0] += region_y_min_index
		index_small[:,1] += region_x_min_index

                print 'index_small.shape', index_small.shape
		for j in np.arange(len(index_small)):
			median_val = 0
			counter = 0
			k = 1
                        while median_val < threshold_min or median_val > threshold_max:
				counter += 1
				y = index_small[j,0]
				x = index_small[j,1]
				img_median = img[y-k:y+k+1,x-k:x+k+1]
					
				#remove middle value
				img_median = np.delete(img_median, (((((2*k)+1))**2)-1)/2, None)
				median_val = np.median(img_median)
				if median_val < threshold_min or median_val > threshold_max:
					print median_val
				img[y,x] = median_val
                                '''
				if counter > 10:
					print 'median value', str(median_val)
					print 'x', str(x)
					print 'y', str(y)
					print 'j', str(j)
					#print 'k', str(k)
					print 'shape', str(img_median.shape)
					if counter > 20:
						a = raw_input('Too many values not within threshold')
                                '''
				k+=1

		hdu = pf.PrimaryHDU(img[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index])
		hdulist = pf.HDUList([hdu])
		hdulist.writeto(directory + '/' + date + '/' + 'cleaned_'+str(i), clobber=True)


def create_reduced_imgs(setnumb):
	setnumb = str(int(setnumb))
	#uncomment the following & comment above if only 1 set is needed.
        #setnumb = raw_input('Enter set number (1,2,3,4, etc.):')

        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
        arr_files = np.arange(start1, end1+1)
	img_sky = pf.open(directory + '/' + date + '/' + 'img_sky_'+ setnumb +'.fits')[0].data
	img_sky = img_sky[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index]

	
	if start1 > 999:
		hdr = pf.open(directory + '/' + date + '/' + 's' + str(start1) + '.fits')[0].header
	else:
		hdr = pf.open(directory + '/' + date + '/' + 's0' + str(start1) + '.fits')[0].header
	filtername = hdr['filt1nam']
        ''' ____DO NOT ERASE!!_____
	if filtername[0] == 'K':
		img_flat = pf.open(directory + '/' + date + '/' + 'cleaned_img_flat_ks.fits')[0].data
		print 'ks'
	elif filtername[0] == 'B':
		img_flat = pf.open(directory + '/' + date + '/' + 'cleaned_img_flat_brg.fits')[0].data
		print 'brg'
	elif filtername[0] == 'J':
		img_flat = pf.open(directory + '/' + date + '/' + 'cleaned_img_flat_j.fits')[0].data
		print 'j'
	else:
		print 'Filter name doesnt start with K, B or J'
        '''
        #####___Change as needed!!!___
        img_flat = pf.open(directory + '/' + date + '/' + 'cleaned_img_flat_ks.fits')[0].data
        #####_________________________

        
	for i in arr_files:
		print i
		print '_____'
		if i > 999:
                        filename = directory + '/' + date + '/' + 's' + str(i) + '.fits'
                        if exists(filename):
                                img = pf.open(filename)[0].data
                        else:
                                print filename, 'doesnt exist'
                                continue
		else:
                        filename = directory + '/' + date + '/' + 's0' + str(i) + '.fits'
                        if exists(filename):
                                img = pf.open(filename)[0].data
                        else:
                                print filename, 'doesnt exist'
                                continue
		img = img[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index]
		img2 = (img-img_sky)/(img_flat)
		hdu = pf.PrimaryHDU(img2)
		hdulist = pf.HDUList([hdu])
		if i > 999:
			hdulist.writeto(directory + '/' + date + '/' + 's' + str(i) + '_reduced_.fits',clobber=True)	
		else:
			hdulist.writeto(directory + '/' + date + '/' + 's0' + str(i) + '_reduced_.fits',clobber=True)

                        
def center_imgs(setnumb):
#create centroid of science img by stacking
        #setnumb = raw_input('Enter set number (1,2,3,4, etc.):')
        setnumb = str(int(setnumb))
	outputname = 'centroid_' + setnumb #raw_input('Enter centroid output name (Dont include .fits):')
        arr_startend = open(directory + '/' + date + '/set_'+ setnumb +'_info.txt', 'rb').read().splitlines()
        start = int(arr_startend[0])
        end = int(arr_startend[1])
	arr_center = []

	for i in np.arange(start, end+1):
                i = 4010
		print '_________'
		print i
		if i > 999:
                        filename = directory + '/' + date + '/s' + str(i) + '_reduced_.fits'
                        if exists(filename):

                                img = pf.open(filename)[0].data
                        else:
                                print filename, 'doesnt exist'
                                continue
		else:
                        filename = directory + '/' + date + '/s0' + str(i) + '_reduced_.fits'
                        if exists(filename):
                                img = pf.open(filename)[0].data
                        else:
                                print filename, 'doesnt exist'
                                continue
		
                #print 'img.shape', img.shape
		index_max = np.argmax(img) ###find index of max value in img
		index_max = np.unravel_index(index_max,img.shape) ###change to 2d dimensional index
                max_val = img[index_max]
		counter = 0
		while img[index_max[0]+1,index_max[1]] < max_val/4. or img[index_max[0]-1,index_max[1]] < max_val/4. or img[index_max[0],index_max[1]-1] < max_val/4. or img[index_max[0],index_max[1]+1] < max_val/4.: ###going through outliers with high values and zeroing them if the surrounding points are much smaller
			img[index_max] = 0
			index_max = np.argmax(img)
			index_max = np.unravel_index(index_max,img.shape)
			print index_max
			counter += 1
			if counter > 20:
				break
		if counter > 20:
			continue
                y = index_max[0] #y index of max value in img
		x = index_max[1] #x index of max value in img
		j = 4 #radius of annulus for which to do center of mass calculation
                arr_chunk = []
                #taking x, y indexes and value of points in circle
                #appending to arr_chunk for center of mass calculations
                for x_circ in np.arange((2*j) + 1) - j:
                        for y_circ in np.arange((2*j) + 1) - j:
                                if distance(x_circ, 0, y_circ, 0) < (j+0.5):
                                        arr_chunk.append([y_circ, x_circ, img[y+y_circ, x+x_circ]])

                arr_y = []
                arr_x = []           
                arr_tot = []
                print 'maximum index y, x', y, ',', x
                for elem in arr_chunk:
                        y_index = elem[0]
                        x_index = elem[1]
                        value_elem = elem[2]
                        arr_y.append(y_index*value_elem)
                        arr_x.append(x_index*value_elem)
                        arr_tot.append(value_elem)
                tot_flux = sum(arr_tot)
                com_y = sum(arr_y)/tot_flux 
                com_x = sum(arr_x)/tot_flux 
                com_x += x #center of mass y index in img
                com_y += y #center of mass x index in img
                dec_com_x = com_x%1
                dec_com_y = com_y%1
                int_com_x = int(com_x)
                int_com_y = int(com_y)

                #fourier transform & shift
                img = np.fft.fft2(img)
                img = fshift(img, [-dec_com_y, -dec_com_x])
                img = np.fft.ifft2(img)
                img = np.real(img)
		k = 80
		output = img[int_com_y-k:int_com_y+k+1,int_com_x-k:int_com_x+k+1]
                if output.shape != (2*k + 1, 2*k +1):
                        print 'star too close to margin, skipping'
                        filename_check1 = directory + '/' + date + '/s' + str(i) + '_reduced_centered.fits'
                        filename_check2 = directory + '/' + date + '/s0' + str(i) + '_reduced_centered.fits'
                        if exists(filename_check1):
                                subprocess.call('rm -rf ' + filename_check1, shell = True) ###
                        if exists(filename_check2):
                                subprocess.call('rm -rf ' + filename_check2, shell = True) ###                                
                        continue
		print 'output.shape', output.shape
                #output.fill(0) ###
                #y_length = output.shape[0] ###
                #x_length = output.shape[1] ###
                #center_index = [(y_length-1)/2, (x_length-1)/2] ###
                #output[center_index[0], center_index[1]] = 0###
                arr_center.append(output)
                hdu = pf.PrimaryHDU(output)
                hdulist = pf.HDUList([hdu])
		if i > 999:
                        hdulist.writeto(directory + '/' + date + '/s' + str(i) + '_reduced_centered.fits', clobber=True)			
#                        subprocess.call('/home/apps/ds9/ds9 ' + directory + '/' + date + '/s' + str(i) + '_reduced_centered.fits', shell = True) ###
		else:
                        hdulist.writeto(directory + '/' + date + '/s0' + str(i) + '_reduced_centered.fits', clobber=True)
#                        subprocess.call('/home/apps/ds9/ds9 ' + directory + '/' + date + '/s0' + str(i) + '_reduced_centered.fits', shell = True) ###
	print 'arr_center shape:', str(np.array(arr_center).shape)
	img_mdn = np.median(np.array(arr_center), axis=0)
	print img_mdn.shape
	hdu = pf.PrimaryHDU(img_mdn)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto(directory + '/' + date + '/' + outputname + '.fits', clobber=True)
        #subprocess.call('/home/apps/ds9/ds9 ' + directory + '/' + date + '/' + outputname + '.fits', shell = True) ###

        
def run_create_sky():
	startset = 1
        arr_setnumbs = np.arange(startset, total_sets+1).astype(int)
        threadnumb = 6
	thepool = ThreadPool(threadnumb)
        arr_imgfinals = thepool.map(create_sky, arr_setnumbs)
        thepool.close()
        thepool.join()

        ''' FOR NO THREADS
	for setnumb in np.arange(startset, total_sets+1):
		print 'set number', setnumb
		create_sky(setnumb)
        '''
        
def run_create_reduced_imgs():

        startset = 1
        #'''
        arr_setnumbs = np.arange(startset, total_sets+1).astype(int)
        threadnumb = 6
	thepool = ThreadPool(threadnumb)
        arr_imgfinals = thepool.map(create_reduced_imgs, arr_setnumbs)
        thepool.close()
        thepool.join()
        #'''
        ''' FOR NO THREADS
	for setnumb in np.arange(startset, total_sets+1):
		print 'set number', setnumb
		create_reduced_imgs(setnumb)
        '''
def run_center_imgs():
	startset = 1
        '''
        arr_setnumbs = np.arange(startset, total_sets+1).astype(int)
        threadnumb = 6
	thepool = ThreadPool(threadnumb)
        arr_imgfinals = thepool.map(center_imgs, arr_setnumbs)
        thepool.close()
        thepool.join()
        '''
	for setnumb in np.arange(startset, total_sets+1):
		print 'set number', setnumb
		center_imgs(setnumb)
        
        
def that_histo():
	startset = 1
        threshold = 35000
	arr_histo = []
	for setnumb in np.arange(startset, total_sets+1):
		arr_startend1 = open(directory + '/' + date + '/set_'+ str(int(setnumb)) +'_info.txt', 'rb').read().splitlines()
		start1 = int(arr_startend1[0])
		end1 = int(arr_startend1[1])
		arr_files = np.arange(start1, end1+1)
		print 'set number', setnumb
		for i in arr_files:
			if i > 999:
				filename = directory + '/' + date + '/' + 's' + str(i) + '_reduced_centered.fits'
				if exists(filename):
					img = pf.open(filename)[0].data
				else:
					print filename, 'does not exist'
					continue
			else:
				filename = directory + '/' + date + '/' + 's0' + str(i) + '_reduced_centered.fits'
				if exists(filename):
					img = pf.open(filename)[0].data
				else:
					print filename, 'does not exist'
					continue

                        '''
                        if np.amax(img) > threshold:
                                print i
                                hdu = pf.PrimaryHDU(img)
                                hdulist = pf.HDUList([hdu])
                                hdulist.writeto(directory + '/' + date + '/test.fits', clobber=True)			
                                subprocess.call('/home/apps/ds9/ds9 ' + directory + '/' + date + '/test.fits', shell = True)
                        '''

                        arr_histo.append(np.amax(img))	
	plt.hist(np.array(arr_histo), bins = 20)
	plt.xlabel('Maximum pixel value')
	plt.show()
	



'''
def mock_binary():
	file1 = raw_input("First centroid? (Don't include .fits):")
	file2 = raw_input("2nd centroid? (Don't include .fits):")
	arr_flux = (np.arange(20.)+1)/100. #try log space, should be better
	arr_rad = np.arange(20)+1
	arr_detect = np.zeros([arr_rad.size, arr_flux.size])
	arr_detect.fill(0)
	#arr_detect = pickle.load(open("mock_binary_detlim.p", "rb" ))
	for index_flux in np.arange(arr_flux.size):
		flux_ratio = arr_flux[index_flux]
		img2 = pf.open('centroid_'+file2+'.fits')[0].data
		img2 /= np.amax(img2)
		img2 *= flux_ratio
		for radius in arr_rad:
			img1 = pf.open('centroid_'+file1+'.fits')[0].data
			img1 /= np.amax(img1)
			img1[:,radius:] += img2[:,:-radius]
			output = img1
			print 'candidate?:', check_binary(output, radius)
			print 'flux ratio', flux_ratio
			print 'radius', radius
			#if flux_ratio > 0.6:
			#	a = raw_input('STOP')
			arr_detect[radius-1, index_flux] = check_binary(output, radius)
			if check_binary(output, radius) > 0.5:
				hdu = pf.PrimaryHDU(output)
				hdulist = pf.HDUList([hdu])
				hdulist.writeto('test.fits', clobber=True)
				print 'radius', str(radius)
				#print 'open -a SAOImage\ DS9 test.fits'
				subprocess.call('open -a SAOImage\ DS9 test.fits', shell = True)

			outputfilename = 'mock_radius'+str(radius)+'_' + 'ratio' + str(flux_ratio)
			hdu = pf.PrimaryHDU(output)
			hdulist = pf.HDUList([hdu])
			hdulist.writeto(outputfilename + '.fits', clobber=True)
			print 'radius', str(radius)
			print 'open -a SAOImage\ DS9 ' + outputfilename + '.fits'
			subprocess.call('open -a SAOImage\ DS9 ' + outputfilename + '.fits', shell = True)
			observable = raw_input('Visible binary? (y or n):')
			if observable == 'y' or observable == 'Y':
				observable = 1
			else:
				observable = 0
			print observable
			arr_detect[radius-1, index_flux] = observable
			print arr_detect
			pickle.dump(arr_detect, open('mock_binary_detlim.p', 'wb'))
	print arr_detect
'''

def psf_subtract_frame():
#perform psf subtraction frame by frame
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])

        #setnumb2 = raw_input('Enter 2nd set number (1,2,3,4, etc.):')
        #arr_startend2 = open('set_' + setnumb2 + '/set_'+ setnumb2 +'_info.txt', 'rb').read().splitlines()
        arr_substar = np.array([])
       
        arr_setnumb = np.array([])
        for setnumb2 in np.arange(1,total_sets+1):
                if setnumb2 != int(setnumb1):
                        print 'adding elements in set no.', setnumb2
                        arr_startend2 = open(directory + '/' + date + '/set_'+ str(setnumb2) +'_info.txt', 'rb').read().splitlines()
                        start2 = int(arr_startend2[0])
                        end2 = int(arr_startend2[1])
                        arr_imgnumbs = np.arange(start2, end2+1)
                        arr_substar = np.append( arr_substar,  arr_imgnumbs)
                        arr_temp = np.zeros(arr_imgnumbs.size)
                        arr_temp.fill(setnumb2)
                        arr_setnumb = np.append(arr_setnumb, arr_temp)

        #define  subtraction ratios to be iterated through
        arr_rad = np.arange(70)+1
        arr_subratio = (np.arange(20)+90)/100.
        #________________________________________________

        #define empty arrays, for psf-subtracted imgs for avging
        arr_img = []
        arr_mdn = np.array([])
        arr_5sd = np.array([])
        #________________________________________________

	for i in np.arange(start1, end1+1): #iterating through images in the 1st set '''start1'''
		print '_________'
		print i
                if i > 999:
                        filename = directory + '/' + date + '/s' + str(int(i)) + '_reduced_centered.fits'
                        if exists(filename):
                                img1 = pf.open(filename)[0].data
                        else:
                                continue
                else:
                        filename = directory + '/' + date + '/s0' + str(int(i)) + '_reduced_centered.fits'
                        if exists(filename):
                                img1 = pf.open(filename)[0].data
                        else:
                                continue
                counter2 = 0
                for j in arr_substar:
                        arr_output = []
                        arr_rss = np.array([])
                        j_index = np.argwhere(arr_substar == j)
                        j_index = j_index[0]
                        j_index = j_index[0]
                        if j > 999:
                                filename = directory + '/' + date + '/s' + str(int(j)) + '_reduced_centered.fits'
                                if exists(filename):
                                        img2 = pf.open(filename)[0].data
                                else:
                                        continue
                        else:
                                filename = directory + '/' + date  + '/s0' + str(int(j)) + '_reduced_centered.fits'
                                if exists(filename):
                                        img2 = pf.open(filename)[0].data
                                else:
                                        continue
                        for subratio in arr_subratio: #iterating through images in the 1st set
                                img1_temp = np.copy(img1)
                                img2_temp = np.copy(img2)
                                img2_temp /= np.amax(img2_temp) #divide by max. brightness
                                img1_temp /= np.amax(img1_temp) #divide by max. brightness
                                center_flux = np.amax(img1_temp)                                
                                img2_temp *= subratio #multiply mock by subtraction ratio
                                img1_temp -= img2_temp #psf subtract
                                arr_rss = np.append(arr_rss, np.sum(img1_temp**2.)) #add residual sum of squares to arr_rss.
                                arr_output.append(img1_temp) #append psf-subtracted img
                        index_min_rss = np.argmin(arr_rss) #index of minimum rss in array of psfsubtracted imgs of diff sub ratios
                        img_best_subtract = arr_output[index_min_rss]
                        #print 'best subtraction ratio', arr_subratio[index_min_rss]
                        if counter2 == 0:
                                rss_best = np.sum(img_best_subtract**2.)
                                img_best = img_best_subtract
                                #jbest = j
                                #print 'best rss', rss_best
                                #print 'best subtraction frame', jbest
                        else:
                                if np.sum(img_best_subtract**2.) < rss_best:
                                        rss_best = np.sum(img_best_subtract**2.)
                                        img_best = img_best_subtract
                                        #jbest = j
                                        #print 'best rss', rss_best
                                        #print 'best subtraction frame', jbest
                        counter2 += 1
                arr_img.append(img_best)
        print 'size of arr_img', len(arr_img)
        img_mdn = np.median(np.array(arr_img), axis = 0)

        for radius in arr_rad:
                struct = check_ring(img_mdn, radius)
                arr_mdn = np.append(arr_mdn, struct['mdn_ring'])
                arr_5sd = np.append(arr_5sd, 5*struct['sd_ring'])

        hdu = pf.PrimaryHDU(img_mdn)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(directory + '/' + date + '/' + 'set' + setnumb1 + '-' + 'others' + '.fits', clobber=True)
        
        pickle.dump(arr_5sd, open(directory + '/' + date + '/' + 'set' + setnumb1 + '-' + 'others' + "_5sdplot.p", 'wb'))
        pickle.dump(arr_mdn, open(directory + '/' + date + '/' + 'set' + setnumb1 + '-' + 'others' + "_mdnplot.p", 'wb'))
        #print 'img_final.size', img_final.size
        #print arr_detect
        #a = raw_input('stopped')

def outdated_Loci():
        #define parameters for loci
        radius_sub = 1 #subtraction radius
        radius_op = 3 #optimization radius                                                                                                             

        #ask for user input regarding which set to take as target psf
        #reads txt file in directory to know which files to use
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])

        arr_substar = np.array([]) #array for adding ref psf filenames
	arr_setnumb = np.array([]) #array for adding target psf filenames
	for setnumb2 in np.arange(1,total_sets+1):
                if setnumb2 != int(setnumb1):
                        print 'adding elements in set no.', setnumb2
                        arr_startend2 = open(directory + '/' + date + '/' + 'set_' + str(setnumb2) +'_info.txt', 'rb').read().splitlines()
                        start2 = int(arr_startend2[0])
                        end2 = int(arr_startend2[1])
                        arr_imgnumbs = np.arange(start2, end2+1)
                        arr_substar = np.append( arr_substar,  arr_imgnumbs) #array of filenames for ref psf
			arr_temp = np.zeros(arr_imgnumbs.size)
                        arr_temp.fill(setnumb2)
                        arr_setnumb = np.append(arr_setnumb, arr_temp) #set numbers

        ##### FOR TESTING #####
	arr_substar = arr_substar[0:20]
        print arr_substar
        arr_setnumb = arr_setnumb[0:20]
        print arr_setnumb
        arr_refpsf = np.arange(start1, start1 + 10)
        ####################### 
        
        #define subtraction ratios to be iterated through                                                                                                       
        arr_rad = np.arange(1, 55, radius_sub)
        arr_rad = arr_rad.astype(int)
        arr_radcheck = np.arange(55)+1
        arr_subratio = (np.arange(20)+90)/100.

        #define empty arrays, for psf-subtracted imgs for avging
        arr_img = []
        arr_mdn = np.array([])
        arr_5sd = np.array([])

        #arr_refpsf = np.arange(start1, end1+1)
        arr_imgfinals = []

        
        threadlist = []
        threadcount = 0
        for i in arr_refpsf: #iterating through images in the 1st set
                print '_________'
                print i
                if i > 999:
                        filename = directory + '/' + date  + '/s' + str(int(i)) + '_reduced_centered.fits'
                        if exists(filename):
                                img1 = pf.open(filename)[0].data
                                img1 /= np.amax(img1)
                        else:
                                continue
                else:
                        filename =  directory + '/' + date + '/s0' + str(int(i)) + '_reduced_centered.fits'
                        if exists(filename):
                                img1 = pf.open(filename)[0].data
                                img1 /= np.amax(img1)
                        else:
                                continue
                img_final = np.zeros(img1.shape)
                for radi in arr_rad:
                        print 'working on radius', radi
                        arr_radsub = np.arange(radi, radi+radius_sub)
                        #print arr_radsub
                        arr_radop = np.arange(radi, radi+radius_op)
                        #print arr_radop
                        img1op = np.array([])
                        arr_indexop = []
                        for rad_op in arr_radop: #creating target optimization annulus
                                arr_ringop, indexop = return_ring(img1, rad_op)
                                for element in indexop:
                                        arr_indexop.append(element)
                                img1op = np.append(img1op, arr_ringop)

                        img1sub = np.array([])
                        arr_indexsub = []
                        for rad_sub in arr_radsub: #creating target subtraction annulus
                                arr_ringsub, indexsub = return_ring(img1, rad_sub)
                                for element in indexsub:
                                        arr_indexsub.append(element)
                                img1sub = np.append(img1sub, arr_ringsub)
                        #arr_refsub = np.array([])
                        arr_refop = np.array([])
                        arr_j = []
                        counter2 = 0
                        for j in arr_substar: #loop through ref imgs
                                arr_output = []
                                arr_rss = np.array([])
                                j_index = np.argwhere(arr_substar == j)
                                j_index = j_index[0]
                                j_index = j_index[0]
                                if j > 999:
                                        filename =  directory + '/' + date  + '/s' + str(int(j)) + '_reduced_centered.fits'
                                        if exists(filename):
                                                img2 = pf.open(filename)[0].data
                                                img2 /= np.amax(img2)
                                        else:
                                                continue
                                else:
                                        filename =  directory + '/' + date  + '/s0' + str(int(j)) + '_reduced_centered.fits'
                                        if exists(filename):
                                                img2 = pf.open(filename)[0].data
                                                img2 /= np.amax(img2)
                                        else:
                                                continue
                                arr_j.append(int(j))

                                img2op = np.array([])
                                for rad_op in arr_radop: #creating reference optimization anunulus
                                        img2op = np.append(img2op, return_ring(img2, rad_op)[0])
                                img2op_temp = img2op.reshape(-1,1)
				
				if counter2 == 0:
                                        arr_refop = img2op_temp
                                        #arr_refsub = img2sub_temp
				else:
			                arr_refop = np.concatenate((arr_refop, img2op_temp), axis = 1)
                                        #arr_refsub = np.concatenate((arr_refsub, img2sub_temp), axis = 1)
                                counter2 += 1
                        Y = img1op
                        print 'Y.shape', Y.shape
                        X = arr_refop
                        print 'X.shape', X.shape
                        alpha = np.dot(np.transpose(X), X)
                        beta = np.dot(np.transpose(X), Y)
                        coeffs = np.dot(np.linalg.inv(alpha), beta)
                        print 'sum of coeffs', np.sum(coeffs)
                        #print 'coeffs.shape', coeffs.shape
                        img_fill_annulus = np.zeros(img1.shape)
                        #print 'arr_j', arr_j
                        for index_imgs in np.arange(len(arr_j)):
                                j = arr_j[index_imgs]
                                j_index = np.argwhere(arr_substar == j)
                                j_index = j_index[0]
                                j_index = j_index[0]
                                if j > 999:
                                        filename =  directory + '/' + date + '/s' + str(int(j)) + '_reduced_centered.fits'
                                        img_temp = pf.open(filename)[0].data
                                else:
                                        filename =  directory + '/' + date + '/s0' + str(int(j)) + '_reduced_centered.fits'
                                        img_temp = pf.open(filename)[0].data
                                img_temp /= np.amax(img_temp)
                                img_fill_annulus += coeffs[index_imgs]*img_temp
                        for coord_annulus in arr_indexsub:
                                img_final[coord_annulus[0], coord_annulus[1]] = img_fill_annulus[coord_annulus[0], coord_annulus[1]]

                img_output = img1 - img_final
                arr_imgfinals.append(img_output)
                #hdu = pf.PrimaryHDU(img_output)
                #hdulist = pf.HDUList([hdu])
                #hdulist.writeto( directory + '/' + date + '/' 'test_loci.fits', clobber=True)
                #subprocess.call('/home/apps/ds9/ds9 test_loci.fits', shell = True)
        img_mdn = np.median(np.array(arr_imgfinals), axis = 0)
        print 'img_mdn.shape', img_mdn.shape
        for radius in arr_radcheck:
                struct = check_ring(img_mdn, radius)
                arr_mdn = np.append(arr_mdn, struct['mdn_ring'])
                arr_5sd = np.append(arr_5sd, 5*struct['sd_ring'])

        hdu = pf.PrimaryHDU(img_mdn)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(directory + '/' + date + 'set' + setnumb1 + '_loci' + '.fits', clobber=True)

        pickle.dump(arr_5sd, open(directory + '/' + date + '/' + 'set' + setnumb1 + '_loci' + "_5sdplot.p", 'wb'))
        pickle.dump(arr_mdn, open(directory + '/' + date + '/' + 'set' + setnumb1 + '_loci' + "_mdnplot.p", 'wb'))

def check_filter(filenumb, date):
        #input: file number
        #output: filter of file number
        if filenumb > 999:
                filename = directory + '/' + date + '/' + 's' + str(int(filenumb)) + '.fits'
        else:
                filename = directory + '/' + date + '/' + 's0' + str(int(filenumb)) + '.fits'
        img_hdr = pf.open(filename)[0].header
        filtername = img_hdr['filt1nam']
	#print filtername
        return filtername


def find_best_rad(mock_factor = 125, radius_mock = 23, halfboxcheck = 40):
	#ask for user input regarding which set to take as target psf
        #reads txt file in directory to know which files to use
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
	arr_targpsf = np.arange(start1, end1+1)

	max_radius_sub = 10
	min_radius_sub = 5
	arr_radius_sub = np.arange(min_radius_sub, max_radius_sub+1)

	max_radius_op = 25
	min_radius_op = 12
	arr_radius_op = np.arange(min_radius_op, max_radius_op+1)
        
	arr_rss = []
	arr_rad = []
	arr_radcheck = np.arange(1, 51)
	for radius_sub in arr_radius_sub:
		for radius_op in arr_radius_op:
			arr_mdn = []
			for i in arr_targpsf:
				filename_finalimg = directory + '/' + date + '/' + str(i) + '_thread_mock_rad' + str(radius_mock) + '_' + 'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_.p'
				if exists(filename_finalimg):
					#print 'rad_sub, rad_op', radius_sub, radius_op
					img_final = pickle.load(open(directory + '/' + date + '/' + str(i) + '_thread_mock_rad' + str(radius_mock) + '_' + 'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_.p', "rb"))
					img_final = img_final[80-halfboxcheck:80+halfboxcheck+1,80-halfboxcheck:80+halfboxcheck+1]
					#print 'img_final.shape', img_final.shape
					arr_rss.append(np.sum(img_final**2.))
					arr_rad.append([radius_sub, radius_op])
					arr_mdn.append(img_final)
				else:
					print filename_finalimg, 'doesnt exist'
			img_final = np.median(np.array(arr_mdn), axis = 0)
                        arr_5sd = np.array([])
			for radius in arr_radcheck:
				struct = check_ring(img_final, radius)
				arr_5sd = np.append(arr_5sd, 5*struct['sd_ring'])			
                        arr_5sd_mag = 2.5*np.log10(1./np.abs(arr_5sd))
                        if all(element > 6 for element in arr_5sd_mag[:10]):
                                plt.plot(arr_radcheck, arr_5sd_mag, 'o-', color = np.random.rand(3,1),  label = '5 s.d. value, ' + 'radop' + str(radius_op) + ', ' + 'radsub' + str(radius_sub))
			index_min = np.argmin(arr_rss)
			'''
			print 'len(arr_rss)', len(arr_rss)
			print 'for img', i
			print 'index_min', index_min
			print 'arr_rad[index_min]', arr_rad[index_min]
			print 'radius subtraction:', arr_rad[index_min][0]
			print 'radius optimization:', arr_rad[index_min][1]
			arr_rss = np.delete(arr_rss, index_min)
			index_min = np.argmin(arr_rss)
	
			print '_____________'
			print '2nd best parameters:'
			print 'len(arr_rss)', len(arr_rss)
			print 'for img', i
			print 'index_min', index_min
			print 'arr_rad[index_min]', arr_rad[index_min]
			print 'radius subtraction:', arr_rad[index_min][0]
			print 'radius optimization:', arr_rad[index_min][1]
			arr_rss = np.delete(arr_rss, index_min)
			index_min = np.argmin(arr_rss)
			print '_____________'
			print '3rd best parameters:'
			print 'len(arr_rss)', len(arr_rss)
			print 'for img', i
			print 'index_min', index_min
			print 'arr_rad[index_min]', arr_rad[index_min]
			print 'radius subtraction:', arr_rad[index_min][0]
			print 'radius optimization:', arr_rad[index_min][1]
			print '______________________________________________'
			print '______________________________________________'
			useless = raw_input('that pause though')
                        '''
        plt.legend(loc = 'upper right')
        plt.xlabel('radius (pixels)')
        plt.ylabel('magnitude')
        plt.gca().invert_yaxis()
        plt.show()

def rename_theta():
        arr_radiusmock = np.arange(15, 60+1, 15).astype(int)
        print 'arr_radiusmock', arr_radiusmock
        arr_mockfactor = np.arange(100, 200+1, 50).astype(int) #np.array([30, 500, 1000]) 
        print 'arr_mockfactor', arr_mockfactor
        arr_theta = np.array([267]).astype(int) #np.array([0]) #np.arange(0, 360, 60).astype(int)
        print 'arr_theta', arr_theta
        #Test parameters_____
        arr_radius_sub = [6]
        arr_radius_op = [19]

        
	#ask for user input regarding which set to take as target psf
        #reads txt file in directory to know which files to use        
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
	arr_targpsf = np.arange(start1, end1+1)

        theta_change = int(raw_input('Change theta to:'))
        
        arr_rss = []
	arr_rad = []
	arr_radcheck = np.arange(1, 61)
        for radius_sub in arr_radius_sub:
                for radius_op in arr_radius_op:
                        for radius_mock in arr_radiusmock:
                                for mock_factor in arr_mockfactor:
                                        for theta_mock in arr_theta:
                                                print 'theta_mock', theta_mock
                                                arr_mdn = []
                                                dx_mock = radius_mock*np.cos(theta_mock)
                                                dx_mock = int(round(dx_mock))
                                                dy_mock = radius_mock*np.sin(theta_mock)
                                                dy_mock = int(round(dy_mock))
                                                print 'dy_mock, dx_mock:', dy_mock, ',', dx_mock
                                                for i in arr_targpsf:
                                                        filename_output = directory + '/' + date + '/' + str(i) + '_thread_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_.p'
                                                        print filename_output
                                                        useless = raw_input('You better wait for input')

                                                        if exists(filename_output):
                                                                filename_change = directory + '/' + date + '/' + str(i) + '_thread_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_change)) + '_.p'
                                                                subprocess.call('mv -f ' + filename_output + ' ' + filename_change,shell = True)

                                                                

        
def plot_mocks():
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 45, 60]).astype(int)#np.arange(15, 60+1, 15).astype(int)
        #arr_radiusmock = arr_radiusmock[::-1]
        print 'arr_radiusmock', arr_radiusmock
        arr_mockfactor = np.array([30, 100, 150, 200, 500, 1000]).astype(int)
        print 'arr_mockfactor', arr_mockfactor
        arr_theta = np.array([0]).astype(int) #np.array([0, 198, 35, 233, 71, 269]).astype(int)
        #arr_theta = np.array([0]) #np.arange(0, 360, 60).astype(int)
        print 'arr_theta', arr_theta
        #Test parameters_____
        arr_radius_sub = [6]
        arr_radius_op = [19]

        halfboxcheck = 45
        
	#ask for user input regarding which set to take as target psf
        #reads txt file in directory to know which files to use        
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
	arr_targpsf = np.arange(start1, end1+1)
        print 'arr_targpsf', arr_targpsf
        
        arr_rss = []
	arr_rad = []
	arr_radcheck = np.arange(1, 61)
        for radius_sub in arr_radius_sub:
                for radius_op in arr_radius_op:
                        for mock_factor in arr_mockfactor:
                                for radius_mock in arr_radiusmock:
                                        for theta_mock in arr_theta:
                                                print 'theta_mock', theta_mock
                                                arr_mdn = []
                                                dx_mock = radius_mock*np.cos(theta_mock*np.pi/180)
                                                dx_mock = int(round(dx_mock))
                                                dy_mock = radius_mock*np.sin(theta_mock*np.pi/180)
                                                dy_mock = int(round(dy_mock))
                                                print 'dy_mock, dx_mock:', dy_mock, ',', dx_mock
                                                for i in arr_targpsf:
                                                        filename_output = directory + '/' + date + '/' + str(i) + '_thread_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_.p'
                                                        filename_output_fits = directory + '/' + date + '/' + str(i) + '_thread_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_.fits'
                                                        print filename_output
                                                        if exists(filename_output):
                                                                #print 'rad_sub, rad_op', radius_sub, radius_op
                                                                pickle.load(open(filename_arrpoolvars, 'rb'))
                                                                file_img = open(filename_output, 'rb')
                                                                img_final = pickle.load(file_img)
                                                                file_img.close()

                                                                        #print 'Pickle error for', filename_output
                                                                        #filename_output = directory + '/' + date + '/' + str(i) + '_thread_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_.fits'
                                                                        #if exists(filename_output):
                                                                        #        img_final = pf.open(filename_output)[0].data
                                                                        #else:
                                                                        #        continue
                                                                #img_final = img_final[80-halfboxcheck:80+halfboxcheck+1,80-halfboxcheck:80+halfboxcheck+1]
                                                                y_length = img_final.shape[0]
                                                                x_length = img_final.shape[1]
                                                                center_index_y, center_index_x = int((y_length-1)/2), int((x_length-1)/2)
                                                                index_mock_y, index_mock_x = center_index_y + dy_mock, center_index_x + dx_mock
                                                                #print 'img_final.shape', img_final.shape
                                                                #arr_rss.append(np.sum(img_final**2.))
                                                                #arr_rad.append([radius_sub, radius_op])
                                                                arr_mdn.append(img_final)
                                                        else:
                                                                print filename_output, 'doesnt exist'
                                                                continue
                                                if arr_mdn:
                                                        print 'numb of imgs for median:', len(arr_mdn)
                                                        img_final = np.median(np.array(arr_mdn), axis = 0)
                                                        flux_mock = img_final[index_mock_y, index_mock_x]
                                                        flux_mock_diff = 2.5*np.log10((1./mock_factor)/flux_mock)
                                                        print 'img_final.shape', img_final.shape
                                                        arr_5sd = np.array([])
                                                        for radius in arr_radcheck:
                                                                struct = check_ring_mock(img_final, radius, theta_mock)
                                                                arr_5sd = np.append(arr_5sd, 5*struct['sd_ring'])			
                                                        arr_5sd_mag = 2.5*np.log10(1./np.abs(arr_5sd))
                                                        #if all(element > 6 for element in arr_5sd_mag[:10]):
                                                        plt.plot(arr_radcheck, arr_5sd_mag, 'o-', color = np.random.rand(3,1),  label = '5 s.d. value, ' + 'theta: ' + str(theta_mock) + ', flux ratio:' + str(mock_factor) + ', mock radius: ' + str(radius_mock) + ', flux lost: ' + str(flux_mock_diff) + 'mag')
                                                        #index_min = np.argmin(arr_rss)
                                plt.legend(loc = 'upper right')
                                plt.xlabel('radius (pixels)')
                                plt.ylabel('magnitude')
                                plt.gca().invert_yaxis()
                                plt.show()
                                        #useless = raw_input('waiting for input')


                        
def create_final_klip():
        arr_radiusmock = np.arange(15, 60+1, 15).astype(int)
        print 'arr_radiusmock', arr_radiusmock
        arr_mockfactor = np.array([30, 500, 1000]) #np.arange(100, 200+1, 50).astype(int)
        print 'arr_mockfactor', arr_mockfactor
        arr_theta = np.array([0]) #np.arange(0, 360, 60).astype(int)
        print 'arr_theta', arr_theta
        '''
        #Test parameters_____
	threadnumb = 6
	thepool = ThreadPool(threadnumb)
        arr_imgfinals = thepool.map(klip, arr_pool_vars)
        thepool.close()
        thepool.join()
        '''

def test_theta(theta_mock):
        radius_mock = 40
        theta_mock *= (2*np.pi)/360
        img_mock = pf.open(directory + '/' + date + '/' + 'centroid_1.fits')[0].data
        img_mock /= np.amax(img_mock)
        
        dx_mock = radius_mock*np.cos(theta_mock)
        dx_mock = int(round(dx_mock))
        dy_mock = radius_mock*np.sin(theta_mock)
        dy_mock = int(round(dy_mock))
        
        img_mock = np.fft.fft2(img_mock) #Fourier shift to construct binary
        img_mock = fshift(img_mock, [dy_mock, dx_mock])
        img_mock = np.fft.ifft2(img_mock)
        img_mock = np.real(img_mock)

        hdu = pf.PrimaryHDU(img_mock)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto('test.fits', clobber = True)
        subprocess.call('/home/apps/ds9/ds9 ' + 'test.fits', shell = True)


def run_loci4():
        #------
        #Define parameters
        #------
        numb_imgs_keep = 200 #number of imgs to keep from correlation matrix
        fwhm_threshold = 10 #maximum FWHM to keep
	maxpix_threshold = 35000 #exclude imgs with max pixel value above this number
	minpix_threshold = 5000 #exclude imgs with min pixel value below this number


        
	#------
        #Define mock-up binary arrays
        #------
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 45, 60]).astype(int) #np.array([10, 15, 20, 25, 30, 45, 60]).astype(int)
        print 'arr_radiusmock', arr_radiusmock
        arr_mockfactor = np.array([30, 100, 150, 200, 500, 1000])#np.arange(100, 200+1, 50).astype(int)
        print 'arr_mockfactor', arr_mockfactor
        arr_theta = np.array([0]) #np.arange(0, 360, 60).astype(int)
        print 'arr_theta', arr_theta


        
        #------
        #Test parameters for subtraction/optimization radius for LOCI
        #------
        arr_radius_sub = [6]
        arr_radius_op = [19]

        

        #------
        #arr_date_ShaneAO: array of date names in telescope directory
        #index_date_targ: index of date in arr_date_ShaneAO
        #------
        print 'arr_date_ShaneAO', arr_date_ShaneAO
        print 'date of target psf', date
        index_date_targ = arr_date_ShaneAO.index(date) #index of date of target psf
        print 'index_date_targ', index_date_targ


        
        #------
        #Ask for user input regarding which set to take as target psf
        #Reads txt file in directory to know which files to use
        #------
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
        arr_targpsf = np.arange(start1, end1+1)
	print 'arr_targpsf', arr_targpsf
        

        
        #------
        #check for filter of 1st img in set
        #------
	filter_of_set = check_filter(start1, date) 
        print 'filter_of_set', filter_of_set



        #------
        # Add reference filenumbs from all dates
        #------
        arr_substar = []
        for index_date in np.arange(len(arr_date_ShaneAO)):
                date_ref = arr_date_ShaneAO[index_date]
                total_sets = struct_ShaneAO_total_sets[date_ref]
                print 'total_sets', total_sets
                arr_substar_temp = np.array([]) #array for adding ref psf filenames
                for setnumb2 in np.arange(1,total_sets+1):
                        if setnumb2 == int(setnumb1) and index_date_targ == index_date:
                                print 'index_date_targ', index_date_targ, 'i.e.', date_ref, 'set', setnumb2, 'skipped'
                                continue
                        else:
                                print 'adding elements in set no.', setnumb2, 'for date', date_ref
                                arr_startend2 = open(directory + '/' + date_ref + '/' + 'set_' + str(setnumb2) +'_info.txt', 'rb').read().splitlines()
                                start2 = int(arr_startend2[0])
                                end2 = int(arr_startend2[1])
                                arr_imgnumbs = np.arange(start2, end2+1)
                                arr_substar_temp = np.append(arr_substar_temp,  arr_imgnumbs) #array of filenames for ref psf
                                arr_temp = np.zeros(arr_imgnumbs.size)
                                arr_temp.fill(setnumb2)
                arr_substar_temp = arr_substar_temp.astype(int)
                arr_substar.append(arr_substar_temp)
        print 'arr_substar', arr_substar


        
        #------
        # load reference imgs one by one
        # Check filter, FWHM, max/min pixel values, img size
        # Include only the ones within parameters set above at beginning of function
        #------
        arr_j = []
        arr_img2 = []
        for index_date in np.arange(len(arr_date_ShaneAO)):
                date_ref = arr_date_ShaneAO[index_date]
                arr_j_temp = []
                arr_img2_temp = []
                arr_substar_temp = arr_substar[index_date]
                print 'arr_substar_temp', arr_substar_temp
                for j in arr_substar_temp: #loop through ref imgs
                        #print 'compiling img', j, 'for subtraction'
                        arr_output = []
                        arr_rss = np.array([])
                        j_index = np.argwhere(arr_substar_temp == j)
                        j_index = j_index[0]
                        j_index = j_index[0]
                        if j > 999:
                                filename =  directory + '/' + date_ref  + '/s' + str(int(j)) + '_reduced_centered.fits'
                                if exists(filename):
                                        img2 = pf.open(filename)[0].data
                                        maxpixval = np.amax(img2)
                                        img2 /= np.amax(img2)
                                else:
                                        print filename, 'doesnt exist'
                                        continue
                        else:
                                filename =  directory + '/' + date_ref  + '/s0' + str(int(j)) + '_reduced_centered.fits'
                                if exists(filename):
                                        img2 = pf.open(filename)[0].data
                                        maxpixval = np.amax(img2)
                                        img2 /= np.amax(img2)
                                else:
                                        print filename, 'doesnt exist'
                                        continue

                        img_filter = check_filter(j, date_ref)
                        
                        if img_filter != filter_of_set: #removing frames of different filter
                                #print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                #print 'img filter: ', img_filter, 'different from set filter:', filter_of_set
                                continue
                        elif maxpixval > maxpix_threshold or maxpixval < minpix_threshold: #removing frames with max pix values not within threshold
                                #print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                #print 'max pixel value of', maxpixval, 'not within allowed region'
                                continue

                        fwhm = get_fwhm(img2)
                        if not fwhm: #remove frames that gaussian func couldn't fit
                                #print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                continue 
                        elif fwhm > fwhm_threshold: #remove frames with large fwhm
                                print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                #print 'FWHM greater than', fwhm_threshold
                                continue
                        arr_j_temp.append(int(j))
                        img2 = img2.flatten() #change to 1D array
                        
                        if img2.size != 25921:
                                print 'img2.shape', img2.shape
                                print j, date_ref
                        arr_img2_temp.append(img2)
                arr_img2_temp = np.array(arr_img2_temp)
                print 'arr_img2_temp.shape', arr_img2_temp.shape
                arr_img2.append(arr_img2_temp)
                arr_j.append(arr_j_temp)


                
        #------        
        #Load array of filenumbs for ShaneAO, correlation matrix, and corresponding filenumbs in matrix
        #------
        filename_structShaneAOfilenumbs = 'struct_ShaneAO_filenumbs.p'
        struct_ShaneAO_filenumbs = pickle.load(open(filename_structShaneAOfilenumbs, 'rb'))
        filename_matrixcorr = 'matrix_corr.p'
        matrix_corr = pickle.load(open(filename_matrixcorr, 'rb'))
        filename_arrfilenumbscorr = 'arr_filenumbs_corr.p'
        arr_filenumbs_corr = pickle.load(open(filename_arrfilenumbscorr, 'rb')).astype(int)
        setnumb1 = int(setnumb1) #set numb for target psf img
        arr_filenumbs = struct_ShaneAO_filenumbs[index_date_targ]


        
        #------
        #Getting rid of frames of the same set in correlation matrix
        #------
        for element in arr_targpsf: 
                index_i_chunk = np.argwhere(arr_filenumbs == element)
                index_i_chunk = index_i_chunk[0][0]
                print 'index_i_chunk', index_i_chunk
                print 'check it', (index_date_targ*max_numb_imgs) + index_i_chunk
                matrix_corr[:, (index_date_targ*max_numb_imgs) + index_i_chunk] = 100



        #------
        #Fill an array with indexes of dates coresponding to correlation matrix
        #------
        numb_dates = len(struct_ShaneAO_filenumbs)
        arr_indexdate_corrs = np.zeros(numb_dates*max_numb_imgs)
        for index_indexdatecorrs in np.arange(numb_dates):
                arr_indexdate_corrs[index_indexdatecorrs*max_numb_imgs:(index_indexdatecorrs+1)*max_numb_imgs] = index_indexdatecorrs
        arr_indexdate_corrs = arr_indexdate_corrs.astype(int)



        #------
        # Create list of structures with info for LOCI function
        # Takes certain number of best imgs by referencing correlation matrix
        #------
        arr_pool_vars = []	
        for index_i in np.arange(arr_targpsf.size):
                arr_j_optimal = []
                arr_img2_optimal = []
                i = arr_targpsf[index_i]
                print i
                index_i = np.argwhere(arr_filenumbs == i)[0][0]
                print 'index_i', index_i
		arr_corrs = matrix_corr[(index_date_targ*max_numb_imgs) + index_i, :]
                index_sort_corr = np.argsort(arr_corrs)
                print 'index_sort_corr', index_sort_corr
                arr_indexdateoptimal = arr_indexdate_corrs[index_sort_corr]
                print 'arr_indexdateoptimal', arr_indexdateoptimal
                arr_filenumbs_optimal = arr_filenumbs_corr[index_sort_corr]
                arr_index_j_optimal = []
                for index_optimal in np.arange(arr_filenumbs_optimal.size).astype(int):
                        index_date_optimal = arr_indexdateoptimal[index_optimal]
                        filenumb_optimal = arr_filenumbs_optimal[index_optimal]
                        index_j_optimal = np.argwhere(arr_j[index_date_optimal] == filenumb_optimal)

                        if index_j_optimal:
                                arr_j_optimal.append(arr_j[index_date_optimal][index_j_optimal[0][0]])
                                arr_img2_optimal.append(arr_img2[index_date_optimal][index_j_optimal[0][0]])
                        if len(arr_j_optimal) >= numb_imgs_keep:
                                break
                arr_img2_optimal = np.array(arr_img2_optimal)
                arr_img2_optimal = np.transpose(arr_img2_optimal)
                print 'arr_img2_optimal.shape', arr_img2_optimal.shape
                print 'len(arr_j_optimal)', len(arr_j_optimal)

                for radius_sub in arr_radius_sub:
                        for radius_op in arr_radius_op:
                                for radius_mock in arr_radiusmock:
                                        for mock_factor in arr_mockfactor:
                                                for theta in arr_theta:
                                                        struct = {'i': i, 'arr_img2': arr_img2_optimal, 'arr_j': arr_j_optimal, 'radius_sub': radius_sub, 'radius_op':radius_op, 'radius_mock':radius_mock, 'mock_factor': mock_factor, 'theta': theta}
                                                        arr_pool_vars.append(struct)


                                                        
        #------
        #Save/load list of structures as needed
        #------
        filename_arrpoolvars = 'loci_arr_pool_vars_' + date + 'set' + str(int(setnumb1)) + '.p'        
        pickle.dump(arr_pool_vars, open(filename_arrpoolvars, 'wb'))
        #arr_pool_vars = pickle.load(open(filename_arrpoolvars, 'rb'))
        
        print '------------------------'
	print 'No of loops to run:', len(arr_pool_vars)
        print '------------------------'
        print 'Now starting threads...'

        

        #------
        #Run LOCI. Comment/uncomment sections for when thread is needed/not needed
        #------
        '''
	for elem in arr_pool_vars:
		loci2(elem)
        '''
	threadnumb = 6
	thepool = ThreadPool(threadnumb)
        arr_pool = thepool.map(loci2, arr_pool_vars)
        thepool.close()
        thepool.join()
        
        

        

def loci2(struct):
        #------
        #load variables from input structure, create arr_rad with appropriate startng subtraction radii
        #------
        i = struct['i']
	radius_sub = struct['radius_sub']
	radius_op = struct['radius_op']
	arr_img2 = struct['arr_img2']
        radius_mock = struct['radius_mock']
        mock_factor = struct['mock_factor']
        theta_mock = struct['theta']
        arr_rad = np.arange(1, 80-radius_op-1, radius_sub)
        arr_rad = arr_rad.astype(int)


        
        #------
        #define filenames for initial img with mock, psf img used for subtraction & final img
        #------
        filename_imgmock = directory + '/' + date + '/' + str(i) + '_mockimg_mock_rad' + str(radius_mock) + '_' + 'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_.fits'
        filename_imgsub = directory + '/' + date + '/' + str(i) + '_subimg_mock_rad' + str(radius_mock) + '_' + 'ratio' + str(mock_factor) +'_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_.fits'
        filename_output = directory + '/' + date + '/' + str(i) + '_thread_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_.fits'
        if exists(filename_imgmock) and exists(filename_imgsub) and exists(filename_output):
                print 'skipping', i, ', radius:', radius_mock, ', flux ratio:', mock_factor, ', theta:', theta_mock
                return 0
        print i, ', radius:', radius_mock, ', flux ratio:', mock_factor, ', theta:', theta_mock



        #------
        #Check pickle files and delete, since we're switching to fits files
        #------
        filename_imgmock_pickle = directory + '/' + date + '/' + str(i) + '_mockimg_mock_rad' + str(radius_mock) + '_' + 'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_.p'
        filename_imgsub_pickle = directory + '/' + date + '/' + str(i) + '_subimg_mock_rad' + str(radius_mock) + '_' + 'ratio' + str(mock_factor) +'_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_.p'
        filename_output_pickle = directory + '/' + date + '/' + str(i) + '_thread_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_.p'
        if exists(filename_imgmock_pickle):
                subprocess.call('rm -rf ' + filename_imgmock_pickle, shell = True)
                print 'deleted', filename_imgmock_pickle
        if exists(filename_imgsub_pickle):
                subprocess.call('rm -rf ' + filename_imgsub_pickle, shell = True)
                print 'deleted', filename_imgsub_pickle                
        if exists(filename_output_pickle):
                subprocess.call('rm -rf ' + filename_output_pickle, shell = True)
                print 'deleted', filename_output_pickle

                
        
        #----------
        #load mock up binary and divide by predetermined factor
        #***Change centroid filename if necessary***
        #----------
        img_mock = pf.open(directory + '/' + date + '/' + 'centroid_1.fits')[0].data
        img_mock /= np.amax(img_mock)
        if mock_factor != 0:
                img_mock /= float(mock_factor)
        else:
                img_mock *= float(mock_factor)


                
        #------
        #load target img, divide by max pixel value
        #------
        if i > 999:
                filename = directory + '/' + date  + '/s' + str(int(i)) + '_reduced_centered.fits'
                if exists(filename):
                        try:
                                img1 = pf.open(filename)[0].data
                                img1 /= np.amax(img1)
                        except:
                                print 'error opening file:', filename
                                return 0
                else:
			#print filename, 'doesnt exist'
                        return 0
        else:
                filename =  directory + '/' + date + '/s0' + str(int(i)) + '_reduced_centered.fits'
                if exists(filename):
                        try:
                                img1 = pf.open(filename)[0].data
                                img1 /= np.amax(img1)
                        except:
                                print 'error opening file:', filename                                
                                return 0
                else:
			#print filename, 'doesnt exist'
                        return 0



        #--------
        #Add mock binary at correct radius & angle, save img as fits file
        #--------
        if radius_mock != 0:
                dx_mock = radius_mock*np.cos(theta_mock*(2*np.pi)/360) #calculating y displacement from center
                dx_mock = int(round(dx_mock))
                dy_mock = radius_mock*np.sin(theta_mock*(2*np.pi)/360) #calculating x displacement from center
                dy_mock = int(round(dy_mock))
        
                img_mock = np.fft.fft2(img_mock) #Fourier shift to construct binary
                img_mock = fshift(img_mock, [dy_mock, dx_mock])
                img_mock = np.fft.ifft2(img_mock)
                img_mock = np.real(img_mock)

                img1 += img_mock
        hdu = pf.PrimaryHDU(img1)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_imgmock, clobber = True)



        #--------|
        #  LOCI  |
        #--------|
        img_final = np.zeros(img1.size)
        for radi in arr_rad:
                arr_radsub = np.arange(radi, radi+radius_sub) #array of starting radi of subtraction region
                arr_radop = np.arange(radi, radi+radius_op) #array of starting radi of optimization region
                img1op = np.array([])
                arr_indexop = []
                for rad_op in arr_radop: #creating target optimization annulus
			arr_ringop, indexop = return_ring(img1, rad_op)			
                        for element in indexop:
                                arr_indexop.append(np.ravel_multi_index(element, img1.shape))
                        img1op = np.append(img1op, arr_ringop)
                img1sub = np.array([])
                arr_indexsub = []
                for rad_sub in arr_radsub: #creating target subtraction annulus
                        arr_ringsub, indexsub = return_ring(img1, rad_sub)
                        for element in indexsub:
                                arr_indexsub.append(np.ravel_multi_index(element, img1.shape))
                        img1sub = np.append(img1sub, arr_ringsub)
		arr_refop = arr_img2[np.array(arr_indexop)]

                #Linear Least Squares fitting
                Y = img1op.reshape([-1,1])
                X = arr_refop
                alpha = np.dot(np.transpose(X), X)
                beta = np.dot(np.transpose(X), Y)
                coeffs = np.dot(np.linalg.inv(alpha), beta)
		coeffs = coeffs.reshape([-1,1])
		img_fill_annulus = np.dot(arr_img2, coeffs)
                for coord_annulus in arr_indexsub:
                        img_final[coord_annulus] = img_fill_annulus[coord_annulus]
	img_final = img_final.reshape(img1.shape)
        img_output = img1 - img_final


        
        #------
        #Save final img, and img used for subtraction as fits files
        #------
        hdu = pf.PrimaryHDU(img_final)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_imgsub, clobber = True)
        
        hdu = pf.PrimaryHDU(img_output)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_output, clobber = True)
        
	return 1


'''loci2 correlation  code: DO NOT ERASE
        
        arr_img2transpose = np.copy(arr_img2.T)
	for img2 in arr_img2transpose:
		for rad_ring in np.arange(1, max_rad_ring+1):
			arr_ringvals, arr_ringindex = return_ring(img1, rad_ring)
			for ringindex in arr_ringindex:
				ringindex = np.ravel_multi_index(ringindex, img1shape)
				img2[ringindex] = 0
				if counter == 0:
					img1_corr[ringindex] = 0
                                        count_points += 1
		counter += 1
		img_res = img1_corr - img2
		sd_imgres = np.std(img_res)
		arr_sd.append(sd_imgres)

	numb_imgs_keep = 100 ###___ VARIABLE___ ###
	index_sort_sd = np.argsort(arr_sd)
	index_sort_sd = index_sort_sd[:numb_imgs_keep]
	#print 'index_sort_sd.shape' + str(index_sort_sd.shape)

	arr_img2 = arr_img2[:, index_sort_sd]
	
	#print 'arr_img2.shape updated' +  str(arr_img2.shape)
	#useless = raw_input('what what')

        
        where0 = np.argwhere(arr_img2[:,0] == 0) #to check if masking worked. Not needed now.
        print 'count_points', count_points
	print 'numb of points where img2 is zero', where0.size
        hdu = pf.PrimaryHDU(img1)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto('test.fits', clobber = True)
        subprocess.call('/home/apps/ds9/ds9 ' + 'test.fits', shell = True)	
'''

def outdated_RunKlip(total_sets = total_sets):
        
        #max_area_klip = 300
        max_rad_img = 80
        rad_annulus = 5
        fwhm_threshold = 10
	maxpix_threshold = 35000
	minpix_threshold = 5000
        #numb_klip_modes = 200
        arr_klipmodes = np.arange(5, 50, 5)
        

        print 'total_sets', total_sets
        print 'arr_date_ShaneAO', arr_date_ShaneAO
        print 'date of target psf', date
        index_date_targ = arr_date_ShaneAO.index(date) #index of date of target psf
        print 'index_date_targ', index_date_targ
        
        #--------------------------------------------------------------
        #Ask for user input regarding which set to take as target psf
        #Reads txt file in directory to know which files to use
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
        arr_targpsf = np.arange(start1, end1+1)
        #-------------------------------------------------------------
        
	filter_of_set = check_filter(start1, date) #check for filter of 1st img in set

        arr_substar = []
        arr_setnumb = []
        print 'setnumb1', int(setnumb1)
        for index_date in np.arange(len(arr_date_ShaneAO)):
                date_ref = arr_date_ShaneAO[index_date]
                total_sets = struct_ShaneAO_total_sets[date_ref]
                print 'total_sets', total_sets
                arr_substar_temp = np.array([]) #array for adding ref psf filenames
                arr_setnumb_temp = np.array([]) #array for adding ref psf set numbers
                for setnumb2 in np.arange(1,total_sets+1):
                        if setnumb2 == int(setnumb1) and index_date_targ == index_date:
                                print 'index_date_targ', index_date_targ, 'i.e.', date_ref, 'set', setnumb2, 'skipped'
                                continue
                        else:
                                print 'adding elements in set no.', setnumb2, 'for date', date_ref
                                arr_startend2 = open(directory + '/' + date_ref + '/' + 'set_' + str(setnumb2) +'_info.txt', 'rb').read().splitlines()
                                start2 = int(arr_startend2[0])
                                end2 = int(arr_startend2[1])
                                arr_imgnumbs = np.arange(start2, end2+1)
                                arr_substar_temp = np.append(arr_substar_temp,  arr_imgnumbs) #array of filenames for ref psf
                                arr_temp = np.zeros(arr_imgnumbs.size)
                                arr_temp.fill(setnumb2)
                                arr_setnumb_temp = np.append(arr_setnumb_temp, arr_temp) #set numbers for ref psfs
                arr_substar_temp = arr_substar_temp.astype(int)
                print 'arr_setnumb_temp', arr_setnumb_temp
                arr_setnumb_temp = arr_setnumb_temp.astype(int)
                arr_substar.append(arr_substar_temp)
                arr_setnumb.append(arr_setnumb_temp)

        print 'arr_setnumb', arr_setnumb
        print 'arr_substar', arr_substar
        #define subtraction ratios to be iterated through
        arr_subratio = (np.arange(20)+90)/100.		

        #define empty arrays, for psf-subtracted imgs for avging
        arr_img = []
        arr_mdn = np.array([])
        arr_5sd = np.array([])

	print 'arr_targpsf', arr_targpsf
        arr_j = []
        arr_img2 = []

        for index_date in np.arange(len(arr_date_ShaneAO)):
                date_ref = arr_date_ShaneAO[index_date]
                #arr_j_temp = []
                #arr_img2_temp = []
                arr_substar_temp = arr_substar[index_date]
                print 'arr_substar_temp', arr_substar_temp
                for j in arr_substar_temp: #loop through ref imgs
                        #print 'compiling img', j, 'for subtraction'
                        arr_output = []
                        arr_rss = np.array([])
                        j_index = np.argwhere(arr_substar_temp == j)
                        j_index = j_index[0]
                        j_index = j_index[0]
                        if j > 999:
                                filename =  directory + '/' + date_ref  + '/s' + str(int(j)) + '_reduced_centered.fits'
                                if exists(filename):
                                        img2 = pf.open(filename)[0].data
                                        maxpixval = np.amax(img2)
                                        img2 /= maxpixval
                                else:
                                        print filename, 'doesnt exist'
                                        continue
                        else:
                                filename =  directory + '/' + date_ref  + '/s0' + str(int(j)) + '_reduced_centered.fits'
                                if exists(filename):
                                        img2 = pf.open(filename)[0].data
                                        maxpixval = np.amax(img2)
                                        img2 /= maxpixval
                                else:
                                        print filename, 'doesnt exist'
                                        continue

                        img_filter = check_filter(j, date_ref)
                        
                        if img_filter != filter_of_set: #removing frames of different filter
                                #print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                #print 'img filter: ', img_filter, 'different from set filter:', filter_of_set
                                continue
                        elif maxpixval > maxpix_threshold or maxpixval < minpix_threshold: #removing frames with max pix values not within threshold
                                #print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                #print 'max pixel value of', maxpixval, 'not within allowed region'
                                continue

                        fwhm = get_fwhm(img2)
                        if not fwhm: #remove frames that gaussian func couldn't fit
                                #print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                continue 
                        elif fwhm > fwhm_threshold: #remove frames with large fwhm
                                print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                #print 'FWHM greater than', fwhm_threshold
                                continue
                        
                        arr_j.append(int(j))

                        img2 = img2.flatten() #change to 1D array
                        
                        if img2.size != 25921:
                                print 'img2 shape isnt looking right'
                                print 'img2.shape', img2.shape
                                print j, date_ref
                        arr_img2.append(img2)
                        
        arr_img2 = np.array(arr_img2)
        print 'arr_img2.shape', arr_img2.shape

        for index_row in np.arange(len(arr_img2)):
                img2 = arr_img2[index_row]
                arr_img2[index_row] -=  np.mean(img2)

        arr_pool_vars = []
        for filenumb_img1 in arr_targpsf:
                for numb_klip_modes in arr_klipmodes:
                        struct = {'filenumb_img1': filenumb_img1, 'arr_img2': arr_img2, 'date_ref': date_ref, 'numb_klip_modes':numb_klip_modes, 'max_rad_img': max_rad_img, 'rad_annulus': rad_annulus}
                        arr_pool_vars.append(struct)

        print  '------------------------' 
        print 'Number of loops to run:', len(arr_pool_vars)
        print '------------------------'
        
        
        filename_arrpoolvars = 'klip_arr_pool_vars_' + date + 'set' + str(int(index_date_targ)) + '.p'
        pickle.dump(arr_pool_vars, open(filename_arrpoolvars, 'wb'))
        
        #arr_pool_vars = pickle.load(open(filename_arrpoolvars, 'rb'))        

        #print  '------------------------' 
        #print 'Number of loops to run:', len(arr_pool_vars)
        #print '------------------------'
        
        '''
        for elem in arr_pool_vars:
                klip(elem)
        '''        
        	
	threadnumb = 6
	thepool = ThreadPool(threadnumb)
        arr_imgfinals = thepool.map(klip, arr_pool_vars)
        thepool.close()
        thepool.join()



def run_klip2():
        #-------
        #Define parameters
        #------
        numb_imgs_keep = 600 #number of imgs to keep from correlation matrix
        numb_pixels_all_imgs = 25921 #number of pixels in all imgs
        max_rad_img = 80 # max radius in imgs
        rad_annulus = 5 # radius of annulus for KLIP search rings
        fwhm_threshold = 10 # upper limit on full-width-half-max
	maxpix_threshold = 35000 # exclude imgs with pixel values above this number
	minpix_threshold = 5000 # exclude imgs with pixel values below this number
        arr_klipmodes = np.arange(5, 50, 5).astype(int) # array with number of klip modes to try


        
	#------
        #Define mock-up binary arrays
        #------
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 45, 60]).astype(int)
        print 'arr_radiusmock', arr_radiusmock
        arr_mockfactor = np.array([30, 100, 150, 200, 500, 1000])
        print 'arr_mockfactor', arr_mockfactor
        arr_theta = np.array([0])
        print 'arr_theta', arr_theta



        #------
        #arr_date_ShaneAO: array of date names in telescope directory
        #index_date_targ: index of date in arr_date_ShaneAO
        #------
        print 'arr_date_ShaneAO', arr_date_ShaneAO
        print 'date of target psf', date
        index_date_targ = arr_date_ShaneAO.index(date) #index of date of target psf
        print 'index_date_targ', index_date_targ


        
        #------
        # Ask for user input regarding which set to take as target psf
        # Reads txt file in directory to know which files to use
        #------
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
        arr_targpsf = np.arange(start1, end1+1)
        print 'arr_targpsf', arr_targpsf
        

        
        #------
        # Add reference filenumbs from all dates
        #------
	filter_of_set = check_filter(start1, date) #check for filter of 1st img in set
        print 'filter_of_set', filter_of_set
        arr_substar = []
        print 'setnumb1', int(setnumb1)
        for index_date in np.arange(len(arr_date_ShaneAO)):
                date_ref = arr_date_ShaneAO[index_date]
                total_sets = struct_ShaneAO_total_sets[date_ref]
                print 'total_sets', total_sets
                arr_substar_temp = np.array([]) #array for adding ref psf filenames
                for setnumb2 in np.arange(1,total_sets+1):
                        if setnumb2 == int(setnumb1) and index_date_targ == index_date:
                                print 'index_date_targ', index_date_targ, 'i.e.', date_ref, 'set', setnumb2, 'skipped'
                                continue
                        else:
                                print 'adding elements in set no.', setnumb2, 'for date', date_ref
                                arr_startend2 = open(directory + '/' + date_ref + '/' + 'set_' + str(setnumb2) +'_info.txt', 'rb').read().splitlines()
                                start2 = int(arr_startend2[0])
                                end2 = int(arr_startend2[1])
                                arr_imgnumbs = np.arange(start2, end2+1)
                                arr_substar_temp = np.append(arr_substar_temp,  arr_imgnumbs) #array of filenames for ref psf
                arr_substar_temp = arr_substar_temp.astype(int)
                print 'arr_setnumb_temp', arr_setnumb_temp
                arr_substar.append(arr_substar_temp)
        print 'arr_substar', arr_substar



        #------
        # load reference imgs one by one
        # Check filter, FWHM, max/min pixel values, img size
        # Include only the ones within parameters set above at beginning of function
        #------
        arr_j = []
        arr_img2 = []
        for index_date in np.arange(len(arr_date_ShaneAO)):
                date_ref = arr_date_ShaneAO[index_date]
                arr_j_temp = []
                arr_img2_temp = []
                arr_substar_temp = arr_substar[index_date]
                print 'arr_substar_temp', arr_substar_temp
                for j in arr_substar_temp: #loop through ref imgs
                        #print 'compiling img', j, 'for subtraction'
                        arr_output = []
                        arr_rss = np.array([])
                        j_index = np.argwhere(arr_substar_temp == j)
                        j_index = j_index[0]
                        j_index = j_index[0]
                        if j > 999:
                                filename =  directory + '/' + date_ref  + '/s' + str(int(j)) + '_reduced_centered.fits'
                                if exists(filename):
                                        img2 = pf.open(filename)[0].data
                                        maxpixval = np.amax(img2)
                                        img2 /= np.amax(img2)
                                else:
                                        print filename, 'doesnt exist'
                                        continue
                        else:
                                filename =  directory + '/' + date_ref  + '/s0' + str(int(j)) + '_reduced_centered.fits'
                                if exists(filename):
                                        img2 = pf.open(filename)[0].data
                                        maxpixval = np.amax(img2)
                                        img2 /= np.amax(img2)
                                else:
                                        print filename, 'doesnt exist'
                                        continue
                                
                        img_filter = check_filter(j, date_ref)
                        if img_filter != filter_of_set: # removing frames of different filter
                                #print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                #print 'img filter: ', img_filter, 'different from set filter:', filter_of_set
                                continue
                        elif maxpixval > maxpix_threshold or maxpixval < minpix_threshold: # removing frames with max pix values not within threshold
                                #print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                #print 'max pixel value of', maxpixval, 'not within allowed region'
                                continue

                        fwhm = get_fwhm(img2)
                        if not fwhm: # remove frames that gaussian func couldn't fit
                                #print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                continue 
                        elif fwhm > fwhm_threshold: # remove frames with large fwhm
                                print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                #print 'FWHM greater than', fwhm_threshold
                                continue
                        
                        arr_j_temp.append(int(j))
                        img2 = img2.flatten() #change to 1D array

                        if img2.size != numb_pixels_all_imgs: # making sure img is the right size 
                                print 'img2.shape', img2.shape
                                print j, date_ref
                                useless = raw_input('IMG NOT CORRECT SIZE. CHECK FILE. Press any key to skip img and continue.')
                                continue
                        arr_img2_temp.append(img2)
                        #print 'len(arr_img2_temp)', len(arr_img2_temp)
                arr_img2_temp = np.array(arr_img2_temp)
                print 'arr_img2_temp.shape', arr_img2_temp.shape
                arr_img2.append(arr_img2_temp)
                arr_j.append(arr_j_temp)



        #------
        #Load array of filenumbs for ShaneAO, correlation matrix, and corresponding filenumbs in matrix
        #------
        filename_structShaneAOfilenumbs = 'struct_ShaneAO_filenumbs.p'
        struct_ShaneAO_filenumbs = pickle.load(open(filename_structShaneAOfilenumbs, 'rb'))
        filename_matrixcorr = 'matrix_corr.p'
        matrix_corr = pickle.load(open(filename_matrixcorr, 'rb')) 
        filename_arrfilenumbscorr = 'arr_filenumbs_corr.p'
        arr_filenumbs_corr = pickle.load(open(filename_arrfilenumbscorr, 'rb')).astype(int) #file numbers corresponding to position in correlation matrix
        setnumb1 = int(setnumb1) #set numb for target psf img
        arr_filenumbs = struct_ShaneAO_filenumbs[index_date_targ] #array of filenumbers in target date



        #------
        # getting rid of frames of the same set in correlation matrix
        #------
        for element in arr_targpsf:
                index_i_chunk = np.argwhere(arr_filenumbs == element)
                index_i_chunk = index_i_chunk[0][0]
                print 'index_i_chunk', index_i_chunk
                print 'check it', (index_date_targ*max_numb_imgs) + index_i_chunk
                matrix_corr[:, (index_date_targ*max_numb_imgs) + index_i_chunk] = 100 #change number if necessary. Now, 100 is bad correlation.
        numb_dates = len(struct_ShaneAO_filenumbs)
        arr_indexdate_corrs = np.zeros(numb_dates*max_numb_imgs) #array to be filled with indexes of dates correspoding to corr matrix
        for index_indexdatecorrs in np.arange(numb_dates):
                arr_indexdate_corrs[index_indexdatecorrs*max_numb_imgs:(index_indexdatecorrs+1)*max_numb_imgs] = index_indexdatecorrs
        arr_indexdate_corrs = arr_indexdate_corrs.astype(int)

        
        
        #------
        # Create list of structures with info for klip
        # Takes certain number of best imgs by referencing correlation matrix
        #------
        arr_pool_vars = []	
        for index_i in np.arange(arr_targpsf.size):
                arr_j_optimal = []
                arr_img2_optimal = []
                i = arr_targpsf[index_i]
                print i
                index_i = np.argwhere(arr_filenumbs == i)[0][0]
                print 'index_i', index_i
		arr_corrs = matrix_corr[(index_date_targ*max_numb_imgs) + index_i, :]
                index_sort_corr = np.argsort(arr_corrs)
                print 'index_sort_corr', index_sort_corr
                arr_indexdateoptimal = arr_indexdate_corrs[index_sort_corr]
                print 'arr_indexdateoptimal', arr_indexdateoptimal
                arr_filenumbs_optimal = arr_filenumbs_corr[index_sort_corr]
                arr_index_j_optimal = []
                for index_optimal in np.arange(arr_filenumbs_optimal.size).astype(int):
                        index_date_optimal = arr_indexdateoptimal[index_optimal]
                        #print 'index_date_optimal', index_date_optimal
                        filenumb_optimal = arr_filenumbs_optimal[index_optimal]
                        #print 'filenumb_optimal', filenumb_optimal
                        #print 'len(arr_substar[index_date_optimal])', len(arr_substar[index_date_optimal])
                        #print 'len(arr_j[index_date_optimal])', len(arr_j[index_date_optimal])
                        index_j_optimal = np.argwhere(arr_j[index_date_optimal] == filenumb_optimal)
                        #print 'index_j_optimal', index_j_optimal
                        if index_j_optimal:
                                arr_j_optimal.append(arr_j[index_date_optimal][index_j_optimal[0][0]])
                                arr_img2_optimal.append(arr_img2[index_date_optimal][index_j_optimal[0][0]])
                        if len(arr_j_optimal) >= numb_imgs_keep:
                                break

                arr_img2_optimal = np.array(arr_img2_optimal)
                #print 'arr_img2_optimal.shape', arr_img2_optimal.shape
                for index_row in np.arange(len(arr_img2_optimal)): # subtract each row by mean of row
                        mean_row = np.mean(arr_img2_optimal[index_row])
                        arr_img2_optimal[index_row] -=  mean_row
                        
                print 'mean of arr_img2_optimal should be zero:'
                print 'np.mean(arr_img2_optimal):', np.mean(arr_img2_optimal)
                for radius_mock in arr_radiusmock: # create list of structures
                        for mock_factor in arr_mockfactor:
                                for theta in arr_theta:
                                        for numb_klip_modes in arr_klipmodes:
                                                struct = {'filenumb_img1': i, 'arr_img2': arr_img2_optimal, 'date_ref': date_ref, 'numb_klip_modes':numb_klip_modes, 'max_rad_img': max_rad_img, 'rad_annulus': rad_annulus, 'radius_mock':radius_mock, 'mock_factor': mock_factor, 'theta': theta}
                                                arr_pool_vars.append(struct)

                                                

        #------
        # save/load list of structures as needed
        #------
        filename_arrpoolvars = 'klip_arr_pool_vars_' + date + 'set' + str(int(setnumb1)) + '.p'
        pickle.dump(arr_pool_vars, open(filename_arrpoolvars, 'wb')) #saves list of structures
        #arr_pool_vars = pickle.load(open(filename_arrpoolvars, 'rb')) #loads list of structures



        #------
        # run KLIP
        #------
        print '------------------------'
	print 'No of loops to run:', len(arr_pool_vars)
        print '------------------------'
        print 'Now starting threads...'

	for elem in arr_pool_vars:
		klip(elem)
        
	'''
	threadnumb = 6
	thepool = ThreadPool(threadnumb)
        arr_pool = thepool.map(klip, arr_pool_vars)
        thepool.close()
        thepool.join()
        '''
        
def create_cube_klip():
        arr_klipmodes = np.arange(5, 50, 5)
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 45, 60]).astype(int) #np.arange(15, 60+1, 15).astype(int)
        print 'arr_radiusmock', arr_radiusmock
        arr_mockfactor = np.array([30, 100, 150, 200, 500, 1000])#np.arange(100, 200+1, 50).astype(int)
        print 'arr_mockfactor', arr_mockfactor
        arr_theta = np.array([0]) #np.arange(0, 360, 60).astype(int)
        print 'arr_theta', arr_theta

        #---------------
        #ask for user input regarding which set to take as target psf
        #reads txt file in directory to know which files to use
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
	arr_targpsf = np.arange(start1, end1+1)


        
        cube = []
        for numb_klip_modes in arr_klipmodes:
                arr_mdn = []
                for filenumb_img1 in arr_targpsf:
                        filename_output = directory + '/' + date + '/' + str(filenumb_img1) + '_klipfinal_' + 'kmode' + str(int(numb_klip_modes)) + '.fits'
                        if exists(filename_output):
                                img = afits.open(filename_output, memmap = False)[0].data
                                arr_mdn.append(img)
                arr_mdn = np.median(np.array(arr_mdn), axis = 0)                
                cube.append(arr_mdn)
        cube = np.array(cube)
        print 'cube.shape', cube.shape
        filename_cube = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_klipcube.fits'
        hdu = pf.PrimaryHDU(cube)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_cube, clobber=True)
                
        
def klip(struct):
        
        filenumb_img1 = struct['filenumb_img1']
        arr_img2 = struct['arr_img2']
        date_ref = struct['date_ref']
        numb_klip_modes = struct['numb_klip_modes']
        max_rad_img = struct['max_rad_img']
        rad_annulus = struct['rad_annulus']
        radius_mock = struct['radius_mock']
        mock_factor = struct['mock_factor']
        theta_mock = struct['theta']
        arr_radstarts = np.arange(1, max_rad_img+1, rad_annulus)


        
        #---------------------
        #creating mock binary to be injected into img
        #---------------------
        try:
                img_mock = pf.open(directory + '/' + date + '/' + 'centroid_1.fits')[0].data
        except:
                return 0
        img_mock /= np.amax(img_mock)
        if mock_factor != 0:
                img_mock /= float(mock_factor)
        else:
                img_mock *= float(mock_factor)



        #---------------------------
        #loading target img, divide by max pixel value
        #---------------------------        
        if filenumb_img1 > 999:
                filename = directory + '/' + date_ref + '/' + 's' + str(filenumb_img1) + '_reduced_centered.fits'
                if exists(filename):
                        img1 = pf.open(filename)[0].data
                        maxpixval = np.amax(img1)
                        img1 /= maxpixval
                else:
                        print filename, 'does not exist'
                        return 0
        else:
                filename = directory + '/' + date_ref + '/' + 's0' + str(filenumb_img1) + '_reduced_centered.fits'
                if exists(filename):
                        img1 = pf.open(filename)[0].data
                        maxpixval = np.amax(img1)
                        img1 /= maxpixval
                else:
                        print filename, 'does not exist'
                        return 0


                
        #------------------------------
        #inject mock-up binary
        #------------------------------
        if radius_mock != 0:
                dx_mock = radius_mock*np.cos(theta_mock*(2*np.pi)/360)
                dx_mock = int(round(dx_mock))
                dy_mock = radius_mock*np.sin(theta_mock*(2*np.pi)/360)
                dy_mock = int(round(dy_mock))
        
                img_mock = np.fft.fft2(img_mock) #Fourier shift to construct binary
                img_mock = fshift(img_mock, [dy_mock, dx_mock])
                img_mock = np.fft.ifft2(img_mock)
                img_mock = np.real(img_mock)
                img1 += img_mock

                

        #-----------------------------------        
        #save initial img with mock up binary added
        #-----------------------------------
        filename_imgmock = directory + '/' + date + '/' + str(filenumb_img1) + '_klipmockimg_' + 'kmode' + str(int(numb_klip_modes)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_theta' + str(int(theta_mock)) + '.fits'
        hdu = pf.PrimaryHDU(img1)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_imgmock, clobber=True) 



        #----------
        #flatten img and save img shape
        #----------
        img1shape = img1.shape
        img1_flat = img1.flatten()
        radi = 0

        

        #-------
        # Define filenames of fits files
        # Exit function if files already exist
        # CHANGE FILENAMES IF NECESSARY        
        #-------
        filename_output = directory + '/' + date + '/' + str(filenumb_img1) + '_klipfinal_' + 'kmode' + str(int(numb_klip_modes)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_theta' + str(int(theta_mock)) + '.fits'
        filename_sub = directory + '/' + date + '/' + str(filenumb_img1) + '_klipsubimg_' + 'kmode' + str(int(numb_klip_modes)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_theta' + str(int(theta_mock)) + '.fits'
        if exists(filename_sub) and exists(filename_output) and exists(filename_imgmock):
                return 0



        #----------|
        #   KLIP   |
        #----------|
        img_sub = np.zeros(img1_flat.shape)
        img_final = np.zeros(img1_flat.shape)
        for index_rad in np.arange(arr_radstarts.size - 1).astype(int):
                arr_ringindex = []
                radi = arr_radstarts[index_rad]
                radi_next = arr_radstarts[index_rad+1]
                while radi < radi_next:
                        ringvals, arr_ringindex_temp = return_ring(img1, radi)
                        for ringindex in arr_ringindex_temp:
                                ringindex = np.ravel_multi_index(ringindex, img1shape)
                                arr_ringindex.append(ringindex)
                        radi += 1
                #print 'len(arr_ringindex)', len(arr_ringindex)
                arr_img2_S = arr_img2[:, arr_ringindex]
                img1_ring = img1_flat[arr_ringindex]

                basis = get_klip_basis(arr_img2_S, numb_klip_modes)
                tZ_dot_Z = np.dot(np.transpose(basis), basis)
                img_sub_ring = np.dot(img1_ring, tZ_dot_Z)
                img_final_ring = img1_ring - img_sub_ring        
                for index_area in np.arange(len(arr_ringindex)).astype(int):
                        index_ring = arr_ringindex[index_area]
                        img_final[index_ring] = img_final_ring[index_area]
                        img_sub[index_ring] = img_sub_ring[index_area]


                        
        #------
        #Save final img, subtraction img as fits files
        #------
        hdu = pf.PrimaryHDU(img_final.reshape(img1shape))
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_output, clobber=True)
        
        hdu = pf.PrimaryHDU(img_sub.reshape(img1shape))
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_sub, clobber=True)
        
        print 'done with', filenumb_img1, 'kmode', int(numb_klip_modes), ', radius', int(radius_mock), ', ratio', int(mock_factor), ', theta', int(theta_mock)
        return 1





def get_klip_basis(R, cutoff):
        #takes matrix R with N rows, M columns. N is number of imgs, M is number of pixels in img
        #cutoff is number of klip modes to use.
        ##

        w, V = np.linalg.eig(np.dot(R, np.transpose(R)))
        sort_ind = np.argsort(w)[::-1] #indices of eigenvals sorted in descending order
        sv = np.sqrt(w[sort_ind]).reshape(-1,1) #column of ranked singular values
        Z = np.dot(1./sv*np.transpose(V[:, sort_ind]), R)
        return Z[0:cutoff, :]



def test_klip():
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])

        arr_for_mdn = []
        for i in np.arange(start1, end1 +1):
                filename = directory + '/' + date + '/' + str(i) + '_klipfinal_test.fits'
                if exists(filename):
                        img = pf.open(filename)[0].data
                        arr_for_mdn.append(img)
        img_test = np.median(np.array(arr_for_mdn), axis = 0)
        hdu = pf.PrimaryHDU(img_test)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto('test.fits', clobber=True)
        subprocess.call('/home/apps/ds9/ds9 ' + 'test.fits', shell = True)        
                


def run_correlation_matrix_thread():
        arr_indexdate1 = np.arange(3)
        #'''
        threadnumb = 3
	thepool = ThreadPool(threadnumb)
        arr_imgfinals = thepool.map(correlation_matrix, arr_indexdate1)
        thepool.close()
        thepool.join()
        #'''
        '''
        for index_date1 in arr_indexdate1:
                correlation_matrix(index_date1)
        '''
        
def merge_mat():
        max_numb_imgs = 5000
        numb_dates = 3
        
        mat0 = pickle.load(open('matrix_corr0.p', 'rb'))
        mat1 = pickle.load(open('matrix_corr1.p', 'rb'))
        mat2 = pickle.load(open('matrix_corr2.p', 'rb'))
        
        mat = np.zeros([numb_dates*max_numb_imgs, numb_dates*max_numb_imgs])
        mat[0:max_numb_imgs,:] = mat0[0:max_numb_imgs,:]
        mat[max_numb_imgs:2*max_numb_imgs,:] = mat1[max_numb_imgs:2*max_numb_imgs,:]
        mat[2*max_numb_imgs:3*max_numb_imgs,:] = mat2[2*max_numb_imgs:3*max_numb_imgs,:]
        
        pickle.dump(mat, open('matrix_corr.p', 'wb'))
        
        
def correlation_matrix(index_date1):
        numb_dates = 3  ### VARIABLE
        fillelement = 100
        matrix_corr = np.zeros([numb_dates*max_numb_imgs, numb_dates*max_numb_imgs])
	imgshape = [161, 161] ###Depends on data set
        matrix_corr.fill(fillelement)
        arr_date_ShaneAO = ['Jun15', 'Nov15', 'Mar16'] #### Why do I have to redefine this?
        #arr_date_ShaneAO = arr_date_ShaneAO ###[:2] ### VARIABLE

        filename_mat = "matrix_corr" + str(int(index_date1)) + ".p"
        
        if exists(filename_mat):
                matrix_corr = pickle.load(open(filename_mat, 'rb'))
        else:
                print filename_mat, 'doesnt exist'
        
        struct_ShaneAO_filenumbs = {}
        for index_date in np.arange(len(arr_date_ShaneAO)):
                date = arr_date_ShaneAO[index_date]
                total_sets = struct_ShaneAO_total_sets[date]
                arr_filenumbs = np.array([])
                for setnumb in np.arange(1, total_sets+1):
                        arr_startend = open(directory + '/' + date + '/set_'+ str(int(setnumb)) +'_info.txt', 'rb').read().splitlines()
                        start = int(arr_startend[0])
                        end = int(arr_startend[1])
                        arr_filenumbs = np.append(arr_filenumbs, np.arange(start, end+1))
                struct_ShaneAO_filenumbs[index_date] = arr_filenumbs

        filename = 'struct_ShaneAO_filenumbs.p'
        pickle.dump(struct_ShaneAO_filenumbs, open(filename, 'wb'))
		
	img_empty = np.zeros(imgshape)
	max_rad_ring = 4
        arr_ringindex = []
        counter_ringindex = 0
        for rad_ring in np.arange(0, max_rad_ring+1):
                ringvals, arr_ringindex_temp = return_ring(img_empty, rad_ring)
                for ringindex in arr_ringindex_temp:
                        counter_ringindex += 1
                        ringindex = np.ravel_multi_index(ringindex, imgshape)
                        arr_ringindex.append(ringindex)
        
        arr_filenumbs_corr = np.zeros(matrix_corr.shape[0])
        for index_datefill in np.arange(len(struct_ShaneAO_filenumbs)):
                arr_filenumbsfill = struct_ShaneAO_filenumbs[index_datefill]
                for index_filenumbsfill in np.arange(arr_filenumbsfill.size):
                        filenumbfill = arr_filenumbsfill[index_filenumbsfill]
                        arr_filenumbs_corr[(index_datefill*max_numb_imgs) + index_filenumbsfill] = filenumbfill
        filename = 'arr_filenumbs_corr.p'
        pickle.dump(arr_filenumbs_corr, open(filename, 'wb'))

        counter = 0
	arr_sd = []
        column_count = 0
        arr_filenumbs1 = struct_ShaneAO_filenumbs[index_date1]
        for index_date2 in np.arange(len(struct_ShaneAO_filenumbs)):
                counter += 1
                arr_filenumbs2 = struct_ShaneAO_filenumbs[index_date2]
                for index_filenumb1 in np.arange(arr_filenumbs1.size):
                        print 100.*index_filenumb1/arr_filenumbs1.size, '%', 'for index_date1:', index_date1, ', column no.', counter
                        filenumb1 = int(arr_filenumbs1[index_filenumb1])
                        if filenumb1 > 999:
                                filename = directory + '/' + arr_date_ShaneAO[index_date1] + '/' + 's' + str(filenumb1) + '_reduced_centered.fits'
                                if exists(filename):
                                        img1 = pf.open(filename)[0].data
                                else:
                                        #print filename, 'does not exist'
                                        continue
                        else:
                                filename = directory + '/' + arr_date_ShaneAO[index_date1] + '/' + 's0' + str(filenumb1) + '_reduced_centered.fits'
                                if exists(filename):
                                        img1 = pf.open(filename)[0].data
                                else:
                                        #print filename, 'does not exist'
                                        continue
                        img1 = img1.flatten()
                        for index_filenumb2 in np.arange(arr_filenumbs2.size):

                                filenumb2 = int(arr_filenumbs2[index_filenumb2])
                                #print filenumb2
                                if matrix_corr[(index_date1*max_numb_imgs) + index_filenumb1, (index_date2*max_numb_imgs) + index_filenumb2]  <= (fillelement-(fillelement/100)):
                                        continue
                                if filenumb2 > 999:
                                        filename = directory + '/' + arr_date_ShaneAO[index_date2] + '/' + 's' + str(filenumb2) + '_reduced_centered.fits'
                                        if exists(filename):
                                                img2 = pf.open(filename)[0].data
                                        else:
                                                #print filename, 'does not exist'
                                                continue
                                else:
                                        filename = directory + '/' + arr_date_ShaneAO[index_date2] + '/' + 's0' + str(filenumb2) + '_reduced_centered.fits'
                                        if exists(filename):
                                                img2 = pf.open(filename)[0].data
                                        else:
                                                #print filename, 'does not exist'
                                                continue
                                #print 'img2.shape', img2.shape
                                img2 = img2.flatten()
                                #print 'img2.shape after flattening', img2.shape
                                for ringindex in arr_ringindex:
                                        img1[ringindex] = 0
                                        img2[ringindex] = 0
                                '''### To view images
                                hdu = pf.PrimaryHDU(img1.reshape(imgshape))
                                hdulist = pf.HDUList([hdu])
                                hdulist.writeto('test.fits', clobber=True)
                                subprocess.call('/home/apps/ds9/ds9 ' + 'test.fits', shell = True)
                                ###
                                hdu = pf.PrimaryHDU(img2.reshape(imgshape))
                                hdulist = pf.HDUList([hdu])
                                hdulist.writeto('test.fits', clobber=True)
                                subprocess.call('/home/apps/ds9/ds9 ' + 'test.fits', shell = True)
                                ###
                                '''
                                img1 /= np.sum(img1)
                                img2 /= np.sum(img2)
                                img_res = img1 - img2
                                sd_img_res = np.std(img_res)
                                matrix_corr[(index_date1*max_numb_imgs) + index_filenumb1, (index_date2*max_numb_imgs) + index_filenumb2] = sd_img_res
                        '''
                        print 'saving matrix_corr'
                        pickle.dump(matrix_corr, open("matrix_corr.p", 'wb'))
                        print 'done saving'
                        '''
                        '''
                        img1 = pf.open(directory + '/' + 'Jun15' + '/' + 's0790' + '_reduced_centered.fits')[0].data
                        img1 /= np.sum(img1)
                        img2 = pf.open(directory + '/' + 'Jun15' + '/' + 's0792' + '_reduced_centered.fits')[0].data
                        img2 /= np.sum(img2)
                        print np.std(img1 - img2)
                        '''
                        #print 'max std value', np.amax(matrix_corr[0,:])
                        #print 'min std value', np.amin(matrix_corr[0,:])
                pickle.dump(matrix_corr, open("matrix_corr" + str(int(index_date1)) + ".p", 'wb'))        
                '''
                if counter == 2:
                print matrix_corr[(index_date1*max_numb_imgs) + index_filenumb1,4995:5010]
                useless = raw_input('press key to continue')
                '''
        pickle.dump(matrix_corr, open("matrix_corr" + str(int(index_date1)) + ".p", 'wb'))
        
        

def view_pickle(radius_mock = 23, mock_factor = 125):
	filenumb = raw_input('Pickle file number (790, 791, etc.):')
	#img_final = pickle.load(open(directory + '/' + date + '/' + filenumb + '_thread_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) +'_.p', "rb"))
	img_final = pickle.load(open(directory + '/' + date + '/' + filenumb + '_thread_mock.p', "rb"))
	hdu = pf.PrimaryHDU(img_final)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto('test_final.fits', clobber=True)

	#img_sub = pickle.load(open(directory + '/' + date + '/' + filenumb + '_subimg_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) +'_.p', "rb"))
        img_sub = pickle.load(open(directory + '/' + date + '/' + filenumb + '_subimg_mock.p', "rb"))
	hdu = pf.PrimaryHDU(img_sub)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto('test_sub.fits', clobber=True)
        '''
	img_mock = pickle.load(open(directory + '/' + date + '/' + filenumb + '_mockimg_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) +'_.p', "rb"))
	hdu = pf.PrimaryHDU(img_mock)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto('test_mock.fits', clobber=True)
        '''
	subprocess.call('/home/apps/ds9/ds9 ' + 'test_final.fits', shell = True)

def plot_pickle(radius_mock = 23, mock_factor = 125):
	
        setnumb1 = raw_input('Enter set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
	
	arr_rad = np.arange(70)+1
        arr_5sd = np.array([])
        arr_mdn = np.array([])
	arr_5sd2 = np.array([])
        arr_mdn2 = np.array([])

	arr_for_mdn = []
	for i in np.arange(start1, end1+1):
		filename = directory + '/' + date + '/' + str(int(i)) + '_thread_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) +'_.p'
		if exists(filename):
			arr_for_mdn.append(pickle.load(open(filename, "rb")))
		else:
			print filename, 'doesnt exist'			
	img_mdn = np.median(np.array(arr_for_mdn), axis = 0)
	#filename = raw_input('Name of pickle file (dont include .p):')
	#img_output = pickle.load(open(filename + ".p", "rb" ))
        hdu = pf.PrimaryHDU(img_mdn)
        hdulist = pf.HDUList([hdu])
	outputfilename = directory + '/' + date + '/' + 'set' + setnumb1 + '_loci.fits'
        print 'writing to ' + outputfilename
        hdulist.writeto(outputfilename, clobber=True)
	for radius in arr_rad:
                struct = check_ring(img_mdn, radius)
                arr_mdn = np.append(arr_mdn, struct['mdn_ring'])
                arr_5sd = np.append(arr_5sd, 5*struct['sd_ring'])
        
        pickle.dump(arr_5sd, open(directory + '/' + date + '/' + 'set' + setnumb1 + '_loci' + "_5sdplot_mock.p", 'wb'))
        pickle.dump(arr_mdn, open(directory + '/' + date + '/' + 'set' + setnumb1 + '_loci' +  "_mdnplot_mock.p", 'wb'))

	arr_for_mdn_mock = []
	for i in np.arange(start1, end1+1):
		filename = directory + '/' + date + '/' + str(int(i)) + '_mockimg_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) +'_.p'
		if exists(filename):
			arr_for_mdn_mock.append(pickle.load(open(filename, "rb")))
		else:
			print filename, 'doesnt exist'			
	img_mdn2 = np.median(np.array(arr_for_mdn_mock), axis = 0)
	#filename = raw_input('Name of pickle file (dont include .p):')
	#img_output = pickle.load(open(filename + ".p", "rb" ))
        hdu = pf.PrimaryHDU(img_mdn2)
        hdulist = pf.HDUList([hdu])
	outputfilename2 = directory + '/' + date + '/' + 'set' + setnumb1 + '_mockimg.fits'
        hdulist.writeto(outputfilename2, clobber=True)
	for radius in arr_rad:
                struct = check_ring(img_mdn2, radius)
                arr_mdn = np.append(arr_mdn2, struct['mdn_ring'])
                arr_5sd = np.append(arr_5sd2, 5*struct['sd_ring'])
        
        pickle.dump(arr_5sd, open(directory + '/' + date + '/' + 'set' + setnumb1 + '_mockimg' + "_5sdplot_mock.p", 'wb'))
        pickle.dump(arr_mdn, open(directory + '/' + date + '/' + 'set' + setnumb1 + '_mockimg' +  "_mdnplot_mock.p", 'wb'))

	subprocess.call('/home/apps/ds9/ds9 ' + outputfilename, shell = True)


def psf_subtract(file1, file2):
        #Subtract different psfs from one another
        #load files
        '''
        file1 = raw_input("First centroid? (Don't include .fits):")
        file2 = raw_input("2nd centroid? (Don't include .fits):")
        '''
	#define flux, radius, subtraction ratio arrays to be iterated through
        #arr_flux = (np.arange(10.)+1)/100.
        arr_flux = np.logspace(-3, -1, num = 10, endpoint = True)
        arr_rad = np.arange(50)+1
        arr_subratio = (np.arange(20)+90)/100
        #________________________________________________

        arr_5sd = np.array([])
        arr_mdn = np.array([])
	#arr_detect = np.empty([arr_rad.size, arr_flux.size])
        #arr_detect.fill(0)
        #arr_detect = pickle.load(open( file1 + "-" + file2 + "_detlim.p", "rb" ))
        
	'''
	####To modify matrix manually####
	arr_detect[y,x] = 1
	pickle.dump(arr_detect, open(file1 + "-" + file2 + "_detlim.p", 'wb'))
	'''
	
        arr_rss = np.array([]) #define empty array of rss values
        arr_output = [] #define empty array of imgs
        for subratio in arr_subratio:
                img2 = pf.open('centroid_'+file2+'.fits')[0].data #open img2, img used for subtraction
                img2 /= np.amax(img2) #divide by max. brightness
                img1 = pf.open('centroid_'+file1+'.fits')[0].data #open img1, the img on which we perform psf subtraction
                img1 /= np.amax(img1) #divide by max. brightness
                center_flux = np.amax(img1)                                
                img2 *= subratio #multiply mock by subtraction ratio
                img1 -= img2 #psf subtract
                arr_rss = np.append(arr_rss, np.sum(img1**2.)) #add residual sum of squares to arr_rss.
                arr_output.append(img1) #append psf-subtracted img
        index_min = np.argmin(arr_rss)
        #print 'best subtraction ratio', arr_subratio[index_min]
        #print 'size arr_rss', arr_rss.size
        #print arr_rss
        img_final = arr_output[index_min]
        #print 'img_final.size', img_final.size
        #print arr_detect
        #a = raw_input('stopped')
        for radius in arr_rad:
                struct = check_ring(img_final, radius)
                arr_mdn = np.append(arr_mdn, struct['mdn_ring'])
                arr_5sd = np.append(arr_5sd, 5*struct['sd_ring'])
        
        pickle.dump(arr_5sd, open(file1 + "-" + file2 + "_5sdplot.p", 'wb'))
        pickle.dump(arr_mdn, open(file1 + "-" + file2 + "_mdnplot.p", 'wb'))
        
def check_ring(img, radi): #inputs are: img array, radius of ring
#check img for possible binaries. 
#Check if any points in a ring are greater than 5sd above the median in a ring
	#file1 = raw_input("Filename? (Don't include .fits):")
	#img = pf.open(file1 + '.fits')[0].data
        #print 'check ring: img.shape', img.shape
	y_length = img.shape[0]
	x_length = img.shape[1]
	center_index = [(y_length-1)/2, (x_length-1)/2]
	#print 'center index', center_index
	center = np.array(center_index) ###_value at center_###
	arr_index = []
	for y_index in np.arange(y_length):
		for x_index in np.arange(x_length):
                        if (radi-0.5) < distance(center_index[0], y_index, center_index[1], x_index) < (radi+0.5): #creating ring of radius: radi
				arr_index.append([y_index, x_index]) #appending indexes of points in ring
	#img *= 0
	#print 'arr_index', arr_index
	arr_ring = np.array([])
	for element in arr_index:
		#print element
		#print 'value', img[element[0]][element[1]]
		arr_ring = np.append(arr_ring, img[element[0]][element[1]])
	arr_candidate = []

        max_ring = np.amax(arr_ring)
        mdn_ring = (np.median(arr_ring))
        sd_ring = np.std(arr_ring)
        
	#print 'arr_ring', arr_ring
	#print 'max', max_ring
	#print 'median', mdn_ring
	#print 'sd', sd_ring
        struct = {'max_ring':max_ring, 'mdn_ring':mdn_ring, 'sd_ring':sd_ring}        
        return struct
	#next_point = True
        #sd_fact = 5.
	#while next_point:
        #if np.amax(arr_ring) > (np.median(arr_ring) + (sd_fact*np.std(arr_ring))):
                #print 'candidate?:', 1
                #return 1
        #arr_candidate.append(np.argwhere(img == np.amax(arr_ring)))
        #index_max = np.argmax(arr_ring)
        #img[index_max[0]][index_max[1]] = 0
        #else:
                #next_point = False
                #print 'candidate?:', 0
                #return 0
	#hdu = pf.PrimaryHDU(img)
	#hdulist = pf.HDUList([hdu])
	#hdulist.writeto('test_ring.fits', clobber=True)
	#subprocess.call('open -a SAOImage\ DS9 ' + 'test_ring.fits', shell = True)


def check_ring_mock(img, radi, theta_mock): #inputs are: img array, radius of ring, angle where binary is positioned
        #check img for possible binaries. 
        #Check if any points in a ring are greater than 5sd above the median in a ring

        theta_pm = 50 #exclude plus-minus this angle around theta_mock
        
        y_length = img.shape[0]
	x_length = img.shape[1]
	center_index = [(y_length-1)/2, (x_length-1)/2]
	#print 'center index', center_index
	center = np.array(center_index) ###_value at center_###
	arr_index = []
	for y_index in np.arange(y_length):
		for x_index in np.arange(x_length):
                        if (radi-0.5) < distance(center_index[0], y_index, center_index[1], x_index) < (radi+0.5): #creating ring of radius: radi
                                if correct_theta(x_index, center_index[1], y_index, center_index[0], theta_mock, theta_pm):
                                        arr_index.append([y_index, x_index]) #appending indexes of points in ring
        '''
        print '--------------------'
        print 'len(arr_index)', len(arr_index)
        print '--------------------'
        '''
        
	arr_ring = np.array([])
	for element in arr_index:
		#print element
		#print 'value', img[element[0]][element[1]]
		arr_ring = np.append(arr_ring, img[element[0]][element[1]])
	arr_candidate = []

        max_ring = np.amax(arr_ring)
        mdn_ring = (np.median(arr_ring))
        sd_ring = np.std(arr_ring)
        
	#print 'arr_ring', arr_ring
	#print 'max', max_ring
	#print 'median', mdn_ring
	#print 'sd', sd_ring
        struct = {'max_ring':max_ring, 'mdn_ring':mdn_ring, 'sd_ring':sd_ring}        
        return struct
	#next_point = True
        #sd_fact = 5.
	#while next_point:
        #if np.amax(arr_ring) > (np.median(arr_ring) + (sd_fact*np.std(arr_ring))):
                #print 'candidate?:', 1
                #return 1
        #arr_candidate.append(np.argwhere(img == np.amax(arr_ring)))
        #index_max = np.argmax(arr_ring)
        #img[index_max[0]][index_max[1]] = 0
        #else:
                #next_point = False
                #print 'candidate?:', 0
                #return 0
	#hdu = pf.PrimaryHDU(img)
	#hdulist = pf.HDUList([hdu])
	#hdulist.writeto('test_ring.fits', clobber=True)
	#subprocess.call('open -a SAOImage\ DS9 ' + 'test_ring.fits', shell = True)

def return_ring(img, radi): #inputs are: img array, radius of ring
#returns array of values in a ring in img of certain radius: radi
#returns array of indexes of ring points in the image

	y_length = img.shape[0]
	x_length = img.shape[1]
	center_index = [(y_length-1)/2, (x_length-1)/2]
        center_index_y = center_index[0]
        center_index_x = center_index[1]
	center = np.array(center_index) ###_value at center_###
	arr_index = []
	for y_index in np.arange(center_index_y - radi, center_index_y + radi +1):
		for x_index in np.arange(center_index_x - radi, center_index_x + radi +1):
                        if (radi-0.5) < distance(center_index[0], y_index, center_index[1], x_index) < (radi+0.5): #creating ring of radius: radi
				arr_index.append([y_index, x_index]) #appending indexes of points in ring
	arr_ring = np.array([])
	for element in arr_index:
		arr_ring = np.append(arr_ring, img[element[0]][element[1]])

        return arr_ring, arr_index

def plotsframe():
        total_sets = 4 ##### subject to change!

        file1 = raw_input("Enter No. of 1st set:")
	minpix = 6



        img_bef = pf.open(directory + '/' + date + '/' + 'centroid_' + file1 + '.fits')[0].data
        img_bef /= np.amax(img_bef) #divide by max. brightness
        center_val = np.amax(img_bef)
        arr_radi = np.arange(70)+1
        arr_colors = ['rb', 'g', 'y']
        arr_centroid = []
        counter = 0
	'''
        for elem in np.arange(total_sets)+1:
                if elem != int(file1):
                        frame1set_5sd = pickle.load(open(directory + '/' + date + '/' + file1 + '-' + str(int(elem)) + '_5sdplot.p', 'rb'))
                        plt.plot(arr_radi[minpix:frame1set_5sd.size], 2.5*np.log10(center_val/np.abs(frame1set_5sd[minpix:])), arr_colors[counter] +'o-', label = '5 s.d. value: set' + str(file1) + '-' + 'set' + str(int(elem)))
                        counter += 1


        arr_5sd_frame = pickle.load(open(directory + '/' + date + '/' + 'set' + file1 + '-' + 'others' + '_5sdplot.p', 'rb'))
        plt.plot(arr_radi[minpix:arr_5sd_frame.size], 2.5*np.log10(center_val/np.abs(arr_5sd_frame[minpix:])), 'bo-', label = '5 s.d. value: set ' + str(file1) + '-' + 'all other sets')   

	'''

        mock_5sd = pickle.load(open(directory + '/' + date + '/' + 'set' + file1 + '_mockimg' + "_5sdplot_mock.p", 'rb'))
        print mock_5sd.size
        plt.plot(arr_radi[minpix:mock_5sd.size], 2.5*np.log10(center_val/np.abs(mock_5sd[minpix:])), 'co-', label = '5 s.d. value: set ' + str(file1) + ' Initial img with mock binary')   

	'''
        loci_5sd = pickle.load(open(directory + '/' + date + '/' + 'set' + file1 + '_loci' + "_5sdplot.p", 'rb'))
        print loci_5sd.size
        plt.plot(arr_radi[minpix:loci_5sd.size], 2.5*np.log10(center_val/np.abs(loci_5sd[minpix:])), 'co-', label = '5 s.d. value: set ' + str(file1) + 'Initial img with mock binary')   
	'''

        loci_mock_5sd = pickle.load(open(directory + '/' + date + '/' + 'set' + file1 + '_loci' + '_5sdplot_mock.p', 'rb'))
        print 'loci_mock_5sd.size', loci_mock_5sd.size
        print 'arr_radi[minpix:loci_mock_5sd.size].shape', arr_radi[minpix:loci_mock_5sd.size].shape
        print 'loci_mock_5sd[minpix:].shape', loci_mock_5sd[minpix:].shape
        plt.plot(arr_radi[minpix:loci_mock_5sd.size], 2.5*np.log10(center_val/np.abs(loci_mock_5sd[minpix:])), 'mo-', label = '5 s.d. value: set ' + str(file1) + '(LOCI) with mock binary at radius 23')   

        plt.legend(loc = 'upper right')
#        plt.yscale('log')
        plt.xlabel('radius (pixels)')
        plt.ylabel('magnitude')
        plt.gca().invert_yaxis()
        plt.show()

def plots1():
        file1 = raw_input("No. of 1st centroid? (Don't include .fits):")
        file2 = raw_input("No. of 2nd centroid? (Don't include .fits):")

        img_bef = pf.open('centroid_' + file1 + '.fits')[0].data
        img_bef /= np.amax(img_bef) #divide by max. brightness
        center_val = np.amax(img_bef)
        arr_mdn_bef = np.array([])
        arr_5sd_bef = np.array([])
        arr_radi = np.arange(50)+1
        for radi in arr_radi:
                struct = check_ring(img_bef, radi)
                arr_mdn_bef = np.append(arr_mdn_bef, struct['mdn_ring'])
                arr_5sd_bef = np.append(arr_5sd_bef, 5*struct['sd_ring'])

        arr_5sd = pickle.load(open(file1 + '-' + file2 + '_5sdplot.p', 'rb'))
        arr_mdn = pickle.load(open(file1 + '-' + file2 + '_mdnplot.p', 'rb'))
        arr_radi = np.arange(50)+1
        '''
        plt.figure()
        plt.plot(arr_radi, 2.5*np.log10(center_val/np.abs(arr_mdn_bef)), 'rx-', label = 'Median value (before subtraction)' )
        plt.plot(arr_radi, 2.5*np.log10(center_val/np.abs(arr_mdn)), 'bx-', label = 'Median value (after subtraction)')
        plt.legend(loc = 'upper right')
#        plt.yscale('log')
        plt.xlabel('radius (pixels)')
        plt.ylabel('magnitude')        
        plt.gca().invert_yaxis()
        '''

        plt.figure()
        plt.plot(arr_radi, 2.5*np.log10(center_val/np.abs(arr_5sd_bef)), 'ro-', label = '5 s.d. value (before subtraction)' )
        plt.plot(arr_radi, 2.5*np.log10(center_val/np.abs(arr_5sd)), 'bo-', label = '5 s.d. value (after subtraction)')        
        plt.legend(loc = 'upper right')
#        plt.yscale('log')
        plt.xlabel('radius (pixels)')
        plt.ylabel('magnitude')
        plt.gca().invert_yaxis()
        plt.show()

def distance(x1, x2, y1, y2):
	return math.sqrt( ((x1 - x2)**2.) + ((y1 - y2)**2.) )

def correct_theta(x1, x_cent, y1, y_cent, theta_mock, theta_pm):
	theta_check = np.arctan2((y1 - y_cent), (x1 - x_cent))
        if theta_check < 0:
                theta_check += 2*np.pi
        theta_check *= (180 / np.pi)
        #print 'theta_check', theta_check
        if np.abs(theta_check - theta_mock) > 180:
                if theta_check < theta_mock:
                        theta_check += 360
                else:
                        theta_mock += 360

        angle_diff =  np.abs(theta_check - theta_mock)
        #print 'angle_diff',angle_diff
        #print 'theta_pm', theta_pm
        if angle_diff > theta_pm:
                return True
        else:
                return False



def gauss_func((y, x), amplitude, yo, xo, sigma_y, sigma_x, theta, offset):#((y, x), sd, offset):
        xo = float(xo)
        yo = float(yo)
        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
        g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
        return g.ravel()
        

def test_gauss_func():
        x = np.linspace(0, 160, 161)
        y = np.linspace(0, 160, 161)
        x_plot = np.copy(x)
        y_plot = np.copy(y)
        x, y = np.meshgrid(x, y)
        x = x.reshape(x.size)
        y = y.reshape(y.size)
        a = gauss_func((y,x), 1000, 100, 80, 15, 5, 10)
        hdu = pf.PrimaryHDU(a.reshape([161,161]))
        hdulist = pf.HDUList([hdu])
        hdulist.writeto('test.fits', clobber=True)
        subprocess.call('/home/apps/ds9/ds9 test.fits', shell = True)

def gauss_fit(start_test, end_test, y_binary = 0, x_binary = 0):
	fwhm_threshold = 10
	maxpix_threshold = 30000
	minpix_threshold = 5000

        x = np.linspace(0, 160, 161)
        y = np.linspace(0, 160, 161)
        x_plot = np.copy(x)
        y_plot = np.copy(y)
        x, y = np.meshgrid(x, y)

        x = x.reshape(x.size)
        y = y.reshape(y.size)
        print 'x.shape', x.shape
        print x
        print 'y.shape', y.shape
        print y

        initial_guess = (10000, y_index_center+y_binary, x_index_center+x_binary, 3, 3, 0, 0)

        arr_setnumbs = np.arange(1, total_sets+1)

        for setnumb in [1]:#arr_setnumbs: ####I
                setnumb = str(int(setnumb))
                arr_startend = open(directory + '/' + date + '/set_'+ setnumb +'_info.txt', 'rb').read().splitlines()
                start = int(arr_startend[0])
                end = int(arr_startend[1])
                arr_files = np.arange(start, end+1)
                
		#-------------------
		arr_files = np.arange(start_test, end_test + 1)
		#___________________

                arr_maxflux = []
                arr_fwhm = []
                for i in arr_files:
                        i = int(i)
                        print i
                        if i > 999:
                                filename = directory + '/' + date + '/' + 's' + str(i) + '_reduced_centered.fits'
                                if exists(filename):
                                        img = pf.open(filename)[0].data
                                else:
                                        print filename, 'does not exist'
                                        continue
                        else:
                                filename = directory + '/' + date + '/' + 's0' + str(i) + '_reduced_centered.fits'
                                if exists(filename):
                                        img = pf.open(filename)[0].data
                                else:
                                        print filename, 'does not exist'
                                        continue
                        #print 'img.shape', img.shape
                        #print 'img', img

                        img_2d = np.copy(img)
                        #print 'img_2d.shape', img_2d.shape

                        img = img.ravel()
                        #print 'img.shape', img.shape
                        try:
                                popt, pcov = curve_fit(gauss_func, np.array([y, x]), img, p0 = initial_guess)
                        except:
                                print 'curve_fit didnt work with file:', i
				subprocess.call('/home/apps/ds9/ds9 ' + filename, shell = True)
                                continue
                        #print 'popt.shape', popt.shape
                        #print 'popt', popt
                        img_fitted = gauss_func((x, y), *popt)
                        #print 'img_fitted.shape', img_fitted.shape
                        #print 'img_fitted', img_fitted
                        sigma_y = abs(popt[3])
                        sigma_x = abs(popt[4])
                        sigma_avg = np.average([sigma_y, sigma_x])
                        arr_maxflux.append(np.amax(img))
			fwhm = 2*np.sqrt(2*np.log(2))*sigma_avg
                        arr_fwhm.append(fwhm)
                        if fwhm > fwhm_threshold or np.amax(img) > maxpix_threshold or np.amax(img) < minpix_threshold:
				print 'FWHM:', fwhm
				subprocess.call('/home/apps/ds9/ds9 ' + filename, shell = True)
			
                        '''
                        hdu = pf.PrimaryHDU(img_fitted.reshape(img_2d.shape))
                        hdulist = pf.HDUList([hdu])
                        hdulist.writeto('test.fits', clobber=True)
                        subprocess.call('/home/apps/ds9/ds9 test.fits', shell = True)
                        #plt.imshow(img.reshape(img_2d.shape), cmap = plt.cm.jet, extent = [np.min(x), np.max(x), np.min(y), np.max(y)])
                        #plt.imshow(img_fitted.reshape(img_2d.shape))
                        #plt.show()
                        '''

                        print '_____________'
                #fig, ax = plt.subplots()
                #pickle.dump(arr_maxflux, open(directory + '/' + date + '/' + 'arr_maxflux.p', 'wb'))
                #pickle.dump(arr_fwhm, open(directory + '/' + date + '/' + 'arr_fwhm.p', 'wb'))
                #plt.plot(arr_maxflux, arr_fwhm, 'o', color = np.random.rand(3,1))
        #plt.xlabel('Maximum pixel value')
        #plt.ylabel('F.W.H.M.')
        #plt.xticks([0, 10000, 20000, 30000, 40000])
        #plt.show()



        #print 'arr_files.size', arr_files.size

def get_fwhm(img, y_binary = 0, x_binary = 0):

	################
	#img = pf.open(directory + '/' + date + '/' + 's' + imgnumbstring + '_reduced_centered.fits')[0].data
	###############
	
	
        halfybox = 4
        halfxbox = 4

        y_length, x_length = img.shape
	#print 'img.shape', img.shape
        y_index_center = int((y_length - 1)/2)
        x_index_center = int((x_length - 1)/2)
        

        img = img[y_index_center-halfybox:y_index_center+halfybox+1, x_index_center-halfxbox:x_index_center+halfxbox+1]
        ybox, xbox = img.shape
        #print 'img box shape', img.shape

	'''
	######################
	hdu = pf.PrimaryHDU(img)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto('test_get_fwhm.fits', clobber=True)
	subprocess.call('/home/apps/ds9/ds9 test_get_fwhm.fits', shell = True)
	##################
	'''

        x = np.arange(xbox)
        y = np.arange(ybox)
        x_plot = np.copy(x)
        y_plot = np.copy(y)
        x, y = np.meshgrid(x, y)

        x = x.reshape(x.size)
        y = y.reshape(y.size)
	'''
        print 'x.shape', x.shape
        print x
        print 'y.shape', y.shape
        print y
	'''


        #initial_guess = (12000, y_index_center+y_binary, x_index_center+x_binary, 4, 4, 0, 0)
	initial_guess = (1, halfybox+y_binary, halfxbox+x_binary, 2, 2, 0, 0)

        arr_setnumbs = np.arange(1, total_sets+1)
	img_2d = np.copy(img)
	#print 'img_2d.shape', img_2d.shape

	img = img.ravel()
	#print 'img.shape', img.shape
	try:
		popt, pcov = curve_fit(gauss_func, np.array([y, x]), img, p0 = initial_guess, maxfev = 3200)
	except:
		print 'curve_fit didnt work with file'
		print '_____________________________'
		print '_____________________________'
		print '_____________________________'
		print '_____________________________'
		return False
	
	#print 'popt.shape', popt.shape
	#print 'popt', popt
	img_fitted = gauss_func((x, y), *popt)
	#print 'img_fitted.shape', img_fitted.shape
	#print 'img_fitted', img_fitted
	sigma_y = abs(popt[3])
	sigma_x = abs(popt[4])
	sigma_avg = np.average([sigma_y, sigma_x])
	fwhm = 2*np.sqrt(2*np.log(2))*sigma_avg
	if fwhm > 10:
		print 'FWHM', fwhm
	return fwhm

        








def view_detlim():
	file1 = raw_input("No. of 1st centroid? (Don't include .fits):")
	file2 = raw_input("No. of 2nd centroid? (Don't include .fits):")
        arr_detect = pickle.load(open(file1 + '-' + file2 + '_detlim.p', 'rb'))
        arr_flux = np.logspace(-3, -1, num = 15, endpoint = True)
        arr_rad = np.arange(50)+1
        plt.imshow(arr_detect, extent = [np.min(arr_flux), np.max(arr_flux), np.max(arr_rad), np.min(arr_rad)], aspect = 'auto', cmap = 'gray')
        plt.xscale('log')
        plt.xlabel('flux ratio between binary and main star')
        plt.ylabel('radius (pixels)')
        plt.show()


def run_detlim():
        arr_cent1 = np.arange(total_sets)+1
        arr_cent2 = np.arange(total_sets)+1
        for file1 in arr_cent1:
                for file2 in arr_cent2:
                        if file1 != file2:
                                print 'centroid_' + str(file1) + '-' + 'centroid_' + str(file2)
                                psf_subtract(str(file1), str(file2))




######################################################################################
#################'''T'was a good effort! (Useless code)'''############################
######################################################################################

def detect_maxima():
#checks image for points of maxima, for potential companions
	file1 = raw_input("Filename? (Don't include .fits):")
	img = pf.open(file1 + '.fits')[0].data
	x_length = img[1].size
	y_length = img[0].size
	print 'x_length', x_length
	print 'y_length', y_length
	median_img = np.median(img)
	arr_bool = np.zeros([y_length, x_length])
	arr_bool.fill(0)
	numpixels = 3

	for i in np.arange(numpixels)+1:
		img_chunk = img[i:-i,i:-i]
		print i
		for index in np.arange(img_chunk.size):
			index = np.unravel_index(index, img_chunk.shape)
			index = np.array(index)
			index += i
			#print 'index', index
			#print 'img.size', img.size
			img_check = img[index[0]-i:index[0]+i+1,index[1]-i:index[1]+i+1]
			#print 'img_check.size', img_check.size
			img_check = np.delete(img_check, (((((2*i)+1))**2)-1)/2, None)
			if all(elem < img[index[0]][index[1]] for elem in img_check):
				if img[index[0]][index[1]] > median_img:
					arr_bool[index[0]][index[1]] += 1
	condition = np.argwhere(arr_bool < numpixels)
	condition2 = np.argwhere(arr_bool == numpixels)
	for item in condition:
		arr_bool[item[0]][item[1]] = 0
	for item2 in condition2:
		arr_bool[item2[0]][item2[1]] = 1
	hdu = pf.PrimaryHDU(arr_bool)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto(file1 + '_test_maxima.fits', clobber=True)


def test_gaussian():
	file1 = raw_input("Centroid? (Don't include .fits):")
	img1 = pf.open('centroid_'+file1+'.fits')[0].data
	x_length = img1.shape[1]
	y_length = img1.shape[0]
	arr = np.zeros([y_length,x_length])
	arr[(y_length-1)/2, (x_length-1)/2] = 1.
	arr_RSS = np.array([])
	n_elem_sd = 50
	range_sd = np.linspace(0.001, 0.4, n_elem_sd)
	for index_sd in range(n_elem_sd):
		arr_gaussian = []
		for y in np.arange(y_length) - ((y_length-1)/2):
			y_index = y + ((y_length-1)/2)
			for x in np.arange(x_length) - ((x_length-1)/2):
				x_index = x +((x_length-1)/2)
				val = gaussian(y,x,range_sd[index_sd])			
				arr_gaussian.append(val)
		arr_gaussian = np.array(arr_gaussian)
		arr_gaussian = np.reshape(arr_gaussian, [y_length, x_length])
		arr_residuals = img1 - arr_gaussian
		arr_sq_res = arr_residuals**2.
		arr_RSS = np.append(arr_RSS, np.sum(arr_sq_res))
	#arr_RSS = np.reshape(arr_RSS, [y_length, x_length])
	print arr_RSS.size
	min_RSS_index = np.argmin(arr_RSS)
	min_sd = range_sd[min_RSS_index]
	print min_RSS_index
        plt.imshow(img1)


def test_check_binary():
        img = np.zeros([161, 161])
        img.fill(0)
        y_length = img.shape[0]
	x_length = img.shape[1]
	center_index = [(y_length-1)/2, (x_length-1)/2]

	center = np.array(center_index) ###_value at center_###
	arr_index = []

        counter = 0
        for element in np.arange(0, 80+1):
                ringvals, ringindexes = return_ring(img, element)
                #print 'len(ringindexes)', len(ringindexes)
                print len(ringindexes) #*(360-80)/360.
                for r_index in ringindexes:
                        img[r_index[0], r_index[1]] = 1
                hdu = pf.PrimaryHDU(img)
                hdulist = pf.HDUList([hdu])
                hdulist.writeto('test_check_binary.fits', clobber=True)
                subprocess.call('/home/apps/ds9/ds9 test_check_binary.fits', shell = True)
                
        '''
	for y_index in np.arange(y_length):
		for x_index in np.arange(x_length):
                        if (radi-0.5) < distance(center_index[0], y_index, center_index[1], x_index) < (radi+0.5): #creating ring of radius: radi
				arr_index.append([y_index, x_index]) #appending indexes of points in ring
        print len(arr_index)
	for element in arr_index:
                img[element[0]][element[1]] = 1

        hdu = pf.PrimaryHDU(img)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto('test_ring.fits', clobber=True)
        subprocess.call('/home/apps/ds9/ds9 test_ring.fits', shell = True)
        '''
        
def open_fits(filename):
        #filename = raw_input('Name of fits file(Do not include .fits):')
        subprocess.call('/home/apps/ds9/ds9 ' + filename, shell = True)

def open_fits_center(filename):
                img = afits.open(filename)[0].data
                img_center = img[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index]
                img_center /= np.median(img_center)
                hdu = afits.PrimaryHDU(img_center)
                hdulist = afits.HDUList([hdu])
                hdulist.writeto('test.fits', clobber = True)
                subprocess.call('/home/apps/ds9/ds9 ' + 'test.fits', shell = True)

                                                        

def openview_pickle(filename):
        img_final = pickle.load(open(filename, "rb"))
        hdu = pf.PrimaryHDU(img_final)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto('test.fits', clobber = True)
        subprocess.call('/home/apps/ds9/ds9 ' + 'test.fits', shell = True)

def open_fits_taurus(filename):
        subprocess.call(ds9_taurus + filename, shell = True)
