import pyfits as pf
import numpy as np
import matplotlib as mpl
import time, datetime
import matplotlib.pyplot as plt
import pickle
import subprocess
import math
import scipy
from scipy.signal import medfilt
from scipy.optimize import curve_fit
from scipy.ndimage.interpolation import shift
from scipy.ndimage.fourier import fourier_shift as fshift
from os.path import exists
import threading
from multiprocessing import Pool
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
date = 'Oct16'
arr_date_ShaneAO = ['Jun15', 'Nov15', 'Mar16', 'Sep16', 'Oct16']
struct_ShaneAO_total_sets = {'Jun15':28, 'Nov15':25, 'Mar16':25, 'Sep16':19, 'Oct16':5}
total_sets = struct_ShaneAO_total_sets[date]

#Showing date when script is run
print '------'
print 'Working with', date, 'files'
print '------'

#Set parameters of region where light hits telescope
arr_region_bounds = open(directory + '/' + date + '/telescope_region.txt', 'rb').read().splitlines()
region_x_min_index = int(arr_region_bounds[0])
region_x_max_index = int(arr_region_bounds[1])
region_y_min_index = int(arr_region_bounds[2])
region_y_max_index = int(arr_region_bounds[3])

#'''
print region_x_min_index
print region_x_max_index
print region_y_min_index
print region_y_max_index
#'''

#------
#***No longer needed, can delete.
#Correlation matrix parameters:
#------
max_numb_imgs = 5000

#------------------
#General guidelines when new data is added
#------------------
#make flat.txt and dark.txt files for dark and flat file numbers
#create_txt() to create txt files with info of different stars
#Change date arrays & dictionary of set numbers above
#Use check_img_center() and change the parameters where light hits telescope
#make telescope_region.txt in date folder e.g. ShaneAO/Jun15/telescope_region.txt
#run_create_sky() creates the sky imgs for all sets for date.
#create_dark()
#print_filters to view available filters for all imgs
#sort_filters create pickle files
#create_flats() to median flats and create master flat
#cleanup_flat() to fix broken pixels
#run_create_reduced_imgs() to subtract imgs by sky img & divide imgs by flat
#run_center_imgs()
#run_create_filts() to perform high-pass filter
#run_create_centroid_filts()
#------------------
#to be revised:
#correlation_matrix2 to create parts ofcorrelation matrix
#merge_mat() to merge
#------------------


def create_info(date = date):
        filename_setinfo = directory + '/' + date + '/' + 'set_' + str(int(setnumb1)) +'_info.txt'
        with open(filename_setinfo, 'rb') as file:
                arr_info = file.read().splitlines()
                for row in arr_info:
                        row = row.rstrip().lstrip()
                        if not row:
                                arr_info.remove(row)
        arr_info = filter(lambda a: a != filtername, arr_info)
        arr_info.append(filtername)
        with open(filename_setinfo, 'w') as newfile:
                for item in arr_info:
                        newfile.write("%s\n" % item)
        print 'updated:', filename_setinfo

def check_img_center():
        date = raw_input('Date:')
        imgnumb = raw_input('Enter 4 digit img number:')
        outputfilename = 'test_imgcenter.fits'
        img = pf.open('ShaneAO/' + date + '/' + 's'+ imgnumb + '.fits')[0].data
        img_center = img[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index]
        print 'img shape', img_center.shape
        hdu = pf.PrimaryHDU(img_center)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(outputfilename ,clobber=True)
        subprocess.call('/home/apps/ds9/ds9 ' + outputfilename, shell = True)

def create_txt(date = date):
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
        print 'dealing with', date
	arr_filternames = []
	#for i in np.arange(1, max_numb_imgs): #CHANGE ACCORDINGLY
        for i in np.arange(119, 218): #CHANGE ACCORDINGLY

                #load img file
                filename = directory + '/' + date + '/' + 's' + str(i).zfill(4) + '.fits'
                if exists(filename):
                        fits = pf.open(filename)
                else:
                        continue
		hdr = fits[0].header

		try:
			filtername = hdr['filt1nam']
		except:
                        print '------'
                        print "img numb", i, "does not have the tag 'FILT1NAM'"
			print hdr.keys()
			print 'number of header tags:', len(hdr.keys())
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

                #------
                #load img files
                #------
                filename = directory + '/' + date + '/' + 's' + str(i).zfill(4) + '.fits'
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
        if arr_ks_flat.size:
                save_fits(arr_ks_flat, directory + '/' + date + '/' + 'arr_ks_flat.fits')
                pickle.dump(arr_ks_flat, open(directory + '/' + date + '/' + 'arr_ks_flat.p', 'wb'))
                print "saved ks flats' img numbers as fits and pickle files"
        if arr_brg_flat.size:
                save_fits(arr_brg_flat, directory + '/' + date + '/' + 'arr_brg_flat.fits')
                pickle.dump(arr_brg_flat, open(directory + '/' + date + '/' + 'arr_brg_flat.p', 'wb'))                
                print "saved brg flats' img numbers"
        if arr_j_flat.size:
                save_fits(arr_j_flat, directory + '/' + date + '/' + 'arr_j_flat.fits')
                pickle.dump(arr_j_flat, open(directory + '/' + date + '/' + 'arr_j_flat.p', 'wb'))
                "saved j flats' img numbers"
                

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
                filename = directory + '/' + date + '/' + 's' + str(i).zfill(4) + '.fits'
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

                #------
                #load img file
                #------
                filename = directory + '/' + date + '/' + 's' + str(i).zfill(4) + '.fits'
                if exists(filename):
                        fits = pf.open(filename, memmap = False)
                        img = fits[0].data
                else:
                        print filename, 'doesnt exist'
                        continue

                img -= img_dark #subtract dark
                img_small = img[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index] #crop img
                mdn_imgsmall = np.median(img_small) #stack
                img /= mdn_imgsmall #normalize by median pixel
                arr_mdn.append(img)
		del img
		del fits
        arr_output = np.median(np.array(arr_mdn), axis = 0)
        return arr_output


def create_flats():

        
        filename_arrksflat_pickle = directory + '/' + date + '/' + "arr_ks_flat.p"
        filename_arrksflat_fits = directory + '/' + date + '/' + "arr_ks_flat.fits"

        filename_arrbrgflat_pickle = directory + '/' + date + '/' + "arr_brg_flat.p"
        filename_arrbrgflat_fits = directory + '/' + date + '/' + "arr_brg_flat.fits"

        filename_arrjflat_pickle = directory + '/' + date + '/' + "arr_h_flat.p"
        filename_arrjflat_fits = directory + '/' + date + '/' + "arr_j_flat.fits"


        arr_ks_flat = np.array([])
        if exists(filename_arrksflat_fits):
                arr_ks_flat = pf.open(filename_arrksflat_fits)[0].data
        elif exists(filename_arrksflat_pickle):
                arr_ks_flat = pickle.load(open(filename_arrksflat_pickle, "rb" ))
        if arr_ks_flat.size:
                try:
                        img_ks_flat = norm_median_subdark(arr_ks_flat)
                        hdu = pf.PrimaryHDU(img_ks_flat)
                        hdulist = pf.HDUList([hdu])
                        hdulist.writeto(directory + '/' + date + '/' + 'img_flat_ks.fits',clobber=True)
                        print 'img_flat_ks.fits created'
                except:
                        print 'error with img_flat_ks'


        arr_brg_flat = np.array([])
        if exists(filename_arrbrgflat_fits):
                arr_brg_flat = pf.open(filename_arrbrgflat_fits)[0].data
        elif exists(filename_arrbrgflat_pickle):
                arr_brg_flat = pickle.load(open(filename_arrbrgflat_pickle, "rb" ))
        if arr_brg_flat.size:
                try:
                        img_brg_flat = norm_median_subdark(arr_brg_flat)
                        hdu = pf.PrimaryHDU(img_brg_flat)
                        hdulist = pf.HDUList([hdu])
                        hdulist.writeto(directory + '/' + date + '/' + 'img_flat_brg.fits',clobber=True)
                        print 'img_flat_brg.fits created'                
                except:
                        print 'error with img_flat_brg'                

        arr_j_flat = np.array([])
        if exists(filename_arrjflat_fits):
                arr_j_flat = pf.open(filename_arrjflat_fits)[0].data
        elif exists(filename_arrjflat_pickle):
                arr_j_flat = pickle.load(open(filename_arrjflat_pickle, "rb" ))
        if arr_j_flat.size:
                try:
                        img_j_flat = norm_median_subdark(arr_j_flat)
                        hdu = pf.PrimaryHDU(img_j_flat)
                        hdulist = pf.HDUList([hdu])
                        hdulist.writeto(directory + '/' + date + '/' + 'img_flat_j.fits',clobber=True)
                        print 'img_flat_j.fits created'                
                except:
                        print 'error with img_flat_j'                


def run_create_filts():
        #Create high-pass filtered images for all recorded stars
        #------

        len_filt_box = 21

        for date_science in arr_date_ShaneAO: #loop through dates for ShaneAO runs
                total_setnumb = struct_ShaneAO_total_sets[date_science] #Total number of "sets" of recorded stars for date
                print 'Date:', date_science
                for setnumb in np.arange(1,total_setnumb+1).astype(int): #loop through sets of stars
                        print 'working with set', setnumb
                        #------
                        # open up txt file for set, add all img numbers for that set into arr_targpsf
                        #------
                        setnumb = str(setnumb)
                        arr_startend = open(directory + '/' + date_science + '/set_'+ setnumb +'_info.txt', 'rb').read().splitlines()
                        start = int(arr_startend[0])
                        end = int(arr_startend[1])
                        arr_targ = np.arange(start, end+1)


                        for img_numb in arr_targ:
                                filename_init =  directory + '/' + date_science  + '/s' + str(int(img_numb)).zfill(4) + '_reduced_centered.fits'
                                filename_output  = directory + '/' + date_science  + '/s' + str(int(img_numb)).zfill(4) + '_RCF.fits'

                                if exists(filename_output):
                                        print filename_output, 'exists'
                                        continue
                                elif exists(filename_init):
                                        img_targ = pf.open(filename_init)[0].data
                                else:
                                        print img_numb, "doesn't exist"
                                        continue
                                        
                                img_targ -= medfilt(img_targ, [len_filt_box, len_filt_box])
                                save_fits(img_targ, filename_output)
                                print 'saved:', filename_output



def run_create_centroid_filts(total_sets = total_sets):
# Create medianed images for all sets up to total_sets
# Runs for LOCI-processed imgs, with high pass filtering
#------
        print 'Date:', date
        for setnumb in np.arange(1,total_sets+1).astype(int): #loop through sets of stars
                print 'working with set', setnumb
                #------
                # open up txt file for set, add all img numbers for that set into arr_targpsf
                #------
                setnumb = str(setnumb)
                arr_startend = open(directory + '/' + date + '/set_'+ setnumb +'_info.txt', 'rb').read().splitlines()
                start = int(arr_startend[0])
                end = int(arr_startend[1])
                arr_targ = np.arange(start, end+1)
                
                filename_output = directory + '/' + date + '/' + 'centroid_filt_' + setnumb + '.fits'
                arr_mdn = median_it_filts(arr_targ)
                
                hdu = pf.PrimaryHDU(arr_mdn)
                hdulist = pf.HDUList([hdu])
                hdulist.writeto(filename_output, clobber=True)
                print 'created:', filename_output
                

def median_it_filts(arr_filenumbs):
        arr_mdn = []
        print arr_filenumbs
        for i in arr_filenumbs:
		#print i
                if i > 999:
			filename = directory + '/' + date + '/' + 's' + str(i) + '_RCF.fits'
			if exists(filename):
				#with pf.open(filename, memmap = False) as fits:
				fits = pf.open(filename, memmap = False)
				img = fits[0].data
			else:
				print filename, 'doesnt exist'
		else:
			filename = directory + '/' + date + '/' + 's0' + str(i) + '_RCF.fits'
			if exists(filename):
				#with pf.open(filename, memmap = False) as fits:
				fits = pf.open(filename, memmap = False)
				img = fits[0].data
			else:
				print filename, 'doesnt exist'
                if 'img' in locals():
                        arr_mdn.append(img)
                        del img
                        del fits

        arr_output = np.median(np.array(arr_mdn), axis = 0)
        return arr_output


def rename_img_subed_stacked_faint(setnumb):
        
        arr_fluxratio = pf.open('set' + str(int(setnumb)) + '_fluxratio_bins_init.fits')[0].data
        mag = -2.5*np.log10(arr_fluxratio)
        print mag
        arr_radius = pf.open(directory + '/' + date + '/' + 'arr_radiusmock.fits')[0].data
        theta = 0
        for index_rad in np.arange(arr_radius.size).astype(int):
                radius = arr_radius[index_rad]
                filename_init = '/home/jjoon/faint_mocks_real/' + 'set'+ str(int(setnumb)) + '_mockrad' + str(int(radius)) + '_theta' + str(int(theta)) + 'locifaintmockfiltfinal.fits'
                filename_out_1 = 'home/jjoon/faint_mocks_real/' + 'set'+ str(int(setnumb)) + '_mockrad' + str(int(radius)) + '_mag'
                filename_out_2 = '%.2f_locifaintmockfiltfinal.fits' % mag[index_rad]
                filename_out = filename_out_1 + filename_out_2
                if exists(filename_init):
                        print 'exists'
                subprocess.call('mv -f ' + filename_init + ' ' + filename_out, shell =True) ###
                print 'created:', filename_out


def create_stacked_mockinit():
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80]).astype(int)
        print 'arr_radiusmock:', arr_radiusmock
        
        arr_theta = np.arange(0, 360, 60).astype(int) 
        print 'arr_theta:', arr_theta

        for setnumb1 in np.arange(1, struct_ShaneAO_total_sets[date]+1).astype(int):
                print '------'
                print 'set:', setnumb1
                print '------'

                #------
                # For set of target imgs,
                # read txt file in directory to know which files to use
                #------
                setnumb1 = str(setnumb1)
                arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
                start1 = int(arr_startend1[0])
                end1 = int(arr_startend1[1])
                arr_targpsf = np.arange(start1, end1+1)

                for radius_mock in arr_radiusmock:
                        for theta_mock in arr_theta:
                                filename_mockinit_output = directory + '/' + date + '/' + str(int(setnumb1)) + '_radmock' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + '_mockinit.fits'
                                arr_imgs = []
                                for imgnumb in arr_targpsf:
                                        filename_mockinit = directory + '/' + date + '/' + str(int(imgnumb)) + '_radmock' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + '_mockinit.fits'
                                        img = pf_loadfits(filename_mockinit, True)
                                        if img.size:
                                                arr_imgs.append(img)
                                img_output = np.median(np.array(arr_imgs), axis = 0)
                                save_fits(img_output, filename_mockinit_output)


def get_mock_shifts(bin_cond = True, faint_cond = True, replace_cond = True):
        # Get x and y shifts of mock binaries before and after
        # Compare for position angles, radius
        #------
        index_date_targ = arr_date_ShaneAO.index(date) #index of date of target psf
        radi_apert = 3
        detlim_mult_factor = 10.
        
        max_rad_img = 80
        len_filt_box = 21 #for high-pass filter, take median of square with this length around points
        
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80]).astype(int)
        print 'arr_radiusmock:', arr_radiusmock
        
        arr_theta = np.arange(0, 360, 60).astype(int) 
        print 'arr_theta:', arr_theta


        for setnumb1 in np.arange(1, struct_ShaneAO_total_sets[date]+1).astype(int):
                print '------'
                print 'set:', setnumb1
                print_timenow()

                #------
                # For set of target imgs,
                # read txt file in directory to know which files to use
                #------
                setnumb1 = str(setnumb1)
                arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
                start1 = int(arr_startend1[0])
                end1 = int(arr_startend1[1])
                arr_targpsf = np.arange(start1, end1+1)
                
                #-----
                # Open stacked, cleaned image for set
                # Record dimensions of array
                # Record total number of elements in array
                #------
                filename_fits = directory + '/' + date + '/' + 'centroid_' + str(int(setnumb1))+'.fits'
                if exists(filename_fits):
                        img_mdn = pf.open(filename_fits)[0].data
                else:
                        print filename_fits, 'doesnt exist. Skipping set'
                        continue

                #------
                #Iterate through radii for mock binaries
                #Check LOCI-subtracted img for 5 sigma value at that radius
                #Create array of corresponding flux ratio
                #------
                print 'getting fluxes of stacked img before LOCI'
                print_timenow()
                arr_fluxratio = []
                for rad_mock in arr_radiusmock:
                        struct_ring = check_ring(img_mdn, rad_mock)
                        sd_ring = struct_ring['sd_ring']
                        arr_fluxratio.append(5*sd_ring*detlim_mult_factor)
                print 'done getting fluxes'
                print_timenow()

                #------
                # Loop through each radius
                #-----
                arr_y_com_init = []
                arr_x_com_init = []
                arr_y_com_filt = []
                arr_x_com_filt = []
                arr_y_com_final = []
                arr_x_com_final = []

                arr_y_com_init_std = []
                arr_x_com_init_std = []
                arr_y_com_filt_std = []
                arr_x_com_filt_std = []
                arr_y_com_final_std = []
                arr_x_com_final_std = []
                for index_rad in np.arange(arr_radiusmock.size).astype(int):
                        
                        
                        radius_mock = arr_radiusmock[index_rad]
                        print '------'
                        print 'radius:', radius_mock
                        print_timenow()
                        arr_y_com_init_thetas = []
                        arr_x_com_init_thetas = []
                        arr_y_com_filt_thetas = []
                        arr_x_com_filt_thetas = []
                        arr_y_com_final_thetas = []
                        arr_x_com_final_thetas = []

                        arr_y_com_init_thetas_std = []
                        arr_x_com_init_thetas_std = []
                        arr_y_com_filt_thetas_std = []
                        arr_x_com_filt_thetas_std = []
                        arr_y_com_final_thetas_std = []
                        arr_x_com_final_thetas_std = []

                        for index_theta in np.arange(arr_theta.size).astype(int):
                                theta_mock = arr_theta[index_theta]
                                arr_y_com_init_frames = np.array([])
                                arr_x_com_init_frames = np.array([])
                                arr_y_com_filt_frames = np.array([])
                                arr_x_com_filt_frames = np.array([])
                                arr_y_com_final_frames = np.array([])
                                arr_x_com_final_frames = np.array([])
                                '''
                                filename_final = directory + '/' + date  + '/' + 'set' + str(int(setnumb1)) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfiltfinal.fits'
                                try:
                                        img_final = pf.open(filename_imgfinal)[0].data
                                except:
                                        continue
                                '''
                                #'''
                                for index_imgnumb in np.arange(arr_targpsf.size).astype(int):
                                        imgnumb = arr_targpsf[index_imgnumb]
                                        percentage_imgnumb = (index_imgnumb+1.)*100/arr_targpsf.size
                                        if not percentage_imgnumb%10:
                                                print 'radius:', radius_mock, ',', 'theta:', theta_mock, ',',
                                                print percentage_imgnumb, '% done with frames in set'
                                        #------
                                        #load target img, divide by max pixel value
                                        #------
                                        filename = directory + '/' + date  + '/s' + str(int(imgnumb)).zfill(4) + '_reduced_centered.fits'
                                        if exists(filename):
                                                try:
                                                        img1 = pf.open(filename)[0].data
                                                        max_pix_primary = np.amax(img1)
                                                        img1 /= max_pix_primary
                                                except:
                                                        print 'error opening file:', filename
                                                        continue
                                        else:
                                                #print filename, 'doesnt exist'
                                                continue
                                        
                                        y_length, x_length = img1.shape
                                        y_index_center = int((y_length - 1)/2)
                                        x_index_center = int((x_length - 1)/2)
                                
                                        #------
                                        #define filename for final img
                                        #------
                                        filename_imgfinal = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfiltfinal.fits'
                                        try:
                                                img_final = pf.open(filename_imgfinal)[0].data
                                        except:
                                                continue

                                        filename_mockinit = directory + '/' + date + '/' + str(int(imgnumb)) + '_radmock' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + '_mockinit.fits'
                                        filename_mockfilt = directory + '/' + date + '/' + str(int(imgnumb)) + '_radmock' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + '_mockfilt.fits'

                                        #----------
                                        #load mock binary and divide by predetermined factor
                                        #***Change centroid filename if necessary***
                                        #----------
                                        img_mock = pf.open(directory + '/' + 'Jun15' + '/' + 'centroid_1.fits')[0].data
                                        img_mock /= np.amax(img_mock)
                                        img_mock *= float(arr_fluxratio[index_rad])
                                        
                                        #--------
                                        #Add mock binary at correct radius & angle, save img as fits file
                                        #--------
                                        dx_mock = radius_mock*np.cos(theta_mock*(2*np.pi)/360) #calculating y displacement from center
                                        dx_mock = int(round(dx_mock))
                                        dy_mock = radius_mock*np.sin(theta_mock*(2*np.pi)/360) #calculating x displacement from center
                                        dy_mock = int(round(dy_mock))
                                        
                                        img_mock = np.fft.fft2(img_mock) #Fourier shift to construct binary for adding to primary
                                        img_mock = fshift(img_mock, [dy_mock, dx_mock])
                                        img_mock = np.fft.ifft2(img_mock)
                                        img_mock = np.real(img_mock)
                                        
                                        # Get center of mass after LOCI
                                        y_com_final, x_com_final = get_com_aperture(img_final, (dy_mock + y_index_center , dx_mock + x_index_center), radi_apert)
                                        arr_y_com_final_frames = np.append(arr_y_com_final_frames, y_com_final)
                                        arr_x_com_final_frames = np.append(arr_x_com_final_frames, x_com_final)
                                        '''
                                        if math.isnan(y_com_final) or math.isnan(x_com_final):
                                                print 'y_com_final', y_com_final
                                                print 'x_com_final', x_com_final
                                        '''
                                        if exists(filename_mockinit):
                                                img1 = pf.open(filename_mockinit)[0].data
                                        else:
                                                img1 += img_mock #Add binary to primary

                                        # Get center of mass before high pass filter
                                        y_com_init, x_com_init = get_com_aperture(img1, (dy_mock + y_index_center , dx_mock + x_index_center), radi_apert)
                                        if not exists(filename_mockinit):
                                                save_fits(img1, filename_mockinit)

                                        arr_y_com_init_frames = np.append(arr_y_com_init_frames, y_com_init)
                                        arr_x_com_init_frames = np.append(arr_x_com_init_frames, x_com_init)
                                        '''
                                        if math.isnan(y_com_init) or math.isnan(x_com_init):
                                                print 'y_com_init', y_com_init
                                                print 'x_com_init', x_com_init
                                        '''

                                        if exists(filename_mockfilt):
                                                img1 = pf.open(filename_mockfilt)[0].data
                                        else:                                
                                                img1 *= max_pix_primary  # Multiply by primary max pixel value from before
                                                img1 -= medfilt(img1, [len_filt_box, len_filt_box]) # Perform high-pass filter on img
                                                img1 /= np.amax(img1) # Normalize img using new max pixel value
                      
                                        # Get center of mass after high pass filter
                                        y_com_filt, x_com_filt = get_com_aperture(img1, (dy_mock + y_index_center , dx_mock + x_index_center), radi_apert)
                                        '''
                                        if math.isnan(y_com_filt) or math.isnan(x_com_filt):
                                                print 'y_com_filt', y_com_filt
                                                print 'x_com_filt', x_com_filt
                                        '''
                                        if not exists(filename_mockfilt):
                                                save_fits(img1, filename_mockfilt)

                                        arr_y_com_filt_frames = np.append(arr_y_com_filt_frames, y_com_filt)
                                        arr_x_com_filt_frames = np.append(arr_x_com_filt_frames, x_com_filt)
          

                                #keep std of frames too
                                #std of mean = std/sqrt{n}
                                arr_y_com_init_thetas_std.append(np.std(arr_y_com_init_frames))
                                arr_x_com_init_thetas_std.append(np.std(arr_x_com_init_frames))
                                arr_y_com_filt_thetas_std.append(np.std(arr_y_com_filt_frames))
                                arr_x_com_filt_thetas_std.append(np.std(arr_x_com_filt_frames))
                                arr_y_com_final_thetas_std.append(np.std(arr_y_com_final_frames))
                                arr_x_com_final_thetas_std.append(np.std(arr_x_com_final_frames))

                                arr_y_com_init_thetas.append(np.mean(arr_y_com_init_frames))
                                arr_x_com_init_thetas.append(np.mean(arr_x_com_init_frames))
                                arr_y_com_filt_thetas.append(np.mean(arr_y_com_filt_frames))
                                arr_x_com_filt_thetas.append(np.mean(arr_x_com_filt_frames))
                                arr_y_com_final_thetas.append(np.mean(arr_y_com_final_frames))
                                arr_x_com_final_thetas.append(np.mean(arr_x_com_final_frames))
                                #'''
                        arr_y_com_init.append(arr_y_com_init_thetas)
                        arr_x_com_init.append(arr_x_com_init_thetas)
                        arr_y_com_filt.append(arr_y_com_filt_thetas)
                        arr_x_com_filt.append(arr_x_com_filt_thetas)
                        arr_y_com_final.append(arr_y_com_final_thetas)
                        arr_x_com_final.append(arr_x_com_final_thetas)

                        arr_y_com_init_std.append(arr_y_com_init_thetas_std)
                        arr_x_com_init_std.append(arr_x_com_init_thetas_std)
                        arr_y_com_filt_std.append(arr_y_com_filt_thetas_std)
                        arr_x_com_filt_std.append(arr_x_com_filt_thetas_std)
                        arr_y_com_final_std.append(arr_y_com_final_thetas_std)
                        arr_x_com_final_std.append(arr_x_com_final_thetas_std)

                '''
                print len(arr_y_com_init)
                print len(arr_x_com_init)
                print len(arr_y_com_filt)
                print len(arr_x_com_filt)
                '''
                arr_y_com_init = np.array(arr_y_com_init)
                arr_x_com_init = np.array(arr_x_com_init)
                arr_y_com_filt = np.array(arr_y_com_filt)
                arr_x_com_filt = np.array(arr_x_com_filt)
                arr_y_com_final = np.array(arr_y_com_final)
                arr_x_com_final = np.array(arr_x_com_final)

                arr_y_com_init_std = np.array(arr_y_com_init_std)
                arr_x_com_init_std = np.array(arr_x_com_init_std)
                arr_y_com_filt_std = np.array(arr_y_com_filt_std)
                arr_x_com_filt_std = np.array(arr_x_com_filt_std)
                arr_y_com_final_std = np.array(arr_y_com_final_std)
                arr_x_com_final_std = np.array(arr_x_com_final_std)
                
                foldername = 'com/'
                filename_y_init = foldername + 'arr_y_com_init_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_y_com_init, filename_y_init)
                filename_x_init = foldername + 'arr_x_com_init_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_x_com_init, filename_x_init)
                filename_y_filt = foldername + 'arr_y_com_filt_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_y_com_filt, filename_y_filt)
                filename_x_filt = foldername + 'arr_x_com_filt_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_x_com_filt, filename_x_filt)
                filename_y_final = foldername + 'arr_y_com_final_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_y_com_final, filename_y_final)
                filename_x_final = foldername + 'arr_x_com_final_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_x_com_final, filename_x_final)

                filename_y_init_std = foldername + 'arr_y_com_init_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_y_com_init_std, filename_y_init_std)
                filename_x_init_std = foldername + 'arr_x_com_init_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_x_com_init_std, filename_x_init_std)
                filename_y_filt_std = foldername + 'arr_y_com_filt_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_y_com_filt_std, filename_y_filt_std)
                filename_x_filt_std = foldername + 'arr_x_com_filt_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_x_com_filt_std, filename_x_filt_std)
                filename_y_final_std = foldername + 'arr_y_com_final_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_y_com_final_std, filename_y_final_std)
                filename_x_final_std = foldername + 'arr_x_com_final_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                save_fits(arr_x_com_final_std, filename_x_final_std)


def plot_histo_com():
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80]).astype(int)
        print 'arr_radiusmock:', arr_radiusmock
        
        arr_theta = np.arange(0, 360, 60).astype(int) 
        print 'arr_theta:', arr_theta

        index_center_img = np.array([80,80])

        for index_rad in np.arange(arr_radiusmock.size).astype(int):
                arr_r_filt = np.array([])
                arr_angle_filt = np.array([])
                arr_r_final = np.array([])
                arr_angle_final = np.array([])
                for setnumb1 in np.arange(1, struct_ShaneAO_total_sets[date]+1).astype(int):
                        setnumb1 = str(int(setnumb1))
                        foldername = 'com/'
                        filename_y_init = foldername + 'arr_y_com_init_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        filename_x_init = foldername + 'arr_x_com_init_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        filename_y_filt = foldername + 'arr_y_com_filt_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        filename_x_filt = foldername + 'arr_x_com_filt_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        filename_y_final = foldername + 'arr_y_com_final_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        filename_x_final = foldername + 'arr_x_com_final_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        filename_y_init_std = foldername + 'arr_y_com_init_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        filename_x_init_std = foldername + 'arr_x_com_init_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        filename_y_filt_std = foldername + 'arr_y_com_filt_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        filename_x_filt_std = foldername + 'arr_x_com_filt_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        filename_y_final_std = foldername + 'arr_y_com_final_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        filename_x_final_std = foldername + 'arr_x_com_final_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'

                        #load files for com
                        if exists(filename_y_init) and exists(filename_x_init):
                                arr_y_init_setnumb = pf.open(filename_y_init)[0].data
                                arr_x_init_setnumb = pf.open(filename_x_init)[0].data
                                #print 'arr_y_init_setnumb', arr_y_init_setnumb
                                #print 'arr_x_init_setnumb', arr_x_init_setnumb
                                #useless = raw_input('press somethign to continue')
                        else:
                                print filename_y_init, 'or', filename_x_init, 'doesnt exist'
                                continue
                        arr_y_init, arr_x_init = arr_y_init_setnumb[index_rad, :], arr_x_init_setnumb[index_rad, :]
                        arr_y_init -= index_center_img[0]
                        arr_x_init -= index_center_img[1]
                        r_init = np.sqrt(arr_y_init**2. + arr_x_init**2.)
                        angle_init = np.arctan2(arr_y_init, arr_x_init)
                        
                        if exists(filename_y_filt) and exists(filename_x_filt):
                                arr_y_filt_setnumb = pf.open(filename_y_filt)[0].data
                                arr_x_filt_setnumb = pf.open(filename_x_filt)[0].data
                                if math.isnan(arr_y_filt_setnumb.flatten()[0]) or math.isnan(arr_x_filt_setnumb.flatten()[0]):
                                        continue
                        else:
                                print filename_y_filt, 'or', filename_x_filt, 'doesnt exist'
                                continue
                        arr_y_filt, arr_x_filt = arr_y_filt_setnumb[index_rad, :], arr_x_filt_setnumb[index_rad, :]
                        arr_y_filt -= index_center_img[0]
                        arr_x_filt -= index_center_img[1]
                        r_filt = np.sqrt(arr_y_filt**2. + arr_x_filt**2.)
                        angle_filt = np.arctan2(arr_y_filt, arr_x_filt)

                        if exists(filename_y_final) and exists(filename_x_final):
                                arr_y_final_setnumb = pf.open(filename_y_final)[0].data
                                arr_x_final_setnumb = pf.open(filename_x_final)[0].data
                                if math.isnan(arr_y_final_setnumb.flatten()[0]) or math.isnan(arr_x_final_setnumb.flatten()[0]):
                                        continue
                        else:
                                print filename_y_final, 'or', filename_x_final, 'doesnt exist'
                                continue
                        arr_y_final, arr_x_final = arr_y_final_setnumb[index_rad, :], arr_x_final_setnumb[index_rad, :]
                        arr_y_final -= index_center_img[0]
                        arr_x_final -= index_center_img[1]
                        r_final = np.sqrt(arr_y_final**2. + arr_x_final**2.)
                        angle_final = np.arctan2(arr_y_final, arr_x_final)

                        #subtract from initial com to get y and x shifts
                        r_shift_filt = r_filt - r_init
                        angle_shift_filt = angle_filt - angle_init
                        r_shift_final = r_final - r_init
                        angle_shift_final = angle_final - angle_init


                        max_rad_shift_histo = 10
                        arr_r_filt = np.append(arr_r_filt, r_shift_filt[np.where(np.abs(r_shift_filt)<max_rad_shift_histo)])
                        arr_r_final = np.append(arr_r_final, r_shift_final[np.where(np.abs(r_shift_final)<max_rad_shift_histo)])
                         
                        max_angle_shift_radians = 0.4
                        arr_angle_filt = np.append(arr_angle_filt, angle_shift_filt[np.where(np.abs(angle_shift_filt)<max_angle_shift_radians)])
                        arr_angle_final = np.append(arr_angle_final, angle_shift_final[np.where(np.abs(angle_shift_final)<max_angle_shift_radians)])

                        '''
                        #load files for std of com
                        arr_y_init_std_setnumb = pf_loadfits(filename_y_init_std, print_cond = True)
                        arr_x_init_std_setnumb = pf_loadfits(filename_x_init_std, print_cond = True)
                        if not arr_y_init_std_setnumb and not arr_x_init_std_setnumb:
                                continue                                
                        arr_y_init_std, arr_x_init_std = arr_y_init_std_setnumb[index_rad, :], arr_x_init_std_setnumb[index_rad, :]
                        r_init_std = np.sqrt(arr_y_init_std**2. + arr_x_init_std**2.)
                        angle_init_std = np.arctan2(arr_y_init_std, arr_x_init_std)

                        arr_y_filt_std_setnumb = pf_loadfits(filename_y_filt_std, print_cond = True)
                        arr_x_filt_std_setnumb = pf_loadfits(filename_x_filt_std, print_cond = True)
                        if not arr_y_filt_std_setnumb and not arr_x_filt_std_setnumb:
                                continue
                        arr_y_filt_std, arr_x_filt_std = arr_y_filt_std_setnumb[index_rad, :], arr_x_filt_std_setnumb[index_rad, :]
                        r_filt_std = np.sqrt(arr_y_filt_std**2. + arr_x_filt_std**2.)
                        angle_filt_std = np.arctan2(arr_y_filt_std, arr_x_filt_std)

                        arr_y_final_std_setnumb = pf_loadfits(filename_y_final_std, print_cond = True)
                        arr_x_final_std_setnumb = pf_loadfits(filename_x_final_std, print_cond = True)
                        if not arr_y_final_std_setnumb and not arr_x_final_std_setnumb:
                                continue
                        arr_y_final_std, arr_x_final_std = arr_y_final_std_setnumb[index_rad, :], arr_x_final_std_setnumb[index_rad, :]
                        r_final_std = np.sqrt(arr_y_final_std**2. + arr_x_final_std**2.)
                        angle_final_std = np.arctan2(arr_y_final_std, arr_x_final_std)
                        '''
                n_bins = 10
                '''
                print 'arr_r_filt', arr_r_filt
                print 'arr_r_final', arr_r_final
                '''
                '''
                plt.hist(arr_r_filt, bins=n_bins)
                plt.xlabel('Shift in radius (After high-pass filtering)')
                plt.show()
                '''
                plt.close()
                plt.hist(arr_r_final, bins = n_bins)
                plt.xlabel('Shift in radius (After LOCI)')
                plt.title('Radial shift for mocks at radius = '+str(int(arr_radiusmock[index_rad])))
                #plt.show()
                plt.savefig("/home/jjoon/rad"+str(int(arr_radiusmock[index_rad]))+".png", dpi = 200)

                plt.close()
                plt.hist(arr_angle_final, bins = n_bins)
                plt.xlabel('Shift in angle (After LOCI)')
                plt.title('Azimuth shift for mocks at radius = '+str(int(arr_radiusmock[index_rad])))
                plt.savefig("/home/jjoon/angle"+str(int(arr_radiusmock[index_rad]))+".png", dpi = 200)
                #plt.show()
                
                #useless = raw_input('enter key to continue')

def pf_loadfits(filename_img, print_cond = False):
        if exists(filename_img):
                try:
                        return pf.open(filename_img)[0].data
                except:
                        return np.array([])
        else:
                if print_cond:
                        print filename_img, 'doesnt exist'
                return np.array([])


def run_loci_mockbins(setnumb1, input_radiusmock = 0, bin_cond = True, faint_cond = False, replace_cond = True):
        # Run LOCI with mock binaries
        # Run with binaries if bin_cond = True
        # Run with replacement if replace_cond = True
        #------
        index_date_targ = arr_date_ShaneAO.index(date) #index of date of target psf
        fillelement = 10 #element for correlation matrix when no correlation is calculated

        #------
        # Open LOCI-subtracted, medianed img for set number
        # Record dimensions of array
        # Record total number of elements in array
        #------
        if bin_cond:
                filename_fits = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_locifiltfinal' + '.fits'
                img_subed_filt = pf.open(filename_fits)[0].data
                img1size = img_subed_filt.size
        else:
                filename_fits = directory + '/' + date + '/' + 'centroid_' + str(int(setnumb1)) + '.fits'
                img_temp = pf.open(filename_fits)[0].data
                img1size = img_temp.size
        

        #------
        # Define LOCI parameters 
        # Define high-pass filter box dimensions(width)
        #------
        if faint_cond:
                detlim_mult_factor = 2. #faint*** ###10. if not faint
        else:
                detlim_mult_factor = 10.
        max_rad_img = 80
        numb_imgs_keep = 200 #number of imgs to keep from correlation matrix
        fwhm_threshold = 10 #maximum FWHM to keep
	maxpix_threshold = 35000 #exclude imgs with max pixel value above this number
	minpix_threshold = 5000 #exclude imgs with min pixel value below this number
        len_filt_box = 21 #for high-pass filter, take median of square with this length around points

        #------
        #Parameters for subtraction/optimization radius for LOCI
        #------
        radius_sub = 6
        radius_op = 19
        
        if bin_cond:
                #------
                #Define mock-up binary arrays
                #------
                if faint_cond:
                        arr_radiusmock = np.array([input_radiusmock]) #faint***
                else:
                        arr_radiusmock = np.array([10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80]).astype(int)
                print 'arr_radiusmock:', arr_radiusmock
                if faint_cond:
                        arr_theta = np.array([0]) #faint***
                else:
                        arr_theta = np.arange(0, 360, 60).astype(int) 
                print 'arr_theta:', arr_theta
        else:
                arr_radiusmock = [0]
                arr_theta = [0]

        #------
        #Iterate through radii for mock binaries
        #Check LOCI-subtracted img for 5 sigma value at that radius
        #Create array of corresponding flux ratio
        #------
        if bin_cond:
                if faint_cond:
                        filename_detlim_final = date + '_set' + str(int(setnumb1)) + '_detlim_final.fits'
                        detlim_final = pf.open(filename_detlim_final)[0].data
                        plt.plot(detlim_final, 'mo-')
                        fluxratio_detlimfinal = (10.**(detlim_final/2.5))**(-1.)
                        arr_fluxratio = fluxratio_detlimfinal*detlim_mult_factor
                        #save_fits(arr_fluxratio, 'set' + str(int(setnumb1)) + '_fluxratio_bins_init.fits')
                else:
                        arr_fluxratio = []
                        for rad_mock in arr_radiusmock:
                                struct_ring = check_ring(img_subed_filt, rad_mock)
                                sd_ring = struct_ring['sd_ring']
                                arr_fluxratio.append(5*sd_ring*detlim_mult_factor)

        #------
        # For set of target imgs,
        # read txt file in directory to know which files to use
        #------
        setnumb1 = str(setnumb1)
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
        arr_targpsf = np.arange(start1, end1+1)
	#print 'arr_targpsf', arr_targpsf
        
        #------
        #check for filter of 1st img in set
        #------
	filter_of_set = ret_filtername(start1, date) 
        print 'filter of set:', filter_of_set

        #------
        #open txt file with list of known binaries
        #------
        arr_bin = open(directory + '/' + directory + '_clearbinaries.txt', 'rb').read().splitlines()
        arr_bin_dates = np.array(arr_bin[0::2])
        arr_bin_setnumbs = np.array(arr_bin[1::2])
        #print 'arr_bin_dates', arr_bin_dates
        #print 'arr_bin_setnumbs', arr_bin_setnumbs
        
        #------
        # Add reference filenumbs from all dates
        #------
        arr_substar = []
        for index_date in np.arange(len(arr_date_ShaneAO)):
                date_ref = arr_date_ShaneAO[index_date]
                total_sets = struct_ShaneAO_total_sets[date_ref]
                print '------'
                print 'total_sets', total_sets
                arr_substar_temp = np.array([]) #array for adding ref psf filenames
                for setnumb2 in np.arange(1,total_sets+1):
                        if setnumb2 == int(setnumb1) and index_date_targ == index_date:
                                print date_ref, 'set', setnumb2, 'skipped'
                                continue
                        elif str(int(setnumb2)) in arr_bin_setnumbs:
                                arr_index_setnumb_bin = np.where(arr_bin_setnumbs == str(int(setnumb2)))[0]
                                #print arr_index_setnumb_bin
                                '''
                                print '------'
                                print 'index_setnumb', index_setnumb
                                print 'set number', setnumb2
                                print 'corresponding date in txt file', arr_bin_dates[index_setnumb]
                                print 'date in loop', date_ref
                                '''
                                for index_setnumb_bin in arr_index_setnumb_bin:
                                        if arr_bin_dates[index_setnumb_bin] == date_ref:
                                                print date_ref, 'set', setnumb2, 'skipped since it is a known binary'
                                                continue
                        else:
                                ####date_ref and setnumb2
                                #print 'adding elements in set no.', setnumb2, 'for date', date_ref
                                arr_startend2 = open(directory + '/' + date_ref + '/' + 'set_' + str(setnumb2) +'_info.txt', 'rb').read().splitlines()
                                start2 = int(arr_startend2[0])
                                end2 = int(arr_startend2[1])
                                arr_imgnumbs = np.arange(start2, end2+1)
                                arr_substar_temp = np.append(arr_substar_temp,  arr_imgnumbs) #array of filenames for ref psf
                                arr_temp = np.zeros(arr_imgnumbs.size)
                                arr_temp.fill(setnumb2)
                arr_substar_temp = arr_substar_temp.astype(int)
                counter_substar_date = arr_substar_temp.size #keep track of imgs added with counter
                print 'number of eligible img files:', counter_substar_date
                arr_substar.append(arr_substar_temp)
        #print 'arr_substar', arr_substar


        #------
        # load reference imgs one by one
        # Check filter, FWHM, max/min pixel values, img size
        # Include only the ones within parameters set above at beginning of function
        #------
        arr_j = []
        arr_img2 = []
        for index_date in np.arange(len(arr_date_ShaneAO)): #FIXFIXFIXFIX
                date_ref = arr_date_ShaneAO[index_date]
                arr_j_temp = []
                arr_img2_temp = []
                arr_substar_temp = arr_substar[index_date]
                #print 'arr_substar_temp', arr_substar_temp
                for j in arr_substar_temp: #loop through ref imgs
                        #print 'compiling img', j, 'for subtraction'
                        arr_output = []
                        arr_rss = np.array([])
                        j_index = np.argwhere(arr_substar_temp == j)
                        j_index = j_index[0]
                        j_index = j_index[0]
                        filename =  directory + '/' + date_ref  + '/s' + str(int(j)).zfill(4) + '_RCF.fits'
                        if exists(filename):
                                img2 = pf.open(filename)[0].data
                                maxpixval = np.amax(img2)
                                img2 /= maxpixval
                        else:
                                #print filename, 'doesnt exist'
                                continue

                        img_filter = ret_filtername(j, date_ref)
                        
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
                                #print 'eliminating', j, 'from subtraction', 'from date', date_ref
                                #print 'FWHM greater than', fwhm_threshold
                                continue
                        arr_j_temp.append(int(j))
                        img2 = img2.flatten() #change to 1D array
                        if img2.size != img1size:
                                print 'img2.shape', img2.shape
                                print j, date_ref
                        arr_img2_temp.append(img2)
                arr_img2_temp = np.array(arr_img2_temp)
                print 'arr_img2_temp.shape', arr_img2_temp.shape, 'for date:', date_ref
                arr_img2.append(arr_img2_temp)
                arr_j.append(arr_j_temp)
        
        '''
        #pickle.dump(arr_j, open('test_arrj.p', 'wb'))
        #pickle.dump(arr_img2, open('test_arrimg2.p', 'wb'))
        arr_j = pickle.load(open('test_arrj.p', 'rb'))
        arr_img2 = pickle.load(open('test_arrimg2.p', 'rb'))
        #print 'saved test_arrj.fits'
        #print 'saved test_arrimg2.fits'
        '''

        #------        
        #Load array of filenumbs for ShaneAO, correlation matrix, and corresponding filenumbs in matrix
        #------
        print 'loading cross-correlation matrix to find best reference imgs'
        filename_structShaneAOfilenumbs = 'struct_ShaneAO_filenumbs.p'
        struct_ShaneAO_filenumbs = pickle.load(open(filename_structShaneAOfilenumbs, 'rb'))
        filename_matrixcorr = 'matrix_corr_filt.fits'
        matrix_corr = pf.open(filename_matrixcorr)[0].data
        filename_arrfilenumbscorr = 'arr_filenumbs_corr.fits'
        arr_filenumbs_corr = pf.open(filename_arrfilenumbscorr)[0].data
        setnumb1 = int(setnumb1) #set numb for target psf img
        arr_filenumbs = struct_ShaneAO_filenumbs[index_date_targ]
        filename_indexdates_corr = 'arr_indexdates_corr.fits'
        arr_indexdates_corr = pf.open(filename_indexdates_corr)[0].data

        #------
        #Getting rid of frames of the same set in correlation matrix
        #------
        for element in arr_targpsf: 
                index_date_in_arr_indexdates_corr = np.where(arr_indexdates_corr == index_date_targ)[0]
                index_imgnumb_in_arr_filenumbs_corr = np.where(arr_filenumbs_corr == element)[0]
                union_imgnumb_and_index = list(set(index_date_in_arr_indexdates_corr) & set(index_imgnumb_in_arr_filenumbs_corr))
                index_final = union_imgnumb_and_index[0]
                matrix_corr[:,index_final] = fillelement

        #------
        # Create list of structures with info for LOCI function
        # Takes certain number of best imgs by referencing correlation matrix
        #------
        print 'creating array of best reference imgs...'
        arr_pool_vars = []
        for index_i in np.arange(arr_targpsf.size):
                print '% of imgs added:', (index_i+1)*100./arr_targpsf.size
                arr_img2_optimal = []
                arr_j_optimal = []
                i = arr_targpsf[index_i] #file number of target img
                print '------'
                print 'target img number', i, date
                print_timenow()
                index_date_in_arr_indexdates_corr = np.where(arr_indexdates_corr == index_date_targ)[0]
                index_imgnumb_in_arr_filenumbs_corr = np.where(arr_filenumbs_corr == i)[0]
                union_imgnumb_and_index = list(set(index_date_in_arr_indexdates_corr) & set(index_imgnumb_in_arr_filenumbs_corr))
                #print 'union_imgnumb_and_index', union_imgnumb_and_index
                index_targ_corr = union_imgnumb_and_index[0]
                #print 'index_targ_corr', index_targ_corr
		arr_corrs = matrix_corr[index_targ_corr, :] #extract row in correlation matrix
                #print 'matrix_corr.shape', matrix_corr.shape
                #print 'arr_corrs[0:6]', arr_corrs[0:6]
                #useless = raw_input('stopped...')
                index_sort_corr = np.argsort(arr_corrs) #sort row for best correlations
                arr_indexdateoptimal = arr_indexdates_corr[index_sort_corr] #corresponding sorted date indexes
                arr_filenumbs_optimal = arr_filenumbs_corr[index_sort_corr] #corresponding sorted reference file numbers
                arr_index_j_optimal = []

                for index_optimal in np.arange(arr_filenumbs_optimal.size).astype(int):
                        index_date_optimal = arr_indexdateoptimal[index_optimal]
                        filenumb_optimal = arr_filenumbs_optimal[index_optimal]
                        if len(arr_j) > index_date_optimal: #make sure it's possible to index date
                                index_j_optimal = np.argwhere(arr_j[index_date_optimal] == filenumb_optimal)
                        if index_j_optimal:
                                img2_optimal = arr_img2[index_date_optimal][index_j_optimal[0][0]] #flattened reference img
                                arr_img2_optimal.append(img2_optimal)
                                index_j_optimal = 0
                        else:
                                continue
                        if len(arr_img2_optimal) >= numb_imgs_keep:
                                break
                arr_img2_optimal = np.array(arr_img2_optimal)
                arr_img2_optimal = np.transpose(arr_img2_optimal)

                #------
                # mock binaries?
                #------
                if bin_cond:
                        for index_rad in np.arange(arr_radiusmock.size): #loop through index of radius array for mocks
                                rad_mock = arr_radiusmock[index_rad]
                                fluxratio = arr_fluxratio[index_rad]

                                #------
                                #loop through position angles
                                #------
                                for theta_mock in arr_theta:
                                        struct = {'i': i, 'arr_img2': arr_img2_optimal, 'radius_sub': radius_sub, 'radius_op':radius_op, 'max_rad_img':max_rad_img, 'theta_mock':theta_mock, 'fluxratio':fluxratio, 'rad_mock':rad_mock, 'len_filt_box':len_filt_box, 'bin_cond':bin_cond, 'replace_cond':replace_cond}
                                        arr_pool_vars.append(struct) #append to array of structures
                else:
                        struct = {'i': i, 'arr_img2': arr_img2_optimal, 'radius_sub': radius_sub, 'radius_op':radius_op, 'max_rad_img':max_rad_img, 'len_filt_box':len_filt_box, 'bin_cond':bin_cond, 'replace_cond':replace_cond}
                        arr_pool_vars.append(struct)  #append to array of structures
                        
                                                        
        #------
        #Save/load list of structures as needed
        #------
        #filename_arrpoolvars = 'loci_arr_pool_vars_' + date + 'set' + str(int(setnumb1)) + '.p'        
        #pickle.dump(arr_pool_vars, open(filename_arrpoolvars, 'wb'))
        #arr_pool_vars = pickle.load(open(filename_arrpoolvars, 'rb'))
        
        print '------------------------'
	print 'No of loops to run:', len(arr_pool_vars)
        print '------------------------'
        print 'Now starting loops...'

        
        #------
        #Run LOCI. Comment/uncomment sections for when thread is needed/not needed
        #------
        max_loop_count = len(arr_pool_vars)
        counter = 0
	for elem in arr_pool_vars:
		theloci(elem)
                counter += 1 
                print counter*100./max_loop_count, '% done'
                print '------'

        '''
	threadnumb = 1
	thepool = ThreadPool(threadnumb)
        arr_pool = thepool.map(theloci, arr_pool_vars)
        thepool.close()
	'''

        
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


def create_reduced_imgs(setnumb, replace_cond = True):
        if not replace_cond:
                print "NOTE: Not replacing files if filenames exist. (replace_cond = False)"
	setnumb = str(int(setnumb))
	#uncomment the following & comment above if only 1 set is needed.
        #setnumb = raw_input('Enter set number (1,2,3,4, etc.):')

        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
        arr_files = np.arange(start1, end1+1)
	img_sky = pf.open(directory + '/' + date + '/' + 'img_sky_'+ setnumb +'.fits')[0].data
	img_sky = img_sky[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index]


	

        hdr = pf.open(directory + '/' + date + '/' + 's' + str(start1).zfill(4) + '.fits')[0].header
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
                filename = directory + '/' + date + '/' + 's' + str(i).zfill(4) + '.fits'
                if exists(filename):
                        img = pf.open(filename)[0].data
                else:
                        print filename, 'doesnt exist'
                        continue

		img = img[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index]
		img2 = (img-img_sky)/(img_flat)
		hdu = pf.PrimaryHDU(img2)
		hdulist = pf.HDUList([hdu])
                filename_output = directory + '/' + date + '/' + 's' + str(i).zfill(4) + '_reduced_.fits'
                if not replace_cond:
                        if exists(filename_output):
                                continue
                                print 'file exists. skipping'
                else:
                        hdulist.writeto(filename_output ,clobber=True)
                        #print 'created file', filename_output
                        
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
		print '_________'
		print i
                
                #------
                #load img file
                #------
                filename = directory + '/' + date + '/s' + str(i).zfill(4) + '_reduced_.fits'
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
		#print 'output.shape', output.shape
                #output.fill(0) ###
                #y_length = output.shape[0] ###
                #x_length = output.shape[1] ###
                #center_index = [(y_length-1)/2, (x_length-1)/2] ###
                #output[center_index[0], center_index[1]] = 0###
                arr_center.append(output)
                hdu = pf.PrimaryHDU(output)
                hdulist = pf.HDUList([hdu])

                hdulist.writeto(directory + '/' + date + '/s' + str(i).zfill(4) + '_reduced_centered.fits', clobber=True)
                #subprocess.call('/home/apps/ds9/ds9 ' + directory + '/' + date + '/s' + str(i).zfill(4) + '_reduced_centered.fits', shell = True) #

	print 'arr_center shape:', str(np.array(arr_center).shape)
	img_mdn = np.median(np.array(arr_center), axis=0)
	print img_mdn.shape
	hdu = pf.PrimaryHDU(img_mdn)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto(directory + '/' + date + '/' + outputname + '.fits', clobber=True)
        print 'created:', outputname + '.fits'
        #subprocess.call('/home/apps/ds9/ds9 ' + directory + '/' + date + '/' + outputname + '.fits', shell = True) ###

        
def run_create_sky():
#run create_sky for all science stars for specific observing date
#create_sky medians science imgs to create sky img

	startset = 1
        arr_setnumbs = np.arange(startset, total_sets+1).astype(int)
        '''
        threadnumb = 6
	thepool = ThreadPool(threadnumb)
        arr_imgfinals = thepool.map(create_sky, arr_setnumbs)
        thepool.close()
        thepool.join()
        '''
        #''' FOR NO THREADS
	for setnumb in np.arange(startset, total_sets+1):
		print 'set number', setnumb
		create_sky(setnumb)
        #'''
        
def run_create_reduced_imgs():
        
        startset = 1
        '''
        arr_setnumbs = np.arange(startset, total_sets+1).astype(int)
        threadnumb = 6
	thepool = ThreadPool(threadnumb)
        arr_imgfinals = thepool.map(create_reduced_imgs, arr_setnumbs)
        thepool.close()
        thepool.join()
        '''

        #''' FOR NO THREADS
	for setnumb in np.arange(startset, total_sets+1):
		print 'set number', setnumb
		create_reduced_imgs(setnumb)
        #'''

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
	

def run_subframe():
        #------
        #ask user for set number on which to perform subtractions
        #------
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')



        #------
        #Arrays for parameters of mock binaries
        #------
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 45, 60])
        print 'arr_radiusmock', arr_radiusmock
        arr_mockfactor = np.array([100, 150, 200, 500])
        print 'arr_mockfactor', arr_mockfactor
        arr_theta = np.array([0])
        print 'arr_theta', arr_theta


        
        #------
        # Append structures with parameters onto list
        # List then acts as input to psf_subtract_frame function
        #------
        arr_struct = []
        for radius_mock in arr_radiusmock:
                for mock_factor in arr_mockfactor:
                        for theta in arr_theta:
                                struct = {'radius_mock':radius_mock, 'mock_factor': mock_factor, 'theta': theta, 'setnumb1':setnumb1}
                                arr_struct.append(struct)


        #------
        # Run psf_subtract_frame
        #------
        print 'Number of subtractions:', len(arr_struct)


        threadnumb = 3
	thepool = ThreadPool(threadnumb)
        arr_pool = thepool.map(psf_subtract_frame, arr_struct)
        thepool.close()
        thepool.join()

        
        '''
        for elem in arr_struct:
                psf_subtract_frame(elem)
        '''


        
        
def psf_subtract_frame(struct):
        #perform psf subtraction frame by frame
        #------
        print struct
        setnumb1 = setnumb1 = struct['setnumb1']
        #print 'Target set number:', setnumb1
        radius_mock = struct['radius_mock']
        mock_factor = struct['mock_factor']
        theta_mock = struct['theta']
        filename_output = directory + '/' + date + '/' + 'set' + str(setnumb1) + '_mockrad' + str(radius_mock) + '_ratio' + str(mock_factor) + '_theta' + str(theta_mock) + '.fits'
        '''
        if exists(filename_output):
                print 'SKIP:', 'mockrad', radius_mock, 'ratio', mock_factor, 'theta', theta_mock
                return 0
        '''
        
        #setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])

        #setnumb2 = raw_input('Enter 2nd set number (1,2,3,4, etc.):')
        #arr_startend2 = open('set_' + setnumb2 + '/set_'+ setnumb2 +'_info.txt', 'rb').read().splitlines()
        arr_substar = np.array([])
       
        arr_setnumb = np.array([])
        for setnumb2 in np.arange(1,total_sets+1):
                if setnumb2 != int(setnumb1):
                        #print 'adding elements in set no.', setnumb2
                        arr_startend2 = open(directory + '/' + date + '/set_'+ str(setnumb2) +'_info.txt', 'rb').read().splitlines()
                        start2 = int(arr_startend2[0])
                        end2 = int(arr_startend2[1])
                        arr_imgnumbs = np.arange(start2, end2+1)
                        arr_substar = np.append( arr_substar,  arr_imgnumbs)
                        arr_temp = np.zeros(arr_imgnumbs.size)
                        arr_temp.fill(setnumb2)
                        arr_setnumb = np.append(arr_setnumb, arr_temp)

        #define  subtraction ratios to be iterated through
        arr_rad = np.arange(80)+1
        arr_subratio = (np.arange(20)+90)/100.
        #________________________________________________

        #define empty arrays, for psf-subtracted imgs for avging
        arr_img = []
        arr_mdn = np.array([])
        arr_5sd = np.array([])
        #________________________________________________

        print 'starting:', 'mockrad', radius_mock, 'ratio', mock_factor, 'theta', theta_mock, datetime.datetime.now().strftime('%m-%d %H:%M:%S')
        
	for i in np.arange(start1, end1+1): #iterating through images in the 1st set
		print '_________'
                print 'i',i
                print 'start1',start1
                print 'end1',end1
                print 100.*(float(i - start1)/float(end1-start1)), '%', 'for mockrad', radius_mock, 'ratio', mock_factor, datetime.datetime.now().strftime('%m-%d %H:%M:%S')
                if i > 999:
                        filename = directory + '/' + date + '/s' + str(int(i)) + '_reduced_centered.fits'
                        if exists(filename):
                                img1 = pf.open(filename)[0].data
                                img1 /= np.amax(img1) #divide by max. brightness
                        else:
                                continue
                else:
                        filename = directory + '/' + date + '/s0' + str(int(i)) + '_reduced_centered.fits'
                        if exists(filename):
                                img1 = pf.open(filename)[0].data
                                img1 /= np.amax(img1) #divide by max. brightness
                        else:
                                continue


                #------
                #Insert mock binary
                #------
                img_mock = pf.open(directory + '/' + date + '/' + 'centroid_1.fits')[0].data
                img_mock /= np.amax(img_mock)
                if mock_factor != 0:
                        img_mock /= float(mock_factor)
                else:
                        img_mock *= float(mock_factor)
                        
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
                        


                #------
                # Loop through reference img numbers
                # Check for best subtraction
                #------
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


                                
                        #------
                        # Iterating through subtraction ratios
                        #------
                        for subratio in arr_subratio: 
                                img1_temp = np.copy(img1)
                                img2_temp = np.copy(img2)
                                img2_temp /= np.amax(img2_temp) #divide by max. brightness
                                img2_temp *= subratio #multiply ref img by subtraction ratio
                                img1_temp -= img2_temp #psf subtract
                                arr_rss = np.append(arr_rss, np.sum(img1_temp**2.)) #add residual sum of squares to arr_rss.
                                arr_output.append(img1_temp) #append psf-subtracted img
                        index_min_rss = np.argmin(arr_rss) #index of min rss in array of psfsubtracted imgs of diff sub-ratios
                        img_best_subtract = arr_output[index_min_rss]
                        #print 'best subtraction ratio', arr_subratio[index_min_rss]
                        if counter2 == 0:
                                rss_best = np.sum(img_best_subtract**2.)
                                img_best = img_best_subtract
                        else:
                                if np.sum(img_best_subtract**2.) < rss_best:
                                        rss_best = np.sum(img_best_subtract**2.)
                                        img_best = img_best_subtract
                        counter2 += 1
                arr_img.append(img_best)
        print 'size of arr_img', len(arr_img)
        img_mdn = np.median(np.array(arr_img), axis = 0)

        print 'creating', filename_output
        hdu = pf.PrimaryHDU(img_mdn)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_output, clobber = True)
        
        return 1


def ret_filtername(filenumb, date):
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


def plot_mocks_klip():
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 45, 60]).astype(int)#np.arange(15, 60+1, 15).astype(int)
        #arr_radiusmock = arr_radiusmock[::-1]
        print 'arr_radiusmock', arr_radiusmock
        arr_mockfactor = np.array([30, 100, 150, 200, 500, 1000]).astype(int)
        print 'arr_mockfactor', arr_mockfactor
        arr_theta = np.array([0]).astype(int) #np.array([0, 198, 35, 233, 71, 269]).astype(int)
        print 'arr_theta', arr_theta
        arr_klipmodes = np.array([15, 25, 35]).astype(int)
        rad_annulus = 10


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


        
        #------
        # plot different graphs for different flux ratio
        #------
        arr_rss = []
	arr_rad = []
	arr_radcheck = np.arange(1, 61)

        for mock_factor in arr_mockfactor:
                for radius_mock in arr_radiusmock:
                        for theta_mock in arr_theta:
                                print 'theta_mock', theta_mock
                                dx_mock = radius_mock*np.cos(theta_mock*np.pi/180)
                                dx_mock = int(round(dx_mock))
                                dy_mock = radius_mock*np.sin(theta_mock*np.pi/180)
                                dy_mock = int(round(dy_mock))
                                print 'dy_mock, dx_mock:', dy_mock, ',', dx_mock
                                for numb_klip_modes in arr_klipmodes:
                                        arr_mdn = []
                                        for i in arr_targpsf:
                                                filename_output = directory + '/' + date + '/' + str(i) + '_klipfinal_' + 'kmode' + str(int(numb_klip_modes)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radannulus' + str(int(rad_annulus)) + '.fits'
                                                #filename_output = directory + '/' + date + '/' + str(i) + '_klipfinal_' + 'kmode' + str(int(numb_klip_modes)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_theta' + str(int(theta_mock)) + '.fits'
                                                if exists(filename_output):
                                                        img_final = pf.open(filename_output)[0].data
                                                else:
                                                        print filename_output, 'doesnt exist'
                                                        continue
                                                y_length = img_final.shape[0]
                                                x_length = img_final.shape[1]
                                                center_index_y, center_index_x = int((y_length-1)/2), int((x_length-1)/2)
                                                index_mock_y, index_mock_x = center_index_y + dy_mock, center_index_x + dx_mock
                                                arr_mdn.append(img_final)
                                        if arr_mdn:
                                                print 'numb of imgs for median:', len(arr_mdn)
                                                img_mdn = np.median(np.array(arr_mdn), axis = 0)
                                                filename_finalfits = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + 'kmode' + str(int(numb_klip_modes)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radannulus' + str(int(rad_annulus)) + '.fits'
                                                hdu = pf.PrimaryHDU(img_mdn)
                                                hdulist = pf.HDUList([hdu])
                                                hdulist.writeto(filename_finalfits, clobber = True)

                                                flux_mock = img_mdn[index_mock_y, index_mock_x]
                                                flux_mock_diff = 2.5*np.log10((1./mock_factor)/flux_mock)
                                                print 'img_mdn.shape', img_mdn.shape
                                                arr_5sd = np.array([])
                                                for radius in arr_radcheck:
                                                        if mock_factor < 50:
                                                                struct = check_ring_mock(img_mdn, radius, theta_mock, 90)
                                                        else:
                                                                struct = check_ring_mock(img_mdn, radius, theta_mock)
                                                        arr_5sd = np.append(arr_5sd, 5*struct['sd_ring'])			
                                                arr_5sd_mag = 2.5*np.log10(1./np.abs(arr_5sd))
                                                #if all(element > 6 for element in arr_5sd_mag[:10]):
                                                #plt.plot(arr_radcheck, arr_5sd_mag, 'o-', color = np.random.rand(3,1),  label = '5 s.d. values: ' + 'klip modes:' + str(int(numb_klip_modes)) + ', flux ratio:' + str(mock_factor) + ', mock radius: ' + str(radius_mock) + ', flux lost: ' + str(flux_mock_diff) + 'mag')
                                                print '------'
                                                        #index_min = np.argmin(arr_rss)
                                #plt.legend(loc = 'upper right')
                                #plt.xlabel('radius (pixels)')
                                #plt.ylabel('magnitude')
                                #plt.gca().invert_yaxis()
                                #plt.show()
                                #useless = raw_input('waiting for input')




        
def plot_mocks_loci():
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 45, 60]).astype(int)#np.arange(15, 60+1, 15).astype(int)
        print 'arr_radiusmock', arr_radiusmock
        arr_mockfactor = np.array([30, 100, 150, 200, 500, 1000]).astype(int)
        print 'arr_mockfactor', arr_mockfactor
        arr_theta = np.array([0]).astype(int) #np.array([0, 198, 35, 233, 71, 269]).astype(int)
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


        
        #------
        # plot different graphs for different flux ratio
        #------
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
                                                        #print filename_output_fits
                                                        if exists(filename_output_fits):
                                                                img_final = pf.open(filename_output_fits)[0].data
                                                                '''
                                                                #print 'Pickle error for', filename_output
                                                                #filename_output = directory + '/' + date + '/' + str(i) + '_thread_mock_rad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_.fits'
                                                                #if exists(filename_output):
                                                                #        img_final = pf.open(filename_output)[0].data
                                                                #else:
                                                                #        continue
                                                                #img_final = img_final[80-halfboxcheck:80+halfboxcheck+1,80-halfboxcheck:80+halfboxcheck+1]'''
                                                        else:
                                                                print 'fits file doesnt exist, checking for pickle file'
                                                                if exists(filename_output):
                                                                        try:
                                                                                with open(filename_output, 'rb') as file_img:
                                                                                        img_final = pickle.load(file_img)
                                                                        except:
                                                                                print filename_output, 'corrupted'
                                                                                continue
                                                                else:
                                                                        print filename_output, 'and', filename_output_fits, 'dont exist'
                                                                        continue
                                                        y_length = img_final.shape[0]
                                                        x_length = img_final.shape[1]
                                                        center_index_y, center_index_x = int((y_length-1)/2), int((x_length-1)/2)
                                                        index_mock_y, index_mock_x = center_index_y + dy_mock, center_index_x + dx_mock
                                                        arr_mdn.append(img_final)
                                                if arr_mdn:
                                                        filename_finalfits = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_loci.fits'
                                                        print 'numb of imgs for median:', len(arr_mdn)
                                                        img_mdn = np.median(np.array(arr_mdn), axis = 0)
                                                        flux_mock = img_mdn[index_mock_y, index_mock_x]
                                                        flux_mock_diff = 2.5*np.log10((1./mock_factor)/flux_mock)
                                                        print 'img_mdn.shape', img_mdn.shape
                                                        arr_5sd = np.array([])

                                                        img_mdn[center_index_y, center_index_x] = 0
                                                        
                                                        hdu = pf.PrimaryHDU(img_mdn)
                                                        hdulist = pf.HDUList([hdu])
                                                        hdulist.writeto(filename_finalfits, clobber = True)
                                                        
                                                        for radius in arr_radcheck:
                                                                if mock_factor < 50:
                                                                        struct = check_ring_mock(img_mdn, radius, theta_mock, 90)
                                                                else:
                                                                        struct = check_ring_mock(img_mdn, radius, theta_mock)
                                                                arr_5sd = np.append(arr_5sd, 5*struct['sd_ring'])			
                                                        arr_5sd_mag = 2.5*np.log10(1./np.abs(arr_5sd))
                                                        #if all(element > 6 for element in arr_5sd_mag[:10]):
                                                        #plt.plot(arr_radcheck, arr_5sd_mag, 'o-', color = np.random.rand(3,1),  label = '5 s.d. value, ' + 'theta: ' + str(theta_mock) + ', flux ratio:' + str(mock_factor) + ', mock radius: ' + str(radius_mock) + ', flux lost: ' + str(flux_mock_diff) + 'mag')
                                                        #index_min = np.argmin(arr_rss)
                                #plt.legend(loc = 'upper right')
                                #plt.xlabel('radius (pixels)')
                                #plt.ylabel('magnitude')
                                #plt.gca().invert_yaxis()
                                #plt.show()
                                        #useless = raw_input('waiting for input')
                        
def compare_plots():
        #------
	#ask for user input regarding which set to take as target psf
        #reads txt file in directory to know which files to use        
        #------
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
	arr_targpsf = np.arange(start1, end1+1)



        #------
        #Change if necessary
        #------
        numb_klip_modes = 25
        ratio = 1000
        radius_sub = 6
        radius_op = 19
        theta_mock = 0
        mock_factor = 30
        arr_x = [10, 15, 20, 25, 30, 45, 60]
        arr_y_klip = []
        arr_y_loci = []
        arr_y_frame = []



        #------
        #Load imgs with weak binary
        #Change filenames if necessary
        #------
        filename_loci_clean = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_mockrad10' + '_' +  'ratio1000' + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta0' + '_loci.fits'
        filename_klip_clean = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + 'kmode25' + '_mockrad10' + '_' +  'ratio1000' + '_theta0' + '_klip.fits'
        filename_frame_clean = directory + '/' + date + '/' + 'set' + str(setnumb1) + '_mockrad10' + '_ratio1000' + '_theta0' + '.fits'
        
        img_loci_clean = pf.open(filename_loci_clean)[0].data
        img_klip_clean = pf.open(filename_klip_clean)[0].data
        img_frame_clean = pf.open(filename_frame_clean)[0].data
                
                
        #------
        # Loop through mock radii,
        # add corrected detection limit to arr_y_klip and arr_y_loci
        #------
        for radius_mock in arr_x:
                filename_loci = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radsub' + str(int(radius_sub)) + '_radop' + str(int(radius_op)) + '_theta' + str(int(theta_mock)) + '_loci.fits'
                filename_klip = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + 'kmode' + str(int(numb_klip_modes)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_theta' + str(int(theta_mock)) + '_klip.fits'
                filename_frame = directory + '/' + date + '/' + 'set' + str(setnumb1) + '_mockrad' + str(radius_mock) + '_ratio' + str(mock_factor) + '_theta' + str(theta_mock) + '.fits'

                

                #------
                #calculate angle
                dx_mock = radius_mock*np.cos(theta_mock*np.pi/180)
                dx_mock = int(round(dx_mock))
                dy_mock = radius_mock*np.sin(theta_mock*np.pi/180)
                dy_mock = int(round(dy_mock))


                
                #------
                #open fits files
                img_loci = pf.open(filename_loci)[0].data
                img_klip = pf.open(filename_klip)[0].data
                img_frame = pf.open(filename_frame)[0].data


                
                #------
                #calculate indexes of mock binary
                y_length = img_loci.shape[0]
                x_length = img_loci.shape[1]
                center_index_y, center_index_x = int((y_length-1)/2), int((x_length-1)/2)
                index_mock_y, index_mock_x = center_index_y + dy_mock, center_index_x + dx_mock


                
                #------
                #Calculate flux losses from klip/loci
                flux_mock_loci = img_loci[index_mock_y, index_mock_x]
                flux_mock_diff_loci = 2.5*np.log10((1./mock_factor)/flux_mock_loci)
                flux_mock_klip = img_klip[index_mock_y, index_mock_x]
                flux_mock_diff_klip = 2.5*np.log10((1./mock_factor)/flux_mock_klip)
                flux_mock_frame = img_frame[index_mock_y, index_mock_x]
                flux_mock_diff_frame = 2.5*np.log10((1./mock_factor)/flux_mock_frame)

                
                print 'radius_mock', radius_mock
                print 'intial mock peak flux', 1./mock_factor
                print 'mock peak flux after', flux_mock_frame
                print 'flux_mock_diff_frame', flux_mock_diff_frame
                print '------'

                
                #------
                #Check clean imgs for 5s.d. values
                #Subtract flux loss from KLIP/LOCI
                #Append to arrays for plotting
                #------
                struct_loci = check_ring_mock(img_loci_clean, radius_mock, theta_mock)
                mag_corrected_loci = 2.5*np.log10(1./np.abs(5*struct_loci['sd_ring'])) - flux_mock_diff_loci
                arr_y_loci.append(mag_corrected_loci)
                        
                struct_klip = check_ring_mock(img_klip_clean, radius_mock, theta_mock)
                mag_corrected_klip = 2.5*np.log10(1./np.abs(5*struct_klip['sd_ring'])) - flux_mock_diff_klip
                arr_y_klip.append(mag_corrected_klip)
                
                struct_frame = check_ring_mock(img_frame_clean, radius_mock, theta_mock)
                mag_corrected_frame = 2.5*np.log10(1./np.abs(5*struct_loci['sd_ring'])) - flux_mock_diff_frame
                arr_y_frame.append(mag_corrected_frame)  
        
        plt.plot(arr_x, arr_y_loci, 'o-', color = 'b',  label = 'LOCI')
        plt.plot(arr_x, arr_y_klip, 'o-', color = 'r',  label = 'KLIP')
        plt.plot(arr_x, arr_y_frame, 'o-', color = 'g', label = 'Frame-by-frame')

        '''
        minpix = 10
        arr_radi = np.arange(70)+1
        arr_5sd_frame = pickle.load(open(directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '-' + 'others' + '_5sdplot.p', 'rb'))
        plt.plot(arr_radi[minpix:arr_5sd_frame.size], 2.5*np.log10(1./np.abs(arr_5sd_frame[minpix:])), '^-', color = 'g', label = 'Frame-by-frame subtraction')
        '''
        
        plt.legend(loc = 'upper right')
        plt.xlabel('Radius (pixels)')
        plt.ylabel('Corrected 5 s.d. values')
        plt.gca().invert_yaxis()
        plt.show()

def ret_filtername_set(setnumb1, date = date):
        #------
        #Ask for user input regarding which set to take as target psf
        #Reads txt file in directory to know which files to use
        #------
        #setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        setnumb1 = str(int(setnumb1))
        print '------'
        #print setnumb1
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        #print arr_startend1
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
        starname = str(arr_startend1[2])
        starname = starname.replace(' ', '')
        #print starname
        arr_targpsf = np.arange(start1, end1+1)

        if start1 > 999:
                filename = directory + '/' + date + '/' + 's' + str(start1) + '.fits'
                if exists(filename):
                        fits = pf.open(filename)
        else:
                filename = directory + '/' + date + '/' + 's0' + str(start1) + '.fits'
                if exists(filename):
                        fits = pf.open(filename)
        hdr = fits[0].header
        filtername = hdr['filt1nam']
        date_taken = hdr['date-obs']
        #print 'date taken', date_taken
        #print 'filter', filtername                
        return filtername


def run_create_final_loci(total_sets = total_sets):
        arr_setnumb = np.arange(1, total_sets+1)
        for setnumb in arr_setnumb:
                create_final_loci(setnumb)

def create_final_loci(setnumb1, filt = True):
        #------
        # Create final LOCI-subtracted imgs
        #------

        
        #------
        #Reads txt file in directory to know which files to use
        #------
        setnumb1 = str(int(setnumb1))
        print '------'
        print setnumb1
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        print arr_startend1
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
        starname = str(arr_startend1[2])
        starname = starname.replace(' ', '')
        print starname
        arr_targpsf = np.arange(start1, end1+1)


        filename = directory + '/' + date + '/' + 's' + str(start1).zfill(4) + '.fits'
        if exists(filename):
                fits = pf.open(filename)

        hdr = fits[0].header
        filtername = hdr['filt1nam']
        date_taken = hdr['date-obs']
        print 'date taken', date_taken
        print 'filter', filtername

        
        #------
        # Load LOCI subtracted images and add to array for medianing
        #------
        arr_for_mdn = []
        for i in arr_targpsf:
                if filt:
                        filename_input = directory + '/' + date + '/' + str(int(i)) + 'locifiltfinal' + '.fits' #change name if filtered/not filtered
                else:
                        filename_input = directory + '/' + date + '/' + str(int(i)) + 'locifinal' + '.fits' #change name if filtered/not filtered

                try:
                        with pf.open(filename_input) as f_temp:
                                img = f_temp[0].data
                except:
                        print filename_input, 'doesnt exist', 'from set:', setnumb1
                        continue
                arr_for_mdn.append(img)
        print 'images for mdn:',len(arr_for_mdn)
        img_final = np.median(np.array(arr_for_mdn), axis = 0)

        
        
        #------
        # Save medianed img as fits file
        #------
        if filt:
                filename_output = directory + '/' + date + '/' + starname +  '__' + filtername + '__' + date_taken + '__locifiltfinal.fits'
                filename_setnumbs = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) +  '_locifiltfinal.fits'
        else:
                filename_output = directory + '/' + date + '/' + starname +  '__' + filtername + '__' + date_taken + '__locifinal.fits'
                filename_setnumbs = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) +  '_locifinal.fits'

        if img_final.size:
                hdu = pf.PrimaryHDU(img_final)
                hdulist = pf.HDUList([hdu])
                hdulist.writeto(filename_output, clobber = True)
                hdulist.writeto(filename_setnumbs, clobber = True)        
                print 'created', filename_output
                print 'created', filename_setnumbs
        else:
                print 'empty array for output. File not created.'
        
        
        

        
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


def theloci(struct):
        #------
        #load variables from input structure, create arr_rad with appropriate starting subtraction radii
        #------
        i = struct['i']
	radius_sub = struct['radius_sub']
	radius_op = struct['radius_op']
	arr_img2 = struct['arr_img2']
        max_rad_img = struct['max_rad_img']
        len_filt_box = struct['len_filt_box']
        arr_rad = np.arange(1, max_rad_img, radius_sub)
        arr_rad = arr_rad.astype(int)

                
        #------
        #define filename for final img
        #------
        if struct['bin_cond']:
                radius_mock = struct['rad_mock']
                mock_factor = struct['fluxratio']
                theta_mock = struct['theta_mock']
                filename_output = directory + '/' + date + '/' + str(i) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfiltfinal.fits' #!!!!!!***filt Change name if not using filt #faint***
                print 'starting', i, 'radius mock:', radius_mock, ',', 'theta', theta_mock, '...', datetime.datetime.now().strftime('%m-%d %H:%M:%S')
        else:
                filename_output = directory + '/' + date + '/' + str(i) + 'locifiltfinal' + '.fits'#!!!!!!***filt Change name if not using filt
                print 'starting', i, '...',
                print_timenow()

        if not struct['replace_cond']:
                if exists(filename_output):
                        print 'skipping', filename_output
                        return 0

        #------
        #load target img, divide by max pixel value
        #------
        filename = directory + '/' + date  + '/s' + str(int(i)).zfill(4) + '_reduced_centered.fits'
        if exists(filename):
                try:
                        img1 = pf.open(filename)[0].data
                        max_pix_primary = np.amax(img1)
                        img1 /= max_pix_primary
                except:
                        print 'error opening file:', filename
                        return 0
        else:
                #print filename, 'doesnt exist'
                return 0
        y_length, x_length = img1.shape
        center_index_y, center_index_x = int((y_length-1)/2), int((x_length-1)/2)

        #------
        # inject mock binary if needed
        #------
        if struct['bin_cond']:                                
                #----------
                #load mock binary and divide by predetermined factor
                #***Change centroid filename if necessary***
                #----------
                img_mock = pf.open(directory + '/' + 'Jun15' + '/' + 'centroid_1.fits')[0].data
                img_mock /= np.amax(img_mock)
                if mock_factor != 0:
                        img_mock *= float(mock_factor)
                        
                #--------
                #Add mock binary at correct radius & angle, save img as fits file
                #--------
                if radius_mock != 0:
                        dx_mock = radius_mock*np.cos(theta_mock*(2*np.pi)/360) #calculating y displacement from center
                        dx_mock = int(round(dx_mock))
                        dy_mock = radius_mock*np.sin(theta_mock*(2*np.pi)/360) #calculating x displacement from center
                        dy_mock = int(round(dy_mock))
        
                        img_mock = np.fft.fft2(img_mock) #Fourier shift to construct binary for adding to primary
                        img_mock = fshift(img_mock, [dy_mock, dx_mock])
                        img_mock = np.fft.ifft2(img_mock)
                        img_mock = np.real(img_mock)

                        img1 += img_mock #Add binary to primary        

        img1 *= max_pix_primary  # Multiply by primary max pixel value from before
        img1 -= medfilt(img1, [len_filt_box, len_filt_box]) # Perform high-pass filter on img
        img1 /= np.amax(img1) # Normalize img using new max pixel value


        #--------|
        #  LOCI  |
        #--------|
        img_final = np.zeros(img1.size)
        for radi in arr_rad:

                #array of starting radi of subtraction region
                if radi+radius_sub > max_rad_img:
                        arr_radsub = np.arange(radi, max_rad_img)
                else:
                        arr_radsub = np.arange(radi, radi+radius_sub) 

                #array of starting radi of optimization region
                if radi+radius_op > max_rad_img:
                        arr_radop = np.arange(radi, max_rad_img)
                else:
                        arr_radop = np.arange(radi, radi+radius_op)
                        
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

        img_output[center_index_y, center_index_x] = 0
        
        #------
        #Save final img as fits file
        #------
        hdu = pf.PrimaryHDU(img_output)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_output, clobber = True)

        print 'done with', filename_output, '...', 
        print_timenow()
	return 1


def loci_gen(img1, radius_sub, radius_op, arr_img2, filename_output):
        #------
        #load variables from input structure, create arr_rad with appropriate starting subtraction radii
        #------

        # perform calculations for max_rad_img***!!!
        arr_rad = np.arange(1, max_rad_img, radius_sub)
        arr_rad = arr_rad.astype(int)

        #------
        #obtain dimensions of img, 2-D index of center pixel
        #------
        y_length, x_length = img1.shape
        center_index_y, center_index_x = int((y_length-1)/2), int((x_length-1)/2)


        #------
        # Perform high-pass filter on img
        #------
        img1 -= medfilt(img1, [len_filt_box, len_filt_box])

        #------
        # Normalize img using new max pixel value
        #------
        img1 /= np.amax(img1)


        #--------|
        #  LOCI  |
        #--------|
        img_final = np.zeros(img1.size)
        for radi in arr_rad:
                #------
                #array of starting radi of subtraction region
                #------
                if radi+radius_sub > max_rad_img:
                        arr_radsub = np.arange(radi, max_rad_img)
                else:
                        arr_radsub = np.arange(radi, radi+radius_sub) 

                #------
                #array of starting radi of optimization region
                #------
                if radi+radius_op > max_rad_img:
                        arr_radop = np.arange(radi, max_rad_img)
                else:
                        arr_radop = np.arange(radi, radi+radius_op)
                        
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


        img_output[center_index_y, center_index_x] = 0
        
        #------
        #Save final img as fits file
        #------
        hdu = pf.PrimaryHDU(img_output)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_output, clobber = True)

        print 'done with', filename_output, '...', datetime.datetime.now().strftime('%m-%d %H:%M:%S')
	return 1




def run_create_final_loci_mockbins(total_sets = total_sets, highpassfilt = True):
        startset = 11 #faint***
        arr_setnumbs = np.arange(startset, total_sets+1).astype(int)
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80]).astype(int)
        print 'arr_radiusmock', arr_radiusmock
        arr_theta = np.arange(0, 360, 60).astype(int)
        print 'arr_theta', arr_theta

        for setnumb1 in arr_setnumbs:
                print '------'
                print setnumb1
                print '------'
                arr_startend1 = open(directory + '/' + date + '/set_'+ str(int(setnumb1)) +'_info.txt', 'rb').read().splitlines()
                start1 = int(arr_startend1[0])
                end1 = int(arr_startend1[1])
                arr_targpsf = np.arange(start1, end1+1)

                for theta_mock in arr_theta:
                        for radius_mock in arr_radiusmock:

                                arr_for_mdn = []
                                if highpassfilt:
                                        filename_output = directory + '/' + date + '/' + 'set'+ str(int(setnumb1)) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locifaintmockfiltfinal.fits' #faint***
                                else:
                                        filename_output = directory + '/' + date + '/' + 'set'+ str(int(setnumb1)) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfinal.fits'


                                for i in arr_targpsf:
                                        if highpassfilt:
                                                filename_input = directory + '/' + date + '/' + str(i) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locifaintmockfiltfinal.fits' #faint***
                                        else:
                                                filename_input = directory + '/' + date + '/' + str(i) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfinal.fits'
                                        try:
                                                with pf.open(filename_input) as f_temp:
                                                        img = f_temp[0].data
                                        except:
                                                print filename_input, 'doesnt exist', 'from set:', setnumb1
                                                continue
                                        arr_for_mdn.append(img)
                                print 'images for mdn:',len(arr_for_mdn)

                                if len(arr_for_mdn) > 1:
                                        img_final = np.median(np.array(arr_for_mdn), axis = 0)

                                        hdu = pf.PrimaryHDU(img_final)
                                        hdulist = pf.HDUList([hdu])
                                        hdulist.writeto(filename_output, clobber=True)
                                        print 'created:', filename_output

def save_fits(arr, filename):
#Inputs are: array for saving to fit & filename for fits.
#------
        hdu = pf.PrimaryHDU(arr)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename, clobber=True)

def print_timenow():
        print datetime.datetime.now().strftime('%m-%d %H:%M:%S')

def plot_detlim_final(filt = True):
        #Plot final detection limits for ALL stars.
        #User correction curves that were saved for specific run and filter
        #------
        global date
        global total_sets
        radi_apert = 3

        #------
        # load array with radii at which we had mock up binaries
        #------
        filename_arrradiusmock = directory + '/' + date + '/' + 'arr_radiusmock.fits'
        arr_radius = pf.open(filename_arrradiusmock)[0].data
        
               
        #------
        #open txt file with list of known binaries
        #------
        arr_bin = open(directory + '/' + directory + '_clearbinaries.txt', 'rb').read().splitlines()
        arr_bin_dates = np.array(arr_bin[0::2])
        arr_bin_setnumbs = np.array(arr_bin[1::2])

         
        #------
        # Obtain detection limits after LOCI for all stars
        #------
        arr_plots = []
        arr_labels = []
        arr_colors = ['r', 'k', 'g', 'm', 'b', 'y']
        for index_date in range(len(arr_date_ShaneAO)):
                date = arr_date_ShaneAO[index_date]
                total_sets = struct_ShaneAO_total_sets[date]
                counter_ks = 0
                counter_brg = 0
                counter = 0
                for setnumb1 in np.arange(total_sets) + 1:


                        #------
                        #Get filter name for this particular set
                        #------
                        filtername = ret_filtername_set(setnumb1, date = date)

                        #------
                        #Uncomment if we need to add filter name to set_setnumb_info.txt (for some reason)
                        #------
                        filename_setinfo = directory + '/' + date + '/' + 'set_' + str(int(setnumb1)) +'_info.txt'
                        with open(filename_setinfo, 'rb') as file:
                                arr_info = file.read().splitlines()
                                for row in arr_info:
                                        row = row.rstrip().lstrip()
                                        if not row:
                                                arr_info.remove(row)
                        arr_info = filter(lambda a: a != filtername, arr_info)
                        arr_info.append(filtername)
                        with open(filename_setinfo, 'w') as newfile:
                                for item in arr_info:
                                        newfile.write("%s\n" % item)
                        print 'updated:', filename_setinfo
                        
                        
                        #------
                        #check if its a known binary. If so, skip.
                        #------
                        cond_bin = False
                        if str(int(setnumb1)) in arr_bin_setnumbs:
                                arr_index_setnumb_bin = np.where(arr_bin_setnumbs == str(int(setnumb1)))[0]
                                for index_setnumb_bin in arr_index_setnumb_bin:
                                        if arr_bin_dates[index_setnumb_bin] == date:
                                                print date, ', set', setnumb1, 'skipped since it is a known binary'
                                                cond_bin = True
                        if cond_bin:
                                continue

                        #filename for array of final detection limit
                        filename_detlim_final = date + '_set' + str(int(setnumb1)) + '_detlim_final.fits'
                        print filename_detlim_final

                        
                        '''
                        with open(filename_setinfo_new, 'rb') as file:
                                arr_info_new = file.read().splitlines()
                                for row in arr_info_new:
                                        row = row.rstrip().lstrip()
                                        if not row:
                                                arr_info_new.remove(row)
                        print 'arr_info_new:', arr_info_new
                        
                        useless = raw_input('press key to continue')
                        '''

                        filename_fits_filt = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_locifiltfinal.fits'

                        filename_fits = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_locifinal.fits'

                        img_subed_filt = pf.open(filename_fits_filt)[0].data
                        img_subed = pf.open(filename_fits)[0].data


                        #------
                        # Open pre-subtraction, medianed img for set number
                        # Obtain flux of aperture around primary
                        #------
                        if filt:
                                filename_avg = directory + '/' + date + '/' + 'centroid_filt_' + str(int(setnumb1)) + '.fits'
                        filename_avg_unfilt = directory + '/' + date + '/' + 'centroid_' + str(int(setnumb1)) + '.fits'
                        if not exists(filename_avg) or not exists(filename_avg_unfilt):
                                print filename_avg, 'or', filename_avg_unfilt, "doesn't exist"
                        img_avg = pf.open(filename_avg)[0].data
                        img_avg /= np.amax(img_avg)
                        img_avg_unfilt = pf.open(filename_avg_unfilt)[0].data
                        img_avg_unfilt /= np.amax(img_avg_unfilt)

                        y_length_avg, x_length_avg = img_avg.shape
                        y_index_center_avg = int((y_length_avg - 1)/2)
                        x_index_center_avg = int((x_length_avg - 1)/2)

                        flux_main_filt = get_flux_aperture(img_avg, [y_index_center_avg, x_index_center_avg], radi_apert)
                        flux_main_unfilt = get_flux_aperture(img_avg_unfilt, [y_index_center_avg, x_index_center_avg], radi_apert)


                        #------
                        # Get 5 sigma values for flux of apertures for pre-subtracted, mdn-ed img
                        #------
                        filename_detlim_init = date + '_set' + str(int(setnumb1)) + '_detlim_init.fits'
                        arr_5sd_pre = []
                        arr_flux_bin_pre = []
                        for rad_mock in arr_radius:
                                arr_index_pre = get_indexap_annulus(rad_mock, radi_apert)
                                arr_5sd_radi_pre = []
                                for elem in arr_index_pre:
                                        y_apert_center_pre = elem[0] + y_index_center_avg
                                        x_apert_center_pre = elem[1] + x_index_center_avg
                                        array_flux_apert_radi_pre = get_flux_aperture(img_avg_unfilt, [y_apert_center_pre, x_apert_center_pre], radi_apert)
                                        arr_5sd_radi_pre.append(array_flux_apert_radi_pre)                        
                                sd_ring_pre = np.std(np.array(arr_5sd_radi_pre))
                                arr_5sd_pre.append(sd_ring_pre*5)
                        arr_mag_5sd_pre = 2.5*np.log10(flux_main_unfilt/np.array(arr_5sd_pre))
                        save_fits(arr_mag_5sd_pre, filename_detlim_init)
                        plt.plot(arr_radius, arr_mag_5sd_pre, 'ko-')
                        print 'created:', filename_detlim_init


                        #------
                        # Get 5 sigma values for flux of apertures for pre-subtracted, mdn-ed, filtered img
                        #------
                        filename_detlim_filt_init = date + '_set' + str(int(setnumb1)) + '_detlim_filt_init.fits'
                        arr_5sd_pre_filt = []
                        arr_flux_bin_pre_filt = []
                        for rad_mock in arr_radius:
                                arr_index_pre_filt = get_indexap_annulus(rad_mock, radi_apert)
                                arr_5sd_radi_pre_filt = []
                                for elem in arr_index_pre_filt:
                                        y_apert_center_pre_filt = elem[0] + y_index_center_avg
                                        x_apert_center_pre_filt = elem[1] + x_index_center_avg
                                        array_flux_apert_radi_pre_filt = get_flux_aperture(img_avg, [y_apert_center_pre_filt, x_apert_center_pre_filt], radi_apert)
                                        arr_5sd_radi_pre_filt.append(array_flux_apert_radi_pre_filt)                        
                                sd_ring_pre_filt = np.std(np.array(arr_5sd_radi_pre_filt))
                                arr_5sd_pre_filt.append(sd_ring_pre_filt*5)
                        arr_mag_5sd_pre_filt = 2.5*np.log10(flux_main_unfilt/np.array(arr_5sd_pre_filt))
                        save_fits(arr_mag_5sd_pre_filt, filename_detlim_filt_init)
                        plt.plot(arr_radius, arr_mag_5sd_pre_filt, 'mo-')
                        print 'created:', filename_detlim_filt_init

                        #------
                        # Iterate through radii for mock binaries
                        # Check medianed, LOCI-subtracted img(without binary) for
                        # 5-sigma value at that radius
                        # Create array of corresponding flux ratio
                        #------
                        arr_5sd = []
                        arr_flux_bin_init = []
                        for rad_mock in arr_radius:
                                arr_index = get_indexap_annulus(rad_mock, radi_apert)
                                arr_flux_radi = []
                                for elem in arr_index:
                                        y_apert_center = elem[0] + y_index_center_avg
                                        x_apert_center = elem[1] + x_index_center_avg
                                        if filt:
                                                array_flux_apert_radi = get_flux_aperture(img_subed_filt, [y_apert_center, x_apert_center], radi_apert)
                                        else:
                                                array_flux_apert_radi = get_flux_aperture(img_subed, [y_apert_center, x_apert_center], radi_apert)
                                        arr_flux_radi.append(array_flux_apert_radi)
                                #print 'arr_flux_radi', arr_flux_radi
                                sd_ring_r = np.std(np.array(arr_flux_radi))
                                arr_5sd.append(sd_ring_r*5)
                        if filt:
                                arr_mag_5sd = 2.5*np.log10(flux_main_filt/arr_5sd)
                        else:
                                arr_mag_5sd = 2.5*np.log10(flux_main_unfilt/arr_5sd)
                        if filtername[0] =='K':
                                filtername_file = 'Ks'
                        elif filtername[0] =='B':
                                filtername_file = 'Brg'
                        else:
                                print 'filtername:', filtername
                                print 'PROBLEM. Filter name doesnt start with K or B'
                        filename_correction = date + '_' + filtername_file + '_final_correction.fits'                 
                        if exists(filename_correction):
                                arr_correction = pf.open(filename_correction)[0].data
                        else:
                                print filename_correction, 'doesnt exist'
                                continue
                        '''
                        if filtername[0] =='K':                        
                                if counter_ks < 1: 
                                        arr_plots += plt.plot(arr_radius, arr_mag_5sd - arr_correction, arr_colors[int(index_date*2)] + 'o-')
                                        arr_labels.append(arr_date_ShaneAO[index_date] + ', ' + 'Ks')
                                else:
                                        plt.plot(arr_radius, arr_mag_5sd - arr_correction, arr_colors[int(index_date*2)] + 'o-')
                                counter_ks +=1 
                        if filtername[0] =='B':                        
                                if counter_brg < 1: 
                                        arr_plots += plt.plot(arr_radius, arr_mag_5sd - arr_correction, arr_colors[int(index_date*2)+1] + 'o-')
                                        arr_labels.append(arr_date_ShaneAO[index_date] + ', ' + 'Brg')
                                else:
                                        plt.plot(arr_radius, arr_mag_5sd - arr_correction, arr_colors[int(index_date*2)+1] + 'o-')
                                counter_brg +=1
                        '''
                        #save final detection limit
                        save_fits(arr_mag_5sd - arr_correction, filename_detlim_final)
                        counter += 1
                        #plt.gca().invert_yaxis()
                        #plt.show()
                #if counter > 7:
        plt.gca().invert_yaxis()
        plt.ylabel('Magnitude difference')
        plt.xlabel('Radius from primary (pixels)')
        plt.legend(arr_plots, arr_labels)
        plt.show()


def plot_correction(date = date):
        #------
        # load array with radii at which we had mock up binaries
        #------
        filename_arrradiusmock = directory + '/' + date + '/' + 'arr_radiusmock.fits'
        arr_radius = pf.open(filename_arrradiusmock)[0].data

        dict_correction = {}
        arr_correction = None
        arr_correction_err = None                
        arr_colors = ['r', 'k', 'g', 'm', 'b', 'y']
        arr_labels = []
        arr_plots = []
        for index_date in range(len(arr_date_ShaneAO)):
                counter = 0
                print 'index_date', index_date
                date = arr_date_ShaneAO[index_date]
                filename_arr_correction = 'arr_correction_' + date + '.fits'
                filename_arr_correction_err = 'arr_correction_err_' + date + '.fits'

                #if arr_correction == None:
                        
                arr_correction = pf.open(filename_arr_correction)[0].data
                arr_correction_err = pf.open(filename_arr_correction_err)[0].data

                #else:
                #        arr_correction = np.concatenate((arr_correction, pf.open(filename_arr_correction)[0].data), axis = 0)
                #        arr_correction_err = np.concatenate((arr_correction_err, pf.open(filename_arr_correction_err)[0].data), axis = 0)

                print 'arr_correction.shape', arr_correction.shape
                print 'arr_correction_err.shape', arr_correction_err.shape

                #------
                # Plot each row in matrix of corrections
                # Will crash if sets arent in order , i.e. 1 - 12 ok, but 1-12 & 14 not ok. ***
                #------
                counter_ks = 0
                counter_brg = 0
                for index_row in range(arr_correction.shape[1]):
                        if ret_filtername(index_row+1, date = arr_date_ShaneAO[index_date])[0] == 'K':
                                if counter_ks < 1: 
                                        arr_plots += plt.plot(arr_radius, arr_correction[index_row], arr_colors[int(index_date*2)] + 'o-')
                                        arr_labels.append(arr_date_ShaneAO[index_date] + ', ' + 'Ks')
                                        dict_correction[arr_date_ShaneAO[index_date] + '_Ks'] = []
                                        dict_correction[arr_date_ShaneAO[index_date] + '_Ks'].append(arr_correction[index_row])
                                else:
                                        plt.plot(arr_radius, arr_correction[index_row], arr_colors[int(index_date*2)] + 'o-')
                                        dict_correction[arr_date_ShaneAO[index_date] + '_Ks'].append(arr_correction[index_row])
                                counter_ks +=1 
                        elif ret_filtername(index_row+1, date = arr_date_ShaneAO[index_date])[0] == 'B':
                                if counter_brg < 1: 
                                        arr_plots += plt.plot(arr_radius, arr_correction[index_row], arr_colors[int(index_date*2) + 1] + 'o-')
                                        arr_labels.append(arr_date_ShaneAO[index_date] + ', ' + 'Brg')
                                        dict_correction[arr_date_ShaneAO[index_date] + '_Brg'] = []
                                        dict_correction[arr_date_ShaneAO[index_date] + '_Brg'].append(arr_correction[index_row])
                                else:
                                        plt.plot(arr_radius, arr_correction[index_row], arr_colors[int(index_date*2) + 1] + 'o-')
                                        dict_correction[arr_date_ShaneAO[index_date] + '_Brg'].append(arr_correction[index_row])
                                counter_brg += 1
                        
        print 'plotting correction curves across runs...'
        plt.legend(arr_plots, arr_labels)
        plt.show()
        plt.close()
        
        
        #------
        # loop through dictionary of each date & filter
        # take mean of correction curves 
        # save as fits file
        #------
        arr_lables = []
        arr_plots = []
        for key_correction in dict_correction:
                filename_arr_correction_final = key_correction + '_final_correction.fits'
                mat_correction = dict_correction[key_correction]
                arr_correction_final = np.mean(mat_correction, axis = 0)
                print 'arr_correction_final.shape', arr_correction_final.shape
                save_fits(arr_correction_final, filename_arr_correction_final)
                print 'created:', filename_arr_correction_final
                arr_plots += plt.plot(arr_radius, arr_correction_final, 'o-')
                arr_labels += key_correction
                plt.legend(arr_plots, arr_labels)
        plt.show()
        
                

        
        print 'shouldve plotted by now...'

        '''
        avg_error = np.mean(arr_correction_err, axis = 0)
        std_correction = np.std(arr_correction, axis = 0)
        plt.plot(arr_radius, avg_error, 'ro-')
        plt.plot(arr_radius, std_correction, 'k^-')
        plt.show()
        
        for index_row in range(arr_correction.shape[1]):
                plt.plot(arr_radius, arr_correction[index_row], 'b')
                #print arr_correction_err[index_row]
                #print 'arr_correction', arr_correction[index_row]
                plt.errorbar(arr_radius, arr_correction[index_row], yerr = arr_correction_err[index_row], ecolor = 'g')


        plt.show()
        '''




def plot_detlim_loci(filt = True, total_sets = total_sets):
        if filt:
                total_sets = 12

        #------
        #Define radius of aperture
        #------
        radi_apert = 3
        startset = 1
        arr_setnumbs = np.arange(startset, total_sets+1).astype(int)
        

        #------
        # define mock binary radii
        # save to fits file
        #------
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80]).astype(int)
        filename_arrradiusmock = directory + '/' + date + '/' + 'arr_radiusmock.fits'
        save_fits(arr_radiusmock, filename_arrradiusmock)
        print 'arr_radiusmock', arr_radiusmock

        arr_theta = np.arange(0, 360, 60).astype(int)
        print 'arr_theta', arr_theta


        #------
        # Create or open stacked correction curves & uncertainties of these curves
        #------
        filename_arr_correction = 'arr_correction_' + date + '.fits'
        filename_arr_correction_err = 'arr_correction_err_' + date + '.fits'
        arr_correction_save = np.empty(shape = (0, arr_radiusmock.size))
        arr_correction_err_save = np.empty(shape = (0, arr_radiusmock.size))
        '''
        if exists(filename_arr_correction):
                arr_correction_save = pf.open(filename_arr_correction)[0].data
        else:
                arr_correction_save = np.empty(shape = (0, arr_radiusmock.size))
        if exists(filename_arr_correction_err):
                arr_correction_err_save = pf.open(filename_arr_correction_err)[0].data
        else:
                arr_correction_err_save = np.empty(shape = (0, arr_radiusmock.size))
        '''


        arr_final_err_save = np.empty(shape = (arr_radiusmock.size))
        arr_5sd_correct = []
        if filt:
                filename_arr5sdcorrectfits = directory + '/' + date + '/' + 'arr_5sd_correct_filt.fits'
        else:
                filename_arr5sdcorrectfits = directory + '/' + date + '/' + 'arr_5sd_correct.fits'
        for setnumb1 in arr_setnumbs:
                print '------'
                print 'set', setnumb1
                print '------'
                

                #------
                # Open LOCI-subtracted, medianed(and perhaps also filtered) img for set number
                # Define img name for plots
                #------
                if filt:
                        filename_plot = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_plotfilt.png'
                        filename_fits_filt = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_locifiltfinal.fits'
                        img_subed_filt = pf.open(filename_fits_filt)[0].data
                else:
                        filename_plot = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_plot.png'



                #------
                #Define and load img after loci without high pass filter
                #------
                filename_fits = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_locifinal.fits'
                if exists(filename_fits):
                        img_subed = pf.open(filename_fits)[0].data
                else:
                        print filename_fits, 'doesnt exist'
                        continue

                

                #------
                # Open pre-subtraction, medianed img for set number
                # Obtain flux of aperture around primary
                #------
                if filt:
                        filename_avg = directory + '/' + date + '/' + 'centroid_filt_' + str(int(setnumb1)) + '.fits'
                else:
                        filename_avg = directory + '/' + date + '/' + 'centroid_' + str(int(setnumb1)) + '.fits'
                if not exists(filename_avg):
                        print "filename_avg doesn't exist"
                img_avg = pf.open(filename_avg)[0].data
                img_avg /= np.amax(img_avg)
                y_length_avg, x_length_avg = img_avg.shape
                y_index_center_avg = int((y_length_avg - 1)/2)
                x_index_center_avg = int((x_length_avg - 1)/2)
                flux_main = get_flux_aperture(img_avg, [y_index_center_avg, x_index_center_avg], radi_apert)
                print 'flux_main', flux_main



                #------
                # Get flux of stacked, medianed primary without high pass filter
                #------
                filename_avg_unfilt = directory + '/' + date + '/' + 'centroid_' + str(int(setnumb1)) + '.fits'
                if not exists(filename_avg_unfilt):
                        print "filename_avg_unfilt doesn't exist"
                img_avg_unfilt = pf.open(filename_avg_unfilt)[0].data
                img_avg_unfilt /= np.amax(img_avg_unfilt)
                y_length_avg_unfilt, x_length_avg_unfilt = img_avg_unfilt.shape
                y_index_center_avg_unfilt = int((y_length_avg_unfilt - 1)/2)
                x_index_center_avg_unfilt = int((x_length_avg_unfilt - 1)/2)
                flux_main_unfilt = get_flux_aperture(img_avg_unfilt, [y_index_center_avg_unfilt, x_index_center_avg_unfilt], radi_apert)



                #------
                # Get 5 sigma values for flux of apertures for pre-subtracted, mdn-ed img
                #------
                arr_5sd_pre = []
                arr_flux_bin_pre = []
                for rad_mock in arr_radiusmock:
                        arr_index_pre = get_indexap_annulus(rad_mock, radi_apert)
                        arr_5sd_radi_pre = []
                        for elem in arr_index_pre:
                                y_apert_center_pre = elem[0] + y_index_center_avg
                                x_apert_center_pre = elem[1] + x_index_center_avg
                                array_flux_apert_radi_pre = get_flux_aperture(img_avg, [y_apert_center_pre, x_apert_center_pre], radi_apert)
                                arr_5sd_radi_pre.append(array_flux_apert_radi_pre)
                        
                        sd_ring_pre = np.std(np.array(arr_5sd_radi_pre))
                        arr_5sd_pre.append(sd_ring_pre*5)
                arr_mag_5sd_pre = 2.5*np.log10(flux_main/np.array(arr_5sd_pre))




                #################
                # TESTING APERTURES AROUND ANNULI
                '''
                for rad_mock in arr_radiusmock:
                        img_test = np.zeros([160, 160])
                        arr_index = get_indexap_annulus(rad_mock, radi_apert)
                        print 'arr_index', arr_index
                        for elem in arr_index:
                                y = elem[0] + y_index_center_avg
                                x = elem[1] + x_index_center_avg
                                img_test = test_get_flux_aperture(img_test, [y,x], radi_apert)
                        hdu = pf.PrimaryHDU(img_test)
                        hdulist = pf.HDUList([hdu])
                        hdulist.writeto('test.fits', clobber=True)
                        subprocess.call('/home/apps/ds9/ds9 ' + 'test.fits', shell = True)    
                '''
                ###################


                #------
                # Run get_flux_aperture to find number of pixels in aperture
                #------
                img_zeros = np.zeros(img_subed.shape)
                arr_empty_apert = get_flux_aperture(img_zeros, [y_index_center_avg, x_index_center_avg], radi_apert, True)
                num_pix_apert = len(arr_empty_apert)



                #------
                #Iterate through radii for mock binaries
                #Check medianed, LOCI-subtracted img(without binary) for 5 sigma value at that radius
                #Create array of corresponding flux ratio
                #------
                arr_5sd = []
                arr_sd = []
                arr_flux_bin_init = []
                for rad_mock in arr_radiusmock:
                        arr_index = get_indexap_annulus(rad_mock, radi_apert)
                        arr_flux_radi = []
                        for elem in arr_index:
                                y_apert_center = elem[0] + y_index_center_avg
                                x_apert_center = elem[1] + x_index_center_avg
                                if filt:
                                        array_flux_apert_radi = get_flux_aperture(img_subed_filt, [y_apert_center, x_apert_center], radi_apert)
                                else:
                                        array_flux_apert_radi = get_flux_aperture(img_subed, [y_apert_center, x_apert_center], radi_apert)
                                arr_flux_radi.append(array_flux_apert_radi)
                        #print 'arr_flux_radi', arr_flux_radi
                        sd_ring_r = np.std(np.array(arr_flux_radi))
                        print 'sd_ring_r', sd_ring_r
                        arr_sd.append(sd_ring_r)
                        arr_5sd.append(sd_ring_r*5)

                        struct_ring = check_ring(img_subed, rad_mock)
                        sd_ring_formocks = struct_ring['sd_ring']
                        mock_factor = sd_ring_formocks*50
                        img_mock = pf.open(directory + '/' + date + '/' + 'centroid_1.fits')[0].data
                        img_mock /= np.amax(img_mock)
                        img_mock *= float(mock_factor)
                        flux_bin_init = get_flux_aperture(img_mock, [y_index_center_avg, x_index_center_avg], radi_apert)
                        arr_flux_bin_init.append(flux_bin_init)
                        
                arr_flux_bin_init = np.array(arr_flux_bin_init)
                arr_fluxratio_init = flux_main_unfilt/arr_flux_bin_init
                arr_mag_5sd = 2.5*np.log10(flux_main/np.array(arr_5sd))
                arr_5sd_err = (flux_main/np.array(arr_sd))/(flux_main/np.array(arr_5sd)) ###Fix this!!!


                counter = 0
                arr_avg_fluxbin = []
                arr_sd = []
                for radius_mock in arr_radiusmock: #iterate through different radii from center
                        arr_delflux = []
                        arr_fluxbin_r = []
                        arr_flux_thetas_r = []
                        for theta_mock in arr_theta: #iterate through different position angles

                                #------
                                # Open image with mock binary after loci subtraction
                                # If file doesn't exist, skip img
                                #------
                                if filt:
                                        filename_final = directory + '/' + date + '/' + 'set'+ str(int(setnumb1)) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfiltfinal.fits'
                                else:
                                        filename_final = directory + '/' + date + '/' + 'set'+ str(int(setnumb1)) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfinal.fits'


                                if exists(filename_final):
                                        try:
                                                img_final = pf.open(filename_final)[0].data
                                                counter += 1
                                        except:
                                                print 'error opening:', filename_final
                                                subprocess.call('rm -rf ' + filename_final, shell = True)
                                                continue
                                else:
                                        print 'File not found:', filename_final
                                        continue



                                #------
                                # find x and y index at center of image
                                # then use radius_mock and theta_mock to get index of mock binary
                                #------
                                y_length, x_length = img_final.shape
                                y_index_center = int((y_length - 1)/2)
                                x_index_center = int((x_length - 1)/2)

                                dx_mock = radius_mock*np.cos(theta_mock*(2*np.pi)/360)
                                dx_mock = int(round(dx_mock))
                                dy_mock = radius_mock*np.sin(theta_mock*(2*np.pi)/360)
                                dy_mock = int(round(dy_mock))

                                y_index_bin = y_index_center + dy_mock
                                x_index_bin = x_index_center + dx_mock
                                flux_bin = get_flux_aperture(img_final, [y_index_bin, x_index_bin], radi_apert)
                                arr_fluxbin_r.append(flux_bin)
                        
                        arr_avg_fluxbin.append(np.mean(arr_fluxbin_r))
                        arr_fluxbin_r = flux_main/arr_fluxbin_r
                        sd_r = np.std(arr_fluxbin_r)
                        arr_sd.append(sd_r)

                print 'arr_sd', arr_sd
                print 'arr_5sd_err', arr_5sd_err
                arr_mag_final_err = 2.5*np.log10(np.sqrt( (np.array(arr_sd)**2.) + (np.array(arr_5sd_err)**2.) ))

                arr_fluxratio_final = flux_main/arr_avg_fluxbin
                arr_mag_init = 2.5*np.log10(arr_fluxratio_init)
                arr_mag_final = 2.5*np.log10(arr_fluxratio_final) 
                arr_mag_diff = arr_mag_final - arr_mag_init #correction curve, diff btwn binary flux before injection and after loci
                print 'arr_mag_final_err', arr_mag_final_err
                print 'arr_mag_diff', arr_mag_diff
                print 'arr_mag_sd', arr_sd
                print '-------'
                
                #------
                # Plot curves
                #------
                
                mag_pre, = plt.plot(arr_radiusmock, arr_mag_5sd_pre, 'ko-') #5 s.d. before subtraction
                mag_init, = plt.plot(arr_radiusmock, arr_mag_init, 'bo-') #Mag. of mock binary before injection
                mag_afterloci, = plt.plot(arr_radiusmock, arr_mag_final, 'ro-') #Mag. of binary after LOCI
                mag_5sd, = plt.plot(arr_radiusmock, arr_mag_5sd,'yo-') #5 s.d. before correction
                mag_5sd_correct, = plt.plot(arr_radiusmock, arr_mag_5sd - arr_mag_diff,'go-') #5 s.d. after correction
                plt.errorbar(arr_radiusmock, arr_mag_5sd - arr_mag_diff, yerr = arr_mag_final_err, ecolor = 'g')
                plt.gca().invert_yaxis()
                plt.legend([mag_pre, mag_init, mag_afterloci, mag_5sd, mag_5sd_correct], ['5 s.d. before subtraction', 'Mag. of mock binary before injection', 'Mag. of binary after LOCI', '5 s.d. before correction', '5 s.d. after correction'])
               
                plt.show()


                #------
                # Append corrected 5sd plots to rows in a table. Save as fits file.
                #-----
                if counter > 50:
                        arr_5sd_correct.append(arr_mag_5sd - arr_mag_diff)
                        arr_correction_save = np.append(arr_correction_save, np.array([arr_mag_diff]), axis = 0)
                        arr_correction_err_save = np.append(arr_correction_err_save, np.array([arr_mag_sd]), axis = 0)
                print 'number of images used:', counter
                #mng = plt.get_current_fig_manager()                                         
                #mng.resize(*mng.window.maxsize())
                #plt.show()
                #plt.savefig(filename_plot, bbox_inches='tight')
                #print 'saved plot as img:', filename_plot
                #plt.close()
                

                hdu = pf.PrimaryHDU(np.array(arr_correction_save))
                hdulist = pf.HDUList([hdu])
                hdulist.writeto(filename_arr_correction, clobber=True)
                print 'created:', filename_arr_correction

                hdu = pf.PrimaryHDU(np.array(arr_correction_err_save))
                hdulist = pf.HDUList([hdu])
                hdulist.writeto(filename_arr_correction_err, clobber=True)
                print 'created:', filename_arr_correction_err

        hdu = pf.PrimaryHDU(np.array(arr_5sd_correct))
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_arr5sdcorrectfits, clobber=True)
        print 'created:', filename_arr5sdcorrectfits        



def plot_all_5sd():
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80]).astype(int)
        arr_5sd = pf.open('ShaneAO/Jun15/arr_5sd_correct.fits')[0].data
        for row_index in np.arange(len(arr_5sd)):
                plt.plot(arr_radiusmock, arr_5sd[row_index, :], 'o-')                        
        plt.gca().invert_yaxis()
        plt.show()



def get_indexap_annulus(radi_annulus, radi_apert):
        #------
        # Input: Distance from center of image at which to center apertures, aperture radius
        #------
        # Output: Array of indices. These indices are the y & x coordinates of aperture centers
        #------


        theta_incr = math.acos(1 - (((radi_apert+0.55)**2.)/(2*((radi_annulus)**2.))))
        theta_incr *= 2
        #print '-----'
        #print 'theta_increase(degrees)', theta_incr*180/math.pi


        arr_index = []
        theta_prev = -1*math.pi
        theta_loop = -1*math.pi + theta_incr
        arr_index.append([0, -1*radi_annulus])
        cond_incircle = True
        while cond_incircle:
                #print 'theta_loop', theta_loop*180/math.pi
                index_y = radi_annulus*math.sin(theta_loop)
                index_y_int = round(index_y)
                index_x = radi_annulus*math.cos(theta_loop)
                index_x_int = round(index_x)
                arr_coord = [(y, x) for y in np.arange(index_y_int - 2, index_y_int + 2 + 1).astype(int) for x in np.arange(index_x_int - 2, index_x_int + 2 + 1).astype(int)]
                #print 'arr_coord', arr_coord
                arr_dist = []
                for coord in arr_coord:
                        y = coord[0]
                        x = coord[1]
                        dist = distance(y, index_y, x, index_x)
                        arr_dist.append(dist)
                cond_safetheta = False
                counter = 0 ####
                while not cond_safetheta:
                        index_min = arr_dist.index(min(arr_dist))
                        min_y = arr_coord[index_min][0]
                        min_x = arr_coord[index_min][1]
                        #print 'min_y', min_y, 'min_x', min_x
                        theta_min = math.atan2(min_y,min_x)
                        #print 'theta_min_index(degrees)', theta_min*180/math.pi
                        #print 'theta_min', theta_min*180/math.pi
                        #print 'theta_prev', theta_prev*180/math.pi                                
                        if theta_min < 0 and theta_prev > 0:
                                break
                        safe_theta = theta_min - theta_prev > theta_incr                        
                        #print 'safe theta? :', safe_theta
                        safe_radi = radi_annulus - 0.5 < distance(min_y, 0, min_x, 0) < radi_annulus + 0.5
                        #print 'safe radi? :', safe_radi
                        if safe_theta and safe_radi:
                                cond_safetheta = True
                        else:
                                'change in angle not large enough, next.'
                                arr_dist[index_min] = 10
                arr_index.append([index_y_int, index_x_int])
                #print 'theta_min', theta_min*180/math.pi
                #print 'theta_prev', theta_prev*180/math.pi                                
                theta_loop = theta_min + theta_incr
                theta_prev = theta_min
                if math.pi - theta_loop < theta_incr:
                        cond_incircle = False
        return arr_index


        





def get_flux_aperture(img, index, radi_ap, ret_array=False):
        #Inputs: image(2d array), indexes(2-element(y,x) array), radius of aperture(number)
        #Returns: total flux in aperture
        #------

        #------
        #deciphering inputs
        #------
        index_y = index[0]
        index_x = index[1]
        j = radi_ap

        #------
        #Knowing index of peak, take aperture of j pixel around this point
        #save total flux of aperture
        #------
        arr_totflux = []
        for x_circ in np.arange((2*j) + 1) - j:
                for y_circ in np.arange((2*j) + 1) - j:
                        if distance(x_circ, 0, y_circ, 0) < (j+0.5):
                                if index_y + y_circ >= img.shape[0] or index_x + x_circ >= img.shape[1]:
                                        continue
                                else:
                                        arr_totflux.append(img[index_y+y_circ, index_x+x_circ])
        tot_flux = sum(arr_totflux)
        if ret_array:
                return arr_totflux
        else:
                return tot_flux


def get_com_aperture(img, index, radi_ap):
        #Returns center of mass in aperture of variable size
        #Inputs: image(2d array), indexes(2-element(y,x) array), radius of aperture(number)
        #Returns: y and x coordinates of center of mass
        #------

        #------
        #deciphering inputs
        #------
        index_y = index[0]
        index_x = index[1]
        j = int(radi_ap)

        #------
        #Knowing index of center, take aperture of j pixels in radius around this point
        #Calculate center of mass pixel, return updated index
        #------
        arr_totflux = []
        arr_com = []
        for x_circ in np.arange((2*j) + 1) - j:
                for y_circ in np.arange((2*j) + 1) - j:
                        if distance(x_circ, 0, y_circ, 0) < (j+0.5):
                                if index_y + y_circ >= img.shape[0] or index_x + x_circ >= img.shape[1]:
                                        continue
                                else:
                                        arr_com.append([y_circ, x_circ, img[index_y+y_circ, index_x+x_circ]])
        arr_y = []
        arr_x = []           
        arr_tot = []
        for elem in arr_com:
                y_index = elem[0]
                x_index = elem[1]
                value_elem = elem[2]
                arr_y.append(y_index*value_elem)
                arr_x.append(x_index*value_elem)
                arr_tot.append(value_elem)
        tot_flux = sum(arr_tot)
        com_y = sum(arr_y)/tot_flux 
        com_x = sum(arr_x)/tot_flux 
        com_y += index_y #center of mass of x index in img
        com_x += index_x #center of mass of y index in img
        return (com_y, com_x)


def fourier_shift(img, (y,x)):
# Fourier shift img by y and x 
# Inputs: img (2D img),  (y,x)  (array/tuple with 2 elements, vertical and horizontal shift)
#------
        img = np.fft.fft2(img)
        img = fshift(img, [y, x])
        img = np.fft.ifft2(img)
        img = np.real(img)
        return img


def test_get_flux_aperture(img, index, radi_ap):
        #Inputs: image(2d array), indexes(2-element(y,x) array), radius of aperture(number)
        #Returns: total flux in aperture
        #------

        #------
        #deciphering inputs
        #------
        index_y = index[0]
        index_x = index[1]
        j = radi_ap

        #------
        #Knowing index of peak, take aperture of j pixel around this point
        #save total flux of aperture
        #------
        arr_totflux = []
        for x_circ in np.arange((2*j) + 1) - j:
                for y_circ in np.arange((2*j) + 1) - j:
                        if distance(x_circ, 0, y_circ, 0) < (j+0.5):
                                if index_y + y_circ >= img.shape[0] or index_x + x_circ >= img.shape[1]:
                                        continue
                                else:
                                        img[index_y+y_circ, index_x+x_circ] += 1

        return img

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
        
	filter_of_set = ret_filtername(start1, date) #check for filter of 1st img in set

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

                        img_filter = ret_filtername(j, date_ref)
                        
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



        

def run_theklip():
        #-------
        #Define parameters
        #------
        numb_imgs_keep = 600 #number of imgs to keep from correlation matrix
        numb_pixels_all_imgs = 25921 #number of pixels in all imgs
        max_rad_img = 80 # max radius in imgs
        rad_annulus = 10 # radius of annulus for KLIP search rings
        fwhm_threshold = 10 # upper limit on full-width-half-max
	maxpix_threshold = 35000 # exclude imgs with pixel values above this number
	minpix_threshold = 5000 # exclude imgs with pixel values below this number
        arr_klipmodes = np.array([15, 25, 35]).astype(int) # array with number of klip modes to try


        

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
        #open txt file with list of known binaries
        #------
        arr_bin = open(directory + '/' + directory + '_clearbinaries.txt', 'rb').read().splitlines()
        arr_bin_dates = np.array(arr_bin[0::2])
        arr_bin_setnumbs = np.array(arr_bin[1::2])
        #print 'arr_bin_dates', arr_bin_dates
        #print 'arr_bin_setnumbs', arr_bin_setnumbs

        #------
        #check for filter of 1st img in set
        #------
	filter_of_set = ret_filtername(start1, date) 
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
                                print date_ref, 'set', setnumb2, 'skipped'
                                continue
                        elif str(int(setnumb2)) in arr_bin_setnumbs:
                                arr_index_setnumb_bin = np.where(arr_bin_setnumbs == str(int(setnumb2)))[0]
                                #print arr_index_setnumb_bin
                                '''
                                print '------'
                                print 'index_setnumb', index_setnumb
                                print 'set number', setnumb2
                                print 'corresponding date in txt file', arr_bin_dates[index_setnumb]
                                print 'date in loop', date_ref
                                '''
                                for index_setnumb_bin in arr_index_setnumb_bin:
                                        if arr_bin_dates[index_setnumb_bin] == date_ref:
                                                print date_ref, 'set', setnumb2, 'skipped since it is a known binary'
                                                continue
                        else:
                                ####date_ref and setnumb2
                                #print 'adding elements in set no.', setnumb2, 'for date', date_ref
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
                                
                        img_filter = ret_filtername(j, date_ref)
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
        print 'Now creating arr_pool_vars, for feeding into KLIP function'
        arr_pool_vars = []	
        for index_i in np.arange(arr_targpsf.size):
                arr_j_optimal = []
                arr_img2_optimal = []
                i = arr_targpsf[index_i]
                print i
                index_i = np.argwhere(arr_filenumbs == i)[0][0]
                #print 'index_i', index_i
		arr_corrs = matrix_corr[(index_date_targ*max_numb_imgs) + index_i, :]
                index_sort_corr = np.argsort(arr_corrs)
                #print 'index_sort_corr', index_sort_corr
                arr_indexdateoptimal = arr_indexdate_corrs[index_sort_corr]
                #print 'arr_indexdateoptimal', arr_indexdateoptimal
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
                        
                #print 'mean of arr_img2_optimal should be nearly zero:'
                #print 'np.mean(arr_img2_optimal):', np.mean(arr_img2_optimal)
                for numb_klip_modes in arr_klipmodes:
                        struct = {'filenumb_img1': i, 'arr_img2': arr_img2_optimal, 'date_ref': date_ref, 'numb_klip_modes':numb_klip_modes, 'max_rad_img': max_rad_img, 'rad_annulus': rad_annulus}
                        arr_pool_vars.append(struct)
                                                
        
                                                
        #------
        # save/load list of structures as needed
        #------
        filename_arrpoolvars = 'klip_arr_pool_vars_' + date + 'set' + str(int(setnumb1)) + '.p'
        #arr_pool_vars = pickle.load(open(filename_arrpoolvars, 'rb')) #loads list of structures

        print '------------------------'
	print 'No of loops to run:', len(arr_pool_vars)
        print '------------------------'
        #pickle.dump(arr_pool_vars, open(filename_arrpoolvars, 'wb')) #saves list of structures
        


        #------
        # run KLIP
        #------
        print 'Now starting threads...'
        
	for elem in arr_pool_vars:
		theklip(elem)
        
        '''
	threadnumb = 6
	thepool = ThreadPool(threadnumb)
        arr_pool = thepool.map(klip, arr_pool_vars)
        thepool.close()
        thepool.join()
        '''

        


        
def theklip(struct):
        #------
        #Extract variables from input structure
        #------
        filenumb_img1 = struct['filenumb_img1']
        arr_img2 = struct['arr_img2']
        date_ref = struct['date_ref']
        numb_klip_modes = struct['numb_klip_modes']
        max_rad_img = struct['max_rad_img']
        rad_annulus = struct['rad_annulus']

        arr_radstarts = np.arange(1, max_rad_img+1, rad_annulus)
        print 'starting:', filenumb_img1, 'kmode', int(numb_klip_modes), datetime.datetime.now().strftime('%m-%d %H:%M:%S')


        
        #-------
        # Define filenames of fits files
        # Exit function if files already exist
        # CHANGE FILENAMES IF NECESSARY        
        #-------
        filename_output = directory + '/' + date + '/' + str(filenumb_img1) + '_klipfinal_' + 'kmode' + str(int(numb_klip_modes)) + '.fits'

        

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


                
        #----------
        #flatten img and save img shape
        #----------
        img1shape = img1.shape
        img1_flat = img1.flatten()
        radi = 0



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
        

        print 'created:', filename_output
        print 'done with:', filenumb_img1, 'kmode', int(numb_klip_modes), datetime.datetime.now().strftime('%m-%d %H:%M:%S')
        return 1



def run_klip2():
        #-------
        #Define parameters
        #------
        numb_imgs_keep = 600 #number of imgs to keep from correlation matrix
        numb_pixels_all_imgs = 25921 #number of pixels in all imgs
        max_rad_img = 80 # max radius in imgs
        rad_annulus = 10 # radius of annulus for KLIP search rings
        fwhm_threshold = 10 # upper limit on full-width-half-max
	maxpix_threshold = 35000 # exclude imgs with pixel values above this number
	minpix_threshold = 5000 # exclude imgs with pixel values below this number
        arr_klipmodes = np.array([15, 25, 35]).astype(int) # array with number of klip modes to try


        
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
        #open txt file with list of known binaries
        #------
        arr_bin = open(directory + '/' + directory + '_clearbinaries.txt', 'rb').read().splitlines()
        arr_bin_dates = np.array(arr_bin[0::2])
        arr_bin_setnumbs = np.array(arr_bin[1::2])
        #print 'arr_bin_dates', arr_bin_dates
        #print 'arr_bin_setnumbs', arr_bin_setnumbs



        #------
        #check for filter of 1st img in set
        #------
	filter_of_set = ret_filtername(start1, date) 
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
                                print date_ref, 'set', setnumb2, 'skipped'
                                continue
                        elif str(int(setnumb2)) in arr_bin_setnumbs:
                                arr_index_setnumb_bin = np.where(arr_bin_setnumbs == str(int(setnumb2)))[0]
                                #print arr_index_setnumb_bin
                                '''
                                print '------'
                                print 'index_setnumb', index_setnumb
                                print 'set number', setnumb2
                                print 'corresponding date in txt file', arr_bin_dates[index_setnumb]
                                print 'date in loop', date_ref
                                '''
                                for index_setnumb_bin in arr_index_setnumb_bin:
                                        if arr_bin_dates[index_setnumb_bin] == date_ref:
                                                print date_ref, 'set', setnumb2, 'skipped since it is a known binary'
                                                continue
                        else:
                                ####date_ref and setnumb2
                                #print 'adding elements in set no.', setnumb2, 'for date', date_ref
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



        '''
        #------
        # Add reference filenumbs from all dates
        #------
	filter_of_set = ret_filtername(start1, date) #check for filter of 1st img in set
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
                arr_substar.append(arr_substar_temp)
        print 'arr_substar', arr_substar
        '''


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
                                
                        img_filter = ret_filtername(j, date_ref)
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
        print 'Now creating arr_pool_vars, for feeding into KLIP function'
        arr_pool_vars = []	
        for index_i in np.arange(arr_targpsf.size):
                arr_j_optimal = []
                arr_img2_optimal = []
                i = arr_targpsf[index_i]
                print i
                index_i = np.argwhere(arr_filenumbs == i)[0][0]
                #print 'index_i', index_i
		arr_corrs = matrix_corr[(index_date_targ*max_numb_imgs) + index_i, :]
                index_sort_corr = np.argsort(arr_corrs)
                #print 'index_sort_corr', index_sort_corr
                arr_indexdateoptimal = arr_indexdate_corrs[index_sort_corr]
                #print 'arr_indexdateoptimal', arr_indexdateoptimal
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
                        
                #print 'mean of arr_img2_optimal should be nearly zero:'
                #print 'np.mean(arr_img2_optimal):', np.mean(arr_img2_optimal)
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
        #arr_pool_vars = pickle.load(open(filename_arrpoolvars, 'rb')) #loads list of structures

        print '------------------------'
	print 'No of loops to run:', len(arr_pool_vars)
        print '------------------------'
        #pickle.dump(arr_pool_vars, open(filename_arrpoolvars, 'wb')) #saves list of structures
        



        #------
        # run KLIP
        #------
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

        
        
def klip(struct):
        #------
        #Extract variables from input structure
        #------
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
        print 'starting:', filenumb_img1, 'kmode', int(numb_klip_modes), ', radius', int(radius_mock), ', ratio', int(mock_factor), ',', datetime.datetime.now().strftime('%m-%d %H:%M:%S')


        
        #-------
        # Define filenames of fits files
        # Exit function if files already exist
        # CHANGE FILENAMES IF NECESSARY        
        #-------
        filename_output = directory + '/' + date + '/' + str(filenumb_img1) + '_klipfinal_' + 'kmode' + str(int(numb_klip_modes)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radannulus' + str(int(rad_annulus)) + '.fits'
        filename_sub = directory + '/' + date + '/' + str(filenumb_img1) + '_klipsubimg_' + 'kmode' + str(int(numb_klip_modes)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radannulus' + str(int(rad_annulus)) + '.fits'
        filename_imgmock = directory + '/' + date + '/' + str(filenumb_img1) + '_klipmockimg_' + 'kmode' + str(int(numb_klip_modes)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_radannulus' + str(int(rad_annulus)) + '.fits'
        if exists(filename_sub) and exists(filename_output) and exists(filename_imgmock):
                print 'skipping', filenumb_img1, 'kmode', int(numb_klip_modes), ', radius', int(radius_mock), ', ratio', int(mock_factor)
                return 0
        

        
        #---------------------
        #creating mock binary to be injected into img
        #---------------------
        img_mock = pf.open(directory + '/' + date + '/' + 'centroid_1.fits')[0].data
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


                
        #----------
        #flatten img and save img shape
        #----------
        img1shape = img1.shape
        img1_flat = img1.flatten()
        radi = 0


        
        #-----------------------------------        
        #save initial img with mock up binary added
        #-----------------------------------
        hdu = pf.PrimaryHDU(img1)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_imgmock, clobber=True) 



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

        print 'done with:', filenumb_img1, 'kmode', int(numb_klip_modes), ', radius', int(radius_mock), ', ratio', int(mock_factor), ',', datetime.datetime.now().strftime('%m-%d %H:%M:%S')
        return 1





def get_klip_basis(R, cutoff):
        #takes matrix R with N rows, M columns. N is number of imgs, M is number of pixels in img
        #cutoff is number of klip modes to use.
        #------

        w, V = np.linalg.eig(np.dot(R, np.transpose(R)))
        sort_ind = np.argsort(w)[::-1] #indices of eigenvals sorted in descending order
        sv = np.sqrt(w[sort_ind]).reshape(-1,1) #column of ranked singular values
        Z = np.dot(1./sv*np.transpose(V[:, sort_ind]), R)
        return Z[0:cutoff, :]


def create_cube_theklip():
        arr_klipmodes = np.arange(15, 35+1, 10)
        

        #---------------
        #ask for user input regarding which set to take as target psf
        #reads txt file in directory to know which files to use
        setnumb1 = raw_input('Enter 1st set number (1,2,3,4, etc.):')
        arr_startend1 = open(directory + '/' + date + '/set_'+ setnumb1 +'_info.txt', 'rb').read().splitlines()
        start1 = int(arr_startend1[0])
        end1 = int(arr_startend1[1])
	arr_targpsf = np.arange(start1, end1+1)



        #------
        # Median imgs for each numb of klip modes for the entire set
        #Now create cube by iterating through klip modes and adding to 'cube' array
        #------
        cube = []
        for numb_klip_modes in arr_klipmodes:
                arr_mdn = []
                for filenumb_img1 in arr_targpsf:
                        filename_output = directory + '/' + date + '/' + str(filenumb_img1) + '_klipfinal_' + 'kmode' + str(int(numb_klip_modes)) + '.fits'
                        if exists(filename_output):
                                img = pf.open(filename_output, memmap = False)[0].data
                                arr_mdn.append(img)
                        else:
                                print filename_output, 'doesnt exist'
                                continue
                arr_mdn = np.median(np.array(arr_mdn), axis = 0)                
                cube.append(arr_mdn)
        cube = np.array(cube)
        print 'cube.shape', cube.shape
        filename_cube = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) +  '_klipcube.fits'
        hdu = pf.PrimaryHDU(cube)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_cube, clobber=True)
        print 'CREATED:', filename_cube
        



def create_cube_klip():
        arr_klipmodes = np.arange(5, 45+1, 10)
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



        #------
        # Create cube by first iterating over each set of parameters
        # Then, median imgs for each numb of klip modes for the entire set
        #------
        for theta_mock in arr_theta:
                for radius_mock in arr_radiusmock:
                        for mock_factor in arr_mockfactor:
                                #------
                                #Now create cube by iterating through klip modes and adding to 'cube' array
                                #------
                                cube = []
                                for numb_klip_modes in arr_klipmodes:
                                        arr_mdn = []
                                        for filenumb_img1 in arr_targpsf:
                                                filename_output = directory + '/' + date + '/' + str(filenumb_img1) + '_klipfinal_' + 'kmode' + str(int(numb_klip_modes)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_theta' + str(int(theta_mock)) + '.fits'
                                                if exists(filename_output):
                                                        img = pf.open(filename_output, memmap = False)[0].data
                                                        arr_mdn.append(img)
                                                else:
                                                        print filename_output, 'doesnt exist'
                                                        continue
                                        arr_mdn = np.median(np.array(arr_mdn), axis = 0)                
                                        cube.append(arr_mdn)
                                cube = np.array(cube)
                                print 'cube.shape', cube.shape
                                filename_cube = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_mockrad' + str(radius_mock) + '_' +  'ratio' + str(mock_factor) + '_theta' + str(int(theta_mock)) + '_klipcube.fits'
                                hdu = pf.PrimaryHDU(cube)
                                hdulist = pf.HDUList([hdu])
                                hdulist.writeto(filename_cube, clobber=True)
                                print 'CREATED:', filename_cube
        


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
        filename_output = 'matrix_corr_filt.fits'
        counter = 0
        for index_chunk_run in np.arange(len(arr_date_ShaneAO)).astype(int):
                filename_matrix = 'matrix_corr_filt' + str(index_chunk_run) + '.fits'
                if not counter:
                        mat_corr = pf.open(filename_matrix)[0].data
                        print mat_corr.shape
                        print 'added', filename_matrix
                else:
                        mat_chunk = pf.open(filename_matrix)[0].data
                        print mat_chunk.shape
                        mat_corr = np.vstack((mat_corr, mat_chunk))
                        print 'added', filename_matrix
                counter+=1
        save_fits(mat_corr, filename_output)


def create_books_corr():
        struct_ShaneAO_filenumbs = {}
        for index_date in np.arange(len(arr_date_ShaneAO)): # loop through INDEX of dates
                date = arr_date_ShaneAO[index_date] #define date for index
                total_sets = struct_ShaneAO_total_sets[date] # define total number of sets for date
                arr_filenumbs = np.array([]) #array for adding img numbers
                for setnumb in np.arange(1, total_sets+1):
                        arr_startend = open(directory + '/' + date + '/set_'+ str(int(setnumb)) +'_info.txt', 'rb').read().splitlines()
                        start = int(arr_startend[0])
                        end = int(arr_startend[1])

                        arr_filenumbs = np.append(arr_filenumbs, np.arange(start, end+1))
                struct_ShaneAO_filenumbs[index_date] = arr_filenumbs

        # Save structure containing img numbers for all dates
        filename = 'struct_ShaneAO_filenumbs.p'
        pickle.dump(struct_ShaneAO_filenumbs, open(filename, 'wb'))        
        print 'created:', filename
        print '^a structure with img number for all dates'

        #------
        # Create and save array of filenumbs with dimensions = length of correlation matrix
        #------
        filename_filenumbs_corr = 'arr_filenumbs_corr' + '.fits'
        filename_indexdates_corr = 'arr_indexdates_corr' + '.fits'
        #arr_filenumbs_corr = np.zeros(matrix_corr.shape[0]*len(arr_date_ShaneAO)) 
        arr_filenumbs_corr = np.array([]) #array of img numbs corresponding to corr matrix
        arr_indexdates_corr = np.array([]) #array of index for dates corresponding to corr matrix

        for index_datefill in np.arange(len(struct_ShaneAO_filenumbs)).astype(int):
                arr_filenumbs_fill = struct_ShaneAO_filenumbs[index_datefill] #array of img numbs for certain date
                arr_indexdates_fill = np.ones(arr_filenumbs_fill.size).astype(int)*index_datefill #array filled with index for certain date
                
                #attach date indexes and img numbs for particular date
                arr_indexdates_corr = np.concatenate((arr_indexdates_corr, arr_indexdates_fill))
                arr_filenumbs_corr = np.concatenate((arr_filenumbs_corr, arr_filenumbs_fill))
        arr_indexdates_corr = arr_indexdates_corr.astype(int)
        #saving arrs for date indexes and img numbs
        save_fits(arr_filenumbs_corr, filename_filenumbs_corr)
        print 'created:', filename_filenumbs_corr
        save_fits(arr_indexdates_corr, filename_indexdates_corr)
        print 'created:', filename_indexdates_corr



def ret_correlation(date1, img1, date2, img2):
#Return correlation of RCF(reduced, centered, filtered) img1 and img2
#date1 and date2 correspond to dates of img1 and img2 respectively
#calculates correlation zeroing aperture around center
#then divides by sum of individual imgs
#then subtract imgs, returns std of subtracted imgs
#------

        img1 = pf.open('ShaneAO/' + date1 + '/s' + str(img1).zfill(4)+ '_RCF.fits')[0].data
        img2 = pf.open('ShaneAO/' + date2 + '/s' + str(img2).zfill(4)+ '_RCF.fits')[0].data

        imgshape = img1.shape

        #------
        # Get 1D indexes of aperture around center of img, with radius = max_rad_img
        # Later in function: will zero out primary before calculating correlation
        #------
	img_empty = np.zeros(imgshape)
	max_rad_ring = 4
        arr_ringindex = []
        for rad_ring in np.arange(0, max_rad_ring+1):
                #------
                #return: img values in aperture, array of 2d indexes for aperture
                #------
                ringvals, arr_ringindex_temp = return_ring(img_empty, rad_ring) 
                
                #------
                #loop through array of 2d indexes, change them to 1d indexes
                #------                
                for ringindex in arr_ringindex_temp:
                        #ringindex = np.ravel_multi_index(ringindex, imgshape)
                        arr_ringindex.append(ringindex) #array of 1d indexes for center aperture of radius max_rad_ring    
 
        for ringindex in arr_ringindex:
                img1[ringindex[0], ringindex[1]] = 0
                img2[ringindex[0], ringindex[1]] = 0

        img1/=np.sum(img1)
        img2/=np.sum(img2)
        return np.std(img1 - img2)

def correlation_matrix2(index_date_run, replace_cond = True): #index chunk is anything from  0 to (N^2)-1, where N is number of dates
        #Correlation matrix for FILTERED imgs
        #------

        numb_dates = len(arr_date_ShaneAO)
        fillelement = 10 #VARIABLE
        imgshape = [161,161]
        len_filt_box = 21
        filename_matrix = 'matrix_corr_filt' + str(int(index_date_run)) + '.fits'

        filename_filenumbs_corr = 'arr_filenumbs_corr' + '.fits'
        filename_indexdates_corr = 'arr_indexdates_corr' + '.fits'
        filename_struct_ShaneAO_filenumbs = 'struct_ShaneAO_filenumbs.p'
        struct_ShaneAO_filenumbs = pickle.load(open(filename_struct_ShaneAO_filenumbs, 'rb'))
        arr_indexdates_corr = pf.open(filename_indexdates_corr)[0].data
        arr_filenumbs_corr = pf.open(filename_filenumbs_corr)[0].data
        if arr_indexdates_corr.size == arr_filenumbs_corr.size:
                len_corr = arr_indexdates_corr.size #length of correlation matrix
        else:
                print '------'
                print 'ERROR!'
                print 'arr_indexdates_corr and arr_filenumbs_corr not the same size'
                print '------'
                return 0
 
        #------
        # Get 1D indexes of aperture around center of img, with radius = max_rad_img
        # Later in function: will zero out primary before calculating correlation
        #------
	img_empty = np.zeros(imgshape)
	max_rad_ring = 4
        arr_ringindex = []
        for rad_ring in np.arange(0, max_rad_ring+1):
                #------
                #return: img values in aperture, array of 2d indexes for aperture
                #------
                ringvals, arr_ringindex_temp = return_ring(img_empty, rad_ring) 
                
                #------
                #loop through array of 2d indexes, change them to 1d indexes
                #------                
                for ringindex in arr_ringindex_temp:
                        #ringindex = np.ravel_multi_index(ringindex, imgshape)
                        arr_ringindex.append(ringindex) #array of 1d indexes for center aperture of radius max_rad_ring


        arr_filenumbs_y = struct_ShaneAO_filenumbs[index_date_run]        
        
        #------
        #check if file already exists. If yes, load it and update as needed
        #------
        if replace_cond:
                matrix_corr = np.ones((arr_filenumbs_y.size,len_corr))*fillelement
        else:
                if exists(filename_matrix):
                        matrix_corr = pf.open(filename_matrix)[0].data
                else:
                        print filename_matrix, 'doesnt exist'
                        matrix_corr = np.ones((arr_filenumbs_y.size,len_corr))*fillelement

        print 'matrix_corr.shape', matrix_corr.shape
        counter = 0
        for index_filenumb_y in np.arange(arr_filenumbs_y.size).astype(int):
                #------
                # open img corresponding to index in y direction of correlation matrix
                #------
                filenumb_y = int(arr_filenumbs_y[index_filenumb_y])
                filename_y = directory + '/' + arr_date_ShaneAO[index_date_run] + '/' + 's' + str(filenumb_y).zfill(4) + '_RCF.fits'
                if exists(filename_y):
                        imgy = pf.open(filename_y)[0].data
                else:
                        print filename_y, 'does not exist'
                        continue
                        #####NOT DONE HERE

                print '------'                                
                print 'imgy', filename_y
                print (index_filenumb_y+1)*100./(arr_filenumbs_y.size), '% done'

                for index_filenumb_x in np.arange(len_corr).astype(int):

                        #------
                        #If correlation matrix has been filled for this, skip
                        #------
                        if matrix_corr[index_filenumb_y, index_filenumb_x] < fillelement:
                                print 'skipping: row already filled'
                                continue

                        #------
                        # open img corresponding to index in x direction of correlatin matrix
                        #------
                        date_x = arr_date_ShaneAO[int(arr_indexdates_corr[index_filenumb_x])]
                        filenumb_x = int(arr_filenumbs_corr[index_filenumb_x])
                        filename_x = directory + '/' + date_x + '/' + 's' + str(filenumb_x).zfill(4) + '_RCF.fits'
                        if exists(filename_x):
                                imgx = pf.open(filename_x)[0].data
                        else:
                                #print filename_x, 'does not exist'
                                continue

                        #------
                        # mask center of img with aperture
                        #------
                        for ringindex in arr_ringindex:
                                imgx[ringindex[0], ringindex[1]] = 0
                                imgy[ringindex[0], ringindex[1]] = 0
                                
                        imgy /= np.sum(imgy)
                        imgx /= np.sum(imgx)
                        img_res = imgy - imgx
                        sd_img_res = np.std(img_res)
                        ####Get rid of middle aperture pixels when taking std
                        matrix_corr[index_filenumb_y, index_filenumb_x] = sd_img_res
                        #print matrix_corr
                if counter%10 == 0:
                        save_fits(matrix_corr, filename_matrix)
                        print_timenow()
                counter += 0

'''
def test_corr_filt():
#For quick test that correlation matrix works
#prints residual value for given images in first 2 lines below

        #------
        # Change as needed***
        #------
        imgx = pf.open(directory + '/' + date + '/' + 's0146_RCF.fits')[0].data
        imgy = pf.open(directory + '/' + 'Jun15' + '/' + 's2715_RCF.fits')[0].data


        imgshape = [161,161]
        len_filt_box = 21

	img_empty = np.zeros(imgshape)
	max_rad_ring = 4
        arr_ringindex = []
        for rad_ring in np.arange(0, max_rad_ring+1):
                #------
                #return: img values in aperture, array of 2d indexes for aperture
                #------
                ringvals, arr_ringindex_temp = return_ring(img_empty, rad_ring) 

                
                #------
                #loop through array of 2d indexes, change them to 1d indexes
                #------                
                for ringindex in arr_ringindex_temp:
                        #ringindex = np.ravel_multi_index(ringindex, imgshape)
                        arr_ringindex.append(ringindex) #array of 1d indexes for center aperture of radius max_rad_ring


        for ringindex in arr_ringindex:
                imgx[ringindex[0], ringindex[1]] = 0
                imgy[ringindex[0], ringindex[1]] = 0

        imgy /= np.sum(imgy)
        imgx /= np.sum(imgx)
        img_res = imgy - imgx
        sd_img_res = np.std(img_res)
        print sd_img_res
'''



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
#returns structure with maximum value in ring, sd of ring, median of ring

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


def check_ring_mock(img, radi, theta_mock, theta_pm = 50): #inputs are: img array, radius of ring, angle where binary is positioned
        #check img for possible binaries. 
        #Check if any points in a ring are greater than 5sd above the median in a ring
        #exclude plus-minus angle theta_pm around theta_mock
        
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
                try:
                        arr_ring = np.append(arr_ring, img[element[0]][element[1]])
                except:
                        print element
                        print radi
                        print arr_index
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
		print '------'
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
	#if fwhm > 10:
	#	print 'FWHM', fwhm
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
#################'''It was a good effort! (Useless code)'''############################
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

