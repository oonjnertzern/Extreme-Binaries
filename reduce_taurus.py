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


#pathname to ds9 in taurus
ds9_taurus = '/Users/gduchene/Applications/Ureka/bin/ds9 '

#change the following when dealing with different date
#or different telescope
#or changes to img files e.g. new stars added
directory = 'ShaneAO'
date = 'Jun15'
arr_date_ShaneAO = ['Jun15', 'Nov15', 'Mar16', 'Sep16', 'Oct16', 'May17']
struct_ShaneAO_total_sets = {'Jun15':28, 'Nov15':25, 'Mar16':25, 'Sep16':19, 'Oct16':5, 'May17':18}
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

#length of radius out to which we keep the img
#our img sizes should extend out to 2*len_half_center+1 pixels for y and x
len_half_center = 80

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
#Change date arrays & dictionary of set numbers above at top of script
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
#run_create_centroid_filts() to create stacked imgs after high-pass filter is applied
#------------------
#to be revised:
#create_books_corr() to create some arrays for bookkeeping with correlation matrix
#correlation_matrix2 to create parts of correlation matrix
#merge_mat() to merge
#------------------

def print_stuff(): ####TEMP
        #prints filters used for all sets for all dates
        for date_targ in arr_date_ShaneAO:
                #print '------'
                #print date_targ
                #print '------'
                print '. '
                print '. '
                for setnumb1 in np.arange(1, struct_ShaneAO_total_sets[date_targ]+1).astype(int):
                        #print setnumb1
                        with open(directory + '/' + date_targ + '/set_'+ str(int(setnumb1)) +'_info.txt', 'rb') as f_temp:
                                arr_startend1 = f_temp.read().splitlines()
                                #print arr_startend1
                        start1 = int(arr_startend1[0])
                        end1 = int(arr_startend1[1])
                        starname = str(arr_startend1[2])
                        '''
                        try:
                                filtername = str(arr_startend1[3])
                                print filtername
                        except:
                                continue #print '.'
                        '''
                        numb_imgs = end1-start1+1
                        filtername = ret_filtername(start1, date_targ)
                        print filtername

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
	for i in np.arange(1, max_numb_imgs): #CHANGE ACCORDINGLY
        #for i in np.arange(119, 218): #CHANGE ACCORDINGLY

                #load img file
                filename = directory + '/' + date + '/' + 's' + str(i).zfill(4) + '.fits'
                if exists(filename):
                        fits = pf.open(filename)
                else:
                        continue
		hdr = fits[0].header

		if 'filt1nam' in hdr:
			filtername = hdr['filt1nam']
                        if filtername not in arr_filternames:
                                print filtername
                                arr_filternames.append(filtername)
		else:
                        print '------'
                        print "img numb", i, "does not have the tag 'FILT1NAM'"
			print hdr.keys()
			print 'number of header tags:', len(hdr.keys())
        print 'filters are:', arr_filternames


def sort_filters():
        arr_ks_flat = []
        arr_brg_flat = []
        arr_j_flat = []
        arr_h_flat = []
        arr_startend = open(directory + '/' + date + '/flat.txt', 'rb').read().splitlines()
        start = int(arr_startend[0])
        end = int(arr_startend[1])
        arr_flats = np.arange(start, end+1)
	for i in arr_flats:
                print '------'
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
                if filtername[0] == 'H':
                        arr_h_flat.append(i)
	arr_ks_flat = np.array(arr_ks_flat)
        print 'arr_ks_flat', arr_ks_flat
	arr_brg_flat = np.array(arr_brg_flat)
        print 'arr_brg_flat', arr_brg_flat
	arr_j_flat = np.array(arr_j_flat)
        print 'arr_j_flat', arr_j_flat
        arr_h_flat = np.array(arr_h_flat)
        print 'arr_h_flat', arr_h_flat
        if arr_ks_flat.size:
                pf_savefits(arr_ks_flat, directory + '/' + date + '/' + 'arr_ks_flat.fits')
                pickle.dump(arr_ks_flat, open(directory + '/' + date + '/' + 'arr_ks_flat.p', 'wb'))
                print "saved ks flats' img numbers as fits and pickle files"
        if arr_brg_flat.size:
                pf_savefits(arr_brg_flat, directory + '/' + date + '/' + 'arr_brg_flat.fits')
                pickle.dump(arr_brg_flat, open(directory + '/' + date + '/' + 'arr_brg_flat.p', 'wb'))                
                print "saved brg flats' img numbers"
        if arr_j_flat.size:
                pf_savefits(arr_j_flat, directory + '/' + date + '/' + 'arr_j_flat.fits')
                pickle.dump(arr_j_flat, open(directory + '/' + date + '/' + 'arr_j_flat.p', 'wb'))
                print "saved j flats' img numbers"
        if arr_h_flat.size:
                pf_savefits(arr_h_flat, directory + '/' + date + '/' + 'arr_h_flat.fits')
                pickle.dump(arr_h_flat, open(directory + '/' + date + '/' + 'arr_h_flat.p', 'wb'))
                print "saved h flats' img numbers"
                

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

def median_it(arr_filenumbs, str_before = directory + '/' + date + '/' + 's', str_after = '.fits'):
        arr_mdn = []
        print 'taking median of imgs:', arr_filenumbs
        for i in arr_filenumbs:
		#print i
                filename = str_before + str(i).zfill(4) + str_after
                if exists(filename):
                        #with pf.open(filename, memmap = False) as fits:
                        fits = pf.open(filename, memmap = False)
                        img = fits[0].data
                else:
                        print filename, 'doesnt exist'
                        continue
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

        filename_arrhflat_fits = directory + '/' + date + '/' + "arr_h_flat.fits"


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

        arr_h_flat = np.array([])
        if exists(filename_arrhflat_fits):
                arr_h_flat = pf.open(filename_arrhflat_fits)[0].data
        elif exists(filename_arrhflat_pickle):
                arr_h_flat = pickle.load(open(filename_arrhflat_pickle, "rb" ))
        if arr_h_flat.size:
                try:
                        img_h_flat = norm_median_subdark(arr_h_flat)
                        hdu = pf.PrimaryHDU(img_h_flat)
                        hdulist = pf.HDUList([hdu])
                        hdulist.writeto(directory + '/' + date + '/' + 'img_flat_h.fits',clobber=True)
                        print 'img_flat_h.fits created'                
                except:
                        print 'error with img_flat_h'


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
                                        #print filename_output, 'exists'
                                        continue
                                elif exists(filename_init):
                                        img_targ = pf.open(filename_init)[0].data
                                else:
                                        print img_numb, "has yet to be reduced & centered"
                                        continue
                                        
                                img_targ -= medfilt(img_targ, [len_filt_box, len_filt_box])
                                pf_savefits(img_targ, filename_output)
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


        filename_subed = directory + '/' + date_targ + '/' + 'set' + str(int(setnumb)) + '_locifiltfinal.fits'
        img_subed = pf_loadfits('')

        #------
        #Iterate through radii for mock binaries
        #Check LOCI-subtracted img for 5 sigma value at that radius
        #Create array of corresponding flux ratio
        #------
        arr_fluxratio = []
        for rad_mock in arr_radiusmock:
                struct_ring = check_ring(img_subed, rad_mock)
                sd_ring = struct_ring['sd_ring']
                arr_fluxratio.append(5*sd_ring*detlim_mult_factor)


        for setnumb1 in np.arange(1, struct_ShaneAO_total_sets[date]-6).astype(int): ###temporary
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
                                filename_old = directory + '/' + date + '/' + str(int(setnumb1)) + '_radmock' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + '_mockinit.fits'
                                filename_mockinit_output = directory + '/' + date + '/set' + str(int(setnumb1)) + '_radmock' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + '_mockinit.fits'

                                if exists(filename_old):
                                        subprocess.call('mv ' + filename_old + ' ' + filename_mockinit_output, shell = True)
                                        continue

                                arr_imgs = []
                                for imgnumb in arr_targpsf:
                                        filename_mockinit = directory + '/' + date + '/' + str(int(imgnumb)) + '_radmock' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + '_mockinit.fits'
                                        img = pf_loadfits(filename_mockinit, True)
                                        if img.size:
                                                arr_imgs.append(img)
                                img_output = np.median(np.array(arr_imgs), axis = 0)
                                pf_savefits(img_output, filename_mockinit_output)

def histo_improve_curves():
        #------
        # Load arrays of difference between 
        # initial detection limit(after high-pass filter)
        # and corrected detection limit. Plot histograms
        #------

        dir_histo_mock = 'Histo_Mock' #directory to save histograms to
        width_bin = 0.1 #width of bins in histograms

        list_filter = []
        list_diff = []
        list_filter_arrdiff = []
        cond_shape = False
        arr_radiusmock = np.array([8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 75]).astype(int)
        for date_targ in arr_date_ShaneAO:
                for setnumb in np.arange(1, struct_ShaneAO_total_sets[date_targ]+1).astype(int):
                        filename_diff_pre_final = directory + '/' + date_targ + '/' + 'arr_diff_pre_final_' + str(int(setnumb)) + '.fits'
                        if exists(filename_diff_pre_final):
                                arr_diff_temp = pf_loadfits(filename_diff_pre_final)
                                hdr_arrdiff = pf.open(filename_diff_pre_final)[0].header
                                if tag_filter in hdr_arrdiff:
                                        filtername_arrdiff = hdr_arrdiff[tag_filter]
                                        list_filter_arrdiff.append(filtername_arrdiff)
                                else:
                                        print '------'
                                        print '**ERROR**'
                                        print tag_filter, 'not in header'
                                        print '------'
                                        useless = raw_input('ERROR WITH HEADER. STOPPED.')
                                #print 'arr shape:', arr_diff_temp.shape
                                if cond_shape:
                                        for row_diff_temp in arr_diff_temp:
                                                if np.array_equal(shape_arrdiff, row_diff_temp.shape):
                                                        if np.isnan(row_diff_temp).any():
                                                                print 'we get nans'                   
                                                                print row_diff_temp
                                                                if sum(np.isnan(row_diff_temp)) > 1:
                                                                        continue
                                                                else:
                                                                        #------
                                                                        # load txt file and save filter used
                                                                        #------
                                                                        filename_info = directory + '/' + date_targ + '/' + 'set_' + str(setnumb) + '_info.txt'
                                                                        arr_startend1 = open(filename_info, 'rb').read().splitlines()
                                                                        imgnumb_start = int(arr_startend1[0])
                                                                        filter_setnumb = ret_filtername(imgnumb_start, date_targ = date_targ)
                                                                        list_filter.append(filter_setnumb)

                                                                        #append array of differences to list
                                                                        row_diff_temp = np.nan_to_num(row_diff_temp)
                                                                        list_diff.append(row_diff_temp)
                                                        else:
                                                                list_diff.append(row_diff_temp)
                                                else:
                                                        print 'array shapes not the same'
                                                        continue
                                else:
                                        shape_arrdiff = arr_diff_temp[0].shape
                                        print shape_arrdiff
                                        cond_shape = True

        print 'list_filter', list_filter
        print 'list_filter_arrdiff', list_filter_arrdiff
        arr_filter = np.array(list_filter)

        #------
        #plot histograms
        #------
        matrix_diff = np.array(list_diff)
        print 'matrix_diff', matrix_diff
        print 'matrix_diff.shape', matrix_diff.shape
        for index_rad in range(matrix_diff.shape[1]):
                arr_hist_temp = matrix_diff[:, index_rad]
                min_bin = math.floor(np.amin(arr_hist_temp)*10)/10.
                if min_bin < -10:
                        min_bin = 0
                loc_bin = np.arange(min_bin, np.amax(arr_hist_temp) + width_bin, width_bin)
                filename_histo_temp = directory + '/' + dir_histo_mock + '/' + 'rad' + str(int(arr_radiusmock[index_rad])) + '.png'
                if exists(filename_histo_temp):
                        continue
                plt.hist(arr_hist_temp, bins = loc_bin)
                plt.xlabel('Shift in magnitude of 5sd (After LOCI)')
                plt.title('5sd difference at radius = '+ str(int(arr_radiusmock[index_rad])))
                plt.savefig(filename_histo_temp, bbox_inches='tight')
                #plt.show()
                plt.close()
        



def compare_clearbin_loci(): #error tag #ccl# #CONTINUE here
        #
        #------
        radi_apert_com = 4
        radi_annul = 12
        radi_apert_flux = 4
        tag_filter = 'FILTER'
        hdr_tag_com_y_guess = 'BINEST-Y' #define header tags for estimated y and x center of mass index
        hdr_tag_com_x_guess = 'BINEST-X'
        hdr_tag_com_y_bef = 'BINCOM-Y' #define header tags for y and x center of mass index
        hdr_tag_com_x_bef = 'BINCOM-X'
        round_precision = 6 #how many decimal places to round com to        

        #------        
        #Load binary dates, set numbers
        #------
        arr_bin_dates, arr_bin_setnumbs = read_clearbin_info()

        for index_clearbin in range(arr_bin_dates.size):
                date_clearbin = arr_bin_dates[index_clearbin]
                setnumb_clearbin = arr_bin_setnumbs[index_clearbin]

                # Save filter used
                filter_setnumb = ret_filtername_set(setnumb_clearbin, date_targ = date_clearbin)

                # Define filenames
                filename_raw = directory + '/' + date_clearbin + '/' + 'centroid_' +  str(int(setnumb_clearbin)) + '.fits'
                filename_before = directory + '/' + date_clearbin + '/' + 'centroid_filt_' +  str(int(setnumb_clearbin)) + '.fits'
                filename_after = directory + '/' + date_clearbin + '/' + 'set' + str(int(setnumb_clearbin)) + '_locifiltfinal.fits'
                
                #------
                # Check that files exist, 
                # that estimated binary location is in header,
                # and that center of mass has been computed for that location.
                #------
                cond_file_exists =  exists(filename_before) and exists(filename_raw)
                hdr_bef = pf.open(filename_before)[0].header
                cond_bin_guess = hdr_tag_com_y_guess in hdr_bef and hdr_tag_com_x_guess in hdr_bef
                cond_bin_guess_com = hdr_tag_com_y_bef in hdr_bef and hdr_tag_com_x_bef in hdr_bef
                if cond_file_exists and cond_bin_guess and cond_bin_guess:
                        #print filename_before

                        y_guess = hdr_bef[hdr_tag_com_y_guess]
                        x_guess = hdr_bef[hdr_tag_com_x_guess]
                        com_y_bef = hdr_bef[hdr_tag_com_y_bef]
                        com_x_bef = hdr_bef[hdr_tag_com_x_bef]

                        print '------'
                        print filename_before
                        print '------'
                        img_bef = pf_loadfits(filename_before) #img after high pass filter but before loci
                        img_bef /= np.amax(img_bef)
                        img_raw = pf_loadfits(filename_raw) #raw stacked img
                        img_raw /= np.amax(img_raw)
                        img_aft = pf_loadfits(filename_after) #img after high pass filter and loci

                        com_y_aft, com_x_aft = get_com_aperture(img_aft, np.array([y_guess, x_guess]), radi_apert_com)
                        index_center_bef= ret_index_center(img_bef)
                        index_center_aft = ret_index_center(img_aft)
                        sep_y_bef, sep_x_bef = com_y_bef - index_center_bef[0], com_x_bef - index_center_bef[1]
                        sep_y_aft, sep_x_aft = com_y_aft - index_center_aft[0], com_x_aft - index_center_aft[1]

                        r_bef, angle_bef = cartesian2polar(sep_x_bef, sep_y_bef)
                        r_aft, angle_aft = cartesian2polar(sep_x_aft, sep_y_aft)

                        print '------'
                        print filter_setnumb
                        print 'binary separation(after LOCI):', r_aft

                        #------
                        # Measuring fluxes for primary,
                        # companions before LOCI(after high-pass filter) and after LOCI
                        #------
                        flux_filt_primary = ret_flux_center_ap(img_bef, radi_apert_flux)
                        flux_primary = ret_flux_center_ap(img_raw, radi_apert_flux)

                        flux_bef = ret_flux_com_aperture(img_bef, np.array([y_guess, x_guess]), radi_apert_com, radi_apert_flux)
                        flux_aft = ret_flux_com_aperture(img_aft, np.array([y_guess, x_guess]), radi_apert_com, radi_apert_flux)
                        mag_bef = 2.5*np.log10(flux_filt_primary/flux_bef) ####ccl#
                        mag_aft_uncorrect = 2.5*np.log10(flux_filt_primary/flux_aft)
                        print 'mag after high-pass filter:', mag_bef
                        print 'mag after LOCI (uncorrected):', mag_aft_uncorrect

                        #------
                        # load all correction curves for this date & filter
                        # Compute avg
                        #------
                        filename_arr_correction_filter = directory + '/' + date_clearbin + '/' + 'arr_correction_' + filter_setnumb[0] + '.fits'
                        if exists(filename_arr_correction_filter):
                                arr_correction_filter = pf_loadfits(filename_arr_correction_filter)
                                arr_correct_mean = np.mean(arr_correction_filter, axis = 0)
                        else:
                                print 'ERROR, file not found:', filename_arr_correction_filter
                                continue

                        #------
                        #load radius array for mock ups
                        #------
                        filename_arrradiusmock = directory + '/' + date_clearbin + '/' + 'arr_radiusmock.fits'
                        arr_radiusmock = pf_loadfits(filename_arrradiusmock)

                        #------
                        #interpolate correction to add to flux_aft
                        #------
                        for index_radiusmock in np.arange(arr_radiusmock.size-1).astype(int):
                                radmock = arr_radiusmock[index_radiusmock]
                                radmock_next = arr_radiusmock[index_radiusmock+1]
                                if radmock < r_aft < radmock_next:
                                        index_rad1 = index_radiusmock
                                        index_rad2 = index_radiusmock+1

                        rad1 = arr_radiusmock[index_rad1]
                        rad2 = arr_radiusmock[index_rad2]
                        frac_r = (r_aft - rad1)/(rad2 - rad1)

                        #------
                        #perform annulus subtraction
                        #------
                        flux_bin_raw_mdnsub = sub_annul(img_raw, [y_guess, x_guess], radi_apert_com, radi_apert_flux, radi_annul)
                        mag_bin_mdnsub = 2.5*np.log10(flux_primary/flux_bin_raw_mdnsub)
                        print 'median subtracted flux:', mag_bin_mdnsub

                        #------
                        # Interpolate correction curve at radius of companion
                        # Correct LOCI subtracted magnitude with interpolated value
                        # Compare with annulus-subtracted magnitude
                        #------
                        ls_mag_diff = []
                        ls_mag_aft_correct = []
                        for index_correct in range(len(arr_correction_filter)):
                                arr_correct = arr_correction_filter[index_correct]
                                mag_interpol_spl = scipy.interpolate.spline(arr_radiusmock, arr_correct, r_aft)
                                arr_radiusmock_fine = np.linspace(np.amin(arr_radiusmock), np.amax(arr_radiusmock), 100, endpoint = True)
                                plt.plot(arr_radiusmock_fine, scipy.interpolate.spline(arr_radiusmock, arr_correct, arr_radiusmock_fine))                             
                                mag_aft_correct = mag_aft_uncorrect - mag_interpol_spl
                                mag_diff = mag_aft_correct - mag_bef
                                #print 'mag_diff', mag_diff
                                ls_mag_diff.append(mag_diff)
                                ls_mag_aft_correct.append(mag_aft_correct)
                        mag_diff_mean = np.mean(ls_mag_diff)
                        mag_diff_sd = np.std(ls_mag_diff)
                                
                        mag_aft_correct_mean = np.mean(ls_mag_aft_correct)
                        print 'magnitude after correction', mag_aft_correct_mean
                        print '% difference in flux:', (mag_aft_correct_mean - mag_bin_mdnsub)*100*np.log(10)/2.5, '%'                     


def compare_brightmock_annul():
        # check bright mock up companions far way from primary
        # compare post-LOCI flux to raw img with annulus subtraction
        #-----
        radi_apert_com = 4
        radi_annul = 12
        radi_apert_flux = 4
        arr_theta = np.arange(0, 360, 60).astype(int)

        #------
        #set numbers of stars that had mock up binaries ran on them
        #-----
        struct_setnumbs = {}
        struct_setnumbs['Jun15'] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 18, 19, 20, 21, 25, 26]
        struct_setnumbs['Nov15'] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19]
        struct_setnumbs['Mar16'] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        struct_setnumbs['Sep16'] = range(1, 19+1)
        struct_setnumbs['Oct16'] = range(1, 5+1)
        struct_setnumbs['May17'] = range(1, 18+1)


        for date_targ in arr_date_ShaneAO:
                arr_setnumbs = struct_setnumbs[date_targ] #change if we only want to generate plots for certain sets

                #load radius array for mock ups
                filename_arrradiusmock = directory + '/' + date_targ + '/' + 'arr_radiusmock.fits'
                arr_radiusmock = pf_loadfits(filename_arrradiusmock)                

                for setnumb in arr_setnumbs:
                        filename_raw = directory + '/' + date_targ + '/' + 'centroid_' + str(int(setnumb)) + '.fits'
                        img_raw = pf_loadfits(filename_raw)
                        img_raw /= np.amax(img_raw)
                        

                        flux_primary = ret_flux_center_ap(img_raw, radi_apert_flux)

                        filename_filt = directory + '/' + date_targ + '/' + 'centroid_filt_' + str(int(setnumb)) + '.fits'
                        img_filt = pf_loadfits(filename_filt)
                        flux_filt_primary = ret_flux_center_ap(img_filt, radi_apert_flux)

                        filename_subed = directory + '/' + date_targ + '/' + 'set' + str(int(setnumb)) + '_locifiltfinal.fits'
                        img_subed = pf_loadfits(filename_subed)


                        for radius_large in arr_radiusmock[np.where(arr_radiusmock > 40)]:
                                for theta_mock in arr_theta:
                                        filename_mockinit = directory + '/' + date_targ + '/' +'set' + str(int(setnumb)) + '_radmock' + str(int(radius_large)) + '_theta' + str(int(theta_mock)) + '_mockinit.fits'

                                        if exists(filename_mockinit):
                                                img_mockinit = pf_loadfits(filename_mockinit)
                                                #pf_savefits(img_mockinit, 'test.fits')
                                                #open_fits('test.fits')
                                        else:
                                                continue
                                        print '------'

                                        print 'r', radius_large, ', theta', theta_mock
                                        shift_xy = cartesian2polar(radius_large, theta_mock, inv = True)
                                        shift_yx = np.array(shift_xy)[::-1]
                                        index_bin = ret_index_center(img_mockinit) + shift_yx
                                        #print 'index center of img', ret_index_center(img_mockinit)
                                        #print 'shift_yx', shift_yx
                                        print 'index_bin', index_bin

                                        flux_annulsub = sub_annul(img_mockinit, index_bin, radi_apert_com, radi_apert_flux, radi_annul)
                                        fluxratio_annulsub = flux_primary/flux_annulsub
                                        #print 'fluxratio_annulsub', fluxratio_annulsub
                                        print 'mag annulus subtraction', fluxratio2mag(fluxratio_annulsub)

                                        struct_ring = check_ring(img_subed, radius_large) #CONTINUE#
                                        sd_ring = struct_ring['sd_ring']
                                        fluxratio = 50*sd_ring
                                        #print 'fluxratio injected', fluxratio
                                        print 'mag injected', fluxratio2mag(1/fluxratio)
                                        
                                        '''
                                        arr_index = get_indexap_annulus(rad_mock, radi_apert)
                                        arr_flux_radi = []
                                        for elem in arr_index:
                                                y_apert_center = elem[0] + y_index_center_avg
                                                x_apert_center = elem[1] + x_index_center_avg
                                                array_flux_apert_radi = ret_flux_aperture(img_subed_filt, [y_apert_center, x_apert_center], radi_apert)
                                                arr_flux_radi.append(array_flux_apert_radi)
                                        #print 'arr_flux_radi', arr_flux_radi
                                        sd_ring_r = np.std(np.array(arr_flux_radi))
                                        print 'sd_ring_r', sd_ring_r
                                        arr_sd.append(sd_ring_r)
                                        arr_5sd.append(sd_ring_r*5)
                                        '''


def sub_annul(img, index_guess, radi_apert_com, radi_apert_flux, radi_annul):
        #Performs annulus median subtraction, returns flux after subtraction
        #------
        #Inputs:
        # 2-d img array, 2-element index of estimated binary position, 
        # radius of aperture for center of mass, radius of aperture for flux measurement,
        # radius of annulus for median calculation
        #------
        #Output:
        # flux of companion after annulus median subtraction
        #------

        #------
        # Do subtraction of median of annulus. Compare to high-pass filtered results
        #------
        com_y_bin, com_x_bin = get_com_aperture(img, index_guess, radi_apert_com) #companion center of mass
        shape_y_img, shape_x_img = img.shape
                        
        #------
        #fourier transform & shift to center binary companion at center of img
        #------
        img_centered = fourier_shift(img, [((shape_y_img-1)/2) - com_y_bin, ((shape_x_img-1)/2) - com_x_bin])
        shape_y_img_centered, shape_x_img_centered = img_centered.shape

        #------
        #Subtract flux of binary by median of annulus
        #------
        arr_ring, arr_index_ring = return_ring(img_centered, radi_annul)
        mdn_ring = np.median(arr_ring)
        img_centered_mdnsub = img_centered - mdn_ring
        flux_bin_mdnsub = ret_flux_center_ap(img_centered_mdnsub, radi_apert_flux)
        #mag_bin_mdnsub = 2.5*np.log10(flux_primary/flux_bin_check_mdnsub)
        return flux_bin_mdnsub


def sort_correction_filter():
        for index_date in range(len(arr_date_ShaneAO)):
                date_targ = arr_date_ShaneAO[index_date]
                struct_filter_setnumbs = {}

                filename_arr_correction = 'arr_correction_' + date_targ + '.fits'
                arr_correction = pf_loadfits(filename_arr_correction)

                filename_filter_correction = 'filter_correction_' + date_targ + '.p'
                ls_filter_correction = pickle.load(open(filename_filter_correction, "rb"))

                for index_correct in range(len(arr_correction)):
                        filter_correction = ls_filter_correction[index_correct]
                        if filter_correction[0] in struct_filter_setnumbs:
                                struct_filter_setnumbs[filter_correction[0]].append(arr_correction[index_correct, :])
                        else:
                                struct_filter_setnumbs[filter_correction[0]] = [arr_correction[index_correct, :]]
                

                print ''
                print struct_filter_setnumbs.keys()
                for filter_acronym in struct_filter_setnumbs:
                        filename_correction_filter = directory + '/' + date_targ + '/' + 'arr_correction_' + filter_acronym + '.fits'
                        arr_correction_filter = np.array(struct_filter_setnumbs[filter_acronym])

                        filename_sd = directory + '/' + date_targ + '/' + 'arr_correction_sd_' + filter_acronym + '.fits'
                        arr_correction_sd = np.std(arr_correction_filter, axis = 0)
                        if arr_correction_filter.shape[0]:
                                pf_savefits(arr_correction_filter, filename_correction_filter)
                                print '------'
                                print 'saved:', filename_correction_filter
                                print 'consisting of', arr_correction_filter.shape[0], 'curves'

                                pf_savefits(arr_correction_sd, filename_sd)

                                #'''
                                #------
                                #extra stuff #TEMP#
                                #------
                                dir_plot = 'Plot_Correction'

                                #load radius array for mock ups
                                filename_arrradiusmock = directory + '/' + date_targ + '/' + 'arr_radiusmock.fits'
                                arr_radiusmock = pf_loadfits(filename_arrradiusmock)

                                plt.plot(arr_radiusmock, np.mean(arr_correction_filter, axis = 0), 'o-')
                                plt.title('Correction curves for ' + date_targ + ', ' + filter_acronym)
                                plt.xlabel('Radius(pixels)')
                                plt.ylabel('Magnitude')
                                filename_plot = directory + '/' + dir_plot + '/' + 'plot_correction_' + date_targ + filter_acronym + '.png'
                                plt.savefig(filename_plot, bbox_inches='tight')
                                plt.close()

                                plt.plot(arr_radiusmock, arr_correction_sd, 'o-')
                                plt.title('S.d. for correction curves for ' + date_targ + ', ' + filter_acronym)
                                plt.xlabel('Radius(pixels)')
                                plt.ylabel('Magnitude')
                                filename_plot = directory + '/' + dir_plot + '/' + 'plot_correction_' + date_targ + '_sd' + filter_acronym + '.png'
                                plt.savefig(filename_plot, bbox_inches='tight')
                                #plt.show()
                                plt.close()
                                #'''



def add_clearbin_loc():
        # Goes through raw images of detected binary companions.
        # If companion is seen, enter approximate index of companion
        # Edit header of **high-pass filtered img** with this estimated position
        #-----

        hdr_tag_com_y_guess = 'BINEST-Y' #define header tags for estimated y and x center of mass index
        hdr_tag_com_x_guess = 'BINEST-X'
        hdr_tag_com_y_bef = 'BINCOM-Y' #define header tags for y and x center of mass index
        hdr_tag_com_x_bef = 'BINCOM-X'
        round_precision = 6

        #------        
        #Load binary dates, set numbers
        #------
        arr_bin_dates, arr_bin_setnumbs = read_clearbin_info()

        for index_clearbin in range(arr_bin_dates.size):
                date_clearbin = arr_bin_dates[index_clearbin]
                setnumb_clearbin = arr_bin_setnumbs[index_clearbin]

                # save filter used
                filter_setnumb = ret_filtername_set(setnumb_clearbin, date_targ = date_clearbin)

                filename_raw = directory + '/' + date_clearbin + '/' + 'centroid_' +  str(int(setnumb_clearbin)) + '.fits'
                filename_hp = directory + '/' + date_clearbin + '/' + 'centroid_filt_' +  str(int(setnumb_clearbin)) + '.fits'
                filename_loci = directory + '/' + date_clearbin + '/' + 'set' + str(int(setnumb_clearbin)) + '_locifiltfinal.fits'
                
                if exists(filename_hp) and exists(filename_raw):
                        #------
                        # Check if center of mass index of binary has already been estimated manually
                        # if yes, load estimated index from header
                        # if no, get user to check if binary is visible in stacked raw img
                        # if binary is visible, user should enter estimated binary's center of mass inde
                        #------
                        hdr_bef = pf.open(filename_hp)[0].header
                        if hdr_tag_com_y_guess in hdr_bef and hdr_tag_com_x_guess in hdr_bef:
                                y_guess = hdr_bef[hdr_tag_com_y_guess]
                                x_guess = hdr_bef[hdr_tag_com_x_guess]
                        else:
                                cond_check = raw_input('Type something then ENTER to view img. Press only ENTER to skip.')
                                if cond_check:
                                        open_fits(filename_raw)
                                #print 'If binary is seen, enter coordinates below. If not, press enter.'
                                yx_guess_str = raw_input('Enter y and x coordinates separated by comma:')
                                if ',' in yx_guess_str:
                                        yx_guess = yx_guess_str.split(',')
                                        try:
                                                y_guess, x_guess = np.array(yx_guess).astype(int)
                                        except:
                                                print 'error with extracting input coordinate'
                                                continue
                                        pf_edit_hdr(filename_hp, hdr_tag_com_y_guess, y_guess)
                                        pf_edit_hdr(filename_hp, hdr_tag_com_x_guess, x_guess)
                                else:
                                        continue

                                #------
                                # Check if binary's center of mass has already been computed
                                # Check is done by searching header
                                # If yes, load C.O.M. from header
                                # If no, calculate using centroid around estimated values
                                #------
                                if hdr_tag_com_y_bef in hdr_bef and hdr_tag_com_x_bef in hdr_bef:
                                        continue
                                else:
                                        com_y_bef, com_x_bef = get_com_aperture(img_bef, np.array([y_guess, x_guess]), radi_apert_com)
                                        pf_edit_hdr(filename_before, hdr_tag_com_y_bef, round(com_y_bef, round_precision))
                                        pf_edit_hdr(filename_before, hdr_tag_com_x_bef, round(com_x_bef, round_precision))


def read_clearbin_info():
        # Loads txt file with info about detected binaries
        # Outputs:
        #  array of dates, array of set numbers (both corresponding to one another)
        #------

        #------
        #open txt file with list of known binaries
        #------
        arr_bin = open(directory + '/' + directory + '_clearbinaries.txt', 'rb').read().splitlines()
        arr_bin_dates = np.array(arr_bin[0::2])
        arr_bin_setnumbs = np.array(arr_bin[1::2])
        
        return arr_bin_dates, arr_bin_setnumbs




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
                        '''
                        arr_y_com_init_thetas_std = []
                        arr_x_com_init_thetas_std = []
                        arr_y_com_filt_thetas_std = []
                        arr_x_com_filt_thetas_std = []
                        arr_y_com_final_thetas_std = []
                        arr_x_com_final_thetas_std = []
                        '''
                        for index_theta in np.arange(arr_theta.size).astype(int):
                                theta_mock = arr_theta[index_theta]

                                arr_y_com_filt_frames = np.array([])
                                arr_x_com_filt_frames = np.array([])

                                
                                filename_final_set = directory + '/' + date  + '/' + 'set' + str(int(setnumb1)) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfiltfinal.fits'
                                img_final_set = pf_loadfits(filename_final_set)
                                if not img_final_set.size:
                                        continue
                                filename_mockinit_set = directory + '/' + date + '/set' + str(int(setnumb1)) + '_radmock' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + '_mockinit.fits'                                                
                                img_init_set = pf_loadfits(filename_mockinit_set)
                                if not img_init_set.size:
                                        continue

                                y_length, x_length = img_final_set.shape
                                y_index_center = int((y_length - 1)/2)
                                x_index_center = int((x_length - 1)/2)

                                #--------
                                #Add mock binary at correct radius & angle, save img as fits file
                                #--------
                                dx_mock = radius_mock*np.cos(theta_mock*(2*np.pi)/360) #calculating y displacement from center
                                dx_mock = int(round(dx_mock))
                                dy_mock = radius_mock*np.sin(theta_mock*(2*np.pi)/360) #calculating x displacement from center
                                dy_mock = int(round(dy_mock))

                                # Get center of mass after LOCI
                                y_com_final, x_com_final = get_com_aperture(img_final_set, (dy_mock + y_index_center , dx_mock + x_index_center), radi_apert)                                
                                # Get center of mass before high-pass filter
                                y_com_init, x_com_init = get_com_aperture(img_init_set, (dy_mock + y_index_center , dx_mock + x_index_center), radi_apert)


                                '''
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

                                        if math.isnan(y_com_final) or math.isnan(x_com_final):
                                                print 'y_com_final', y_com_final
                                                print 'x_com_final', x_com_final

                                        if exists(filename_mockinit):
                                                img1 = pf.open(filename_mockinit)[0].data
                                        else:
                                                img1 += img_mock #Add binary to primary

                                        # Get center of mass before high pass filter
                                        y_com_init, x_com_init = get_com_aperture(img1, (dy_mock + y_index_center , dx_mock + x_index_center), radi_apert)
                                        if not exists(filename_mockinit):
                                                pf_savefits(img1, filename_mockinit)

                                        arr_y_com_init_frames = np.append(arr_y_com_init_frames, y_com_init)
                                        arr_x_com_init_frames = np.append(arr_x_com_init_frames, x_com_init)

                                        if math.isnan(y_com_init) or math.isnan(x_com_init):
                                                print 'y_com_init', y_com_init
                                                print 'x_com_init', x_com_init


                                        if exists(filename_mockfilt):
                                                img1 = pf.open(filename_mockfilt)[0].data
                                        else:                                
                                                img1 *= max_pix_primary  # Multiply by primary max pixel value from before
                                                img1 -= medfilt(img1, [len_filt_box, len_filt_box]) # Perform high-pass filter on img
                                                img1 /= np.amax(img1) # Normalize img using new max pixel value
                      
                                        # Get center of mass after high pass filter
                                        y_com_filt, x_com_filt = get_com_aperture(img1, (dy_mock + y_index_center , dx_mock + x_index_center), radi_apert)

                                        if math.isnan(y_com_filt) or math.isnan(x_com_filt):
                                                print 'y_com_filt', y_com_filt
                                                print 'x_com_filt', x_com_filt
                                                
                                        if not exists(filename_mockfilt):
                                                pf_savefits(img1, filename_mockfilt)

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

                                arr_y_com_filt_thetas.append(np.mean(arr_y_com_filt_frames))
                                arr_x_com_filt_thetas.append(np.mean(arr_x_com_filt_frames))
                                '''
                                arr_y_com_init_thetas.append(y_com_init)
                                arr_x_com_init_thetas.append(x_com_init)
                                
                                arr_y_com_final_thetas.append(y_com_final)
                                arr_x_com_final_thetas.append(x_com_final)
                               
                        
                        arr_y_com_init.append(arr_y_com_init_thetas)
                        arr_x_com_init.append(arr_x_com_init_thetas)

                        arr_y_com_final.append(arr_y_com_final_thetas)
                        arr_x_com_final.append(arr_x_com_final_thetas)
                        '''
                        arr_y_com_init_std.append(arr_y_com_init_thetas_std)
                        arr_x_com_init_std.append(arr_x_com_init_thetas_std)
                        arr_y_com_filt_std.append(arr_y_com_filt_thetas_std)
                        arr_x_com_filt_std.append(arr_x_com_filt_thetas_std)
                        arr_y_com_final_std.append(arr_y_com_final_thetas_std)
                        arr_x_com_final_std.append(arr_x_com_final_thetas_std)
                        '''
                '''
                print len(arr_y_com_init)
                print len(arr_x_com_init)
                print len(arr_y_com_filt)
                print len(arr_x_com_filt)
                '''
                
                arr_y_com_init = np.array(arr_y_com_init)
                print 'arr_y_com_init', arr_y_com_init
                arr_x_com_init = np.array(arr_x_com_init)
                print 'arr_x_com_init', arr_x_com_init
                arr_y_com_final = np.array(arr_y_com_final)
                print 'arr_y_com_final', arr_y_com_final
                arr_x_com_final = np.array(arr_x_com_final)
                print 'arr_x_com_final', arr_x_com_final
                '''
                arr_y_com_filt = np.array(arr_y_com_filt)
                arr_x_com_filt = np.array(arr_x_com_filt)

                arr_y_com_init_std = np.array(arr_y_com_init_std)
                arr_x_com_init_std = np.array(arr_x_com_init_std)
                arr_y_com_filt_std = np.array(arr_y_com_filt_std)
                arr_x_com_filt_std = np.array(arr_x_com_filt_std)
                arr_y_com_final_std = np.array(arr_y_com_final_std)
                arr_x_com_final_std = np.array(arr_x_com_final_std)
                '''
                try:
                        foldername = 'com/'
                        filename_y_init = foldername + 'arr_y_com_init_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        pf_savefits(arr_y_com_init, filename_y_init)
                        filename_x_init = foldername + 'arr_x_com_init_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        pf_savefits(arr_x_com_init, filename_x_init)
                        filename_y_final = foldername + 'arr_y_com_final_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        pf_savefits(arr_y_com_final, filename_y_final)
                        filename_x_final = foldername + 'arr_x_com_final_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                        pf_savefits(arr_x_com_final, filename_x_final) 
                except:
                        print 'something went wrong with saving files'

                '''

                filename_y_filt = foldername + 'arr_y_com_filt_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                pf_savefits(arr_y_com_filt, filename_y_filt)
                filename_x_filt = foldername + 'arr_x_com_filt_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                pf_savefits(arr_x_com_filt, filename_x_filt)
                filename_y_final = foldername + 'arr_y_com_final_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                pf_savefits(arr_y_com_final, filename_y_final)
                filename_x_final = foldername + 'arr_x_com_final_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                pf_savefits(arr_x_com_final, filename_x_final) 
                '''
                '''
                filename_y_init_std = foldername + 'arr_y_com_init_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                pf_savefits(arr_y_com_init_std, filename_y_init_std)
                filename_x_init_std = foldername + 'arr_x_com_init_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                pf_savefits(arr_x_com_init_std, filename_x_init_std)
                filename_y_filt_std = foldername + 'arr_y_com_filt_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                pf_savefits(arr_y_com_filt_std, filename_y_filt_std)
                filename_x_filt_std = foldername + 'arr_x_com_filt_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                pf_savefits(arr_x_com_filt_std, filename_x_filt_std)
                filename_y_final_std = foldername + 'arr_y_com_final_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                pf_savefits(arr_y_com_final_std, filename_y_final_std)
                filename_x_final_std = foldername + 'arr_x_com_final_std_' + '_' + 'setnumb' + str(int(setnumb1)) + '_' + date + '.fits'
                pf_savefits(arr_x_com_final_std, filename_x_final_std)
                '''

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
                if len(arr_r_final):
                        plt.hist(arr_r_final, bins = n_bins)
                        plt.xlabel('Shift in radius (After LOCI)')
                        plt.title('Radial shift for mocks at radius = '+str(int(arr_radiusmock[index_rad])))
                        #plt.show()
                        plt.savefig("/home/jjoon/rad"+str(int(arr_radiusmock[index_rad]))+"sets.png", dpi = 200)

                        plt.close()
                        plt.hist(arr_angle_final, bins = n_bins)
                        plt.xlabel('Shift in angle (After LOCI)')
                        plt.title('Azimuth shift for mocks at radius = '+str(int(arr_radiusmock[index_rad])))
                        plt.savefig("/home/jjoon/angle"+str(int(arr_radiusmock[index_rad]))+"sets.png", dpi = 200)
                        #plt.show()
                        
                        #useless = raw_input('enter key to continue')


def run_loci_date():
        for setnumb in np.arange(1, struct_ShaneAO_total_sets[date] +1 ).astype(int):
                run_loci_mockbins(setnumb, bin_cond = False, faint_cond = False, replace_cond = True)


def run_loci_mockbins(setnumb1, input_radiusmock = [0], bin_cond = True, faint_cond = False, replace_cond = True):
        # Run LOCI with mock binaries
        # Run with binaries if bin_cond = True
        # Run with replacement if replace_cond = True
        # input_radiusmock must be array with radii at which to place binaries
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
                if input_radiusmock[0]:
                        arr_radiusmock = np.array(input_radiusmock)
                else:
                        arr_radiusmock = np.array([8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 75]).astype(int)
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
                        #pf_savefits(arr_fluxratio, 'set' + str(int(setnumb1)) + '_fluxratio_bins_init.fits')
                else:
                        arr_fluxratio = []
                        for rad_mock in arr_radiusmock:
                                struct_ring = check_ring(img_subed_filt, rad_mock)
                                sd_ring = struct_ring['sd_ring']
                                arr_fluxratio.append(5*sd_ring*detlim_mult_factor)
                                print 'r:', rad_mock, ',', 'flux ratio of binary:', 5*sd_ring*detlim_mult_factor
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
        print 'adding reference img numbers from all dates...'
        arr_substar = []
        for index_date in np.arange(len(arr_date_ShaneAO)):
                date_ref = arr_date_ShaneAO[index_date]
                total_sets = struct_ShaneAO_total_sets[date_ref]
                print '------'
                print 'date:', date_ref, 'total_sets', total_sets
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
	array_flatname = ['img_flat_']*2

        array_flatname[0] += 'ks'
        array_flatname[1] += 'h'
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
                      


  
def center_imgs(setnumb): #error tag: #cim#
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

                output = center_star(img)

                if output.shape != (2*len_half_center + 1, 2*len_half_center +1): #cim#
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


def center_star(img, len_half_center = len_half_center): #error tag: #cst#
        #input is 2d numpy array
        #finds pixel with highest count, excludes outliers
        #centers using centroid
        #returns centered img of size k x k (k is defined here)
        #------

        index_max = np.argmax(img) ###find index of max value in img
        index_max = np.unravel_index(index_max,img.shape) ###change to 2d dimensional index
        max_val = img[index_max] #max pixel count
        counter = 0
        while img[index_max[0]+1,index_max[1]] < max_val/4. or img[index_max[0]-1,index_max[1]] < max_val/4. or img[index_max[0],index_max[1]-1] < max_val/4. or img[index_max[0],index_max[1]+1] < max_val/4.: ###going through outliers with high values and zeroing them if the surrounding points are much smaller
                img[index_max] = 0
                index_max = np.argmax(img)
                index_max = np.unravel_index(index_max, img.shape)
                print index_max
                counter += 1
                if counter > 20:
                        break
        if counter > 20:
                print '------'
                print 'too many outliers'
                print '------'
                return np.array([])
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
        print 'maximum index y, x:', y, ',', x
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
        com_x += x #center of mass x index in img
        com_y += y #center of mass y index in img
        dec_com_x = com_x%1
        dec_com_y = com_y%1
        int_com_x = int(com_x)
        int_com_y = int(com_y)

        #fourier transform & shift
        img = np.fft.fft2(img)
        img = fshift(img, [-dec_com_y, -dec_com_x])
        img = np.fft.ifft2(img)
        img = np.real(img)
        shape_y_img, shape_x_img = img.shape

        img_template = np.zeros([2*len_half_center+1, 2*len_half_center+1])
        y_shape_template, x_shape_template = img_template.shape
        index_y_center_template, index_x_center_template = len_half_center, len_half_center
        index_y_start_template_offset, index_x_start_template_offset = -len_half_center, -len_half_center
        index_y_end_template_offset, index_x_end_template_offset = len_half_center+1, len_half_center+1

        #for cropping img
        index_y_start = int_com_y-len_half_center
        index_y_end = int_com_y+len_half_center+1
        index_x_start = int_com_x-len_half_center
        index_x_end = int_com_x+len_half_center+1
        
        if index_y_start < 0:
                index_y_start = 0
                index_y_start_template_offset = -int_com_y

        if index_x_start < 0:
                index_x_start = 0
                index_x_start_template_offset = -int_com_x

        if index_y_end > shape_y_img:
                index_y_end_template_offset = shape_y_img - int_com_y 
        if index_x_end > shape_x_img:
                index_x_end_template_offset = shape_x_img - int_com_x 
        
        img_cropped = img[index_y_start:index_y_end,index_x_start:index_x_end]
        print '------'
        print 'img_cropped.shape'
        print img_cropped.shape
        print '------'
        print img_template[index_y_center_template + index_y_start_template_offset:index_y_center_template + index_y_end_template_offset, index_x_center_template + index_x_start_template_offset:index_x_center_template + index_x_end_template_offset].shape
        print '------'
        img_template[index_y_center_template + index_y_start_template_offset:index_y_center_template + index_y_end_template_offset, index_x_center_template + index_x_start_template_offset:index_x_center_template + index_x_end_template_offset] += img_cropped
        output = img_template
        return output
        
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

def print_calbin(): #TEMP#
        arr_dates = []
        for index_date in np.arange(len(arr_date_ShaneAO)).astype(int): #looping through each date
                date_targ = arr_date_ShaneAO[index_date]
                #print date_targ

                #------
                #open txt file with calibration binaries' img numbers
                #------
                filename_txt = directory + '/' + date_targ + '/' + 'calibration_binaries_' + date_targ + '.txt'
                if exists(filename_txt): 

                        #load txt file as list, with each row as one element of list
                        with open(filename_txt, 'rb') as f: 
                                arr_bin = f.read().splitlines()

                        #------
                        # create array of structures
                        # each structure contains calibration binary names,
                        # and img file numbers for each position
                        #------
                        arr_structbin = []
                        struct_temp = {}
                        counter_star = 1
                        for index_row in np.arange(len(arr_bin)).astype(int):
                                row = arr_bin[index_row]
                                if row[:4] == 'NAME':
                                        struct_temp['name'] = row[6:]
                                        #print row[6:]
                                elif not row:
                                        arr_structbin.append(struct_temp)
                                        struct_temp = {}
                                        counter_star = 1
                                else:
                                        try:
                                                arr_imgnumbs = [int(imgnumb) for imgnumb in arr_bin[index_row].split(' ')]
                                                struct_temp[counter_star] = arr_imgnumbs
                                                if counter_star == 1:
                                                        imgnumb1 = arr_imgnumbs[0]
                                                        filename_imgnumb1 = directory + '/' + date_targ + '/' + 's' + str(imgnumb1).zfill(4) + '.fits'
                                                        img_hdr = pf.open(filename_imgnumb1)[0].header
                                                        #for key in img_hdr:
                                                        #        print key
                                                        print img_hdr['DATE-BEG'][:10]
                                                        #print time_obs
                                                        #print img_hdr.keys()
                                                        
                                                        #useless = raw_input('enter key to continue')
                                                counter_star += 1
                                        except:
                                                print '------'
                                                print 'problem with row', int(index_row+1), 'in', filename_txt
                                                print '------'




def ret_separation_calbins(): #error tag: #rsc#
        # Prints separation in y,x and radius,separation angle
        # Do so for all calibration binaries

        radi_apert = 6
        radi_apert_mask = 15
        const_mask = 0

        counter_calbin = 1
        cond_calbin_numb = True
        while cond_calbin_numb:
                counter_calbin_pos = 1
                cond_calbin_posnumb = True
                arr_y_sep_avg = []
                arr_x_sep_avg = []
                arr_angle_sep_avg = []
                arr_rad_sep_avg = []
                while cond_calbin_posnumb:
                        filename_calbin_mdn_centered = filename_centered_calbin = directory + '/' + 'Calibration_Binaries' + '/' + 'calbin_' + str(counter_calbin) + '_RCpos' + str(counter_calbin_pos) + '.fits'
                        if exists(filename_calbin_mdn_centered):
                                '''
                                cond_ds9 = raw_input('Type something to open DS9, press ENTER to skip.')
                                if cond_ds9:
                                        subprocess.call('/home/apps/ds9/ds9 ' + filename_calbin_mdn_centered, shell = True)
                                '''
                                '''
                                str_y_guess = raw_input('Estimated y of binary center: ')
                                str_x_guess = raw_input('Estimated x of binary center: ')
                                y_guess, x_guess = int(str_y_guess), int(str_x_guess)
                                img = pf_loadfits(filename_calbin_mdn_centered)
                                com_y, com_x = get_com_aperture(img, (y_guess, x_guess), radi_apert)
                                index_center_y, index_center_x = int((img.shape[0]-1)/2), int((img.shape[1]-1)/2)
                                sep_y, sep_x = com_y -index_center_y, com_x - index_center_x
                                print 'separation in y,x:', sep_y, sep_x
                                '''
                                img = pf_loadfits(filename_calbin_mdn_centered)
                                index_center = int((img.shape[0]-1)/2), int((img.shape[1]-1)/2)
                                img = mask_aperture(img, index_center, radi_apert_mask, const_mask) #mask center star, leaving only companion
                                index_maxpixval_1d = np.argmax(img) #index of max pixel value after masking
                                index_maxpixval_2d = np.unravel_index(index_maxpixval_1d, img.shape)
                                #print '------'
                                #print 'max pixel index(y,x:)', index_maxpixval_2d
                                com_y, com_x = get_com_aperture(img, index_maxpixval_2d, radi_apert)
                                #print 'binary companion index(y,x):', com_y, ',', com_x
                                index_center_y, index_center_x = index_center
                                sep_y, sep_x = com_y - index_center_y, com_x - index_center_x
                                arr_y_sep_avg.append(sep_y)
                                arr_x_sep_avg.append(sep_x)
                                angle_sep = np.arctan2(sep_y, sep_x) #separation angle
                                rad_sep = np.sqrt((sep_y**2.) + (sep_x**2.)) #separation distance
                                arr_angle_sep_avg.append(angle_sep)
                                arr_rad_sep_avg.append(rad_sep)
                                #subprocess.call('/home/apps/ds9/ds9 ' + filename_calbin_mdn_centered, shell = True)
                                counter_calbin_pos += 1
                        else:
                                cond_calbin_posnumb = False
                                counter_calbin += 1
                #print np.mean(arr_x_sep_avg)
                print np.mean(arr_angle_sep_avg)*360./(2*math.pi) #rsc# #TEMP# Change depending on wanted output
                filename_pos1 = filename_centered_calbin = directory + '/' + 'Calibration_Binaries' + '/' + 'calbin_' + str(counter_calbin) + '_RCpos1' + '.fits'
                if not exists(filename_pos1):
                        print filename_pos1, 'doesnt exist' 
                        cond_calbin_numb = False


def process_calbins(): #error tag: #pcl#
        #performs sky subtraction, flat-fielding, centering of calibration binaries
        #imgs must be larger in order to contain wide binaries
        #define len_half_enter_calbin below

        len_half_center_calbin = 250 #radius away from primary to contain within img.

        arr_dates = []
        counter_calbin = 1
        for index_date in np.arange(len(arr_date_ShaneAO)).astype(int): #looping through each date
                date_targ = arr_date_ShaneAO[index_date]
                #------
                #open txt file with calibration binaries' img numbers
                #------
                filename_txt = directory + '/' + date_targ + '/' + 'calibration_binaries_' + date_targ + '.txt'
                if exists(filename_txt): 

                        #load txt file as list, with each row as one element of list
                        with open(filename_txt, 'rb') as f: 
                                arr_bin = f.read().splitlines()
                        
                        #------
                        # create array of structures
                        # each structure contains calibration binary names,
                        # and img file numbers for each position
                        # arr_structbin only containts info for one date at a time.
                        #------
                        arr_structbin = []
                        struct_temp = {}
                        counter_star = 1
                        for index_row in np.arange(len(arr_bin)).astype(int):
                                row = arr_bin[index_row]
                                if not row:
                                        print struct_temp['name']
                                        arr_structbin.append(struct_temp)
                                        struct_temp = {}
                                        counter_star = 1
                                else:
                                        try:
                                                arr_imgnumbs = [int(imgnumb) for imgnumb in arr_bin[index_row].split(' ')]
                                                struct_temp[counter_star] = arr_imgnumbs
                                                str_imgnumb1 = str(arr_imgnumbs[0])
                                                filename_img1 = directory + '/' + date_targ + '/' + 's' + str_imgnumb1.zfill(4) + '.fits'
                                                hdr_img1 = pf.open(filename_img1)[0].header
                                                struct_temp['name'] = hdr_img1['OBJECT']
                                                counter_star += 1
                                        except:
                                                continue
                                                #print '------'
                                                #print 'problem with row', int(index_row+1), 'in', filename_txt
                                                #print '------'
                        if len(struct_temp):
                                print struct_temp['name']
                                arr_structbin.append(struct_temp)
                                struct_temp = {}
                        #Don't delete
                        #------
                        # Create sky imgs for all positions
                        # There are typically 2 positions per binary pair
                        # Save sky imgs to structure with key 'sky#',
                        # '#' represents position number(typically only 1 & 2)
                        #------
                        for index_structbin in np.arange(len(arr_structbin)).astype(int):
                                struct_bin = arr_structbin[index_structbin]
                                counter_sky = 1
                                end_cond = False
                                while not end_cond:
                                        if counter_sky in struct_bin:
                                                arr_imgnumbs = struct_bin[counter_sky]
                                                str_before = directory + '/' + date_targ + '/' + 's'
                                                img_sky = median_it(arr_imgnumbs, str_before = str_before)
                                                arr_structbin[index_structbin]['sky'+str(int(counter_sky))] = img_sky
                                                counter_sky += 1
                                        else:
                                                end_cond = True
                        

                        for index_structbin in np.arange(len(arr_structbin)).astype(int):
                                struct_bin = arr_structbin[index_structbin]

                                # creating array with number of positions
                                # for N positions, arr_posnumb = [1, 2, ..., N]
                                arr_posnumb = []
                                counter_pos = 1
                                cond_pos = True
                                while cond_pos:
                                        if counter_pos in struct_bin:
                                                arr_posnumb.append(counter_pos)
                                                counter_pos+=1
                                        else:
                                                cond_pos = False

                                for index_countersky in np.arange(len(arr_posnumb)).astype(int):
                                        counter_sky = arr_posnumb[index_countersky]
                                        img_sky = struct_bin['sky' + str(int(arr_posnumb[index_countersky-1]))]
                                        img_sky = img_sky[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index] #crop sky img
                                        filename_flat = directory + '/' + date_targ + '/' + 'cleaned_img_flat_ks.fits'#pcl# #TEMP#, may need to edit filename of flat
                                        img_flat = pf_loadfits(filename_flat)

                                        arr_imgnumbs = struct_bin[counter_sky]
                                        for imgnumb in arr_imgnumbs:
                                                filename = directory + '/' + date_targ + '/' + 's' + str(imgnumb).zfill(4) + '.fits'
                                                if exists(filename):
                                                        img = pf.open(filename)[0].data
                                                else:
                                                        print filename, 'doesnt exist'
                                                        continue

                                                img = img[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index]
                                                #print 'np.amax(img_sky)', np.amax(img_sky)
                                                img2 = (img-img_sky)/(img_flat)
                                                filename_reduced = directory + '/' + date_targ + '/' + 's' + str(imgnumb).zfill(4) + '_reduced_.fits'
                                                pf_savefits(img2, filename_reduced)
                                                img_centered = center_star(img2, len_half_center_calbin)
                                                #print img_centered
                                                
                                                if img_centered.shape != (2*len_half_center_calbin + 1, 2*len_half_center_calbin +1): 
                                                        #print 'img shape aint right. CHECK'
                                                        #print 'cropped img shape:', img_centered.shape
                                                        filename_test = 'test.fits'
                                                        pf_savefits(img_centered, filename_test)
                                                        subprocess.call('/home/apps/ds9/ds9 ' + filename_reduced, shell = True)
                                                        #subprocess.call('/home/apps/ds9/ds9 ' + filename_test, shell = True)
                                                        useless = raw_input('enter key to continue')
                                                        continue
                                                
                                                filename_centered = directory + '/' + date_targ + '/' + 's' + str(imgnumb).zfill(4) + '_reduced_centered.fits'
                                                pf_savefits(img_centered, filename_centered)
                                        str_before = directory + '/' + date_targ + '/' + 's'
                                        str_after = '_reduced_centered.fits'
                                        filename_calbin_mdn_centered = directory + '/' + 'Calibration_Binaries' + '/' + 'calbin_' + str(counter_calbin) + '_RCpos' + str(counter_sky) + '.fits' # '_' + struct_bin['name'] + 
                                        img_mdn_centered = median_it(arr_imgnumbs, str_before, str_after)
                                        pf_savefits(img_mdn_centered, filename_calbin_mdn_centered)
                                        file_img = pf.open(filename_calbin_mdn_centered, mode = 'update')
                                        img_hdr = file_img[0].header
                                        img_hdr['OBJECT'] = struct_bin['name']
                                        file_img.flush()
                                        file_img.close()
                                        #hdulist.writeto(filename_calbin_mdn_centered)
                                counter_calbin += 1
                                        #subprocess.call('/home/apps/ds9/ds9 ' + filename_centered, shell = True)
        print counter_calbin-1


        '''
                        # #TEMP# Delete when done
                        for struct_bin in arr_structbin:
                                for key in struct_bin:
                                        print key
        '''

        '''
                        arr_bin_startnumb = filter(None, np.array(arr_bin[0::3]))
                        arr_bin_endnumb = filter(None, np.array(arr_bin[1::3]))
                        arr_bin_starname = filter(None, np.array(arr_bin[2::3]))
                        print len(arr_bin_startnumb), len(arr_bin_endnumb), len(arr_bin_starname)

                        for index_calbin in np.arange(len(arr_bin_startnumb)).astype(int):
                                start = int(arr_bin_startnumb[index_calbin])
                                end = int(arr_bin_endnumb[index_calbin])
                                arr_imgnumb = np.arange(start, end+1).astype(int)
                                starname = arr_bin_starname[index_calbin]
                                arr_sky = np.arange(start, end)
                                img_sky = median_it(arr_sky)
                                filename_output = directory + '/' + date_targ + '/' + 'img_sky_calbin' + str(index_calbin + 1) +'.fits'
                                pf_savefits(img_sky, filename_output)
                                print 'created sky img for', starname
                                img_sky = img_sky[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index] #crop sky img
                                filename_flat = directory + '/' + date_targ + '/' + 'cleaned_img_flat_ks.fits'#pcl# #TEMP#, may need to edit filename of flat
                                img_flat = pf_loadfits(filename_flat)
                                
                                for imgnumb in arr_imgnumb:
                                        print '------'
                                        print imgnumb
                                        filename = directory + '/' + date_targ + '/' + 's' + str(imgnumb).zfill(4) + '.fits'
                                        if exists(filename):
                                                img = pf.open(filename)[0].data
                                        else:
                                                print filename, 'doesnt exist'
                                                continue

                                        img = img[region_y_min_index:region_y_max_index,region_x_min_index:region_x_max_index]
                                        print 'np.amax(img_sky)', np.amax(img_sky)
                                        img2 = (img-img_sky)/(img_flat)
                                        filename_reduced = directory + '/' + date_targ + '/' + 's' + str(imgnumb).zfill(4) + '_reduced_.fits'
                                        pf_savefits(img2, filename_reduced)
                                        img_centered = center_star(img2, len_half_center_calbin)
                                        if img_centered.shape != (2*len_half_center_calbin + 1, 2*len_half_center_calbin +1): 
                                                print 'img shape aint right. CHECK'
                                                print img_centered.shape
                                                filename_test = 'test.fits'
                                                pf_savefits(img_centered, filename_test)
                                                subprocess.call('/home/apps/ds9/ds9 ' + filename_reduced, shell = True)
                                                #subprocess.call('/home/apps/ds9/ds9 ' + filename_test, shell = True)
                                                useless = raw_input('enter key to continue')
                                                continue
                                        filename_centered = directory + '/' + date_targ + '/' + 's' + str(imgnumb).zfill(4) + '_reduced_centered.fits'
                                        pf_savefits(img_centered, filename_centered)
                                        #subprocess.call('/home/apps/ds9/ds9 ' + filename_centered, shell = True)
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


def fluxratio2mag(fluxratio):
        # Convert flux ratio to apparent magnitude
        # fluxratio = (reference flux/observed flux)
        #------
        return 2.5*np.log10(fluxratio)


def ret_filtername_set(setnumb_targ, date_targ):
        #Returns filter used for **SET** of frames
        #------
        # Input: set number(int), date(str)
        # Output: filter name of file number(str)
        #------

        # load txt file and save filter used
        filename_info = directory + '/' + date_targ + '/' + 'set_' + str(setnumb_targ) + '_info.txt'
        arr_startend1 = open(filename_info, 'rb').read().splitlines()
        imgnumb_start = int(arr_startend1[0])
        filter_setnumb = ret_filtername(imgnumb_start, date_targ)
        return filter_setnumb

        
def ret_filtername(filenumb, date_targ = date):
        #Returns filter used for img frame
        #------
        #input: file number
        #output: filter of file number
        #------

        filename = directory + '/' + date_targ + '/' + 's' + str(int(filenumb)).zfill(4) + '.fits'
        img_hdr = pf.open(filename)[0].header
        filtername = img_hdr['filt1nam']
	#print filtername
        return filtername


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

        #------
        #Test parameters
        #------
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
        
        
        

        
def test_theta(theta_mock): ##TEMP#
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


def theloci(struct): #error tag: #tlo#
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
                filename_output = directory + '/' + date + '/' + str(i) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfiltfinal.fits' #!!!!!!***filt Change name if not using filt #faint*** #tlo#
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

                        filename_mockinit = directory + '/' + date + '/' + str(int(i)) + '_radmock' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + '_mockinit.fits'
                        pf_savefits(img1, filename_mockinit)
                        pf_edit_hdr(filename_mockinit, 'FLRATIO', mock_factor) #tlc#

        ##img1 *= max_pix_primary  # Multiply by primary max pixel value from before ##not necessary
        img1 -= medfilt(img1, [len_filt_box, len_filt_box]) # Perform high-pass filter on img
        img1 /= np.amax(img1) # Normalize img using new max pixel value #tlo#


        # save file after high-pass filter as fits fils
        if struct['bin_cond']:
                filename_mockfilt = directory + '/' + date + '/' + str(int(i)) + '_radmock' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + '_mockfilt.fits'
                pf_savefits(img1, filename_mockfilt)


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


def run_create_final_loci_mockbins(total_sets = total_sets, highpassfilt = True): #error tag: #rcflm#
        #create stacked loci mock up binary imgs
        #------

        startset = 1
        arr_setnumbs = np.arange(startset, total_sets+1).astype(int)
        arr_radiusmock = np.array([8, 75]).astype(int)#np.array([8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 75]).astype(int) 
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
                                        filename_output = directory + '/' + date + '/' + 'set'+ str(int(setnumb1)) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfiltfinal.fits' #faint***
                                else:
                                        filename_output = directory + '/' + date + '/' + 'set'+ str(int(setnumb1)) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfinal.fits'

                                for i in arr_targpsf:
                                        if highpassfilt:
                                                filename_input = directory + '/' + date + '/' + str(i) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfiltfinal.fits' #faint***
                                        else:
                                                filename_input = directory + '/' + date + '/' + str(i) + '_mockrad' + str(int(radius_mock)) + '_theta' + str(int(theta_mock)) + 'locimockfinal.fits'
                                        if exists(filename_input):
                                                try:
                                                        with pf.open(filename_input) as f_temp:
                                                                img = f_temp[0].data
                                                except:
                                                        print 'error OPENING img file, but it exists:', filename_input
                                                        continue
                                        else:
                                                #print filename_input, 'doesnt exist', 'from set:', setnumb1
                                                continue
                                        arr_for_mdn.append(img)
                                print 'images for mdn:',len(arr_for_mdn)

                                if len(arr_for_mdn) > 5:
                                        img_final = np.median(np.array(arr_for_mdn), axis = 0)
                                        pf_savefits(img_final, filename_output)
                                        print 'created:', filename_output

def plot_detlim_final(filt = True):
        #Plot final detection limits for ALL stars.
        #User correction curves that were saved for specific run and filter
        #------
        global date
        global total_sets
        radi_apert = 4

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
                        filtername = ret_filtername_set(setnumb1, date)

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

                        flux_main_filt = ret_flux_aperture(img_avg, [y_index_center_avg, x_index_center_avg], radi_apert)
                        flux_main_unfilt = ret_flux_aperture(img_avg_unfilt, [y_index_center_avg, x_index_center_avg], radi_apert)


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
                                        array_flux_apert_radi_pre = ret_flux_aperture(img_avg_unfilt, [y_apert_center_pre, x_apert_center_pre], radi_apert)
                                        arr_5sd_radi_pre.append(array_flux_apert_radi_pre)                        
                                sd_ring_pre = np.std(np.array(arr_5sd_radi_pre))
                                arr_5sd_pre.append(sd_ring_pre*5)
                        arr_mag_5sd_pre = 2.5*np.log10(flux_main_unfilt/np.array(arr_5sd_pre))
                        pf_savefits(arr_mag_5sd_pre, filename_detlim_init)
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
                                        array_flux_apert_radi_pre_filt = ret_flux_aperture(img_avg, [y_apert_center_pre_filt, x_apert_center_pre_filt], radi_apert)
                                        arr_5sd_radi_pre_filt.append(array_flux_apert_radi_pre_filt)                        
                                sd_ring_pre_filt = np.std(np.array(arr_5sd_radi_pre_filt))
                                arr_5sd_pre_filt.append(sd_ring_pre_filt*5)
                        arr_mag_5sd_pre_filt = 2.5*np.log10(flux_main_unfilt/np.array(arr_5sd_pre_filt))
                        pf_savefits(arr_mag_5sd_pre_filt, filename_detlim_filt_init)
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
                                                array_flux_apert_radi = ret_flux_aperture(img_subed_filt, [y_apert_center, x_apert_center], radi_apert)
                                        else:
                                                array_flux_apert_radi = ret_flux_aperture(img_subed, [y_apert_center, x_apert_center], radi_apert)
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
                        pf_savefits(arr_mag_5sd - arr_correction, filename_detlim_final)
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
                        if ret_filtername(index_row+1, date_targ = arr_date_ShaneAO[index_date])[0] == 'K':
                                if counter_ks < 1: 
                                        arr_plots += plt.plot(arr_radius, arr_correction[index_row], arr_colors[int(index_date*2)] + 'o-')
                                        arr_labels.append(arr_date_ShaneAO[index_date] + ', ' + 'Ks')
                                        dict_correction[arr_date_ShaneAO[index_date] + '_Ks'] = []
                                        dict_correction[arr_date_ShaneAO[index_date] + '_Ks'].append(arr_correction[index_row])
                                else:
                                        plt.plot(arr_radius, arr_correction[index_row], arr_colors[int(index_date*2)] + 'o-')
                                        dict_correction[arr_date_ShaneAO[index_date] + '_Ks'].append(arr_correction[index_row])
                                counter_ks +=1 
                        elif ret_filtername(index_row+1, date_targ = arr_date_ShaneAO[index_date])[0] == 'B':
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
                pf_savefits(arr_correction_final, filename_arr_correction_final)
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


def ret_flux_com_aperture(img, index, radi_com, radi_ap):
        #------
        # Return flux in aperture around center of mass
        #Inputs:
        # 2d img array, 2d index, 
        # radius of center of mass aperture(int), radius of aperture for flux calculation(int)
        #Output:
        # total flux in aperture
        #------
        index_com = get_com_aperture(img, index, radi_com)
        com_y, com_x = index_com
        dec_com_x = com_x%1
        dec_com_y = com_y%1
        int_com_x = int(com_x)
        int_com_y = int(com_y)
        img = fourier_shift(img, [-dec_com_y, -dec_com_x])
        flux_ap = ret_flux_aperture(img, [int_com_y, int_com_x], radi_ap, ret_array=False)
        return flux_ap


def run_detlim_loci():
        #Run plot_detlim_loci for all dates
        #------
        
        for date in arr_date_ShaneAO:
                plot_detlim_loci(date)


def plot_detlim_loci(date=date): #error tag: #pdl#
        filt = True
        startset = 1 #pdl#
        tag_filter = "FILTER"


        #------
        #Define radius of apertures
        #For center of mass calculations and flux measurements
        #------
        radi_apert = 4
        radi_com = 4
        

        #------
        #set numbers of stars that had mock up binaries ran on them
        #-----
        struct_setnumbs = {}
        struct_setnumbs['Jun15'] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 18, 19, 20, 21, 25, 26]
        struct_setnumbs['Nov15'] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19]
        struct_setnumbs['Mar16'] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        struct_setnumbs['Sep16'] = range(1, 19+1)
        struct_setnumbs['Oct16'] = range(1, 5+1)
        struct_setnumbs['May17'] = range(1, 18+1)
        arr_setnumbs = struct_setnumbs[date] #change if we only want to generate plots for certain sets
        

        #------
        # define mock binary radii
        # save to fits file
        #------
        arr_radiusmock = np.array([8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 75]).astype(int)
        filename_arrradiusmock = directory + '/' + date + '/' + 'arr_radiusmock.fits'
        pf_savefits(arr_radiusmock, filename_arrradiusmock)
        print 'arr_radiusmock', arr_radiusmock

        arr_theta = np.arange(0, 360, 60).astype(int)
        print 'arr_theta', arr_theta


        #------
        # Create or open stacked correction curves & uncertainties of these curves
        #-----
        filename_filtername_correction_date = 'filter_correction_' + date + '.p'
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

        arr_final_err_save = np.zeros(shape = (arr_radiusmock.size))
        arr_5sd_correct = []

        # Name of filename for 5sd values after LOCI correction
        if filt:
                filename_arr5sdcorrectfits = directory + '/' + date + '/' + 'arr_5sd_correct_filt.fits'
        else:
                filename_arr5sdcorrectfits = directory + '/' + date + '/' + 'arr_5sd_correct.fits'

        
        #loop through sets of stars for this date.
        ls_filtername_set = []
        for setnumb1 in arr_setnumbs:
                print '------'
                print 'set', setnumb1
                print '------'
                
                filtername_set = ret_filtername_set(setnumb1, date)
                
                #------
                # Open LOCI-subtracted, medianed(and perhaps also filtered) img for set number
                # Define img name for plots
                #------
                if filt:
                        filename_plot = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_plotfilt.png'
                        filename_fits_filt = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_locifiltfinal.fits'
                        img_subed_filt = pf_loadfits(filename_fits_filt)
                else:
                        filename_plot = directory + '/' + date + '/' + 'set' + str(int(setnumb1)) + '_plot.png'
                
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
                        continue
                img_avg = pf.open(filename_avg)[0].data
                img_avg /= np.amax(img_avg)
                y_length_avg, x_length_avg = img_avg.shape
                y_index_center_avg = int((y_length_avg - 1)/2)
                x_index_center_avg = int((x_length_avg - 1)/2)
                flux_main = ret_flux_aperture(img_avg, [y_index_center_avg, x_index_center_avg], radi_apert)

                print 'flux_main', flux_main


                #------
                # Get flux of stacked, medianed primary without high pass filter
                #------
                filename_avg_unfilt = directory + '/' + date + '/' + 'centroid_' + str(int(setnumb1)) + '.fits'
                if not exists(filename_avg_unfilt):
                        print "filename_avg_unfilt doesn't exist"
                        continue
                img_avg_unfilt = pf_loadfits(filename_avg_unfilt)
                img_avg_unfilt /= np.amax(img_avg_unfilt)
                y_length_avg_unfilt, x_length_avg_unfilt = img_avg_unfilt.shape
                y_index_center_avg_unfilt = int((y_length_avg_unfilt - 1)/2)
                x_index_center_avg_unfilt = int((x_length_avg_unfilt - 1)/2)
                flux_main_unfilt = ret_flux_aperture(img_avg_unfilt, [y_index_center_avg_unfilt, x_index_center_avg_unfilt], radi_apert)


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
                                array_flux_apert_radi_pre = ret_flux_aperture(img_avg, [y_apert_center_pre, x_apert_center_pre], radi_apert)
                                arr_5sd_radi_pre.append(array_flux_apert_radi_pre)
                        
                        sd_ring_pre = np.std(np.array(arr_5sd_radi_pre))
                        arr_5sd_pre.append(sd_ring_pre*5)
                arr_mag_5sd_pre = 2.5*np.log10(flux_main/np.array(arr_5sd_pre))


                #------
                # Run get_flux_aperture to find number of pixels in aperture
                #------
                img_zeros = np.zeros(img_subed_filt.shape)
                arr_empty_apert = ret_flux_aperture(img_zeros, [y_index_center_avg, x_index_center_avg], radi_apert, True)
                num_pix_apert = len(arr_empty_apert)


                #------
                # Iterate through radii for mock binaries
                # Check medianed, LOCI-subtracted img(without binary) for 5 sigma value at that radius
                # Create array of corresponding flux ratio
                #------
                arr_5sd = []
                arr_sd = []
                arr_flux_bin_init = []
                for rad_mock in arr_radiusmock:
                        '''
                        arr_index = get_indexap_annulus(rad_mock, radi_apert)
                        arr_flux_radi = []
                        for elem in arr_index:
                                y_apert_center = elem[0] + y_index_center_avg
                                x_apert_center = elem[1] + x_index_center_avg
                                array_flux_apert_radi = ret_flux_aperture(img_subed_filt, [y_apert_center, x_apert_center], radi_apert)
                                arr_flux_radi.append(array_flux_apert_radi)
                        #print 'arr_flux_radi', arr_flux_radi
                        sd_ring_r = np.std(np.array(arr_flux_radi))
                        print 'sd_ring_r', sd_ring_r
                        '''
                        sd_ring_r = ret_sd_fluxap(img, rad_mock, radi_apert_flux)
                        arr_sd.append(sd_ring_r)
                        arr_5sd.append(sd_ring_r*5)

                        struct_ring = check_ring(img_subed_filt, rad_mock)
                        sd_ring_formocks = struct_ring['sd_ring']
                        mock_factor = sd_ring_formocks*50
                        img_mock = pf.open(directory + '/' + 'Jun15' + '/' + 'centroid_1.fits')[0].data
                        img_mock /= np.amax(img_mock)
                        img_mock *= float(mock_factor)

                        #filename_mockinit = directory + '/' + date_targ + '/' + str(int(setnumb)) + '_radmock' + str(int(radius_large)) + '_theta' + str(int(theta_mock)) + '_mockinit.fits'
                        #pf_savefits(img_mock)

                        flux_bin_init = ret_flux_center_ap(img_mock, radi_apert_flux)
                        arr_flux_bin_init.append(flux_bin_init)
                        
                arr_flux_bin_init = np.array(arr_flux_bin_init)
                arr_fluxratio_init = flux_main_unfilt/arr_flux_bin_init
                arr_mag_5sd = 2.5*np.log10(flux_main/np.array(arr_5sd)) #uncorrected 5sd values after loci
                arr_5sd_err = (flux_main/np.array(arr_sd))/(flux_main/np.array(arr_5sd)) ###Fix this!!!

                counter = 0
                arr_avg_fluxbin = []
                arr_sd = []
                arr_mag_alltheta_final_trans = []
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
                                        #print 'File not found:', filename_final
                                        continue

                                #------
                                # find x and y index at center of image
                                #------
                                y_length, x_length = img_final.shape
                                y_index_center = int((y_length - 1)/2)
                                x_index_center = int((x_length - 1)/2)

                                #------
                                # Use radius_mock and theta_mock to get index of mock binary
                                #------
                                dx_mock = radius_mock*np.cos(theta_mock*(2*np.pi)/360)
                                dx_mock = int(round(dx_mock))
                                dy_mock = radius_mock*np.sin(theta_mock*(2*np.pi)/360)
                                dy_mock = int(round(dy_mock))

                                y_index_bin = y_index_center + dy_mock
                                x_index_bin = x_index_center + dx_mock
                                
                                #------
                                # obtain flux of mock binary
                                #------
                                flux_bin = ret_flux_com_aperture(img_final, [y_index_bin, x_index_bin], radi_com, radi_apert)
                                arr_fluxbin_r.append(flux_bin)
                                #arr_fluxbin_r is flux after

                        arr_mag_alltheta_final_trans.append(2.5*np.log10(flux_main/np.array(arr_fluxbin_r))) #magnitudes for binaries for each position angle after loci
                        arr_avg_fluxbin.append(np.mean(arr_fluxbin_r))
                        arr_fluxbin_r = flux_main/arr_fluxbin_r
                        sd_r = np.std(arr_fluxbin_r)
                        arr_sd.append(sd_r)


                #print 'arr_sd', arr_sd
                arr_mag_sd = 2.5*np.log10(np.sqrt( (np.array(arr_sd)**2.) + (np.array(arr_5sd_err)**2.) ))
                #print 'arr_5sd_err', arr_5sd_err
                arr_mag_final_err = 2.5*np.log10(np.sqrt( (np.array(arr_sd)**-2.) + (np.array(arr_5sd_err)**2.) )) # exponent for arr_sd is -2, since we need to invert flux ratio

                arr_fluxratio_final = flux_main_unfilt/arr_avg_fluxbin #pdl#
                arr_mag_init = 2.5*np.log10(arr_fluxratio_init)
                arr_mag_final = 2.5*np.log10(arr_fluxratio_final)
                arr_mag_diff = arr_mag_final - arr_mag_init #correction curve, diff btwn binary flux before injection and after loci
                #print 'arr_mag_final_err', arr_mag_final_err
                #print 'arr_mag_diff', arr_mag_diff
                #print 'arr_mag_sd', arr_sd
                #print '-------'
                
                arr_mag_alltheta_final = np.transpose(np.array(arr_mag_alltheta_final_trans))
                print 'arr_mag_alltheta_final.shape', arr_mag_alltheta_final.shape

                arr_mag_diff_pre_final = []
                if counter > 50:
                        for index_row in range(len(arr_mag_alltheta_final)):
                                arr_mag_theta_temp = arr_mag_alltheta_final[index_row]
                                arr_mag_diff_pre_final_temp = (arr_mag_5sd - (arr_mag_theta_temp - arr_mag_init)) - arr_mag_5sd_pre
                                arr_mag_diff_pre_final.append(arr_mag_diff_pre_final_temp)
                        filename_diff_pre_final = directory + '/' + date + '/' + 'arr_diff_pre_final_' + str(int(setnumb1)) +  '.fits'
                        arr_mag_diff_pre_final = np.array(arr_mag_diff_pre_final)
                        print 'arr_mag_diff_pre_final', arr_mag_diff_pre_final
                        print 'arr_mag_diff_pre_final.shape', arr_mag_diff_pre_final.shape
                        pf_savefits(arr_mag_diff_pre_final, filename_diff_pre_final)
                        pf_edit_hdr(filename_diff_pre_final, tag_filter, filtername_set )
                        print '------'
                        print 'Created:', filename_diff_pre_final
                        print filename_diff_pre_final
                        print '------'

                
                #Do not delete
                #------
                # Plot curves
                #------
                plt.close()
                mag_pre, = plt.plot(arr_radiusmock, arr_mag_5sd_pre, 'ko-') #5 s.d. before subtraction
                mag_init, = plt.plot(arr_radiusmock, arr_mag_init, 'bo-') #Mag. of mock binary before injection
                mag_afterloci, = plt.plot(arr_radiusmock, arr_mag_final, 'ro-') #Mag. of binary after LOCI
                mag_5sd, = plt.plot(arr_radiusmock, arr_mag_5sd,'yo-') #5 s.d. before correction
                mag_5sd_correct, = plt.plot(arr_radiusmock, arr_mag_5sd - arr_mag_diff,'go-') #5 s.d. after correction
                print 'mag_5sd_correct', arr_mag_5sd - arr_mag_diff
                #useless = raw_input('stopped...')
                #plt.errorbar(arr_radiusmock, arr_mag_5sd - arr_mag_diff, yerr = arr_mag_final_err, ecolor = 'g', color ='g') ### FIX ######
                plt.gca().invert_yaxis()
                plt.legend([mag_pre, mag_init, mag_afterloci, mag_5sd, mag_5sd_correct], ['5 s.d. before subtraction', 'Mag. of mock binary before injection', 'Mag. of binary after LOCI', '5 s.d. before correction', '5 s.d. after correction'])
               
                #plt.show() # #pdl#
                #plt.close()
                #------
                # Append corrected 5sd plots to rows in a table. Save as fits file.
                #-----
                if counter > 50:
                        arr_5sd_correct.append(arr_mag_5sd - arr_mag_diff)
                        arr_correction_save = np.append(arr_correction_save, np.array([arr_mag_diff]), axis = 0)
                        ls_filtername_set.append(filtername_set)
                        #arr_correction_err_save = np.append(arr_correction_err_save, np.array([arr_mag_sd]), axis = 0)
                pf_savefits(np.array(arr_correction_save), filename_arr_correction)
                pickle.dump(ls_filtername_set, open(filename_filtername_correction_date, 'wb')) #saves list filters
                print 'number of images used:', counter
                print 'created:', filename_arr_correction

                mng = plt.get_current_fig_manager()                                         
                mng.resize(*mng.window.maxsize())
                if exists(filename_plot):
                        subprocess.call('rm -rf ' + filename_plot, shell = True)  #delete graph if exists
                if not np.isnan(arr_mag_5sd - arr_mag_diff).any():
                        plt.savefig(filename_plot)#, bbox_inches='tight')
                        print 'saved plot as img:', filename_plot
                        #plt.show()
                

                plt.close()
                
                
                hdu = pf.PrimaryHDU(np.array(arr_correction_err_save))
                hdulist = pf.HDUList([hdu])
                hdulist.writeto(filename_arr_correction_err, clobber=True)
                print 'created:', filename_arr_correction_err
                #'''

        ''' #Do not delete
        hdu = pf.PrimaryHDU(np.array(arr_5sd_correct))
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename_arr5sdcorrectfits, clobber=True)
        print 'created:', filename_arr5sdcorrectfits        
        '''

def test_apertures_ring(): ###TEMP#
#for radius from center of img, highlight apertures for obtaining 5sigma values
        radi_apert = 3
        arr_radiusmock = np.arange(2, 10, 1)##np.array([10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80]).astype(int)
        y_index_center_avg, x_index_center_avg = 80, 80
        for rad_mock in arr_radiusmock:
                img_test = np.zeros([160, 160])
                arr_index = get_indexap_annulus(rad_mock, radi_apert)
                print 'arr_index', arr_index
                for elem in arr_index:
                        y = elem[0] + y_index_center_avg
                        x = elem[1] + x_index_center_avg
                        img_test = test_ret_flux_aperture(img_test, [y,x], radi_apert)
                hdu = pf.PrimaryHDU(img_test)
                hdulist = pf.HDUList([hdu])
                hdulist.writeto('test.fits', clobber=True)
                subprocess.call('/home/apps/ds9/ds9 ' + 'test.fits', shell = True)    
                useless = raw_input('enter any key to move on')


def plot_all_5sd():
        arr_radiusmock = np.array([10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80]).astype(int)
        arr_5sd = pf.open('ShaneAO/Jun15/arr_5sd_correct.fits')[0].data
        for row_index in np.arange(len(arr_5sd)):
                plt.plot(arr_radiusmock, arr_5sd[row_index, :], 'o-')                        
        plt.gca().invert_yaxis()
        plt.show()


def ret_sd_fluxap(img, rad_mock, radi_apert_flux):
        #perform aperture flux measurements at certain radius
        #calculate standard deviation of flux measurements
        #------

        index_imgcenter = ret_index_center(img)
        arr_index = get_indexap_annulus(rad_mock, radi_apert_flux)
        arr_flux_radi = []
        for index2d in arr_index:
                index_apcenter = np.array(index2d) + index_imgcenter
                flux_apert_radi = ret_flux_aperture(img, index_apcenter, radi_apert_flux)
                arr_flux_radi.append(flux_apert_radi)

        sd_ring_r = np.std(np.array(arr_flux_radi))
        return sd_ring_r


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


def ret_index_center(img):
        # Returns index of center of img
        #------
        # Input: img(N-dimensional array)
        # Output: N-element array with index of center
        #------
        if not all(np.array(img.shape)%2 == 1):
                print '------'
                print 'ERROR, img shape not odd'
                print '------'

        imgshape = img.shape
        return (np.array(img.shape).astype(int)-1)/2
        

def ret_flux_center_ap(img, radi_ap, ret_array=False):
        # Inputs: image(2d array(NxM), both N and M must be odd), radius of aperture(number)
        # Returns: total flux in aperture located at **CENTER OF IMG**(or array of counts)
        #------

        #------
        #check if img shape is odd in both dimensions
        #------
        if not all(np.array(img.shape)%2 == 1):
                print '------'
                print 'ERROR, img shape not odd'
                print '------'

        #------
        #deciphering inputs
        #------
        index_y = int(img.shape[0]-1)/2
        index_x = int(img.shape[1]-1)/2
        j = radi_ap

        #------
        #Knowing index of peak, take aperture of j pixel around this point
        #save total flux of aperture
        #------
        arr_totflux = []
        for x_circ in np.arange((2*j) + 1) - j:
                for y_circ in np.arange((2*j) + 1) - j:
                        if distance(x_circ, 0, y_circ, 0) < (j+0.5):
                                if index_y + y_circ >= img.shape[0] or index_y + y_circ < 0  or index_x + x_circ >= img.shape[1] or index_x + x_circ < 0:
                                        continue
                                else:
                                        arr_totflux.append(img[index_y+y_circ, index_x+x_circ])
        tot_flux = sum(arr_totflux)
        if ret_array:
                return arr_totflux
        else:
                return tot_flux
        


def ret_flux_aperture(img, index, radi_ap, ret_array=False):
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
                                if index_y + y_circ >= img.shape[0] or index_y + y_circ < 0  or index_x + x_circ >= img.shape[1] or index_x + x_circ < 0:
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


def mask_aperture(img, index, radi_ap, const_mask):
        #Inputs: 
        # image(2d array), indexes(2-element(y,x) array), 
        # radius of aperture(number), constant used to replace pixel values(number)
        #Returns: 
        # updated image with masked aperture(2d array)
        #------


        #deciphering inputs
        index_y = index[0]
        index_x = index[1]
        j = radi_ap

        #------
        #Knowing index of aperture center, take aperture of radius j pixel around this point
        #save total flux of aperture
        #------
        arr_totflux = []
        for x_circ in np.arange((2*j) + 1) - j:
                for y_circ in np.arange((2*j) + 1) - j:
                        if distance(x_circ, 0, y_circ, 0) < (j+0.5):
                                if index_y + y_circ >= img.shape[0] or index_x + x_circ >= img.shape[1]:
                                        continue
                                else:
                                        img[index_y+y_circ, index_x+x_circ] = const_mask

        return img

def pf_loadfits(filename_img, print_cond = False):
        #Inputs:
        # path to img file(str), print_cond set to True if you want to be informed if file doesn't exist(bool)
        #Output:
        # array stored within fits file
        #------

        if exists(filename_img):
                try:
                        return pf.open(filename_img)[0].data
                except:
                        return np.array([])
        else:
                if print_cond:
                        print filename_img, 'doesnt exist'
                return np.array([])


def pf_savefits(arr, filename):
        #Inputs are: array for saving to fit & filename for fits.
        #------

        hdu = pf.PrimaryHDU(arr)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto(filename, clobber=True)


def pf_edit_hdr(filename, hdr_tag, hdr_elem):
        #Edits header of fits file. Can be used to replace or add tag
        #Inputs:
        # filename of fits file(str), tag of header element(str), info to add as header element
        #------

        file_fits = pf.open(filename, mode = 'update')
        img_hdr = file_fits[0].header
        img_hdr[hdr_tag] = hdr_elem
        file_fits.flush()
        file_fits.close()

def print_timenow():
        print datetime.datetime.now().strftime('%m-%d %H:%M:%S')



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
        pf_savefits(mat_corr, filename_output)


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
        pf_savefits(arr_filenumbs_corr, filename_filenumbs_corr)
        print 'created:', filename_filenumbs_corr
        pf_savefits(arr_indexdates_corr, filename_indexdates_corr)
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
        
        max_rad_ring = 4 #radius of masking aperture for primary

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
                        pf_savefits(matrix_corr, filename_matrix)
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

def distance(x1, x2, y1, y2):
	return math.sqrt( ((x1 - x2)**2.) + ((y1 - y2)**2.) )

def cartesian2polar(var1, var2, inv = False):
        if inv:
                r = var1 
                angle = var2 #in degrees initially
                x = r*np.cos(angle*(2*np.pi)/360)
                y = r*np.sin(angle*(2*np.pi)/360)
                return x, y
        else:
                x = var1
                y = var2
                r = np.sqrt(x**2. + y**2.)
                angle = np.arctan2(y, x)
                return r, angle*360/(2*np.pi) #returns angle in degrees

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

        for setnumb in [1]: #arr_setnumbs:
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
        #calculate fwhm of star centered at middle of img, given 


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

        
def open_fits(filename):
        #filename = raw_input('Name of fits file(Do not include .fits):')
        subprocess.call('/home/apps/ds9/ds9 ' + filename, shell = True)


def openview_pickle(filename):
        img_final = pickle.load(open(filename, "rb"))
        hdu = pf.PrimaryHDU(img_final)
        hdulist = pf.HDUList([hdu])
        hdulist.writeto('test.fits', clobber = True)
        subprocess.call('/home/apps/ds9/ds9 ' + 'test.fits', shell = True)

def open_fits_taurus(filename):
        subprocess.call(ds9_taurus + filename, shell = True)
