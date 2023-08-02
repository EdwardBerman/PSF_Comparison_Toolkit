import glob
from astropy.io import fits
from astropy.table import Table, vstack
from catalogaugmenter import catalog, psf
from catalogaugmenter import webb_psf, epsfex, shopt, piff_psf 
#from catalogplotter import ResidPlots
import os
import re 
from datetime import datetime, timedelta

#Make augmented catalogs with columns for each psf fitter than use the plotter with these new catalogs

ims = glob.glob('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/jw0*cal.fits')
print(len(ims))
f115_a1_list=[]
f115_a2_list=[]
f115_a3_list=[]
f115_a4_list=[]
f115_b1_list=[]
f115_b2_list=[]
f115_b3_list=[]
f115_b4_list=[]
f150_a1_list=[]
f150_a2_list=[]
f150_a3_list=[]
f150_a4_list=[]
f150_b1_list=[]
f150_b2_list=[]
f150_b3_list=[]
f150_b4_list=[]
f277_a5_list=[]
f277_b5_list=[]
f444_a5_list=[]
f444_b5_list=[]

for im in ims:
    imhead = fits.getheader(im, ext=0)
    filt = imhead['FILTER']
    det = imhead['DETECTOR'].replace('LONG', '5')
    if (det =='NRCA1') & (filt == 'F115W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f115_a1_list.append(catalog_obj)
        #f115_a1_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCA1') & (filt == 'F150W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f150_a1_list.append(catalog_obj)
        #f150_a1_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCA2') & (filt == 'F115W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f115_a2_list.append(catalog_obj)
        #f115_a2_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCA2') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_starcat.psf'
            print(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_a2_list.append(catalog_obj)
        except:
            pass
        #f150_a2_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCA3') & (filt == 'F115W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f115_a3_list.append(catalog_obj)
        #f115_a3_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCA3') & (filt == 'F150W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f150_a3_list.append(catalog_obj)
        #f150_a3_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCA4') & (filt == 'F115W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f115_a4_list.append(catalog_obj)
        #f115_a4_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCA4') & (filt == 'F150W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f150_a4_list.append(catalog_obj)
        #f150_a4_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCA5') & (filt == 'F277W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f277_a5_list.append(catalog_obj)
        #f277_a5_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCA5') & (filt == 'F444W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f444_a5_list.append(catalog_obj)
        #f444_a5_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCB1') & (filt == 'F115W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f115_b1_list.append(catalog_obj)
        #f115_b1_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCB1') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_b1_list.append(catalog_obj)
            #f150_b1_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB2') & (filt == 'F115W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f115_b2_list.append(catalog_obj)
            #f115_b2_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB2') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_b2_list.append(catalog_obj)
            #f150_b2_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB3') & (filt == 'F115W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f115_b3_list.append(catalog_obj)
        #f115_b3_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCB3') & (filt == 'F150W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f150_b3_list.append(catalog_obj)
        #f150_b3_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCB4') & (filt == 'F115W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f115_b4_list.append(catalog_obj)
        except:
            pass
        #f115_b4_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCB4') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_b4_list.append(catalog_obj)
            #f150_b4_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB5') & (filt == 'F277W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f277_b5_list.append(catalog_obj)
        #f277_b5_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCB5') & (filt == 'F444W'):
        catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_starcat.fits' ))
        catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
        string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
        string2 = string+'_starcat.psf'
        catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
        f444_b5_list.append(catalog_obj)
        #f444_b5_list.append(fits.getheader(im, ext=0)['DATE'])
'''
cat = f277_b5_list[0]
cat.concatenate_catalogs(f277_b5_list[1])
cat.concatenate_catalogs(f277_b5_list[2])
cat.concatenate_catalogs(f277_b5_list[3])
cat.save_new('f277_b5_augmented_corrected_psf_starcat.fits')

cat = f444_a5_list[0]
cat.concatenate_catalogs(f444_a5_list[1])
cat.concatenate_catalogs(f444_a5_list[2])
cat.concatenate_catalogs(f444_a5_list[3])
cat.save_new('f444_a5_augmented_corrected_psf_starcat.fits')

cat = f277_a5_list[0]
cat.concatenate_catalogs(f277_a5_list[1])
cat.concatenate_catalogs(f277_a5_list[2])
cat.concatenate_catalogs(f277_a5_list[3])
cat.save_new('f277_a5_augmented_corrected_psf_starcat.fits')

cat = f150_b4_list[0]
cat.concatenate_catalogs(f150_b4_list[1])
cat.concatenate_catalogs(f150_b4_list[2])
#cat.concatenate_catalogs(f150_b4_list[3])
cat.save_new('f150_b4_augmented_corrected_psf_starcat.fits')

cat = f115_b4_list[0]
cat.concatenate_catalogs(f115_b4_list[1])
#cat.concatenate_catalogs(f115_b4_list[2])
#cat.concatenate_catalogs(f115_b4_list[3])
cat.save_new('f115_b4_augmented_corrected_psf_starcat.fits')

cat = f150_b3_list[0]
cat.concatenate_catalogs(f150_b3_list[1])
cat.concatenate_catalogs(f150_b3_list[2])
cat.concatenate_catalogs(f150_b3_list[3])
cat.save_new('f150_b3_augmented_corrected_psf_starcat.fits')

cat = f115_b3_list[0]
cat.concatenate_catalogs(f115_b3_list[1])
cat.concatenate_catalogs(f115_b3_list[2])
cat.concatenate_catalogs(f115_b3_list[3])
cat.save_new('f115_b3_augmented_corrected_psf_starcat.fits')

cat = f150_b2_list[0]
#cat.concatenate_catalogs(f150_b2_list[1])
#cat.concatenate_catalogs(f150_b2_list[2])
#cat.concatenate_catalogs(f150_b2_list[3])
cat.save_new('f150_b2_augmented_corrected_psf_starcat.fits')

cat = f115_b2_list[0]
#cat.concatenate_catalogs(f115_b2_list[1])
#cat.concatenate_catalogs(f115_b2_list[2])
#cat.concatenate_catalogs(f115_b2_list[3])
cat.save_new('f115_b2_augmented_corrected_psf_starcat.fits')

cat = f150_b1_list[0]
cat.concatenate_catalogs(f150_b1_list[1])
cat.concatenate_catalogs(f150_b1_list[2])
#cat.concatenate_catalogs(f150_b1_list[3])
cat.save_new('f150_b1_augmented_corrected_psf_starcat.fits')

cat = f115_b1_list[0]
cat.concatenate_catalogs(f115_b1_list[1])
cat.concatenate_catalogs(f115_b1_list[2])
cat.concatenate_catalogs(f115_b1_list[3])
cat.save_new('f115_b1_augmented_corrected_psf_starcat.fits')

cat = f150_a4_list[0]
cat.concatenate_catalogs(f150_a4_list[1])
cat.concatenate_catalogs(f150_a4_list[2])
cat.concatenate_catalogs(f150_a4_list[3])
cat.save_new('f150_a4_augmented_corrected_psf_starcat.fits')

cat = f115_a4_list[0]
cat.concatenate_catalogs(f115_a4_list[1])
cat.concatenate_catalogs(f115_a4_list[2])
cat.concatenate_catalogs(f115_a4_list[3])
cat.save_new('f115_a4_augmented_corrected_psf_starcat.fits')

cat = f150_a3_list[0]
cat.concatenate_catalogs(f150_a3_list[1])
cat.concatenate_catalogs(f150_a3_list[2])
cat.concatenate_catalogs(f150_a3_list[3])
cat.save_new('f150_a3_augmented_corrected_psf_starcat.fits')

cat = f115_a3_list[0]
cat.concatenate_catalogs(f115_a3_list[1])
cat.concatenate_catalogs(f115_a3_list[2])
cat.concatenate_catalogs(f115_a3_list[3])
cat.save_new('f115_a3_augmented_corrected_psf_starcat.fits')

cat = f150_a2_list[0]
cat.concatenate_catalogs(f150_a2_list[1])
#cat.concatenate_catalogs(f150_a2_list[2])
#cat.concatenate_catalogs(f150_a2_list[3])
cat.save_new('f150_a2_augmented_corrected_psf_starcat.fits')

cat = f115_a2_list[0]
cat.concatenate_catalogs(f115_a2_list[1])
cat.concatenate_catalogs(f115_a2_list[2])
cat.concatenate_catalogs(f115_a2_list[3])
cat.save_new('f115_a2_augmented_corrected_psf_starcat.fits')

cat = f150_a1_list[0]
cat.concatenate_catalogs(f150_a1_list[1])
cat.concatenate_catalogs(f150_a1_list[2])
cat.concatenate_catalogs(f150_a1_list[3])
cat.save_new('f150_a1_augmented_corrected_psf_starcat.fits')
'''
cat = f115_a1_list[0]
cat.concatenate_catalogs(f115_a1_list[1])
cat.concatenate_catalogs(f115_a1_list[2])
cat.concatenate_catalogs(f115_a1_list[3])
cat.save_new('/home/eddieberman/research/mcclearygroup/cweb_psf/augmented_starcatalogs/f115_a1_augmented_corrected_psf_starcat.fits')

