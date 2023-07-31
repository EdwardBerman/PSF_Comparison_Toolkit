import unittest
test = unittest.TestCase
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc,rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
import matplotlib.pyplot as plt
import ipdb, pdb
from astropy.io import fits
from scipy.stats import chi2
import galsim
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import mean_squared_error, normalized_root_mse
from catalogutils import AttrDict, set_rc_params
from catalog_hsm_fitter import do_hsm_fit
import catalogaugmenter as ca
import glob
import re
from astropy.table import Table, vstack, hstack, Column
import pandas as pd

class plot():
    def __init__(self, catalog_object, psf_object):
        self.catalog_obj = catalog_object
        self.psf_obj = psf_object
        catalog_object.add_noise_flux([psf_object])
        self.stars = catalog_object.data['VIGNET'] #formerly starmaker
        self.psfs  = catalog_object.data[psf_object.nameColumn()] #formerly psfmaker
        self.err_stamps = catalog_object.data['ERR_VIGNET']
        catalog_object.crop([psf_object], vignet_size=75)
        catalog_object.save_new(outname=catalog_object.catalog)
        try:
            self.cropped_stars = catalog_object.data['VIGNET_CROPPED']
        except:
            self.cropped_stars = catalog_object.data['VIGNET']
        try:
            self.cropped_psfs = catalog_object.data[psf_object.nameColumn()+'_CROPPED']
        except:
            self.cropped_psfs = catalog_object.data[psf_object.nameColumn()]
        try:
            self.cropped_err_stamps = catalog_object.data['ERR_VIGNET_CROPPED']
        except:
            self.cropped_err_stamps = catalog_object.data['ERR_VIGNET']


class resid_plot(plot):
    def __init__(self, catalog_object, psf_object):
        super().__init__(catalog_object, psf_object)
        self.avg_star = np.nanmean(self.cropped_stars, axis=0)
        self.avg_psf = np.nanmean(self.cropped_psfs, axis=0)
        self.avg_residual = np.zeros(self.avg_star.shape)
        set_rc_params(fontsize=16)
        self.titles = []
        self.sum_residuals = 0.0
        self.fwhm_model = []
        self.fwhm_data = []

    def preprocessing(self):
        stars = self.cropped_stars
        for i in range(len(stars)):
            for j in range(stars[0].shape[0]):
                for k in range(stars[0].shape[1]):
                    if stars[i][j, k] < -1000:
                        stars[i][j, k] = np.nan
        for i in range(len(stars)):
            stars[i] = stars[i]/np.nansum(stars[i])
        self.cropped_stars = stars
        
        psfs = self.cropped_psfs
        for i in range(len(psfs)):
            for j in range(psfs[0].shape[0]):
                for k in range(psfs[0].shape[1]):
                    if psfs[i][j, k] < -1000:
                        psfs[i][j, k] = np.nan
        for i in range(len(psfs)):
            psfs[i] = psfs[i]/np.nansum(psfs[i] + 1e-10)
        self.cropped_psfs = psfs

    def return_residuals_sum(self):
        return self.sum_residuals

    def set_residuals_sum(self):
        self.sum_residuals = np.nansum(self.avg_residual)

    def calc_fwhm(self, pix_scale=0.03):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        fwhm_star = []
        fwhm_psf = []
        
        for i in range(len(stars)):
            if type(stars[i]) is galsim.image.Image:
                gs_object = stars[i]
            else:
                # HSM fits fail if there are too many negative pixels
                gs_object = galsim.Image(stars[i], wcs=galsim.PixelScale(pix_scale))    
            fwhm_star.append(gs_object.calculateFWHM())

        for i in range(len(psfs)):
            if type(psfs[i]) is galsim.image.Image:
                gs_object = psfs[i]
            else:
                gs_object = galsim.Image(psfs[i], wcs=galsim.PixelScale(pix_scale))
            fwhm_psf.append(gs_object.calculateFWHM())
        
        self.fwhm_model = fwhm_star
        self.fwhm_data = fwhm_psf

    def return_fwhm_star(self):
        return self.fwhm_model

    def return_fwhm_psf(self):
        return self.fwhm_data

    def set_titles(self, titles):
        self.titles = titles

    def save_figure(self, outname):
        vmin = -1e-2 #np.nanmin(self.avg_psf)
        vmax = np.maximum(np.nanmax(self.avg_star), np.nanmax(self.avg_psf))
        
        vmin2 = np.nanmin(self.avg_residual)
        vmax2 = np.nanmax(self.avg_residual)

        norm = colors.SymLogNorm(vmin=vmin, vmax=vmax, linthresh=1e-6)
        norm2 = colors.SymLogNorm(vmin=vmin2, vmax=vmax2, linthresh=1e-4)
        cmap=plt.cm.turbo
        cmap2=plt.cm.bwr_r
        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=[15,7], tight_layout=True)
        
        im = axs[0].imshow(self.avg_star, norm=norm, cmap=cmap)
        axs[0].set_title(self.titles[0])
        divider = make_axes_locatable(axs[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[0].axvline((self.avg_star.shape[0]-1)*0.5,color='black')
        axs[0].axhline((self.avg_star.shape[1]-1)*0.5,color='black')
        
        im = axs[1].imshow(self.avg_psf, norm=norm, cmap=cmap)
        axs[1].set_title(self.titles[1])
        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[1].axvline((self.avg_psf.shape[0]-1)*0.5,color='black')
        axs[1].axhline((self.avg_psf.shape[1]-1)*0.5,color='black')
        
        im = axs[2].imshow(self.avg_residual, norm=norm2, cmap=cmap2)
        axs[2].set_title(self.titles[2])
        divider = make_axes_locatable(axs[2])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[2].axvline((self.avg_residual.shape[0]-1)*0.5,color='black')
        axs[2].axhline((self.avg_residual.shape[1]-1)*0.5,color='black')

        fig.tight_layout()

        plt.savefig(outname)

class mean_relative_error_plot(resid_plot):
    def __init__(self, catalog_object, psf_object):
        super().__init__(catalog_object, psf_object)

    def set_residuals(self):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        resids = []
        for i in range(len(stars)):
            relative_error = np.zeros((stars[0].shape[0], stars[0].shape[1]))
            for j in range(stars[0].shape[0]):
                for k in range(stars[0].shape[1]):
                    relative_error[j,k] = (stars[i][j,k] - psfs[i][j,k])/(stars[i][j,k] + 1.0e-6)
            resids.append(relative_error)
        self.avg_residual = np.nanmean(resids, axis=0)
        self.set_residuals_sum()

class mean_absolute_error_plot(resid_plot):
    def __init__(self, catalog_object, psf_object):
        super().__init__(catalog_object, psf_object)

    def set_residuals(self):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        resids = []
        for i in range(len(stars)):
            absolute_error = np.zeros((stars[0].shape[0], stars[0].shape[1]))
            for j in range(stars[0].shape[0]):
                for k in range(stars[0].shape[1]):
                    absolute_error[j,k] = np.abs(stars[i][j,k] - psfs[i][j,k])
            resids.append(absolute_error)
        self.avg_residual = np.nanmean(resids, axis=0)
        self.set_residuals_sum()

class mean_error_plot(resid_plot):
    def __init__(self, catalog_object, psf_object):
        super().__init__(catalog_object, psf_object)

    def set_residuals(self):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        resids = []
        for i in range(len(stars)):
            error = np.zeros((stars[0].shape[0], stars[0].shape[1]))
            for j in range(stars[0].shape[0]):
                for k in range(stars[0].shape[1]):
                    error[j,k] = (stars[i][j,k] - psfs[i][j,k]) 
            resids.append(error)
        self.avg_residual = np.nanmean(resids, axis=0)
        self.set_residuals_sum()

class chi_2_error_plot(resid_plot):
    def __init__(self, catalog_object, psf_object):
        super().__init__(catalog_object, psf_object)
        self.chi2_vals = []

    def set_residuals(self, polydim=2, polydeg=1):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        noise_map = self.cropped_err_stamps

        resids = []
        npix = stars[0].shape[0] * stars[0].shape[1]
        nparams = np.prod([polydeg+i+1 for i in range(polydim)])/polydim
        dof = npix - nparams
        chi2_vals = []

        for i in range(len(stars)):
            resids.append(np.divide(np.square(stars[i] - psfs[i]), np.square(noise_map[i] + 1.0e-6)))
            chi2_finite = np.divide(np.square(stars[i] - psfs[i]), np.square(noise_map[i] + 1.0e-6))[np.isfinite(np.divide(np.square(stars[i] - psfs[i]), np.square(noise_map[i] + 1.0e-6)))]
            ddof = chi2_finite.size-nparams
            chi2_vals.append(chi2_finite.sum()/ddof)
        self.avg_residual = np.nanmean(resids, axis=0)
        self.set_residuals_sum()
        self.chi2_vals = chi2_vals

    def return_chi2_vals(self):
        return self.chi2_vals


class ssim_error_plot(resid_plot):
    def __init__(self, catalog_object, psf_object):
        super().__init__(catalog_object, psf_object)

    def set_residuals(self):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        
        ssims = []
        ssim_val = []
        nrmse = []
        
        for i in range(len(stars)):
            ssim_res = ssim(psfs[i], stars[i], full=True, win_size=3, data_range=psfs.max()-psfs.min())
            ssim_val.append(ssim_res[0])
            ssims.append(ssim_res[1] - np.median(ssim_res[1]))
            nrmse.append(normalized_root_mse(stars[i], psfs[i]))
        
        self.avg_residual = np.nanmean(ssims, axis=0)
        self.set_residuals_sum()


directory_path = "/home/eddieberman/research/mcclearygroup/cweb_psf/augmented_starcatalogs"
file_list = glob.glob(directory_path + "/*augmented_corrected_psf_starcat.fits")

def extract_3_numbers(filename):
    pattern = r'\d{3}'
    matches = re.findall(pattern, filename)
    return matches

#columns = ['Detector', 'Filter', 'PSFex: SSIM ', 'PSFex: Chi-Square', 'PSFex: Absolute Error', 'PSFex: Relative Error','WebbPSF: SSIM ', 'WebbPSF: Chi-Square', 'WebbPSF: Absolute Error', 'WebbPSF: Relative Error']
columns = ['Detector', 'Filter', 'PSFex: Absolute Error', 'WebbPSF: Absolute Error', 'PSFex: Chi-Square', 'WebbPSF: Chi-Square','PSFex: Relative Error', 'WebbPSF: Relative Error', 'PSFex: SSIM ', 'WebbPSF: SSIM ', 'PSFex: FWHM Residual', 'WebbPSF: FWHM Residual', 'PSFex: Reduced Chi Square', 'WebbPSF: Reduced Chi Square']
data_table = Table(names=columns, dtype=['S10', 'S10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'])

for file in file_list:
    try:
        catalog_object_name = file
        epsfex_object_name = '/home/eddieberman/research/mcclearygroup/cweb_psf/augmented_starcatalogs/f150w_b2_starcat.psf'
        
        chi2_name = file.replace("augmented_corrected_psf_starcat.fits", "psfex_chi2.png")
        mre_name = file.replace("augmented_corrected_psf_starcat.fits", "psfex_mre_resid.png")
        abs_name = file.replace("augmented_corrected_psf_starcat.fits", "psfex_abs_resid.png")
        ssim_name = file.replace("augmented_corrected_psf_starcat.fits", "psfex_ssim.png")

        mean_absolute_error_plot_psfex = mean_absolute_error_plot(ca.catalog(catalog_object_name), ca.epsfex(epsfex_object_name))
        mean_absolute_error_plot_psfex.preprocessing()
        mean_absolute_error_plot_psfex.set_residuals()
        sum_residuals_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
        mean_absolute_error_plot_psfex.calc_fwhm()
        avg_star_fwhm = np.nanmean(mean_absolute_error_plot_psfex.return_fwhm_star())
        avg_psf_fwhm = np.nanmean(mean_absolute_error_plot_psfex.return_fwhm_psf())
        mean_absolute_error_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Absolute Error \n Sum Absolute Error = {sum_residuals_abs_psfex}'])
        mean_absolute_error_plot_psfex.save_figure(outname=abs_name)
        psfex_fwhm_residual = np.nanmean(np.array(mean_absolute_error_plot_psfex.return_fwhm_star()) - np.array(mean_absolute_error_plot_psfex.return_fwhm_psf()))

        mean_chi2_plot_psfex = chi_2_error_plot(ca.catalog(catalog_object_name), ca.epsfex(epsfex_object_name))
        mean_chi2_plot_psfex.preprocessing()
        mean_chi2_plot_psfex.set_residuals()
        sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
        mean_chi2_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Chi2 Error \n Sum Chi2 Error = {sum_residuals_chi2_psfex}'])
        mean_chi2_plot_psfex.save_figure(outname=chi2_name)
        reduced_chi2_psfex = np.nanmean(mean_chi2_plot_psfex.return_chi2_vals())

        mean_mre_plot_psfex = mean_relative_error_plot(ca.catalog(catalog_object_name), ca.epsfex(epsfex_object_name))
        mean_mre_plot_psfex.preprocessing()
        mean_mre_plot_psfex.set_residuals()
        sum_residuals_mre_psfex = mean_mre_plot_psfex.return_residuals_sum()
        mean_mre_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average MRE Error \n Sum MRE Error = {sum_residuals_mre_psfex}'])
        mean_mre_plot_psfex.save_figure(outname=mre_name)

        ssim_plot_psfex = ssim_error_plot(ca.catalog(catalog_object_name), ca.epsfex(epsfex_object_name))
        ssim_plot_psfex.preprocessing()
        ssim_plot_psfex.set_residuals()
        sum_residuals_ssim_psfex = ssim_plot_psfex.return_residuals_sum()
        ssim_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average SSIM Error \n Sum SSIM Error = {sum_residuals_ssim_psfex}'])
        ssim_plot_psfex.save_figure(outname=ssim_name)

        webbpsf_object_name = '/home/eddieberman/research/mcclearygroup/cweb_psf/augmented_starcatalogs/f150w_b2_WebbPSF.fits'
        chi2_name = file.replace("augmented_corrected_psf_starcat.fits", "webb_chi2.png")
        mre_name = file.replace("augmented_corrected_psf_starcat.fits", "webb_mre_resid.png")
        abs_name = file.replace("augmented_corrected_psf_starcat.fits", "webb_abs_resid.png")
        ssim_name = file.replace("augmented_corrected_psf_starcat.fits", "webb_ssim.png")

        mean_absolute_error_plot_webbpsf = mean_absolute_error_plot(ca.catalog(catalog_object_name), ca.webb_psf(webbpsf_object_name))
        mean_absolute_error_plot_webbpsf.preprocessing()
        mean_absolute_error_plot_webbpsf.set_residuals()
        sum_residuals_abs_webb = mean_absolute_error_plot_webbpsf.return_residuals_sum()
        mean_absolute_error_plot_webbpsf.calc_fwhm()
        avg_star_fwhm = np.nanmean(mean_absolute_error_plot_webbpsf.return_fwhm_star())
        avg_psf_fwhm = np.nanmean(mean_absolute_error_plot_webbpsf.return_fwhm_psf())
        mean_absolute_error_plot_webbpsf.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Absolute Error \n Sum Absolute Error = {sum_residuals_abs_webb}'])
        mean_absolute_error_plot_webbpsf.save_figure(outname=abs_name)
        webbpsf_fwhm_residual = np.nanmean(np.array(mean_absolute_error_plot_webbpsf.return_fwhm_star()) - np.array(mean_absolute_error_plot_webbpsf.return_fwhm_psf()))

        mean_chi2_plot_webbpsf = chi_2_error_plot(ca.catalog(catalog_object_name), ca.webb_psf(webbpsf_object_name))
        mean_chi2_plot_webbpsf.preprocessing()
        mean_chi2_plot_webbpsf.set_residuals()
        sum_residuals_chi2_webb = mean_chi2_plot_webbpsf.return_residuals_sum()
        mean_chi2_plot_webbpsf.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Chi2 Error \n Sum Chi2 Error = {sum_residuals_chi2_webb}'])
        mean_chi2_plot_webbpsf.save_figure(outname=chi2_name)
        reduced_chi2_webb = np.nanmean(mean_chi2_plot_webbpsf.return_chi2_vals())

        mean_mre_plot_webbpsf = mean_relative_error_plot(ca.catalog(catalog_object_name), ca.webb_psf(webbpsf_object_name))
        mean_mre_plot_webbpsf.preprocessing()
        mean_mre_plot_webbpsf.set_residuals()
        sum_residuals_mre_webb = mean_mre_plot_webbpsf.return_residuals_sum()
        mean_mre_plot_webbpsf.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average MRE Error \n Sum MRE Error = {sum_residuals_mre_webb}'])
        mean_mre_plot_webbpsf.save_figure(outname=mre_name)

        ssim_plot_webbpsf = ssim_error_plot(ca.catalog(catalog_object_name), ca.webb_psf(webbpsf_object_name))
        ssim_plot_webbpsf.preprocessing()
        ssim_plot_webbpsf.set_residuals()
        sum_residuals_ssim_webb = ssim_plot_webbpsf.return_residuals_sum()
        ssim_plot_webbpsf.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average SSIM Error \n Sum SSIM Error = {sum_residuals_ssim_webb}'])
        ssim_plot_webbpsf.save_figure(outname=ssim_name)
        
        letter_number_codes = re.findall(r'[A-Za-z]\d', file)
        print('\n', 'F'+extract_3_numbers(file)[0]+'W_', letter_number_codes[1], ' Sum Absolute Error PSFex = ', sum_residuals_abs_psfex, ' Sum Absolute Error WebbPSF = ', sum_residuals_abs_webb, ' Sum Chi2 Error PSFex = ', sum_residuals_chi2_psfex, ' Sum Chi2 Error WebbPSF = ', sum_residuals_chi2_webb, ' Sum MRE Error PSFex = ', sum_residuals_mre_psfex, ' Sum MRE Error WebbPSF = ', sum_residuals_mre_webb, ' Sum SSIM Error PSFex = ', sum_residuals_ssim_psfex, ' Sum SSIM Error WebbPSF = ', sum_residuals_ssim_webb, 'FWHM Residual PSFex', psfex_fwhm_residual, 'FWHM Residual WebbPSF', webbpsf_fwhm_residual, 'Reduced Chi Square PSFex', reduced_chi2_psfex, 'Reduced Chi Square WebbPSF', reduced_chi2_webb)
        data_table.add_row(['F'+extract_3_numbers(file)[0]+'W', letter_number_codes[1], sum_residuals_abs_psfex, sum_residuals_abs_webb, sum_residuals_chi2_psfex, sum_residuals_chi2_webb, sum_residuals_mre_psfex, sum_residuals_mre_webb, sum_residuals_ssim_psfex, sum_residuals_ssim_webb, psfex_fwhm_residual, webbpsf_fwhm_residual, reduced_chi2_psfex, reduced_chi2_webb])
    except Exception as e:
        print("An error occurred: ", e)

df = data_table.to_pandas()
sorted_df = df.sort_values(by=['Filter', 'Detector'])
pd.set_option("display.max_columns", None)  # Display all columns
pd.set_option("display.width", None)  # Disable column width truncation

# Print the DataFrame (which will show all columns)
print(sorted_df)

# Optionally, convert the DataFrame back to an Astropy table
table = Table.from_pandas(sorted_df)
