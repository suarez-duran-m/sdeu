"""
Created on Wed Mar  2 14:13:12 2022

@author: katarina
"""

import numpy as np
import matplotlib.pyplot as plt
import ROOT
from array import array

plt.rcParams.update({'font.size': 12})
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import os
from scipy.optimize import curve_fit
# from astropy.time import Time
import glob
from scipy import stats


# import csv, pickle

class CCH:
    # GPStime = unixTime - 315964782; cchs are in gpstime, ochs in unix
    uni_gps = 315964782
    pm = 3  # max. allowed discrepancy between event (och) and cch files [s]
    faulty_pmt = 0.06  # PMT excluded if maximum within fty_pmt*len(histo x-axis)

    def __init__(self):
        self.uberdict = {}
        self.dirs = {}
        self.plot_files = {}
        self.delta_mean = {}
        self.cols = ['r', 'orange', 'royalblue']
        return

    def load_data(self, histo_cq, histo_q, offLine_cqpk, offLine_cqpkErr, \
                  offLine_qpk, offLine_qpkErr, vemFromHisto, signal, ldfR, spDist, \
                  energy, energyErr, zenith, key_st, save_cch_list=False, bad_pmts={}, och_diff=True, \
                  keep_individual=False):
        """
        Loads CCH and CH. Combines them into averages for different SDs and
        stores them in the uberdict under uberdict[SD_id]['data_cch'] and
        uberdict[SD_id]['data_och'] as a matrix with 0th column the FADC count,
        and 1st-3rd column corresponding to the count in each PMT.
        Automatically ommits SDs with ID<=100.
        Plots all the unpaired CCH.
        If specified, saves a list of all the CCH with their timestamps.
        If specified in bad_pmts, plots all the raw CCH and their paired OCH.

        Parameters
        ----------
        cch_file : string
            Path to the .txt file with CCH data.
        och_dir : string
            Path to the directory with OCH data in .dat format.
            Has to end with '/'.
        output_dir : string
            Path to the directory where all the results will be saved.
        save_cch_list : boolean
            If true, saves 'CCH_list.dat' in the output directory.
            This is a list of all the CCH in the form SD_ID GPS_time UTC_time.
            Times are those of CCH, not events.
            Default : False.
        bad_pmts : python dictionary
            The form is key: SD ID, element: list of PMT numbers.
            Plots all the raw CCH and their paired OCH (before averaging).
        och_diff : boolean
            Default is True. Has to be True if counts of OCH are not divided
            by bin width.
        keep_individual : boolean
            Default is False. If True, the data structure will keep indiidual
            CCH and OCH as well as their averages.

        Returns
        -------
        None.

        """
        output_dir = "results/"
        self.count = [0, 0, 0]
        self.dirs['out'] = output_dir
        self.dirs['plot'] = output_dir + "plots/"
        if not os.path.isdir(self.dirs['out']):
            os.makedirs(self.dirs['out'])
        if not os.path.isdir(self.dirs['plot']):
            os.makedirs(self.dirs['plot'])

        # ==============
        # loading coincidence histograms data
        # ==============
        ch_split = 403  # number of channel :403 from where the tail is extended
        scale_bulk = 8

        n = [0, 0, 0]
        stId = str(key_st.split()[0])
        self.uberdict[stId] = {}
        gps_times = []
        charge = np.zeros(histo_cq.shape)
        charge[:, 0] = histo_cq[:, 0]
        charge_o = np.zeros(histo_q.shape)
        charge_o[:, 0] = histo_q[:, 0]
        # summing
        gps_times.append(float(key_st.split()[1]))
        for pmt in range(1, 4):
            # if not faulty pmt
            if np.average(histo_cq[:, pmt]) > 0.05:
                charge[:, pmt] = np.add(charge[:, pmt], histo_cq[:, pmt])
                charge_o[:, pmt] = np.add(charge_o[:, pmt], histo_q[:, pmt])
                n[pmt - 1] += 1

        # calculating average and error:
        charge_err = np.zeros(charge.shape)
        charge_o_err = np.zeros(charge_o.shape)
        for pmt in range(1, 4):
            if n[pmt - 1] != 0:
                charge_err[:, pmt] = np.sqrt(charge[:, pmt]) / n[pmt - 1]
                charge[:, pmt] = charge[:, pmt] / n[pmt - 1]
                # due to bin width
                charge_err[:ch_split, pmt] = charge_err[:ch_split, pmt] \
                                             / np.sqrt(8)
                charge_err[ch_split:, pmt] = charge_err[ch_split:, pmt] \
                                             / np.sqrt(32)
                charge_o_err[:, pmt] = np.sqrt(charge_o[:, pmt]) / n[pmt - 1]
                charge_o[:, pmt] = charge_o[:, pmt] / n[pmt - 1]
                # due to bin width
                charge_o_err[:401, pmt] = charge_o_err[:401, pmt] \
                                          / np.sqrt(8)
                charge_o_err[401:, pmt] = charge_o_err[401:, pmt] \
                                          / np.sqrt(32)
        #
        # storing the results
        self.uberdict[stId]['histo_count'] = n
        self.uberdict[stId]['gps_times'] = gps_times
        self.uberdict[stId]['cch'] = charge
        self.uberdict[stId]['cch_err'] = charge_err
        self.uberdict[stId]['och'] = charge_o
        self.uberdict[stId]['och_err'] = charge_o_err
        self.uberdict[stId]['offLine_qpk'] = offLine_qpk
        self.uberdict[stId]['offLine_qpkErr'] = offLine_qpkErr
        self.uberdict[stId]['offLine_cqpk'] = offLine_cqpk
        self.uberdict[stId]['offLine_cqpkErr'] = offLine_cqpkErr
        self.uberdict[stId]['vemFromHistogram'] = vemFromHisto
        self.uberdict[stId]['signal'] = signal
        self.uberdict[stId]['ldfR'] = ldfR
        self.uberdict[stId]['spDist'] = spDist
        self.uberdict[stId]['Energy'] = energy
        self.uberdict[stId]['EnergyErr'] = energyErr
        self.uberdict[stId]['Zenith'] = zenith
        #
        return

    def count_avg_cch(self):
        self.count = [0, 0, 0]
        for sd in self.uberdict.keys():
            for i in range(1, 4):
                if np.average(self.uberdict[sd]['cch'][:, i]) > 0:
                    self.count[i - 1] += 1
        return self.count

    def number_cch_fits(self):
        self.fits_no = [0, 0, 0]
        for sd in self.uberdict.keys():
            for i in range(1, 4):
                if not (self.uberdict[sd][i]['nofit_cch']):
                    self.fits_no[i - 1] += 1
        #
        return self.fits_no

    def drop_data(self, sd, pmt):
        self.uberdict[sd]['cch'][:, pmt] = 0
        self.uberdict[sd]['och'][:, pmt] = 0
        return

    def drop_data_sd(self, sd):
        del self.uberdict[sd]
        return

    def plot_nofit(self, sd, pmt=[1, 2, 3], histo_type='cch', err=True, mat=[], \
                   save=[], title='default'):
        """
        Plots a single PMT

        Parameters
        ----------
        sd : int
            SD ID.
        pmt : array of int
            PMT number(s).
        histo_type : str, optional
            The type of histogram out of 'cch', 'och'. The default is cch.
        err : bool, optional
            If errors should be plotted. The default is True.

        Returns
        -------
        None.

        """
        h, w = 12, 9

        plt.figure(figsize=(h, w))
        for n, i in enumerate(pmt):
            if len(mat) > 0 and not err:
                plt.errorbar(mat[:, 0], mat[:, i], marker='.', label='PMT %i' % i, \
                             c=self.cols[i - 1], linestyle='-')
                continue
            elif len(mat) > 0 and err and histo_type == 'cch':
                plt.errorbar(mat[:, 0], mat[:, i], \
                             yerr=np.append(np.sqrt(mat[:403, i]) / np.sqrt(8), \
                                            np.sqrt(mat[403:, i]) / np.sqrt(32)), capsize=2, \
                             marker='.', label=r'PMT %i' % i, c=self.cols[i - 1], linestyle='')
                plt.ylim(-1, 20)
                continue
            elif len(mat) > 0 and err and histo_type == 'och':
                plt.errorbar(mat[:, 0], mat[:, i], \
                             yerr=np.append(np.sqrt(mat[:401, i]) / np.sqrt(8), \
                                            np.sqrt(mat[401:, i]) / np.sqrt(32)), capsize=2, \
                             marker='.', label=r'PMT %i' % i, c=self.cols[i - 1], linestyle='')
                plt.ylim(-1, 150)
                continue
            if (self.uberdict[sd][histo_type][10, i] == 0) or err == False:
                plt.errorbar(self.uberdict[sd][histo_type][:, 0], \
                             self.uberdict[sd][histo_type][:, i], \
                             marker='.', label=r'PMT %i' % i, \
                             c=self.cols[i - 1], linestyle='')
                continue
            plt.errorbar(self.uberdict[sd][histo_type][:, 0], \
                         self.uberdict[sd][histo_type][:, i], \
                         yerr=self.uberdict[sd][histo_type + '_err'][:, i], \
                         marker='.', label='PMT %i' % i, c=self.cols[i - 1], \
                         linestyle='', capsize=2)
        plt.xlabel(r'Charge [FADC counts]')
        plt.ylabel('$\mathrm{d}N/\mathrm{d}Q$')
        plt.legend()
        plt.grid()
        # if title=='default':
        # plt.title('SD '+str(sd)+' '+histo_type.upper()+' '+\
        #     str(np.round(self.uberdict[list(self.uberdict.keys())[0]]\
        #                  ['gps_times'][0])) ) 
        # else:
        #     plt.title(title)
        if len(save) == 0:
            plt.show()
        else:
            save[0].savefig(bbox_inches='tight')
        plt.close('all')
        return

    def plot_nofit_all(self):
        for key in self.uberdict.keys():
            self.plot_nofit(key)
            self.plot_nofit(key, histo_type='och')
        return

    def rebin(self, data_h, cut_off=600, smooth_chan=8):
        new_dath = [[], [], [], []]
        for k in range(int(cut_off / smooth_chan)):
            new_dath[0].append( \
                np.sum(data_h[smooth_chan * k:(smooth_chan * (k + 1)), 0]) / smooth_chan)
            for j in range(1, 4):
                new_dath[j].append( \
                    np.sum(data_h[smooth_chan * k:(smooth_chan * (k + 1)), j]))
        return np.array(new_dath).T

    def differentiate(self, smooth_h):
        deri = [[], [], [], []]
        deri[0].append(smooth_h[0, 0])
        for j in range(1, 4):
            deri[j].append((smooth_h[1, j] - smooth_h[0, j]) / (smooth_h[1, 0] - \
                                                                smooth_h[0, 0]))
        for i, ch in enumerate(smooth_h[1:(len(smooth_h) - 1), :]):
            deri[0].append(ch[0])
            for j in range(1, 4):
                deri[j].append((smooth_h[i + 1, j] - smooth_h[i - 1, j]) / \
                               (smooth_h[i + 1, 0] - smooth_h[i - 1, 0]))
        return np.array(deri).T

    def find_bounds_cch(self, deri, smooth_h, data_h, pmtno):
        cut_start = 5
        h_fac = 0.25
        deri_min = np.argmin(deri[:, pmtno][cut_start:]) + cut_start
        zeroos = np.where(np.diff(np.sign(deri[:, pmtno][cut_start:])))[0] \
                 + cut_start
        try:
            hump = zeroos[zeroos < deri_min][-1]
        except:
            #print('No zeros b4 deri_min (l 429). CCH fit unsuccessful.')
            return 1
        left = hump - 1
        right = hump + 1

        hump_max = np.max(smooth_h[:, pmtno])
        while smooth_h[left, pmtno] > (1 - h_fac) * hump_max:
            left -= 1
        while smooth_h[right, pmtno] > (1 - h_fac) * hump_max:
            right += 1
        #
        # changing to data index
        left = np.abs(data_h[:, 0] - smooth_h[left, 0]).argmin()
        right = np.abs(data_h[:, 0] - smooth_h[right, 0]).argmin()
        #
        return [left, right], \
               [deri[hump, 0], (deri[hump + 1, 0] - deri[hump, 0]) / np.sqrt(12)]

    def find_bounds_och(self, derivat, smooth_h, data_h, pmtno):
        cut_start = 3
        broader_fit_right = -0.4
        val_fac_l = 0.1
        h_fac = 0.25
        low_vh_thresh = 0.85
        left_v_height = 1.3

        zeroos = np.where(np.diff(np.sign(derivat[:, pmtno][cut_start:])))[0] + \
                 cut_start
        derimin = np.argmin(derivat[:, pmtno][cut_start:]) + cut_start
        second_max = np.argmax(derivat[:, pmtno][derimin:]) + derimin
        second_min = np.argmin(derivat[:, pmtno][second_max:]) + second_max
        dist = second_max - derimin

        try:
            zero_b4_smax = zeroos[zeroos < second_max][-1]
        except:
            #print('Problem with zero_b4_smax (l 464).')
            return 1
        try:
            zero_after_smax = zeroos[zeroos > second_max][0]
        except:
            #print('Problem with zero_after_smax (l 467).')
            return 1

        left_v = derimin + int(np.round(val_fac_l * dist))
        hv = (zero_after_smax + zero_b4_smax) // 2

        right_h = second_min + int(broader_fit_right * (second_min - zero_after_smax))

        # verify the right&left hump border is below half of the hump         
        val = zero_b4_smax
        hump_max = np.max(smooth_h[val:, pmtno])
        val_min = np.min([smooth_h[left_v:hv, pmtno]])
        #
        # MSD, Checking if right limits is on fat bins
        #
        if smooth_h[right_h][0] > 3200:
            return 1
        while smooth_h[right_h, pmtno] > (1 - h_fac) * hump_max:
            right_h += 1

        # adjusting valley left to be at left_v_height*valley of the smoothed        
        if smooth_h[left_v, pmtno] < left_v_height * val_min:
            while True:
                if smooth_h[left_v - 1, pmtno] > left_v_height * val_min:
                    break
                left_v -= 1
        else:
            while True:
                if smooth_h[left_v + 1, pmtno] < left_v_height * val_min:
                    break
                left_v += 1

        if val_min / hump_max > low_vh_thresh:
            left_h = val
            right_v = zero_after_smax
            # converting to original x-indices
            right_h = np.abs(data_h[:, 0] - smooth_h[right_h, 0]).argmin()
            left_h = np.abs(data_h[:, 0] - smooth_h[left_h, 0]).argmin()
            left_v = np.abs(data_h[:, 0] - smooth_h[left_v, 0]).argmin()
            right_v = np.abs(data_h[:, 0] - smooth_h[right_v, 0]).argmin()
            _ = 0
            return [left_h, right_h], [left_v, right_v], _
        else:
            while (smooth_h[hv, pmtno] - val_min) > (hump_max - val_min) / 2:
                hv -= 1
            # converting to original x-indices
            right_h = np.abs(data_h[:, 0] - smooth_h[right_h, 0]).argmin()
            hv = np.abs(data_h[:, 0] - smooth_h[hv, 0]).argmin()
            left_v = np.abs(data_h[:, 0] - smooth_h[left_v, 0]).argmin()
            _ = 0
            return [hv, right_h], [left_v, hv], _

    def quadrat(self, x, a, b, c):
        return a * x ** 2 + b * x + c

    def quadrat_err(self, x, cov):
        ea, eb, ec = np.sqrt(np.diag(cov))
        eps_sq = x ** 4 * ea ** 2 + x ** 2 * eb ** 2 + ec ** 2 + 2 * x ** 3 * cov[0, 1] + \
                 2 * x ** 2 * cov[0, 2] + 2 * x * cov[1, 2]
        return np.sqrt(eps_sq)
    #
    # MSD: LogNormal function, from Kata's
    def lin_lognorm(self, x, m, k, x0, mu, sigma):
        A = (m * x0 + k) / np.exp(-(np.log(x0) - mu) ** 2 / (2 * sigma ** 2))
        y = np.piecewise(x, [x < x0, x >= x0], \
            [lambda x: m * x + k, lambda x: A * np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2))])
        return y

    def quad_fit(self, sd, pmt, x_arr, y_arr, y_err, cch='cch', fitrange=[]):
        smooth_var = 15
        smooth_var_och = 10
        #
        # MSD, ignoring wrong charge histograms, for instances
        #
        # St.: 850, pmt 3, GPS: 1343601187
        if self.uberdict[sd]['offLine_qpk'][pmt - 1] < 500:
            self.uberdict[sd][pmt]['nofit_och'] = True
            return
        # finding fitting range
        if cch == 'cch':
            #
            if len(fitrange) == 0:
                smooth_h = self.rebin(self.uberdict[sd]['cch'], \
                                      smooth_chan=smooth_var)  # COULD pass deri instead of doing it 3x
                deri = self.differentiate(smooth_h=smooth_h)
                try:
                    [xleft, xright], qpk2 = self.find_bounds_cch(deri=deri, \
                                                                 smooth_h=smooth_h, \
                                                                 data_h=self.uberdict[sd]['cch'], pmtno=pmt)
                except ValueError:
                    #print("Unable to find CCH bounds for ", sd, pmt)
                    self.uberdict[sd][pmt]['nofit_cch'] = True
                    return
                except TypeError:
                    #print("Unable to find CCH bounds for ", sd, pmt)
                    self.uberdict[sd][pmt]['nofit_cch'] = True
                    return
            else:
                xleft, xright = np.abs(x_arr - fitrange[0]).argmin(), \
                                np.abs(x_arr - fitrange[1]).argmin()
            h_type = '_cch'
        elif cch == 'och':
            h_type = '_och'
            if len(fitrange) == 0:
                smooth_h = self.rebin(self.uberdict[sd]['och'], \
                                      smooth_chan=smooth_var_och)  # COULD pass deri instead of doing it 3x
                deri = self.differentiate(smooth_h=smooth_h)
                try:
                    [xleft, xright], _, qpk2 \
                        = self.find_bounds_och(deri, \
                                               smooth_h=smooth_h, \
                                               data_h=self.uberdict[sd]['och'], pmtno=pmt)
                except ValueError:
                    #print("Unable to find OCH bounds for ", sd, pmt)
                    self.uberdict[sd][pmt]['nofit_och'] = True
                    return
                except TypeError:
                    #print("Unable to find OCH bounds for ", sd, pmt)
                    self.uberdict[sd][pmt]['nofit_och'] = True
                    return
            else:
                xleft, xright = np.abs(x_arr - fitrange[0]).argmin(), \
                                np.abs(x_arr - fitrange[1]).argmin()
        elif cch == 'valley':
            h_type = '_valley'
            if len(fitrange) == 0:
                smooth_h = self.rebin(self.uberdict[sd]['och'], \
                                      smooth_chan=smooth_var_och)  # COULD pass deri instead of doing it 3x
                deri = self.differentiate(smooth_h=smooth_h)
                try:
                    _, [xleft, xright], _ = self.find_bounds_och(deri, smooth_h, \
                                                                 data_h=self.uberdict[sd]['och'], pmtno=pmt)
                except ValueError:
                    #print("Unable to find OCH bounds for ", sd, pmt)
                    self.uberdict[sd][pmt]['nofit_och'] = True
                    return
                except TypeError:
                    #print("Unable to find OCH bounds for ", sd, pmt)
                    self.uberdict[sd][pmt]['nofit_och'] = True
                    return
            else:
                xleft, xright = np.abs(x_arr - fitrange[0]).argmin(), \
                                np.abs(x_arr - fitrange[1]).argmin()
        #
        # fitting
        if cch == 'cch':
            #
            # Fitting Pol1 + logNormal
            #pars0 = [-3e-3, 7.5, 1700, 7.22, 0.4]
            pars0 = [4.6, 0.0017, 1045, 7.32, 0.32]
            #bounds0 = [[-0.5, 1., 400, 0, 1e-3], [0.0, 20, 2000, 100, 1]]
            bounds0 = [[-0.1, -1e-3, 100, 0, 1e-3], [10.0, 1.0, 3200, 20, 1]]
            tmpLeft = np.where(x_arr == 604)[0][0]
            tmpRight = np.where(x_arr == 4016)[0][0]
            #
            try:
                pars_poLogNorm, pcov_poLogNorm = curve_fit(self.lin_lognorm, x_arr[tmpLeft:tmpRight], \
                                       y_arr[tmpLeft:tmpRight], sigma=y_err[tmpLeft:tmpRight], absolute_sigma=True, \
                                       p0=pars0, bounds=bounds0)
            except:
                print('l. 628: cch fit unsuccessful')
                return
            #
            x_max_poLogNorm = np.exp(pars_poLogNorm[3])
            mu, sigma = pars_poLogNorm[3], pars_poLogNorm[4]
            emu, esigma = np.sqrt(pcov_poLogNorm[3, 3]), np.sqrt(pcov_poLogNorm[4, 4])
            x_max_err_poLogNorm = np.exp(mu) * emu
            ndf_poLogNorm = xright - xleft - 5
            max_count_poLogNorm = self.lin_lognorm(x_max_poLogNorm, *pars_poLogNorm)
            max_count_err_poLogNorm = 0
            funvals = self.lin_lognorm(x_arr[xleft:xright], *pars_poLogNorm)
            chisq_arr = np.divide(np.subtract(y_arr[xleft:xright], funvals), \
                                  y_err[xleft:xright]) ** 2
            p_val_poLogNorm = stats.chi2.sf(np.sum(chisq_arr), ndf_poLogNorm)
            chisq_per_dof_poLogNorm = np.sum(chisq_arr) / ndf_poLogNorm
            #
            pars0 = [-5.5e-05, 1.9e-01, -20.0]
            bounds0 = [[-1, 0, -np.inf], [0, np.inf, 0]]
            #
            # Fitting pol2
            try:
                quad_pars, pcov = curve_fit(self.quadrat, x_arr[xleft:xright], \
                                            y_arr[xleft:xright], sigma=y_err[xleft:xright], \
                                            absolute_sigma=True, p0=pars0, bounds=bounds0)
            except:
                #print(cch + ' fit unsuccessful in l. 575.')
                self.uberdict[sd][pmt]['nofit_och'] = True
                return
            #
            # calculating the parameters, qpk, qpk_err, max_count
            a, b, c = quad_pars
            ea, eb, ec = np.sqrt(np.diag(pcov))
            x_max = -b / (2 * a)
            x_max_err = np.sqrt(eb ** 2 / 4 / a ** 2 + b ** 2 * ea ** 2 / 4 / a ** 4 \
                                - b / 2 / a ** 3 * pcov[0, 1])
            #
            max_count = self.quadrat(x_max, *quad_pars)
            max_count_err = self.quadrat_err(x_max, cov=pcov)

            funvals = self.quadrat(x_arr[xleft:xright], *quad_pars)
            chisq_arr = np.divide(np.subtract(y_arr[xleft:xright], funvals), \
                                  y_err[xleft:xright]) ** 2
            p_val = stats.chi2.sf(np.sum(chisq_arr), (xright - xleft - 3))
            chisq_per_dof = np.sum(chisq_arr) / (xright - xleft - 3)
        #
        else:
            if cch == 'valley':
                pars0 = [5.5e-05, -1.9e-01, 20.0]
                bounds0 = [[0, -100, 0], [1, 0, np.inf]]
            else:
                pars0 = [-5.5e-05, 1.9e-01, -20.0]
                bounds0 = [[-1, 0, -np.inf], [0, np.inf, 0]]
            #
            try:
                quad_pars, pcov = curve_fit(self.quadrat, x_arr[xleft:xright], \
                                            y_arr[xleft:xright], sigma=y_err[xleft:xright], \
                                            absolute_sigma=True, p0=pars0, bounds=bounds0)
            except:
                #print(cch + ' fit unsuccessful in l. 575.')
                self.uberdict[sd][pmt]['nofit_och'] = True
                return
            # calculating the parameters, qpk, qpk_err, max_count
            a, b, c = quad_pars
            ea, eb, ec = np.sqrt(np.diag(pcov))
            x_max = -b / (2 * a)
            x_max_err = np.sqrt(eb ** 2 / 4 / a ** 2 + b ** 2 * ea ** 2 / 4 / a ** 4 \
                                - b / 2 / a ** 3 * pcov[0, 1])
            #
            max_count = self.quadrat(x_max, *quad_pars)
            max_count_err = self.quadrat_err(x_max, cov=pcov)

            funvals = self.quadrat(x_arr[xleft:xright], *quad_pars)
            chisq_arr = np.divide(np.subtract(y_arr[xleft:xright], funvals), \
                                  y_err[xleft:xright]) ** 2
            p_val = stats.chi2.sf(np.sum(chisq_arr), (xright - xleft - 3))
            chisq_per_dof = np.sum(chisq_arr) / (xright - xleft - 3)
        #
        # MSD
        if x_max_err / x_max > 0.1:
            #print(cch + ' fit unsuccessful; Higher error:', \
                #x_max, x_max_err, x_max_err / x_max)
            self.uberdict[sd][pmt]['nofit_och'] = True
            return                                                 
        #
        # storing the results
        self.uberdict[sd][pmt]['fitpars' + h_type] = quad_pars
        self.uberdict[sd][pmt]['cov_matrix' + h_type] = pcov
        self.uberdict[sd][pmt]['fitpars_err' + h_type] = [ea, eb, ec]
        self.uberdict[sd][pmt]['fitrange' + h_type] = [xleft, xright]
        self.uberdict[sd][pmt]['ndf' + h_type] = xright - xleft - 3
        self.uberdict[sd][pmt]['qpk' + h_type] = [x_max, x_max_err]
        self.uberdict[sd][pmt]['qpk_val' + h_type] = [max_count, max_count_err]
        self.uberdict[sd][pmt]['p_val' + h_type] = p_val
        self.uberdict[sd][pmt]['chisq/ndf' + h_type] = chisq_per_dof
        #
        if h_type == '_cch':
            self.uberdict[sd][pmt]['ndf_poLogNorm' + h_type] = ndf_poLogNorm
            self.uberdict[sd][pmt]['cqpk' + h_type + '_poLogNorm'] = [x_max_poLogNorm, x_max_err_poLogNorm]
            self.uberdict[sd][pmt]['cqpk_val' + h_type + '_poLogNorm'] = [max_count_poLogNorm, max_count_err_poLogNorm]
            self.uberdict[sd][pmt]['p_val' + h_type + '_poLogNorm'] = p_val_poLogNorm
            self.uberdict[sd][pmt]['chisq/ndf' + h_type + '_poLogNorm'] = chisq_per_dof_poLogNorm
        #        
        if (h_type == '_cch' or h_type == '_och') and len(fitrange) == 0:
            self.uberdict[sd][pmt]['qpk2' + h_type] = qpk2
        #
        return

    def fit_one(self, sd, pmt, h_type='cch', fitrange=[]):
        if not (pmt in list(self.uberdict[sd].keys())):
            self.uberdict[sd][pmt] = {}
        # if self.uberdict[sd]['cch'][30, pmt] < 1.:
        #
        # MSD
        if np.sum(self.uberdict[sd]['cch'][:, pmt]) < 1.:
            self.uberdict[sd][pmt]['nofit_' + h_type] = True
            return
        self.uberdict[sd][pmt]['nofit_' + h_type] = False
        if h_type == 'cch':
            if len(fitrange) > 0:
                self.quad_fit(sd, pmt, self.uberdict[sd][h_type][:, 0], \
                              self.uberdict[sd][h_type][:, pmt], \
                              self.uberdict[sd][h_type + '_err'][:, pmt],
                              fitrange=fitrange)
            else:
                self.quad_fit(sd, pmt, self.uberdict[sd][h_type][:, 0], \
                              self.uberdict[sd][h_type][:, pmt], \
                              self.uberdict[sd][h_type + '_err'][:, pmt])
            # check if fit worked
            #if self.uberdict[sd][pmt]['nofit_cch']:
                #print('cch fit unsuccessful for ' + str(sd) + ', ' + str(pmt) \
                #      + ', ' + str(self.uberdict[sd]['gps_times'][0]))
                #print(':((((((((((((')
            return
        else:
            if len(fitrange) > 0 and h_type == 'och':
                self.quad_fit(sd, pmt, self.uberdict[sd][h_type][:, 0], \
                              self.uberdict[sd][h_type][:, pmt], \
                              self.uberdict[sd][h_type + '_err'][:, pmt], cch='och', \
                              fitrange=fitrange)
            elif len(fitrange) > 0 and h_type == 'valley':
                self.quad_fit(sd, pmt, self.uberdict[sd]['och'][:, 0], \
                              self.uberdict[sd]['och'][:, pmt], \
                              self.uberdict[sd]['och' + '_err'][:, pmt], cch='valley', \
                              fitrange=fitrange)
            else:
                self.quad_fit(sd, pmt, self.uberdict[sd][h_type][:, 0], \
                              self.uberdict[sd][h_type][:, pmt], \
                              self.uberdict[sd][h_type + '_err'][:, pmt], cch='och')
                self.quad_fit(sd, pmt, self.uberdict[sd][h_type][:, 0], \
                              self.uberdict[sd][h_type][:, pmt], \
                              self.uberdict[sd][h_type + '_err'][:, pmt], cch='valley')

            if self.uberdict[sd][pmt]['nofit_och']:
                #print('och fit unsuccessful for ' + str(sd) + ', ' + str(pmt) \
                #      + ', ' + str(self.uberdict[sd]['gps_times'][0]))
                #print(':((((((((((((')
                #
                return
            h = self.uberdict[sd][pmt]['qpk_val_och'][0]
            v = self.uberdict[sd][pmt]['qpk_val_valley'][0]
            eh, ev = self.uberdict[sd][pmt]['qpk_val_och'][1], \
                     self.uberdict[sd][pmt]['qpk_val_valley'][1]
            self.uberdict[sd][pmt]['v/h'] = [v / h, np.sqrt(ev ** 2 / h ** 2 + \
                                                            v ** 2 * eh ** 2 / h ** 4)]
            return

    def fit_all(self):
        for sd_id in self.uberdict.keys():
            for pmtno in range(1, 4):
                self.fit_one(sd_id, pmtno)
                self.fit_one(sd_id, pmtno, h_type='och')
        return

    def repair_fit(self, sd_id, pmtno, fitrange, h_type='cch'):
        if h_type == 'cch':
            self.fit_one(sd_id, pmtno, fitrange=fitrange)
        else:
            self.fit_one(sd_id, pmtno, h_type=h_type, fitrange=fitrange)
        return

    def plot_fit(self, sd, pmt, zoomed=False, save=[], fit_errs=False,
                 qpk2=False, contours=False):
        plt.rcParams['legend.fontsize'] = 12  # 16
        if self.uberdict[sd][pmt]['nofit_cch']:
            return
        if zoomed:
            h, w = 5, 9
        else:
            h, w = 8, 6
        qpk_cch = np.round(self.uberdict[sd][pmt]['qpk_cch'])
        qpk_och = np.round(self.uberdict[sd][pmt]['qpk_och'])
        maxno_cch = self.uberdict[sd][pmt]['qpk_val_cch'][0]
        max_err_cch = self.uberdict[sd][pmt]['qpk_val_cch'][1]
        maxno_och = self.uberdict[sd][pmt]['qpk_val_och'][0]
        max_err_och = self.uberdict[sd][pmt]['qpk_val_och'][1]

        xcch, xoch = self.uberdict[sd]['cch'][:, 0], \
                     self.uberdict[sd]['och'][:, 0]
        ycch, yoch = self.uberdict[sd]['cch'][:, pmt] / maxno_cch, \
                     self.uberdict[sd]['och'][:, pmt] / maxno_och
        yerrcch = self.uberdict[sd]['cch_err'][:, pmt] / maxno_cch
        yerroch = self.uberdict[sd]['och_err'][:, pmt] / maxno_och
        ############ plottin' ####################################    
        fig, ax = plt.subplots(figsize=(h, w))

        ax.errorbar(xcch, ycch, yerr=yerrcch, c='dimgrey', zorder=2, \
                    linestyle='', capsize=2)  # , label='coincidence charge histogram')
        if zoomed:
            ax.set_xlim(-5, 3000)
        ax2 = ax.twinx()
        ax2.errorbar(xoch, yoch, yerr=yerroch, c='darkgrey', zorder=1, \
                     linestyle='', capsize=2)  # , label='charge histogram')
        ax.set_ylim(ax2.get_ylim())
        ax3 = ax.twinx()
        ax3.set_ylim(ax.get_ylim())

        x_pts = xcch[self.uberdict[sd][pmt]['fitrange_cch'][0]: \
                     self.uberdict[sd][pmt]['fitrange_cch'][1]]
        y_pts = self.quadrat(x_pts, \
                             *self.uberdict[sd][pmt]['fitpars_cch']) / maxno_cch
        if zoomed:
            ax3.plot(x_pts, y_pts, c='orange', zorder=3, linewidth=2, \
                     label='CCH fit:$Q_\mathrm{pk}$=%i$\pm$%i' \
                           % (qpk_cch[0], qpk_cch[1]))
        else:
            ax3.plot(x_pts, y_pts, c='orange', zorder=3, linewidth=2, \
                     label='CCH fit:\n$Q_\mathrm{pk}$=%i$\pm$%i\np-value=%.3f\n'\
                           +'$\chi^2$/ndf=%.2f\nndf=%i' \
                           % (*qpk_cch, self.uberdict[sd][pmt]['p_val_cch'], \
                              self.uberdict[sd][pmt]['chisq/ndf_cch'], \
                              self.uberdict[sd][pmt]['ndf_cch']) )
            # ax3.plot(x_pts,y_pts,c='orange',zorder=3,linewidth=2,\
            #     label='CCH fit:\n$dN/dQ_{max}$=%.2f$\pm$%.2f\np-value=%.3f\n$\chi^2$/ndf=%.2f\nndf=%i'\
            #     %(maxno_cch, max_err_cch,self.uberdict[sd][pmt]['p_val_cch'],\
            #       self.uberdict[sd][pmt]['chisq/ndf_cch'],\
            #           self.uberdict[sd][pmt]['ndf_cch']))
            if contours:
                pcov = self.uberdict[sd][pmt]['cov_matrix_cch']
                err = 2 * self.quadrat_err(x_pts, pcov) / maxno_cch
                ax3.plot(x_pts, np.add(y_pts, err), c='orange', \
                         zorder=1, linewidth=1, linestyle='-')  # (0, (1, 10)))
                ax3.plot(x_pts, np.subtract(y_pts, err), c='orange', \
                         zorder=3, linewidth=.5, linestyle='-')  # (0, (1, 10)))
        if fit_errs:
            max_y = self.uberdict[sd][pmt]['qpk_val_cch'][0]
            ax3.errorbar(self.uberdict[sd][pmt]['qpk_cch'][0], \
                         1, xerr=self.uberdict[sd][pmt]['qpk_cch'][1], \
                         yerr=self.uberdict[sd][pmt]['qpk_val_cch'][1] / max_y, \
                         c='orange', marker='', capsize=2)
        if qpk2:
            q, q_err = self.uberdict[sd][pmt]['qpk2_cch']
            ax3.axvspan(q - q_err, q + q_err, alpha=0.2, color='orange', \
                        label='Rebinned $Q_\mathrm{pk, CCH}$=%i$\pm$%i' \
                              % (int(q), int(np.round(q_err))))
        ax3.set_yticks([])
        ax.set_xlabel(r'charge [FADC counts]')
        ax.set_ylabel(r'normalised $\mathrm{d}N/\mathrm{d}Q$', c='k')

        x_pts = xoch[self.uberdict[sd][pmt]['fitrange_och'][0]: \
                     self.uberdict[sd][pmt]['fitrange_och'][1]]
        y_pts = self.quadrat(x_pts, \
                             *self.uberdict[sd][pmt]['fitpars_och']) / maxno_och
        x_pts_val = xoch[self.uberdict[sd][pmt]['fitrange_valley'][0]: \
                         self.uberdict[sd][pmt]['fitrange_valley'][1]]
        y_pts_val = self.quadrat(x_pts_val, \
                                 *self.uberdict[sd][pmt]['fitpars_valley']) / maxno_och
        if zoomed:
            ax2.plot(x_pts, y_pts, c='k', zorder=4, linewidth=2, \
                     label='CH valley fit')
            ax2.plot(x_pts, y_pts, c='r', zorder=4, linewidth=2, \
                     label='CH fit:$Q_\mathrm{pk}$=%i$\pm$%i' \
                           % (qpk_och[0], qpk_och[1]))
        else:
            ax2.plot(x_pts, y_pts, c='r', zorder=4, linewidth=2, \
                     label='CH fit:\n$Q_\mathrm{pk}$=%i$\pm$%i\np-value=%.3f\n$\chi^2$/ndf=%.2f\nndf=%i' \
                           % (*qpk_och, \
                              self.uberdict[sd][pmt]['p_val_och'], \
                              self.uberdict[sd][pmt]['chisq/ndf_och'], \
                              self.uberdict[sd][pmt]['ndf_och']))
            # ax2.plot(x_pts,y_pts,c='r',zorder=4,linewidth=2,\
            #       label='CH fit:\n$dN/dQ_{max}$=%.2f$\pm$%.2f\np-value=%.3f\n$\chi^2$/ndf=%.2f\nndf=%i'\
            #               %(maxno_och, max_err_och,\
            #                 self.uberdict[sd][pmt]['p_val_och'],\
            #                 self.uberdict[sd][pmt]['chisq/ndf_och'],\
            #                     self.uberdict[sd][pmt]['ndf_och']))
            ax2.plot(x_pts_val, y_pts_val, c='k', zorder=4, linewidth=2, \
                     label='CH valley fit:\np-value=%.3f\n$\chi^2$/ndf=%.2f\nndf=%i' \
                           % (self.uberdict[sd][pmt]['p_val_valley'], \
                              self.uberdict[sd][pmt]['chisq/ndf_valley'], \
                              self.uberdict[sd][pmt]['ndf_valley']) + \
                           '\n\n$\Delta=$%.2f$\pm$%.2f%%\n$v/h=$%.3f$\pm$%.3f%%' \
                           % (self.uberdict[sd][pmt]['Delta'][0] * 100, \
                              self.uberdict[sd][pmt]['Delta'][1] * 100, \
                              *self.uberdict[sd][pmt]['v/h']))

        if fit_errs:
            max_y = self.uberdict[sd][pmt]['qpk_val_och'][0]
            ax2.errorbar(self.uberdict[sd][pmt]['qpk_och'][0], \
                         1, xerr=self.uberdict[sd][pmt]['qpk_och'][1], \
                         yerr=self.uberdict[sd][pmt]['qpk_val_och'][1] / max_y, \
                         c='r', marker='', capsize=2)
            ax2.errorbar(self.uberdict[sd][pmt]['qpk_valley'][0], \
                         self.uberdict[sd][pmt]['qpk_val_valley'][0] / max_y, \
                         xerr=self.uberdict[sd][pmt]['qpk_valley'][1], \
                         yerr=self.uberdict[sd][pmt]['qpk_val_valley'][1] / max_y, \
                         c='k', marker='', capsize=2)
        if qpk2:
            q, q_err = self.uberdict[sd][pmt]['qpk2_och']
            ax3.axvspan(q - q_err, q + q_err, alpha=0.2, color='red', \
                        label='Rebinned $Q_\mathrm{pk, CH}$=%i$\pm$%i' \
                              % (int(q), int(np.round(q_err))))
        #
        a, b = 0.0, 1.5
        # a,b=0.5,1.1
        c, d = 500, 4000
        ax2.set_ylim(a, b)
        # ax2.set_xlim(c,d)
        ax.set_ylim(a, b)
        # ax.set_xlim(c,d)
        ax3.set_ylim(a, b)
        # ax3.set_xlim(c,d)

        ax2.axes.yaxis.set_ticklabels([])
        ax3.axes.yaxis.set_ticklabels([])
        ax2.axes.yaxis.set_ticks([])
        ax3.axes.yaxis.set_ticks([])

        # legend
        lines, labels = ax3.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines3, labels3 = ax.get_legend_handles_labels()
        ax2.legend(lines + lines2 + lines3, labels + labels2 + labels3, loc=0)
        # ax.grid()
        # if 'Delta' in self.uberdict[sd][pmt].keys():
        #     plt.title('SD %s, PMT %i, $\Delta$=%.1f$\pm$%.1f%%, v/h=%.2f$\pm$%.2f'\
        #         %(sd, pmt,*(100*np.array(self.uberdict[sd][pmt]['Delta'])),\
        #           *self.uberdict[sd][pmt]['v/h']))
        # else:
        plt.title('SD %s, PMT %i, time %i' % (sd, pmt, \
                                              np.round(self.uberdict[sd]['gps_times'][0])))

        if len(save) > 0:
            if zoomed:
                self.plot_files['fit_zoomed'].savefig(bbox_inches='tight')
            else:
                save[0].savefig(bbox_inches='tight')
        else:
            plt.show()
        plt.close('all')
        return

    def save_fits(self, sd_pmt, fit_errs=False, qpk2=False, \
                  filename='fit.pdf'):
        # self.plot_files['fit_zoomed'] = PdfPages(self.dirs['plot']+\
        #                                          'fit_zoomed.pdf')
        self.plot_files['fit'] = PdfPages(self.dirs['plot'] + filename)
        if sd_pmt == 'all':
            for sd in sorted(self.uberdict.keys()):
                for pmt in range(1, 4):
                    if not self.uberdict[sd][pmt]['nofit_cch']:
                        self.plot_fit(sd, pmt, save=True, fit_errs=fit_errs, \
                                      qpk2=qpk2)
        else:
            for i in sd_pmt:
                sd, pmt = i
                if sd in self.uberdict.keys():
                    if not self.uberdict[sd][pmt]['nofit_cch']:
                        self.plot_fit(sd, pmt, save=True, fit_errs=fit_errs, \
                                      qpk2=qpk2)
        self.plot_files['fit'].close()
        # self.plot_files['fit_zoomed'].close()
        return

    def plot_Delta_outlier(self, sd, pmt, vh, vhErr, qpkOffline, qpkErrOffline):
        xoch = np.copy(self.uberdict[sd]['och'][:, 0])
        yoch = np.copy(self.uberdict[sd]['och'][:, pmt])
        ycch = np.copy(self.uberdict[sd]['cch'][:, pmt])
        #
        outRootFile = ROOT.TFile('results/rootFiles2/fittedHisto_delta_' + \
                                 str(int(self.uberdict[sd]['gps_times'][0])) + '_' + \
                                 sd + '_' + str(pmt) + '.root', 'recreate')
        tree = ROOT.TTree("T", "")
        #
        gpsTime = array('d', [0.])
        stid = array('d', [0.])
        pmtid = array('d', [0.])
        qpk = array('d', [0.])
        qpkErr = array('d', [0.])
        qpkOff = array('d', [0.])
        qpkErrOff= array('d', [0.])
        cqpk = array('d', [0.])
        cqpkErr = array('d', [0.])
        vhVal = array('d', [0.])
        vhValErr = array('d', [0.])
        sgnl = array('d', [0.])
        ldfr = array('d', [0.])
        dist = array('d', [0.])
        ener = array('d', [0.])
        enerErr = array('d', [0.])
        angle = array('d', [0.])
        #
        cqpk_offline = array('d', [0.])
        cqpkErr_offline = array('d', [0.])
        cqpk_poLogNomr = array('d', [0.])
        cqpkErr_poLogNomr = array('d', [0.])
        #
        tree.Branch("gpsTime", gpsTime, 'gpsTime/D')
        tree.Branch("sdId", stid, 'sdId/D')
        tree.Branch("pmtId", pmtid, 'pmtId/D')
        tree.Branch("qpk", qpk, 'qpk/D')
        tree.Branch("qpkErr", qpkErr, 'qpkErr/D')
        tree.Branch("qpkOffLine", qpkOff, 'qpkOffLine/D')
        tree.Branch("qpkErrOffLine", qpkErrOff, 'qpkErrOffLine/D')
        tree.Branch("cqpkPy", cqpk, 'cqpkPy/D')
        tree.Branch("cqpkErrPy", cqpkErr, 'cqpkErrPy/D')
        tree.Branch("vh", vhVal, 'vh/D')
        tree.Branch("vhErr", vhValErr, 'vhErr/D')
        tree.Branch("signal", sgnl, 'signal/D')
        tree.Branch("LDF", ldfr, 'LDF/D')
        tree.Branch("spDist", dist, 'spDist/D')
        tree.Branch("Energy", ener, 'Energy/D')
        tree.Branch("EnergyErr", enerErr, 'EnergyErr/D')
        tree.Branch("Zenith", angle, 'Zenith/D')
        #
        tree.Branch("cqpkPoLogNorm", cqpk_poLogNomr, 'cqpkPoLogNorm/D')
        tree.Branch("cqpkErrPoLogNorm", cqpkErr_poLogNomr, 'cqpkErrPoLogNorm/D')
        tree.Branch("cqpkOffLine", cqpk_offline, 'cqpkOffLine/D')
        tree.Branch("cqpkErrOffLine", cqpkErr_offline, 'cqpkErrOffLine/D')
        #
        gpsTime[0] = 1. * int(self.uberdict[sd]['gps_times'][0])
        stid[0] = int(sd)
        pmtid[0] = pmt
        qpk[0] = self.uberdict[sd][pmt]['qpk_och'][0]
        qpkErr[0] = self.uberdict[sd][pmt]['qpk_och'][1]
        qpkOff[0] = qpkOffline
        qpkErrOff[0] = qpkErrOffline
        cqpk[0] = self.uberdict[sd][pmt]['qpk_cch'][0]
        cqpkErr[0] = self.uberdict[sd][pmt]['qpk_cch'][1]
        vhVal[0] = vh
        vhValErr[0] = vhErr
        sgnl[0] = self.uberdict[sd]['signal'][pmt-1]
        ldfr[0] = self.uberdict[sd]['ldfR']
        dist[0] = self.uberdict[sd]['spDist']
        ener[0] = self.uberdict[sd]['Energy']
        enerErr[0] = self.uberdict[sd]['EnergyErr']
        angle[0] = self.uberdict[sd]['Zenith']
        #
        cqpk_poLogNomr[0] = self.uberdict[sd][pmt]['cqpk_cch_poLogNorm'][0]
        cqpkErr_poLogNomr[0] = self.uberdict[sd][pmt]['cqpk_cch_poLogNorm'][1]
        cqpk_offline[0] = self.uberdict[sd]['offLine_cqpk'][pmt-1]
        cqpkErr_offline[0] = self.uberdict[sd]['offLine_cqpkErr'][pmt-1]
        #
        rch = ROOT.TH1D('rch', 'rH St. '+str(sd)+' PMT '+str(pmt), len(xoch)-1, xoch)
        cch = ROOT.TH1D('cch', 'cH St. '+str(sd)+' PMT '+str(pmt), len(xoch)-1, xoch)
        #
        for bin_i in range(len(xoch)):
            rch.SetBinContent(bin_i, yoch[bin_i])
            cch.SetBinContent(bin_i, ycch[bin_i])
        #
        tree.Fill()
        tree.Write("", ROOT.TObject.kOverwrite)
        outRootFile.Write()
        outRootFile.Close()
        #
        return

    def calculate_delta(self, sd, pmt):
        if self.uberdict[sd][pmt]['nofit_cch'] or \
                self.uberdict[sd][pmt]['nofit_och']:
            return np.inf, np.inf
        # 
        # MSD, for comparison with Qpk from Offline
        tmpQpk = self.uberdict[sd]['offLine_qpk'][pmt-1]
        tmpQpkErr = self.uberdict[sd]['offLine_qpkErr'][pmt-1]
        #
        Delta = (self.uberdict[sd][pmt]['qpk_cch'][0] \
                 - self.uberdict[sd][pmt]['qpk_och'][0]) \
                / self.uberdict[sd][pmt]['qpk_och'][0]
        Delta_err = np.sqrt( \
            self.uberdict[sd][pmt]['qpk_cch'][1] ** 2 / \
            self.uberdict[sd][pmt]['qpk_och'][0] ** 2 + \
            self.uberdict[sd][pmt]['qpk_cch'][0] ** 2 * \
            self.uberdict[sd][pmt]['qpk_och'][1] ** 2 / \
            self.uberdict[sd][pmt]['qpk_och'][0] ** 4 )
        #
        # MSD
        self.plot_Delta_outlier(sd, pmt, self.uberdict[sd][pmt]['v/h'][0], \
                                self.uberdict[sd][pmt]['v/h'][1], \
                                tmpQpk, tmpQpkErr)
        #
        self.uberdict[sd][pmt]['Delta'] = [Delta, Delta_err]
        return self.uberdict[sd][pmt]['Delta']

    def calculate_all_deltas(self):
        num = np.zeros(3)
        den = np.zeros(3)
        d_mean = np.zeros(3)
        d_err = np.zeros(3)
        for sd_id in self.uberdict.keys():
            for pmtno in range(1, 4):
                D, De = self.calculate_delta(sd_id, pmtno)
                if D != np.inf:
                    num[pmtno - 1] += D / De ** 2
                    den[pmtno - 1] += 1 / De ** 2
        for i in range(3):
            if num[i] > 0. and den[i] > 0.:
                d_mean[i] = num[i] / den[i]
                d_err[i] = np.sqrt(1. / den[i])
            else:
                d_mean[i] = np.inf
                d_err[i] = np.inf
        for i in range(1, 4):
            self.delta_mean[i] = [d_mean[i - 1], d_err[i - 1]]
        return

    def extract_deltas(self, pmtno):
        Delta, Delta_err, sd_id = [], [], []
        for sd in self.uberdict.keys():
            if not self.uberdict[sd][pmtno]['nofit_cch']:
                sd_id.append(sd)
                Delta.append(self.uberdict[sd][pmtno]['Delta'][0])
                Delta_err.append(self.uberdict[sd][pmtno]['Delta'][1])
        return np.array([sd_id, Delta, Delta_err]).T

    def plot_delta_histo(self, save=False, fit_distrib=False, \
                         only_peak_fit=False):
        fig, axs = plt.subplots(1, 3, figsize=(21, 7), sharex=True, sharey=True)

        for i in range(3):
            buf = self.extract_deltas(i + 1)[:, 1] * 100
            bin_num = int(np.round(3 * (1 + 3.322 * np.log10(len(buf)))))
            n, bins, patches = axs[i].hist(buf, bin_num, \
                                           density=False, facecolor=self.cols[i], alpha=0.75)
            if fit_distrib:
                (mu, sigma) = stats.norm.fit(buf)
                if only_peak_fit:
                    # first empty bins from the peak
                    peak = np.argmax(n)
                    zer = np.where(n == 0)[0]
                    # first zero bin from right
                    dif = zer - peak
                    try:
                        r = bins[np.min(dif[np.where(dif > 0)]) + peak]
                    except ValueError:
                        r = bins[-1]
                    try:
                        l = bins[np.max(dif[np.where(dif < 0)]) + peak]
                    except ValueError:
                        l = bins[0]
                    buf = buf[np.where(buf < r)]
                    buf = buf[np.where(buf > l)]
                    cut = [l, r]
                    (mu, sigma) = stats.norm.fit(buf)
                pdf_fitted = stats.norm.pdf(np.linspace(bins[0], bins[-1], 100), \
                                            mu, sigma)
                # Get bin edges
                xh = [0.5 * (bins[r] + bins[r + 1]) for r in range(len(bins) - 1)]
                # Get bin width from this
                binwidth = (max(xh) - min(xh)) / len(bins)
                pdf_fitted = pdf_fitted * (len(buf) * binwidth)
                if only_peak_fit:
                    axs[i].plot(np.linspace(bins[0], bins[-1], 100), pdf_fitted, \
                                'k--', linewidth=2, label='Gaussian fit peak only,\n' + \
                                                          '$\mu$=%.1f, $\sigma$=%.1f\nrange: [%.1f,%.1f]' \
                                                          % (mu, sigma, *cut))
                else:
                    axs[i].plot(np.linspace(bins[0], bins[-1], 100), pdf_fitted, \
                                'k--', linewidth=2, label= \
                                    'Gaussian fit, $\mu$=%.1f, $\sigma$=%.1f' % (mu, sigma))

            axs[i].set_title("PMT %i, $\Delta_\mathrm{mean}$=%.1f$\pm$%.1f%%" \
                             % (i + 1, self.delta_mean[i + 1][0] * 100, self.delta_mean[i + 1][1] * 100))
            axs[i].set_xlabel('$\Delta$ [%]')
            axs[i].set_ylabel('#')
            axs[i].legend()
            axs[i].xaxis.set_tick_params(labelbottom=True)
            axs[i].yaxis.set_tick_params(labelbottom=True)

        plt.tight_layout()
        if save:
            ciel1 = PdfPages(self.dirs['plot'] + 'Deltas_histo.pdf')
            ciel1.savefig(bbox_inches='tight')
            ciel1.close()
        else:
            plt.show()

    def plot_delta_bands(self, save=False):
        stations = np.array(sorted(list(self.uberdict.keys())))
        wid = int(0.28 * len(stations))
        split_pg = int(np.round(wid / 20))
        wid = np.floor(wid / split_pg)
        split_st = np.array_split(stations, split_pg)
        if save:
            ciel1 = PdfPages(self.dirs['plot'] + 'Deltas_band.pdf')
        for page in range(split_pg):
            stations = split_st[page]
            fig, axs = plt.subplots(3, 1, figsize=(wid, 12), \
                                    sharex=True, sharey=True)
            for i in range(3):
                delta_res = self.extract_deltas(i + 1)
                have_delta = (delta_res[:, 0]).astype(int)
                y_ax = np.array(delta_res)[:, 1] * 100
                y_er = np.array(delta_res)[:, 2] * 100
                x_tick = stations
                x_ax = np.linspace(0, 1, len(x_tick))

                delta, no_delta = [], []
                for n, j in enumerate(stations):
                    ind = np.where(have_delta == j)[0]
                    if len(ind) == 0:
                        no_delta.append([x_ax[n], \
                                         self.delta_mean[i + 1][0] * 100, 0])
                    else:
                        delta.append([x_ax[n], y_ax[ind[0]], y_er[ind[0]]])
                delta, no_delta = np.array(delta), np.array(no_delta)
                axs[i].errorbar(delta[:, 0], delta[:, 1], yerr=delta[:, 2], \
                                marker='_', capsize=2, linestyle='', c=self.cols[i])
                axs[i].errorbar(no_delta[:, 0], no_delta[:, 1], c='darkred', \
                                marker='x', linestyle='')

                axs[i].axhspan(100 * (self.delta_mean[i + 1][0] - \
                                      self.delta_mean[i + 1][1]), \
                               100 * (self.delta_mean[i + 1][0] + \
                                      self.delta_mean[i + 1][1]), alpha=0.2, color='red', \
                               label='$\Delta_\mathrm{mean}$=%.1f$\pm$%.1f%%' \
                                     % (self.delta_mean[i + 1][0] * 100, self.delta_mean[i + 1][1] * 100))
                axs[i].set_title("PMT %i" % (i + 1))
                axs[i].set_ylabel('$\Delta$ [%]')
                axs[i].set_xlabel('SD ID')
                axs[i].set_xticks(x_ax)
                axs[i].set_xticklabels(x_tick, rotation=90)
                axs[i].legend(loc=1, fancybox=True, framealpha=0.5)
                axs[i].xaxis.set_tick_params(labelbottom=True)
                axs[i].yaxis.set_tick_params(labelbottom=True)

            plt.tight_layout()
            if save:
                ciel1.savefig(bbox_inches='tight')
            else:
                plt.show()

        if save:
            ciel1.close()
        return

    def drop_fit_data(self, sd, pmt):
        self.uberdict[sd][pmt] = {'nofit_cch': True}
        return

    def cut_start(self, sd, pmt, Q, histo_type='och'):
        index = np.argmin(np.abs(self.uberdict[sd][histo_type][:, 0] - Q)) + 1
        self.uberdict[sd][histo_type][:index, pmt] = 0
        self.uberdict[sd][histo_type + '_err'][:index, pmt] = 0
        return

    def plot_och_fit(self, sd, pmt=[1, 2, 3], err=True, save=[], title='default', \
                     fit_errs=False, fit_contour=[]):
        if not (sd in self.uberdict.keys()):
            return
        plt.rcParams.update({'font.size': 12})
        h, w = 15, 9
        histo_type = 'och'
        self.grey_c = ['dimgrey', 'darkgrey', 'lightgrey']

        fig, ax = plt.subplots(figsize=(h, w))
        for n, i in enumerate(pmt):
            if self.uberdict[sd][i]['nofit_och']:
                continue
            qpk_och = np.round(self.uberdict[sd][i]['qpk_och'])
            maxno_och = self.uberdict[sd][i]['qpk_val_och'][0]
            xoch = self.uberdict[sd]['och'][:, 0]
            yoch = self.uberdict[sd]['och'][:, i] / maxno_och
            yerroch = self.uberdict[sd]['och_err'][:, i] / maxno_och
            if (self.uberdict[sd][histo_type][10, i] == 0):
                ax.errorbar(self.uberdict[sd][histo_type][:, 0], \
                            self.uberdict[sd][histo_type][:, i], \
                            marker='.', label='PMT %i' % i, \
                            c=self.grey_c[n], linestyle='')
                continue

            ax.errorbar(xoch, yoch, yerr=yerroch, marker='_', \
                        label='PMT %i, v/h=%.2f$\pm$%.2f' % \
                              (i, *self.uberdict[sd][i]['v/h']), c=self.grey_c[n], \
                        linestyle='', capsize=2)

            ax.set_xlabel('charge [FADC counts]')
            ax.set_ylabel('normalised dN/dQ', c='k')

            x_pts = xoch[self.uberdict[sd][i]['fitrange_och'][0]: \
                         self.uberdict[sd][i]['fitrange_och'][1]]
            y_pts = self.quadrat(x_pts, \
                                 *self.uberdict[sd][i]['fitpars_och']) / maxno_och
            x_pts_val = xoch[self.uberdict[sd][i]['fitrange_valley'][0]: \
                             self.uberdict[sd][i]['fitrange_valley'][1]]
            y_pts_val = self.quadrat(x_pts_val, \
                                     *self.uberdict[sd][i]['fitpars_valley']) / maxno_och

            ax.plot(x_pts_val, y_pts_val, c=self.cols[n], zorder=4, linewidth=0.1, \
                    linestyle='-', \
                    label=r'$\bf{PMT\:' + str(i) + '\:valley\:fit:}$' + '\np-value=%.3f\n$\chi^2$/ndf=%.2f\nndf=%i' \
                          % (self.uberdict[sd][i]['p_val_valley'], \
                             self.uberdict[sd][i]['chisq/ndf_valley'], \
                             self.uberdict[sd][i]['ndf_valley']))
            ax.plot(x_pts, y_pts, c=self.cols[n], zorder=4, linewidth=0.1, \
                    linestyle='-', \
                    label=r'$\bf{PMT\:' + str(
                        i) + '\:hump\:fit:}$' + '\n$Q_\mathrm{pk}$=%i$\pm$%i\np-value=%.3f\n$\chi^2$/ndf=%.2f\nndf=%i' \
                          % (*qpk_och, self.uberdict[sd][i]['p_val_och'], \
                             self.uberdict[sd][i]['chisq/ndf_och'], \
                             self.uberdict[sd][i]['ndf_och']))

            if i in fit_contour:
                fit_name, xpt_name = ['_och', '_valley'], [x_pts, x_pts_val]
                for m, j in enumerate(fit_name):
                    pars = self.uberdict[sd][i]['fitpars' + j]
                    pcov = self.uberdict[sd][i]['cov_matrix' + j]
                    err = self.quadrat_err(xpt_name[m], pcov) / maxno_och
                    y_pts = self.quadrat(xpt_name[m], *pars) / maxno_och
                    ax.plot(xpt_name[m], np.add(y_pts, err), c=self.cols[n], \
                            zorder=6, linewidth=.1, linestyle='-')  # (0, (1, 10)))
                    ax.plot(xpt_name[m], np.subtract(y_pts, err), c=self.cols[n], \
                            zorder=6, linewidth=.1, linestyle='-')  # (0, (1, 10)))
            if fit_errs:
                max_y = self.uberdict[sd][i]['qpk_val_och'][0]
                ax.errorbar(self.uberdict[sd][i]['qpk_och'][0], \
                            1, xerr=self.uberdict[sd][i]['qpk_och'][1], \
                            yerr=self.uberdict[sd][i]['qpk_val_och'][1] / max_y, \
                            c=self.cols[n], marker='', capsize=2, zorder=6)
                ax.errorbar(self.uberdict[sd][i]['qpk_valley'][0], \
                            self.uberdict[sd][i]['qpk_val_valley'][0] / max_y, \
                            xerr=self.uberdict[sd][i]['qpk_valley'][1], \
                            yerr=self.uberdict[sd][i]['qpk_val_valley'][1] / max_y, \
                            c=self.cols[n], marker='', capsize=2, zorder=6)

        ax.legend()
        ax.grid()
        # if 'Delta' in self.uberdict[sd][pmt].keys():
        #     plt.title('SD %s, PMT %i, $\Delta$=%.1f$\pm$%.1f%%, v/h=%.2f$\pm$%.2f'\
        #         %(sd, pmt,*(100*np.array(self.uberdict[sd][pmt]['Delta'])),\
        #           *self.uberdict[sd][pmt]['v/h']))
        # else:
        plt.title('SD %i, OCH fit' % (sd))

        if save:
            self.plot_files['fit'].savefig(bbox_inches='tight')
        else:
            plt.show()
        plt.close('all')
        return

    def save_och_fits(self, sd_ids, fit_errs=False, fit_contour=[]):
        # self.plot_files['fit_zoomed'] = PdfPages(self.dirs['plot']+\
        #                                          'fit_zoomed.pdf')
        self.plot_files['fit'] = PdfPages(self.dirs['plot'] + 'fit_och.pdf')
        for sd in sd_ids:
            if sd in self.uberdict.keys():
                self.plot_och_fit(sd, save=True, fit_errs=fit_errs, \
                                  fit_contour=fit_contour)
        self.plot_files['fit'].close()
        return
    #
    # MSD
    #    
    def comparison_offline(self):
        success_fit = 0
        for sd_id in self.uberdict.keys():
            for pmt_i in range(1, 4):
                try:
                    ifNotCQpkfit = self.uberdict[sd_id][pmt_i]['nofit_cch']
                    isVEMhisto = self.uberdict[sd_id]['vemFromHistogram'][pmt_i - 1]
                    if not ifNotCQpkfit and not isVEMhisto:
                        success_fit += 1
                        self.plot_comparison_offline(sd_id, pmt_i)
                except:
                    print('No key for', pmt_i, sd_id, \
                          self.uberdict[sd_id]['gps_times'][0])
        return success_fit

    def plot_comparison_offline(self, sd, pmt):
        plt.rcParams.update({'font.size': 24})
        h = 15
        w = 9
        name = str(self.uberdict[sd]['gps_times'][0]) + sd + '_' + str(pmt) + '.pdf'
        ciel1 = PdfPages(self.dirs['plot'] + 'noQpkYesCQpk_' + name)
        self.grey_c = ['dimgrey', 'darkgrey', 'lightgrey']
        self.grey_cc = ['lightgrey', 'dimgrey', 'darkgrey']
        fig, ax = plt.subplots(figsize=(h, w))
        xoch = np.copy(self.uberdict[sd]['och'][:, 0])
        yoch = np.copy(self.uberdict[sd]['och'][:, pmt])
        yerroch = np.copy(self.uberdict[sd]['och_err'][:, pmt])

        qpk = self.uberdict[sd][pmt]['qpk_och'][0]
        errQpk = self.uberdict[sd][pmt]['qpk_och'][1]
        cqpk = self.uberdict[sd][pmt]['qpk_cch'][0]
        errCqpk = self.uberdict[sd][pmt]['qpk_cch'][1]

        Ycqpk = self.uberdict[sd][pmt]['qpk_val_cch'][0]
        errYcqpk = np.sqrt(Ycqpk)
        Yqpk = self.uberdict[sd][pmt]['qpk_val_och'][0]
        errYqpk = np.sqrt(Yqpk)

        tmp0 = 0.
        tmp1 = 0.
        tmp2 = 0.
        norm_Yerr = []
        for i in range(len(yerroch)):
            tmp0 = (yerroch[i] / Yqpk) ** 2
            tmp1 = (yoch[i] / (Yqpk ** 2) * errYqpk) ** 2
            tmp2 = np.sqrt(tmp0 + tmp1)
            norm_Yerr.append(tmp2)
        yoch /= Yqpk
        yerroch = np.array(norm_Yerr)
        ax.errorbar(xoch, yoch, yerr=yerroch, marker='_', \
                    label='St %s, PMT %i' % (sd, pmt), \
                    c=self.grey_c[pmt - 1], linestyle='', capsize=2)
        ycch = np.copy(self.uberdict[sd]['cch'][:, pmt])
        yerrcch = np.copy(self.uberdict[sd]['cch_err'][:, pmt])
        norm_Yerr = []
        for i in range(len(yerrcch)):
            tmp0 = (yerrcch[i] / Ycqpk) ** 2
            tmp1 = (ycch[i] / (Ycqpk ** 2) * errYcqpk) ** 2
            tmp2 = np.sqrt(tmp0 + tmp1)
            norm_Yerr.append(tmp2)
        ycch /= Ycqpk
        yerrcch = np.array(norm_Yerr)
        ax.errorbar(xoch, ycch, yerr=yerrcch, marker='_', \
                    label='CQpk: %.2f $\pm$ %.2f' % (qpk, errQpk), \
                    c=self.grey_cc[pmt - 1], linestyle='', capsize=2)
        ax.set_xlabel('charge [FADC counts]')
        ax.set_ylabel('Counts/FADC', c='k')

        plt.legend()
        plt.plot()
        ciel1.savefig(bbox_inches='tight')
        ciel1.close()
        plt.close(fig)
