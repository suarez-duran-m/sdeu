"""
Created on Tue Mar  29 13:39:15 2022

@author: katarina
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
from matplotlib.backends.backend_pdf import PdfPages
import os
from scipy.optimize import curve_fit
# from astropy.time import Time
import glob
from scipy import stats
from matplotlib.transforms import Bbox
# import csv, pickle

class COMBINE_RESULTS:
    def __init__(self, obj, DIR_OUT):
        self.obj = obj
        self.results = {}
        self.cols=['r', 'orange', 'royalblue']
        self.dirs={'out':DIR_OUT, 'plot':DIR_OUT}
        self.signal = {}
        self.stationsWithCoinci = []
        self.nDeltasPerPMT = 1
        
        if not os.path.isdir(self.dirs['out']):
            os.makedirs(self.dirs['out'])
        if not os.path.isdir(self.dirs['plot']):
            os.makedirs(self.dirs['plot'])
        return
    
    def print_objects(self):
        for n,i in enumerate(self.obj):
            print('Analysis object '+str(n+1)+':', i)
        return
    
    def average_objects(self):
        # creating container for all SDs
        for obj in self.obj:
            for sd_id in obj.uberdict.keys():
                self.results[sd_id] = {}
                for pmtno in range(1, 4):
                    self.results[sd_id][pmtno] = {'deltas':[], 'delta_errs':[],\
                                              'v/hs':[], 'v/h_errs':[]}
        # taking results from all the CCH objects
        aveSignal = []
        aveCSignal = []
        signalPmt = [[], [], []]
        cSignalPmt = [[], [], []]
        for obj in self.obj:
            for sd_id in obj.uberdict.keys():
                self.stationsWithCoinci.append( sd_id )
                sgnlPmt = [0., 0., 0.]
                cSgnlPmt = [0., 0., 0.]
                nPmtsWithSgnl = 0
                tmpAveSgnl = 0.
                tmpAveCSgnl = 0.
                for pmtno in range(1, 4):
                    if 'Delta' not in obj.uberdict[sd_id][pmtno]:
                        continue
                    if obj.uberdict[sd_id][pmtno]['nofit_cch']:
                        continue
                    if obj.uberdict[sd_id][pmtno]['nofit_och']:
                        continue
                    self.results[sd_id][pmtno]['deltas'].append(\
                            obj.uberdict[sd_id][pmtno]['Delta'][0])
                    self.results[sd_id][pmtno]['delta_errs'].append(\
                            obj.uberdict[sd_id][pmtno]['Delta'][1])
                    self.results[sd_id][pmtno]['v/hs'].append(\
                            obj.uberdict[sd_id][pmtno]['v/h'][0])
                    self.results[sd_id][pmtno]['v/h_errs'].append(\
                            obj.uberdict[sd_id][pmtno]['v/h'][1])
                    #
                    # MSD, signals
                    if obj.uberdict[sd_id][pmtno]['Delta'][0] < 0:
                        continue
                    if obj.uberdict[sd_id][pmtno]['Delta'][0] > 0.1:
                        continue
                    if obj.uberdict[sd_id][pmtno]['v/h'][0] < 0.65:
                        continue
                    if obj.uberdict[sd_id][pmtno]['v/h'][0] > 0.8:
                        continue
                    crrSgnl = obj.uberdict[sd_id]['signal'][pmtno-1]
                    if crrSgnl < 1:
                        continue
                    #
                    # MSD, from VEM_Qpk to VEM_cQpk
                    crrCSgnl = ( 1. \
                            + obj.uberdict[sd_id][pmtno]['Delta'][0] )\
                            / obj.uberdict[sd_id][pmtno]['qpk_cch'][0]
                    crrCSgnl *= crrSgnl*obj.uberdict[sd_id]['offLine_qpk'][pmtno-1]
                    #
                    tmpAveSgnl += crrSgnl
                    sgnlPmt[pmtno-1] = crrSgnl
                    tmpAveCSgnl += crrCSgnl
                    cSgnlPmt[pmtno-1] = crrCSgnl
                    nPmtsWithSgnl += 1
                #
                if nPmtsWithSgnl < 3:
                    continue
                tmpAveSgnl /= nPmtsWithSgnl
                tmpAveCSgnl /= nPmtsWithSgnl
                aveSignal.append( tmpAveSgnl )
                aveCSignal.append( tmpAveCSgnl )
                for pmt_i in range(3):
                    if sgnlPmt[pmt_i] == 0:
                        signalPmt[pmt_i].append(0)
                    if cSgnlPmt[pmt_i] == 0:
                        cSignalPmt[pmt_i].append(0)
                    signalPmt[pmt_i].append( sgnlPmt[pmt_i] )
                    cSignalPmt[pmt_i].append( cSgnlPmt[pmt_i] )
        #
        self.signal['aveSgnl'] = aveSignal
        self.signal['aveSgnlCoinci'] = aveCSignal
        self.signal['sgnlPmt'] = signalPmt
        self.signal['sgnlPmtCoinci'] = cSignalPmt
        #
        # calculating average Delta and h/v
        for sd_id in self.results.keys():

            deltasPMT1 = len(self.results[sd_id][1]['deltas'])
            deltasPMT2 = len(self.results[sd_id][2]['deltas'])
            deltasPMT3 = len(self.results[sd_id][3]['deltas'])
            #
            # MSD, Filter per number of deltas for each PMT
            if deltasPMT1 < self.nDeltasPerPMT or deltasPMT2 < self.nDeltasPerPMT \
                    or deltasPMT3 < self.nDeltasPerPMT:
                continue
            #
            if deltasPMT1 != deltasPMT2 or deltasPMT1 != deltasPMT3:
                continue
            #
            '''
            if np.average(self.results[sd_id][1]['deltas']) < 0 or \
                    np.average(self.results[sd_id][2]['deltas']) < 0 or \
                    np.average(self.results[sd_id][3]['deltas']) < 0:
                continue
            '''
            #
            for pmtno in range(1, 4):
                num, den = 0, 0
                if len(self.results[sd_id][pmtno]['deltas']) == 0:
                    continue
                for n, i in enumerate(self.results[sd_id][pmtno]['deltas']):
                    num += i / self.results[sd_id][pmtno]['delta_errs'][n]**2
                    den += 1 / self.results[sd_id][pmtno]['delta_errs'][n]**2
                self.results[sd_id][pmtno]['Delta'] = [num/den, np.sqrt(1/den)]
                #
                num, den = 0, 0
                for n,i in enumerate(self.results[sd_id][pmtno]['v/hs']):
                    num += i/self.results[sd_id][pmtno]['v/h_errs'][n]**2
                    den += 1/self.results[sd_id][pmtno]['v/h_errs'][n]**2
                self.results[sd_id][pmtno]['v/h'] = [num/den, np.sqrt(1/den)]
                #
        #
        return

    def plot_EventsPerStation(self):
        plt.rcParams.update({'font.size': 24})
        h = 15
        w = 9
        #
        ciel1 = PdfPages('results/distEventPerStation.pdf')
        self.grey_c = ['dimgrey', 'darkgrey', 'lightgrey']
        self.grey_cc = ['lightgrey', 'dimgrey', 'darkgrey']
        fig, ax = plt.subplots(figsize=(h, w))
        #
        counter = np.zeros(2000)
        for sd_id in self.stationsWithCoinci:
            counter[int(sd_id)] += 1
        bins = np.linspace(1, 1000, 1000)
        #
        # Plotting
        n, bins, patches = ax.hist(counter, bins)
        #
        # Fitting
        def Gauss(x, a, mean, std):
            return a * np.exp(-(x - mean)**2/(2 * std**2))
        #x = np.linspace(1, bins[np.where(n > 0)[-1][-1]]+1, 100)
        popt, pcov = curve_fit(Gauss, bins[:-1], n)
        #
        ax.plot(bins, Gauss(bins, popt[0], popt[1], popt[2]), label = \
                'mean: %.2f $\pm$ %.2f\nStd: %.2f $\pm$ %.2f' \
                % (popt[1], popt[2]/np.sqrt(len(n)), popt[2], pcov[2][2]))
        ax.set_xlabel('Events/Station')
        ax.set_ylabel('Counts [au]')
        ax.set_xlim(0, bins[np.where(n > 0)[-1][-1]]+1)
        ax.legend()
        #
        ciel1.savefig(bbox_inches='tight')
        ciel1.close()
        plt.close()

        return
    
    def count_deltas(self, vh_thresh=2, Delta_thresh=[-1, 1]):
        c=[0, 0, 0]
        for sd in self.results.keys():
            for pmt in range(1, 4):
                if len(self.results[sd][pmt]['deltas']) > 0:
                    if (self.results[sd][pmt]['v/h'][0] < vh_thresh) \
                            and self.results[sd][pmt]['Delta'][0] > Delta_thresh[0] \
                            and self.results[sd][pmt]['Delta'][0] < Delta_thresh[1]:
                        c[pmt-1]+=1
        return c
    
    def make_matrix(self, pmtno, exclude):
        to_plot = []
        for sd_id in self.results.keys():
            if sd_id in exclude:
                continue
            if 'Delta' not in self.results[sd_id][pmtno].keys():
                continue
            if 'v/h' not in self.results[sd_id][pmtno].keys():
                continue
            if self.results[sd_id][pmtno]['v/h'][0] < 0.8: #\
                # and self.results[sd_id][pmtno]['Delta'][0] < 0.1:
                to_plot.append([sd_id, *self.results[sd_id][pmtno]['Delta'], \
                                *self.results[sd_id][pmtno]['v/h']])
                '''
                if self.results[sd_id][pmtno]['v/h'][0] < 0.8 \
                        and self.results[sd_id][pmtno]['Delta'][0] > 0.04 \
                        and self.results[sd_id][pmtno]['Delta'][0] < 0.05:
                    
                    print("MSD Delta-Ave:", sd_id, pmtno, \
                          self.results[sd_id][pmtno]['Delta'][0], \
                          self.results[sd_id][pmtno]['v/h'][0] )
                '''
        #
        to_plot = np.array(to_plot)
        to_plot_float = to_plot.astype(np.float)
        return to_plot_float
    
    def calculate_deltas_mean(self):
        self.delta_mean = {}
        num, den = np.zeros(3), np.zeros(3)
        for sd_id in self.results.keys():
            for pmtno in range(1, 4):
                if 'Delta' in self.results[sd_id][pmtno].keys():
                    D, De = self.results[sd_id][pmtno]['Delta']
                    if D != np.inf:
                        num[pmtno-1] += D/De**2
                        den[pmtno-1] += 1/De**2
        d_mean, d_err = np.divide(num,den), np.sqrt(1/den)
        for i in range(1, 4):
            self.delta_mean[i] = [d_mean[i-1], d_err[i-1]]
        return
    
    def set_delta_mean(self, dictio):
        self.delta_mean = dictio
        return
    #
    # MSD
    #
    def rebin(self, data_h, cut_off, smooth_chan):
             new_dath = [[],[]]
             for k in range(int(cut_off/smooth_chan)):
                 new_dath[0].append(\
                         np.average(data_h[smooth_chan*k:(smooth_chan*(k+1)), 0])) 
                         #np.sum(data_h[smooth_chan*k:(smooth_chan*(k+1)),0])/smooth_chan)
                 new_dath[1].append(\
                         np.average(data_h[smooth_chan*k:(smooth_chan*(k+1)), 1]))
                         #np.sum(data_h[smooth_chan*k:(smooth_chan*(k+1)),1]))
             return np.array(new_dath).T
    
    def plot_Delta_vs_vh(self, save=False, errs=True, exclude={1:[],2:[],3:[]}):
        plt.rcParams.update({'font.size': 18})
        fig, axs = plt.subplots(1, 3, figsize=(15,5), sharex=True, sharey=True)
        for i in range(3):
            to_plot = self.make_matrix(i+1, exclude[i+1])
            '''
            # 
            # MSD, for fatbins
            fatBin = 12
            tmp_to_plot = np.array([to_plot[:, 3], to_plot[:, 1]])
            tmp_to_plot = tmp_to_plot.T
            new_to_plot = self.rebin(tmp_to_plot, len(tmp_to_plot), fatBin)
            ys = new_to_plot[:, 1]*100
            xs = new_to_plot[:, 0]
            #            
            tmp_to_plot = np.array([to_plot[:, 4], to_plot[:, 2]])
            tmp_to_plot = tmp_to_plot.T
            new_to_plot = self.rebin(tmp_to_plot, len(tmp_to_plot), fatBin)
            y_ers = new_to_plot[:, 1]*100
            x_ers = new_to_plot[:, 0]
            '''
            ys = to_plot[:, 1]*100
            y_ers = to_plot[:, 2]*100
            xs = to_plot[:, 3]
            x_ers = to_plot[:, 4]
            bins = np.where(xs < 0.8)
            '''
            ys = ys[bins] # for profile
            y_ers = y_ers[bins]
            xs = xs[bins]
            x_ers = x_ers[bins]
            '''
            # rsq = stats.pearsonr(xs, ys)
            # print('116: ', rsq[0])
            xmean, ymean = np.mean(xs)*np.ones(len(xs)),\
                            np.mean(ys)*np.ones(len(ys))
            pear = np.sum((xs-xmean) * (ys-ymean)) / np.sqrt(np.sum((xs-xmean)**2))\
                                               / np.sqrt(np.sum((ys-ymean)**2))
            err_r = (1-pear**2) / np.sqrt(len(xs))
            #
            if errs:
                axs[i].errorbar(xs, ys, yerr=y_ers, xerr=x_ers,capsize=2,\
                            c=self.cols[i], marker='', linestyle='',
                            label="Pearson's r=%.2f$\pm$%.2f"\
                                             %(pear, err_r))
                # axs[i].plot(xs, ys, c='k', marker='o', linestyle='')
            else:
                axs[i].errorbar(xs, ys,\
                            c=self.cols[i], marker='o', linestyle='')
            
            # axs[i].set_title("PMT %i, Pearson's r=%.2f$\pm$%.2f"\
            #                  %(i+1,pear,err_r))
            axs[i].set_ylabel('$\Delta$ [%]')
            axs[i].set_xlabel(r'valley/hump')
            axs[i].grid()
            axs[i].legend(fontsize=18)
            axs[i].xaxis.set_tick_params(labelbottom=True)
            axs[i].yaxis.set_tick_params(labelbottom=True)
        plt.tight_layout()
        if save:
            if errs:
                #ciel1= PdfPages(self.dirs['plot']+'Delta_vs_vh_profile.pdf')
                ciel1= PdfPages(self.dirs['plot']+'Delta_vs_vh.pdf')
            else:
                ciel1= PdfPages(self.dirs['plot']+'Delta_vs_vh_noerrs.pdf')
            for i in range(3):
                extent = axs[i].get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
                ciel1.savefig(bbox_inches=extent)
            ciel1.close()
        else:
            plt.show()
        return
    
    def plot_delta_bands(self, save=False, y_lim=[]):
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({'font.size': 18})
        stations = np.array(sorted(list(self.results.keys())))
        wid = int(0.28*len(stations))
        split_pg = int(np.round(wid/10))
        wid = np.floor(wid/split_pg)
        split_st = np.array_split(stations,split_pg)
        # finding limits of y axis for all pages the same
        delta_res_all = self.make_matrix(1, [])
        delta_res_all = np.concatenate((delta_res_all, self.make_matrix(2, [])))
        delta_res_all = np.concatenate((delta_res_all, self.make_matrix(3, [])))
        indmax = np.argmax(delta_res_all[:, 1] + delta_res_all[:, 2])
        indmin = np.argmin(delta_res_all[:, 1] - delta_res_all[:, 2])
        ymax, ymin = (delta_res_all[:, 1] + delta_res_all[:, 2])[indmax]*100,\
            (delta_res_all[:, 1] - delta_res_all[:, 2])[indmin]*100
        Delta_ind = ['1,3', '2', '1,3']
        
        if save:
            ciel1= PdfPages(self.dirs['plot']+'Deltas_band.pdf')
        for page in range(split_pg):
            stations=split_st[page]
            fig, axs = plt.subplots(3, 1, figsize=(wid, 12),\
                                    sharex=True, sharey=True)
            for i in range(3):
                delta_res = self.make_matrix(i+1, [])
                have_delta = (delta_res[:, 0]).astype(int)
                y_ax = np.array(delta_res)[:, 1]*100
                y_er = np.array(delta_res)[:, 2]*100
                x_tick = stations
                x_ax = np.linspace(0, 1, len(x_tick))
               
                delta, no_delta = [], []
                for n, j in enumerate(stations):
                    ind = np.where(have_delta==float(j))[0]                    
                    if len(ind)==0:
                        no_delta.append([x_ax[n], self.delta_mean[i+1][0]*100, 0])
                    else:
                        delta.append([x_ax[n], y_ax[ind[0]], y_er[ind[0]]])
                delta, no_delta = np.array(delta), np.array(no_delta)
                axs[i].errorbar(delta[:, 0], delta[:, 1], yerr=delta[:, 2], \
                                c=self.cols[i], marker='_', capsize=2, linestyle='')
                if no_delta.shape[0]>0:
                    axs[i].errorbar(no_delta[:, 0], no_delta[:, 1], c='darkred',\
                                    marker='x', linestyle='')
                #
                # MSD
                num = 0.
                den = 0.
                for jj in range(len(delta[:, 1])):
                    num += delta[jj, 1]/(delta[jj, 2])**2
                    den += 1./(delta[jj, 2])**2
                tmpDeltaMean = np.divide(num, den)
                tmpDeltaMeanErr = np.sqrt(1./den)
                axs[i].axhspan(tmpDeltaMean\
                        - tmpDeltaMeanErr,\
                        tmpDeltaMean\
                        + tmpDeltaMeanErr,\
                        alpha=0.2, color='red',\
                        label='$<\Delta$'+'$_{%s}$'%Delta_ind[i]\
                        +'$>$=%.2f$\pm$%.2f%%'\
                        %(tmpDeltaMean,\
                        tmpDeltaMeanErr))
                '''
                axs[i].axhspan(100*(self.delta_mean[i+1][0]-\
                                    self.delta_mean[i+1][1]),\
                               100*(self.delta_mean[i+1][0]+\
                                self.delta_mean[i+1][1]), alpha=0.2, color='red',\
                        label='$<\Delta$'+'$_{%s}$'%Delta_ind[i]+'$>$=%.2f$\pm$%.2f%%'\
                       %(self.delta_mean[i+1][0]*100, \
                         self.delta_mean[i+1][1]*100))
                '''
                axs[i].set_title("PMT"+" %i"%(i+1))
                axs[i].set_ylabel('$\Delta$ [\%]')
                axs[i].set_xlabel('SD ID')
                axs[i].set_xticks(x_ax)
                axs[i].set_ylim(ymin-1,ymax+1)
                axs[i].set_xticklabels(x_tick, rotation=90)
                axs[i].legend(loc=1, fancybox=True, framealpha=0.5)
                axs[i].xaxis.set_tick_params(labelbottom=True)
                axs[i].yaxis.set_tick_params(labelbottom=True)
            if len(y_lim)>0:
                plt.ylim((y_lim[0],y_lim[-1]))
            plt.tight_layout()
            if save:
                ciel1.savefig(bbox_inches='tight')
            else:
                plt.show()
            
        if save:   
            ciel1.close()
        return
    
    def plot_delta_histo(self, save=False, fit_distrib=False,\
                         only_peak_fit=False):
        plt.rcParams.update({'font.size': 18})
        fig, axs = plt.subplots(1, 3, figsize=(21,7), sharex=True, sharey=True)
        #
        for pmt_i in range(3):
            buf = self.make_matrix(pmt_i+1, [])
            buf = buf[:, 1]*100
            #
            bin_num = int(np.round( 3*(1 + 3.322*np.log10(len(buf))) ) )
            n, bins, patches = axs[pmt_i].hist(buf, bin_num,\
                             density=False, facecolor=self.cols[pmt_i], alpha=0.75)
            print("MSD bins ", pmt_i, bins)
            print("MSD n ", pmt_i, n)
            if fit_distrib:
                #(mu, sigma) = stats.norm.fit(buf)
                if only_peak_fit:
                    # MSD
                    #                    
                    #first empty bins from the peak
                    peak = np.argmax(n)
                    zer = np.where(n==0)[0]
                    # first zero bin from right
                    dif = zer - peak
                    try:
                        r = bins[np.min(dif[np.where(dif>0)]) + peak]
                    except ValueError:
                        r = bins[-1]
                    try:
                        l = bins[np.max(dif[np.where(dif<0)]) + peak]
                    except ValueError:
                        l = bins[0]
                    '''
                    if pmt_i == 0 or pmt_i == 2:
                        r = 7.0
                    
                    if pmt_i == 1:
                        l = 2.5
                    '''
                    buf = buf[np.where(buf<r)]
                    buf = buf[np.where(buf>l)]
                    cut = [l, r]
                    (shape, loc, scale) = stats.lognorm.fit(buf)
                    #(mu, sigma) = stats.norm.fit(buf)
                    mode = scale * np.exp(-shape**2) + loc
                    #
                    sumShape = 0.
                    sumScale = 0.
                    sumLoc1 = 0.
                    sumLoc2 = 0.
                    for i in buf:
                        sumShape += np.log((i-loc)/scale) * np.log((i-loc)/scale)
                        sumScale += 1. + np.log((i-loc)/scale)
                        sumLoc1 += 1. / (i-loc)**2
                        sumLoc2 += np.log((i-loc)/scale) / (i-loc)**2 - 1./(i-loc)**2
                    sumShape *= 3./(shape**4)
                    varShape = len(buf)/shape**2 - sumShape
                    varShape = 1. / varShape
                    #
                    sumScale += -len(buf)/scale**2 - sumScale/(shape**2 * scale**2)
                    varScale = 1. / sumScale
                    #
                    varLoc = sumLoc1 + sumLoc2/shape**2
                    varLoc = 1. / varLoc
                    #
                    ErrMode = np.sqrt(\
                            (np.exp(-1.*shape**2)*varScale)**2 + varLoc**2 \
                            + (2.*shape*scale*np.exp(-1.*shape**2)*varShape)**2 )
                #
                #pdf_fitted = stats.norm.pdf(np.linspace(bins[0], bins[-1], 100),\
                #                           mu, sigma)
                pdf_fitted = stats.lognorm.pdf(np.linspace(bins[0], bins[-1], 100),\
                        shape, loc, scale)
                #Get bin edges
                xh = [0.5 * (bins[r] + bins[r+1]) for r in range(len(bins)-1)]
                #Get bin width from this
                binwidth = (max(xh) - min(xh)) / len(bins)
                pdf_fitted = pdf_fitted * (len(buf) * binwidth)
                if only_peak_fit:
                    axs[pmt_i].plot(np.linspace(bins[0],bins[-1],100), pdf_fitted,\
                            'k--', linewidth=2, label='Gaussian fit,\n'+\
                            #'$\mu$=%.2f$\pm$%.2f\n'%(mu, sigma/len(buf))\
                            #+'$\sigma$=%.2f$\pm$%.2f\n'\
                            #%(sigma, sigma/np.sqrt(2*len(buf)))
                            'LogNorm fit\n$s$=%.1f loc=%.1f scale=%.1f\n'\
                                    %(shape, loc, scale)\
                            +'mode=%.2f $\pm$ %.2f\n'%(mode, ErrMode)\
                            +'range: [%.1f, %.1f]'%(*cut,))
                            #'k--', linewidth=2, label='Gaussian fit peak only,\n'+\
                            #'$\mu$=%.1f, $\sigma$=%.1f\nrange: [%.1f, %.1f]'\
                            #%(mu, sigma, *cut))
                else:
                    axs[pmt_i].plot(np.linspace(bins[0],bins[-1],100), pdf_fitted,\
                            'k--', linewidth=2, label=\
                            'LogNormal fit\nmode=%.2f $\pm$ %.2f'%( mode, ErrMode ) )
                            #'Gaussian fit\n$\mu$=%.1f$\pm$%.1f\n'
                            #%(mu,sigma/len(buf))\
                            #+'$\sigma$=%.1f$\pm$%.1f' \
                            #%(sigma,sigma/np.sqrt(2*len(buf))))
            axs[pmt_i].set_title("PMT %i"%(pmt_i+1))
            axs[pmt_i].set_xlabel('$\Delta$ [\%]')
            axs[pmt_i].set_ylabel('\#')
            axs[pmt_i].legend()
            axs[pmt_i].xaxis.set_tick_params(labelbottom=True)
            axs[pmt_i].yaxis.set_tick_params(labelbottom=True)
            print("MSD total entries", pmt_i, np.sum(n))

        plt.tight_layout()
        if save:
            ciel1 = PdfPages(self.dirs['plot']+'Deltas_histo.pdf')
            #ciel1= PdfPages(self.dirs['plot']+'Deltas_histo_profile.pdf')
            ciel1.savefig(bbox_inches='tight')
            ciel1.close()
        else:
            plt.show()
        return

    # MSD
    #                    
    def plot_deltaMix_histo(self, save=False):
        plt.rcParams.update({'font.size': 18})
        fig, axs = plt.subplots(1, 2, figsize=(21,7), sharex=True, sharey=True)

        for pmt_i in range(2):            
            buf_pmt = self.make_matrix(pmt_i+1, [])
            buf_pmt = buf_pmt[:, 1]*100
            #
            # Mixing PMT1 with PMT3
            if pmt_i == 0:
                buf_tmp = self.make_matrix(3, [])
                buf_tmp = buf_tmp[:, 1]*100
                buf_pmt = np.append(buf_pmt, buf_tmp)

            bin_num = int(np.round( 3*(1 + 3.322*np.log10(len(buf_pmt))) ) )
            n, bins, patches = axs[pmt_i].hist(buf_pmt, bin_num, density=False,\
                    facecolor=self.cols[pmt_i], alpha=0.75)
            #
            #first empty bins from the peak
            peak = np.argmax(n)
            zer = np.where(n==0)[0]
            # first zero bin from right
            dif = zer - peak
            try:
                r = bins[np.min(dif[np.where(dif>0)])+peak]
            except ValueError:
                r = bins[-1]
            try:
                l = bins[np.max(dif[np.where(dif<0)])+peak]
            except ValueError:
                l = bins[0]
            '''
            if pmt_i == 1:
                l = 2.5
            '''
            buf_pmt = buf_pmt[np.where(buf_pmt < r)]
            buf_pmt = buf_pmt[np.where(buf_pmt > l)]
            cut = [l, r]
            (shape, loc, scale) = stats.lognorm.fit(buf_pmt)
            mode = scale * np.exp(-shape**2) + loc
            sumShape = 0.
            sumScale = 0.
            sumLoc1 = 0.
            sumLoc2 = 0.
            for i in buf_pmt:
                sumShape += np.log( (i-loc)/scale ) * np.log( (i-loc)/scale )
                sumScale += 1. + np.log( (i-loc)/scale )
                sumLoc1 += 1. / (i-loc)**2
                sumLoc2 += np.log( (i-loc)/scale ) / (i-loc)**2 - 1./(i-loc)**2
                sumShape *= 3./(shape**4)
            varShape = len(buf_pmt)/shape**2 - sumShape
            varShape = 1. / varShape
            
            sumScale += -len(buf_pmt)/scale**2 - sumScale/(shape**2 * scale**2)
            varScale = 1. / sumScale
            
            varLoc = sumLoc1 + sumLoc2/shape**2
            varLoc = 1. / varLoc
            
            ErrMode = np.sqrt(\
                    (np.exp(-1.*shape**2)*varScale)**2 + varLoc**2 \
                    + (2.*shape*scale*np.exp(-1.*shape**2)*varShape)**2 )
            #
            pdf_fitted = stats.lognorm.pdf(np.linspace(bins[0], bins[-1], 100),\
                    shape, loc, scale)            
            #Get bin edges
            xh = [0.5 * (bins[r] + bins[r+1]) for r in range(len(bins)-1)]
            #Get bin width from this
            binwidth = (max(xh) - min(xh)) / len(bins)
            pdf_fitted = pdf_fitted * (len(buf_pmt) * binwidth)
            #
            # Plotting
            axs[pmt_i].plot(np.linspace(bins[0],bins[-1],100), pdf_fitted,\
                    'k--', linewidth=2, label='Gaussian fit,\n'+\
                    'LogNorm fit\n$s$=%.1f loc=%.1f scale=%.1f\n'\
                    %(shape, loc, scale)\
                    +'mode=%.2f $\pm$ %.2f\n'%(mode, ErrMode)\
                    +'range: [%.1f, %.1f]'%(*cut,))
            if pmt_i == 0:
                axs[pmt_i].set_title("PMT %i + PMT %i"%(1, 3))
            else:
                axs[pmt_i].set_title("PMT %i"%(pmt_i+1))
            axs[pmt_i].set_xlabel('$\Delta$ [\%]')
            axs[pmt_i].set_ylabel('\#')
            axs[pmt_i].legend()
            axs[pmt_i].xaxis.set_tick_params(labelbottom=True)
            axs[pmt_i].yaxis.set_tick_params(labelbottom=True)

        plt.tight_layout()
        if save:
            ciel1 = PdfPages(self.dirs['plot']+'DeltasMix_histo.pdf')
            ciel1.savefig(bbox_inches='tight')
            ciel1.close()
        else:
            plt.show()
        return

    def plot_logP(self, save=False, log_counts=False, log_bins=30, log=True):
        cols=['orange', 'red', 'black']
        plt.rcParams.update({'font.size': 18})
        fig, axs = plt.subplots(1, 3, figsize=(15,5),sharex=True, sharey=True)
        # getting logP
        self.logP = {'cch':[], 'och':[], 'valley':[]}
        for obj in self.obj:
            for key in obj.uberdict.keys():
                for pmtno in range(1,4):
                    if not obj.uberdict[key][pmtno]['nofit_cch'] and\
                        not obj.uberdict[key][pmtno]['nofit_och']:
                        # self.logP[pmtno].append([\
                        #    np.log10(obj.uberdict[key][pmtno]['p_val_cch']),\
                        #    np.log10(obj.uberdict[key][pmtno]['p_val_och']),\
                        #    np.log10(obj.uberdict[key][pmtno]['p_val_valley'])])
                        self.logP['cch'].append(obj.uberdict[key][pmtno]['p_val_cch'])
                        self.logP['och'].append(obj.uberdict[key][pmtno]['p_val_och'])
                        self.logP['valley'].append(\
                                    obj.uberdict[key][pmtno]['p_val_valley'])
        mins,maxs=[],[]
        for key in self.logP.keys():
            self.logP[key]=np.array(self.logP[key])
            mins.append(np.min(self.logP[key]))
            maxs.append(np.max(self.logP[key]))
        # finding range of histogram
        if log:
            buf_min, buf_max = np.min(mins), np.max(maxs)
            tot_width=np.log10(buf_max)-np.log10(buf_min)
            plot_range=(np.log10(buf_min)-0.1*tot_width, np.log10(buf_max)+0.1*tot_width)
        else:
            plot_range=(-0.01,1.01)
        for i,key in enumerate(self.logP.keys()):
            if log:
                buf = np.log10(self.logP[key])
            else:
                buf = self.logP[key]
            # log_bins =10**(np.arange(np.floor(np.min(lg10)),\
            #             np.ceil(np.max(lg10))+1 ))
            if log_counts:
                n, bins, patches = axs[i].hist(buf,log_bins,log=True,\
                    range=plot_range, density=False, facecolor=cols[i], alpha=0.75)
            else:
                n, bins, patches = axs[i].hist(buf,log_bins,\
                   range=plot_range, density=False, facecolor=cols[i], alpha=0.75)
            # axs[i,j-1].set_title("PMT %i, %s fits"\
            #     %(j,titlee[i]))
            if log:
                axs[i].set_xlabel('log(p-value)')
            else:
                axs[i].set_xlabel('p-value')
            axs[i].set_ylabel('\#')
            axs[i].xaxis.set_tick_params(labelbottom=True)
            axs[i].xaxis.set_ticks([0.0,0.2,0.4,0.6,0.8,1.0])
            axs[i].yaxis.set_tick_params(labelbottom=True)
        plt.tight_layout()
        if save:
            if log:
                ciel1= PdfPages(self.dirs['plot']+'logP.pdf')
            else:
                ciel1= PdfPages(self.dirs['plot']+'p_vals.pdf')
            for i in range(3):
                extent = axs[i].get_tightbbox(fig.canvas.get_renderer()).\
                    transformed(fig.dpi_scale_trans.inverted())
                ciel1.savefig(bbox_inches=extent)
            ciel1.close()
        else:
            plt.show()
        return
    
    def vh_check(self):
        outliers=[]
        pmt1_2, pmt2_3, pmt3_1 = [], [], []
        for key in self.results.keys():
            #
            # MSD
            '''
            deltasPMT1 = len(self.results[key][1]['deltas'])
            deltasPMT2 = len(self.results[key][2]['deltas'])
            deltasPMT3 = len(self.results[key][3]['deltas'])
            #
            if deltasPMT1 < self.nDeltasPerPMT or deltasPMT2 < self.nDeltasPerPMT \
                    or deltasPMT3 < self.nDeltasPerPMT:
                continue
            #
            if deltasPMT1 != deltasPMT2 or deltasPMT1 != deltasPMT3:
                continue
            '''
            #
            '''
            if np.average(self.results[key][1]['deltas']) < 0 or\
                    np.average(self.results[key][2]['deltas']) < 0 or\
                    np.average(self.results[key][3]['deltas']) < 0:
                continue
            '''
            #
            if len(self.results[key][1]['v/hs'])>0 and\
                len(self.results[key][2]['v/hs'])>0:
                pmt1_2.append([key, *self.results[key][1]['v/h'],\
                               *self.results[key][2]['v/h']])
            if len(self.results[key][2]['v/hs'])>0 and\
                len(self.results[key][3]['v/hs'])>0:
                pmt2_3.append([key, *self.results[key][2]['v/h'],\
                               *self.results[key][3]['v/h']])
            if len(self.results[key][3]['v/hs'])>0\
                and len(self.results[key][1]['v/hs'])>0:
                pmt3_1.append([key, *self.results[key][3]['v/h'],\
                               *self.results[key][1]['v/h']])
    
        comp_vh = [pmt1_2, pmt2_3, pmt3_1]
        for n,i in enumerate(comp_vh):
            comp_vh[n] = np.array(i).astype(np.float)

        #comp_vh_float = comp_vh.astype(np.float)
        errs=False
        ax_labels = [[1,2],[2,3],[3,1]] 
        grad=[0.1, 0.11, 0.1]
        fig, axs = plt.subplots(1, 3, figsize=(15,5),sharex=True, sharey=True)
        for i in range(3):
            to_plot = comp_vh[i]
            #
            xs = to_plot[:, 1]
            x_ers = to_plot[:, 2]
            ys = to_plot[:, 3]
            y_ers = to_plot[:, 4]
            #
            if errs:
                axs[i].errorbar(xs, ys, yerr=y_ers, xerr=x_ers,\
                            c='k', marker='', linestyle='')
            else:
                axs[i].errorbar(xs, ys,\
                            c='k', marker='.', linestyle='')
            for n, j in enumerate(xs):
                if np.abs(ys[n]/j-1)>grad[i]:
                    if int(to_plot[:,0][n])==914:
                           continue
                    axs[i].annotate(int(to_plot[:,0][n]), xy =(j,ys[n]),\
                                         xytext =(j,ys[n]))
                    outliers.append(int(to_plot[:,0][n]))
            axs[i].set_ylabel('valley/hump PMT %i'%ax_labels[i][1])
            axs[i].set_xlabel('valley/hump PMT %i'%ax_labels[i][0])
            axs[i].xaxis.set_tick_params(labelbottom=True)
            axs[i].yaxis.set_tick_params(labelbottom=True)
        # plotting y=x
        xs = axs[0].get_xlim()
        ys = axs[0].get_ylim()
        for i in range(3):
            axs[i].plot(xs, ys, c='k', marker='', linestyle='-')
            axs[i].set_xlim(xs[0], 0.88)
            axs[i].set_ylim(ys[0], 0.85)
        plt.tight_layout()
        plt.savefig(self.dirs['plot']+'vh_pmt_check.pdf',bbox_inches='tight')
        plt.show()
        return sorted(list(set(outliers)))
    
    def deltas_check(self):
        pmt1_2, pmt2_3, pmt3_1 = [], [], []
        Dm = self.delta_mean[1][0]*100,self.delta_mean[2][0]*100,\
            self.delta_mean[3][0]*100
        Dm_err = self.delta_mean[1][1]*100,self.delta_mean[2][1]*100,\
            self.delta_mean[3][1]*100
        for key in self.results.keys():
            if len(self.results[key][1]['deltas'])>0 and\
                len(self.results[key][2]['deltas'])>0:
                D1, eD1 = self.results[key][1]['Delta']
                D2, eD2 = self.results[key][2]['Delta']
                # pmt1_2.append([key, D1-Dm1, eD1, D2-Dm2, eD2])
                pmt1_2.append([key, D1, eD1, D2, eD2])
            if len(self.results[key][2]['deltas'])>0 and\
                len(self.results[key][3]['deltas'])>0:
                D2, eD2 = self.results[key][2]['Delta']
                D3, eD3 = self.results[key][3]['Delta']
                # pmt2_3.append([key, D2-Dm2, eD2, D3-Dm3, eD3])
                pmt2_3.append([key, D2, eD2, D3, eD3])
            if len(self.results[key][3]['deltas'])>0\
                and len(self.results[key][1]['deltas'])>0:
                D1, eD1 = self.results[key][1]['Delta']
                D3, eD3 = self.results[key][3]['Delta']
                # pmt3_1.append([key, D3-Dm3, eD3, D1-Dm1, eD1])
                pmt3_1.append([key, D3, eD3, D1, eD1])
    
        comp_D = [pmt1_2, pmt2_3, pmt3_1]
        for n,i in enumerate(comp_D):
            comp_D[n] = np.array(i)
    
        errs=False
        ax_labels = [[1,2],[2,3],[3,1]]
        out_thr = {1:7, 2:11, 3:10}
        d_max=1000
        fig, axs = plt.subplots(1, 3, figsize=(15,5),sharex=True, sharey=True)
        for i in range(3):
            to_plot=comp_D[i]
    
            xs = to_plot[:,1]*100
            x_ers = to_plot[:,2]*100
            ys = to_plot[:,3]*100
            y_ers = to_plot[:,4]*100
            
            if errs:
                axs[i].errorbar(xs, ys, yerr=y_ers, xerr=x_ers,\
                            c='k', marker='', linestyle='')
            else:
                axs[i].errorbar(xs, ys,\
                            c='k', marker='.', linestyle='')
            for n, j in enumerate(xs):
                # diff_x, diff_y = j-Dm[ax_labels[i][0]-1],\
                #                 ys[n]-Dm[ax_labels[i][1]-1]
                # if np.abs(diff_x)>d_max*Dm_err[ax_labels[i][0]-1] or\
                #     np.abs(diff_y)>d_max*Dm_err[ax_labels[i][1]-1]:
                #     axs[i].annotate(int(to_plot[:,0][n]), xy =(j,ys[n]),\
                #                          xytext =(j,ys[n]), fontsize=16)
                if j>out_thr[ax_labels[i][0]] or ys[n]>out_thr[ax_labels[i][1]]:
                    axs[i].annotate(int(to_plot[:,0][n]), xy =(j,ys[n]),\
                                         xytext =(j,ys[n]), fontsize=16)
            # axs[i].set_ylabel('$\Delta-\overline{\Delta}$  PMT %i [%%]'\
            #                   %ax_labels[i][1])
            # axs[i].set_xlabel('$\Delta-\overline{\Delta}$ PMT %i [%%]'\
            #                   %ax_labels[i][0])
            axs[i].set_ylabel('$\Delta$  PMT %i'%ax_labels[i][1]+' $[\%]$')
            axs[i].set_xlabel('$\Delta$ PMT %i'%ax_labels[i][0]+'$[\%]$')
            axs[i].xaxis.set_tick_params(labelbottom=True)
            axs[i].yaxis.set_tick_params(labelbottom=True)
        # plotting y=x
        xs = axs[0].get_xlim()
        ys = axs[0].get_ylim()
        # for i in range(3):
        #     axs[i].plot(xs, ys, c='k', marker='', linestyle='-')
        #     axs[i].set_xlim(xs)
        #     axs[i].set_ylim(ys)
        plt.tight_layout()
        plt.savefig(self.dirs['plot']+'Deltas_pmt_check.pdf',\
                    bbox_inches='tight')
        plt.show()
        return

    def plot_signalsRelDiffPmt(self, save=False):
        colors=['r', 'orange', 'royalblue', 'violet']
        fig, axs = plt.subplots(2, 3, figsize=(19.2, 10.8), sharex=True, sharey='row')
        nBins = 10
        minSgnl = 0.5
        maxSgnl = 1.5
        xBins = [minSgnl + (i*0.125) for i in range(nBins)]
        yDiff = np.zeros( (3, nBins) )
        yErr = np.zeros( (3, nBins) )
        yDiff2 = np.zeros( (3, nBins) )
        cntsBin = np.zeros( (3, nBins) )
        sgnl_i = 0
        tmpSgnlPmt = 0
        for i in range(len(self.signal['aveSgnl'])):
            sgnl_i =  self.signal['aveSgnl'][i]
            if np.log10(sgnl_i) < minSgnl or np.log10(sgnl_i) > maxSgnl:
                continue
            tmpBin = int( (np.log10(sgnl_i) - minSgnl)*10 )
            for pmt_i in range(3):
                if self.signal['sgnlPmt'][pmt_i][i] == 0:
                    continue
                tmpSgnlPmt = self.signal['sgnlPmt'][pmt_i][i] / sgnl_i - 1.
                tmpSgnlPmt *= 100.
                yDiff[pmt_i][tmpBin] += tmpSgnlPmt
                yDiff2[pmt_i][tmpBin] += tmpSgnlPmt**2
                cntsBin[pmt_i][tmpBin] += 1
        yDiffCoinci = np.zeros( (3, nBins) )
        yErrCoinci = np.zeros( (3, nBins) )
        yDiff2Coinci = np.zeros( (3, nBins) )
        cntsBinCoinci = np.zeros( (3, nBins) )
        sgnl_i = 0
        tmpSgnlPmt = 0
        for i in range(len(self.signal['aveSgnlCoinci'])):
            sgnl_i = self.signal['aveSgnlCoinci'][i]
            if np.log10(sgnl_i) < minSgnl or np.log10(sgnl_i) > maxSgnl:
                continue
            tmpBin = int( (np.log10(sgnl_i) - minSgnl)*10 )
            for pmt_i in range(3):
                if self.signal['sgnlPmtCoinci'][pmt_i][i] == 0:
                    continue
                tmpSgnlPmt = self.signal['sgnlPmtCoinci'][pmt_i][i] / sgnl_i - 1.
                tmpSgnlPmt *= 100.
                yDiffCoinci[pmt_i][tmpBin] += tmpSgnlPmt
                yDiff2Coinci[pmt_i][tmpBin] += tmpSgnlPmt**2
                cntsBinCoinci[pmt_i][tmpBin] += 1
        for pmt_i in range(3):
            yDiff[pmt_i] /= cntsBin[pmt_i]
            yErr[pmt_i] = yDiff2[pmt_i] - cntsBin[pmt_i] * yDiff[pmt_i]**2
            yErr[pmt_i] /= cntsBin[pmt_i]
            yErr[pmt_i] = np.sqrt( yErr[pmt_i] ) / np.sqrt( cntsBin[pmt_i] )
            axs[0][pmt_i].errorbar(xBins, yDiff[pmt_i], yerr=yErr[pmt_i], xerr=0.125, \
                    linestyle='', marker='o', capsize=2, color=colors[0], \
                    label="Signal")
            #
            yDiffCoinci[pmt_i] /= cntsBinCoinci[pmt_i]
            yErrCoinci[pmt_i] = yDiff2Coinci[pmt_i] \
                    - cntsBinCoinci[pmt_i] * yDiffCoinci[pmt_i]**2
            yErrCoinci /= cntsBinCoinci[pmt_i]
            yErrCoinci[pmt_i] = np.sqrt( yErrCoinci[pmt_i] ) \
                    / np.sqrt( cntsBinCoinci[pmt_i] )
            axs[0][pmt_i].errorbar(xBins, yDiffCoinci[pmt_i], yerr=yErrCoinci[pmt_i], \
                    xerr=0.125, linestyle='', marker='*', capsize=2, \
                    color=colors[1], label="Signal (with $\Delta$)")
            #
            axs[0][pmt_i].set_title("PMT %i" % (pmt_i+1))
            axs[0][pmt_i].set_ylabel('$S_{%i}/'%(pmt_i+1)+r'\left\langle S_{1,2,3}\right\rangle-1\,\mathrm{[\%]}$')
            axs[0][pmt_i].xaxis.set_tick_params(labelbottom=True)
            axs[0][pmt_i].yaxis.set_tick_params(labelbottom=True)
            axs[0][pmt_i].legend()
            #
            residual = yDiff[pmt_i] - yDiffCoinci[pmt_i]
            residualErr = np.sqrt( yErr[pmt_i]**2 + yErrCoinci[pmt_i]**2 )
            axs[1][pmt_i].errorbar( xBins, residual, yerr=residualErr, \
                    xerr=0.125, linestyle='', marker='*', capsize=2, \
                    color=colors[2], label="Average: %.2f"\
                    %(np.sum(residual)/len(residual)) )
            axs[1][pmt_i].set_xlabel(r'$log_{10}\left(\left\langle S_{1,2,3}\right\rangle/\mathrm{VEM}\right)$')
            axs[1][pmt_i].set_ylabel("Residuals")
            axs[1][pmt_i].xaxis.set_tick_params(labelbottom=True)
            axs[1][pmt_i].yaxis.set_tick_params(labelbottom=True)
            axs[1][pmt_i].legend()
        #
        if save:
            fig.tight_layout()
            ciel1 = PdfPages(self.dirs['plot']+'signalsPmt_vs_average.pdf')
            ciel1.savefig(bbox_inches='tight')
            ciel1.close()
        return
    
    def getPrtlRMS_DiffSgnlPairPMT(self, readKey, keyPmt, minSgnl, maxSgnl, nBins):
        retDiffSgnlPmt2 = np.zeros( (3, nBins) )
        retMeanDiffSgnlPmt = np.zeros( (3, nBins) )
        retCntsPerBin = np.zeros( (3, nBins) )
        for i in range(len(self.signal[readKey])):
            sgnl_i = self.signal[readKey][i]
            if np.log10(sgnl_i) < minSgnl or np.log10(sgnl_i) > maxSgnl:
                continue
            tmpBin = int( (np.log10(sgnl_i) - minSgnl)*10 )
            for pmt_i in range(3):
                tmpSgnlPmt_i = self.signal[keyPmt][pmt_i][i]
                pmt_j = (pmt_i+1) if ( pmt_i == 0 or pmt_i == 1 ) else 0
                tmpSgnlPmt_j = self.signal[keyPmt][pmt_j][i]
                if tmpSgnlPmt_i < 1. and tmpSgnlPmt_j < 1.:
                    continue
                tmpDiffPmt = np.sqrt(2) * (tmpSgnlPmt_i - tmpSgnlPmt_j)
                tmpDiffPmt /= (tmpSgnlPmt_i + tmpSgnlPmt_j)
                #
                retDiffSgnlPmt2[pmt_i][tmpBin] += tmpDiffPmt**2
                retMeanDiffSgnlPmt[pmt_i][tmpBin] += tmpDiffPmt
                retCntsPerBin[pmt_i][tmpBin] += 1
        #
        return retDiffSgnlPmt2, retMeanDiffSgnlPmt, retCntsPerBin

    def plot_signalsResolution(self, save=False):
        colors=['r', 'orange', 'royalblue', 'violet']
        fig, axs = plt.subplots(2, 3, figsize=(19.2, 10.8), sharex=True, sharey='row')
        nBins = 10
        minSgnl = 0.5
        maxSgnl = 1.5
        xBins = [minSgnl + (i*0.125) for i in range(nBins)]
        yErr = np.zeros( (3, nBins) )
        cntsBin = np.zeros( (3, nBins) )
        #
        # Doing average per pair of PMT
        diffSgnlPmt2, meanDiffSgnlPmt, cntsPerBin \
                = self.getPrtlRMS_DiffSgnlPairPMT('aveSgnl', 'sgnlPmt', \
                minSgnl, maxSgnl, nBins)
        #
        diffCSgnlPmt2, meanDiffCSgnlPmt, cntsPerBinCoinci \
                = self.getPrtlRMS_DiffSgnlPairPMT('aveSgnlCoinci', 'sgnlPmtCoinci', \
                minSgnl, maxSgnl, nBins)
        #
        for pmt_i in range(3):
            tmpAve = meanDiffSgnlPmt[pmt_i] / cntsPerBin[pmt_i]
            rmsDiffSgnlPmt = diffSgnlPmt2[pmt_i] \
                    - ( cntsPerBin[pmt_i] * tmpAve**2 )
            rmsDiffSgnlPmt = np.sqrt( rmsDiffSgnlPmt / cntsPerBin[pmt_i] )
            errRMS = rmsDiffSgnlPmt / np.sqrt(2*cntsPerBin[pmt_i])
            #
            tmpAveCoinci = meanDiffCSgnlPmt[pmt_i] / cntsPerBinCoinci[pmt_i]
            rmsDiffCSgnlPmt = diffCSgnlPmt2[pmt_i] \
                    - ( cntsPerBinCoinci[pmt_i] * tmpAveCoinci**2 )
            rmsDiffCSgnlPmt = np.sqrt( rmsDiffCSgnlPmt / cntsPerBinCoinci[pmt_i] )
            errRMSconci = rmsDiffCSgnlPmt / np.sqrt(2*cntsPerBinCoinci[pmt_i])
            #
            axs[0][pmt_i].errorbar(xBins, rmsDiffSgnlPmt, xerr=0.125, \
                    yerr=errRMS, linestyle='', marker='o', capsize=2, \
                    color=colors[0], label="Signal")
            axs[0][pmt_i].errorbar(xBins, rmsDiffCSgnlPmt, xerr=0.125, \
                    yerr=errRMSconci, linestyle='', marker='o', capsize=2, \
                    color=colors[1], label="Signal (With Delta)")
            #
            axs[0][pmt_i].set_title("PMT %i" % (pmt_i+1))
            pmt_j = pmt_i+2 if pmt_i != 2 else 1
            axs[0][pmt_i].set_ylabel(r"RMS $\left(S_"+"{%i}"%(pmt_i+1)
            +r"\,-\,S_"+"{%i}"%(pmt_j)+r"\right) / \left(\langle S_{"
            +"%i,%i}"%(pmt_i+1, pmt_j)+r"\rangle\sqrt(2)\right)$")
           
            axs[0][pmt_i].xaxis.set_tick_params(labelbottom=True)
            axs[0][pmt_i].yaxis.set_tick_params(labelbottom=True)
            axs[0][pmt_i].legend()
            #
            residual = rmsDiffSgnlPmt - rmsDiffCSgnlPmt
            residualErr = np.sqrt( errRMS**2 + errRMSconci**2 )
            axs[1][pmt_i].errorbar( xBins, residual, yerr=residualErr, \
                    xerr=0.125, linestyle='', marker='*', capsize=2, \
                    color=colors[2], label="Average: %.4f"\
                    %(np.sum(residual)/len(residual)) )
            axs[1][pmt_i].set_xlabel(r'$log_{10}\left(\left\langle S_{1,2,3}\right\rangle\,\mathrm{VEM}\right)$')
            axs[1][pmt_i].set_ylabel("Residuals")
            axs[1][pmt_i].label_outer()
            axs[1][pmt_i].xaxis.set_tick_params(labelbottom=True)
            axs[1][pmt_i].yaxis.set_tick_params(labelbottom=True)
            axs[1][pmt_i].legend()
        #
        if save:
            fig.tight_layout()
            ciel1 = PdfPages(self.dirs['plot']+'signalsDiff_Pmt.pdf')
            ciel1.savefig(bbox_inches='tight')
            ciel1.close()
        #
        #
        return


