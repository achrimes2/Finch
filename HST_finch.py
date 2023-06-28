#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 14:37:17 2021

@author: ashley
"""


# FAQs
# https://www.stsci.edu/scientific-community/software/drizzlepac/frequently-asked-questions

#from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.wcs import WCS
import glob
import os

angscale = 3.80 #kpc/arcsec


drz_or_phot = input('Enter 1 for drizzling, 2 for phot, 3 for PSF, 4 for offsets, 5 for Pch and 6 for offsets compare: ')

if drz_or_phot == '1':
    from drizzlepac import tweakreg
    from drizzlepac import astrodrizzle

    folders = ['EPOCH1/if2901010']
    
    for names in folders:
        WCSkey = ['IDC_w3m18525i-FIT_REL_GAIADR2'] 
        TYPE = [1]     

        keys = []
        init = 0
        
        FILES = glob.glob('/home/ashley/Desktop/FBOTS_2023/HST/'+names+'/*_flc.fits') #Use FLC for UVIS - CTE corrected
        print(FILES)

        for file in FILES:
            input_image = fits.open(file)
            hdul = input_image[1].header
            print('FLT image = ',file)
            print('Date = ',hdul['DATE'])  
            print('Exptime = ',hdul['SAMPTIME'])
            print('WCSkey = ',hdul['WCSNAME'])
            if hdul['WCSNAME'] != init:
                keys.append(hdul['WCSNAME'])
            init = hdul['WCSNAME']
            input_image.close()
    
        astrodrizzle.AstroDrizzle(FILES, output='output.fits', wcskey=keys[0], driz_sep_bits='64,512', driz_sep_rot=0.0, combine_type='minmed', final_wcs=True, final_pixfrac=0.8, final_scale=0.025, driz_cr_corr=True, driz_cr_scale="1.2 0.7",  driz_cr_snr="3.5 3.0", final_wht_type='EXP', final_rot=0.0)
        print('Drizzle complete for '+names)

    
    print('Finished.')







##########################################################################################################

if drz_or_phot == '2':

    #### IMPORT ####    
    from astropy.stats import sigma_clipped_stats
    from astropy.stats import SigmaClip
    from astropy.coordinates import SkyCoord   
    from photutils.aperture import aperture_photometry
    from photutils.aperture import CircularAperture
    from photutils.aperture import CircularAnnulus
    from photutils.background import Background2D, MedianBackground
    from photutils.utils import calc_total_error

    ee_corr = np.loadtxt('encircled_correction.txt')
    folders = ['EPOCH1/if2901010','EPOCH1/if2901020']
    gain = 1.5

    NN = 0
    for names in folders: #first 555W then 814W
        os.chdir('/home/ashley/Desktop/FBOTS_2023/HST/'+names)
        FILE = 'output_sci.fits'
 
        print(names)
        
        plt.figure(1)
        
        #### THE DATA ####  
        # Loading in image data and header information        
        ##### 010=F555W, 020=F814W ######
        hstimageB = fits.open(FILE)
        hduB = hstimageB #[1]
        dataB = hduB[0].data
        hdulB = hstimageB[0].header
        LB = len(dataB[0,:]) #(dataB) #(dataB[0,:])
        w = WCS(hdulB)
        meanexpt = hdulB['EXPTIME']
        data =  dataB 
        
        #### ZERO POINTS ####  
        #Zero point in AB magnitudes
        ZP = -2.5*np.log10(hdulB['PHOTFLAM']) - 21.10 - 5*np.log10(hdulB['PHOTPLAM']) + 18.692
        print('ZP is ',ZP)
        # See also https://www.stsci.edu/hst/instrumentation/wfc3/data-analysis/photometric-calibration/uvis-photometric-calibration

        sky = SkyCoord('10h08m03.815s', '21d04m26.85', frame='icrs') 
        X, Y = w.world_to_pixel(sky)   
        theta = np.linspace(0,np.pi*2,100)
    
        if NN == 0: #F555W
            VMIN,VMAX = -0.0329089,0.0302982
            coord_finch = (1111,1126) #x,y
            FWHM = 0.0835 #in arcsec, uses figure(2) below, 3.34*0.025 arcsec F555W 0.0835"
        elif NN == 1: #F814W
            VMIN,VMAX = -0.0300408,0.0260148
            coord_finch = (1112,3215) #x,y
            FWHM = 0.0835  #in arcsec, uses figure(2) below, 3.35*0.025 arcsec F814W 0.0835"

        ax1 = plt.subplot(1,2,NN+1,projection=w) 
        ra = ax1.coords[0]
        dec = ax1.coords[1]
        ra.set_major_formatter('dd:mm:ss.s') 
        dec.set_major_formatter('dd:mm:ss.s')
        lon = ax1.coords[0]
        lat = ax1.coords[1]
        lon.set_ticklabel(size=10)
        lat.set_ticklabel(size=10) # normally 5
        ax1.set_xlabel('Right Ascension',fontsize=13)
        ax1.set_ylabel('Declination',fontsize=13)
        
        #Aperture selection and plotting
        arcsecsize = 4 
        pixscale = 0.025   #Final scale for WFC3/UVIS after drizzling
        imgsize = int(arcsecsize/pixscale)
        
        ax1.imshow(data, origin='lower', cmap='Greys', vmin=VMIN, vmax=VMAX,interpolation='None')
        
        aper=0.4 #0.4 arsecond radius with enc correction. Add number in arcsec, converted to pixels.
        R = aper/pixscale 
        eenum = 3+NN #3 or 4, f555w or f814w
        aperture_correction = ee_corr[5,eenum] #2 is 0.2, 5 is 0.4. eenum is F555W or F814W
        
        ax1.plot(coord_finch[0]+R*np.sin(theta),coord_finch[1]+R*np.cos(theta),'-b')
        ax1.set_xlim([X-imgsize,X+imgsize])
        ax1.set_ylim([Y-imgsize,Y+imgsize])
        

        ############# Photometry ############# 
        positions = [(coord_finch[0], coord_finch[1])]
        inner,outer=1.5,4
        aperture = CircularAperture(positions, r=R)
        annulus_aperture = CircularAnnulus(positions, r_in=inner*R, r_out=outer*R)
        apers = [aperture, annulus_aperture]
        
        bkg_estimator = MedianBackground()   
        BOX = 50
        PHOTSIG = 3.0
        bkg = Background2D(data, (BOX, BOX), filter_size=(3, 3),bkg_estimator=bkg_estimator, sigma_clip = SigmaClip(sigma=PHOTSIG, sigma_lower=None, sigma_upper=None, maxiters=5, cenfunc='median', stdfunc='std'), edge_method='pad')
        #Image is blurred 3x3 and the bkg for each pixel determined by 25x25 box centered on it.
        effective_gain = meanexpt/gain 
        error = calc_total_error(data - bkg.background, bkg.background_rms, effective_gain)
    
        # Do the source photometry
        phot_table = aperture_photometry(data, apers, error=error)
        for col in phot_table.colnames:
            phot_table[col].info.format = '%.8g'
            
        N = int(np.pi*(R**2)) #Number of pixels in aperture

        for bkgtype in ['bkg','annulus']: #two background methods
            if bkgtype == 'annulus':
                annulus_masks = annulus_aperture.to_mask(method='center')
                annulus_data = annulus_masks[0].multiply(data)
                mask = annulus_masks[0].data
                annulus_data_1d = annulus_data[mask > 0]
                mean_sigclip, median_sigclip, std_sigclip = sigma_clipped_stats(annulus_data_1d,sigma=PHOTSIG)
                bkg_sum = median_sigclip * aperture.area   #sigma clipped background summed 
                
                
                aacorr = (outer**2) / (inner**2) #Correction to error in bkg sum due to size difference between aperture and annulus.
                #Uncertainties: we have from phot_table aperture_sum_0 and aperture_sum_1
                bkg_sum_err = phot_table['aperture_sum_err_1']/aacorr  
                aper_sum_err = phot_table['aperture_sum_err_0'] 
                tot_aper_err = np.sqrt(bkg_sum_err**2 + aper_sum_err**2)

                final_sum = (phot_table['aperture_sum_0'] - bkg_sum)
                phot_table['residual_aperture_sum'] = final_sum
                phot_table['residual_aperture_sum'].info.format = '%.8g'
                source_flux = phot_table['residual_aperture_sum']
                
            elif bkgtype == 'bkg': 
                bkgaperlist = np.loadtxt('bkgaper.txt')
                extraX,extraY = bkgaperlist[:,0],bkgaperlist[:,1]
                listXpos = np.ndarray.tolist(X+extraX)
                listYpos = np.ndarray.tolist(Y+extraY)
                merged_list = tuple(zip(listXpos, listYpos))
                bkgaper = CircularAperture(positions=merged_list, r=aper/pixscale)
                bkg_table = aperture_photometry(data, bkgaper, error=error)
                bkg_sum = np.median(bkg_table['aperture_sum']) #median over 4 apertures.
                if aper != 0.4:
                    bkg_sum = np.median(bkg_table['aperture_sum'] / ((0.4**2)/(aper**2)) )
                aper_masks = bkgaper.to_mask(method='center')
                aper_data1,aper_data2,aper_data3,aper_data4 = aper_masks[0].multiply(data),aper_masks[1].multiply(data),aper_masks[2].multiply(data),aper_masks[3].multiply(data)
                maskaper = aper_masks[0].data
                aper_data_1d = np.append(aper_data4[maskaper > 0],np.append(aper_data3[maskaper > 0],np.append(aper_data1[maskaper > 0],aper_data2[maskaper > 0])))
                #Redefine the following, overwriting annulus stats
                mean_sigclip, median_sigclip, std_sigclip = sigma_clipped_stats(aper_data_1d,sigma=PHOTSIG)
                bkg_sum_err =  std_sigclip*np.sqrt(N) 
                tot_aper_err = np.sqrt(bkg_sum_err**2 + aper_sum_err**2)
                
                final_sum = (phot_table['aperture_sum_0'] - bkg_sum)
                phot_table['residual_aperture_sum'] = final_sum
                phot_table['residual_aperture_sum'].info.format = '%.8g'
                source_flux = phot_table['residual_aperture_sum']

    
    

            # Signal to noise ratio.
            # Signal in aperture divided by root(N)*rms of background
            SNR = source_flux/(np.sqrt(N)*std_sigclip)
            
            print('Aperture photometry (ECF corrected) with '+bkgtype+' estimator used')
            if SNR >= 3:
                print('SNR is ',SNR[0])
                print('Expected magnitude errors is ~ ',1/SNR[0])
                #(This is encircled energy corrected.)
                STMAG_source = -2.5*np.log10((source_flux/aperture_correction)*hdulB['PHOTFLAM']) - 21.10
                ABMAG_source = STMAG_source - 5*np.log10(hdulB['PHOTPLAM']) + 18.692
                print('Counterpart magnitude = '+str(ABMAG_source[0]))
                
                #Error up 
                STMAG_source_eu = -2.5*np.log10((np.abs(source_flux-tot_aper_err)/aperture_correction)*hdulB['PHOTFLAM']) - 21.10
                ABMAG_source_eu = (STMAG_source_eu - 5*np.log10(hdulB['PHOTPLAM']) + 18.692) - ABMAG_source
                #Error dn /aperture_correction
                STMAG_source_ed = -2.5*np.log10((np.abs(source_flux+tot_aper_err)/aperture_correction)*hdulB['PHOTFLAM']) - 21.10
                ABMAG_source_ed = ABMAG_source - (STMAG_source_ed - 5*np.log10(hdulB['PHOTPLAM']) + 18.692)
                
                print('Err_up = ',ABMAG_source_eu[0])
                print('Err_dn = ',ABMAG_source_ed[0])
                
                plt.figure(2)
                plt.subplot(1,2,NN+1)
                plt.title('FWHM calculate')
                plt.plot(data[coord_finch[1]-15:coord_finch[1]+15,coord_finch[0]]-np.mean(bkg.background_rms),'-k')
                halfmaxi = np.max(data[coord_finch[1]-15:coord_finch[1]+15,coord_finch[0]]-np.mean(bkg.background_rms))/2
                plt.plot([0,30],[halfmaxi,halfmaxi])
                plt.xlabel('Pixels')
                
            elif SNR < 3:
                STMAG_source = -2.5*np.log10((3*np.sqrt(N)*std_sigclip*hdulB['PHOTFLAM'])/aperture_correction) - 21.10
                ABMAG_source = STMAG_source - 5*np.log10(hdulB['PHOTPLAM']) + 18.692
                print('3 sigma limit at location (annulus bkg) = '+str(ABMAG_source[0]))
            
            pos_err_finch = FWHM/(2.35*SNR.value[0])
            print('Positional error of Finch is ',pos_err_finch)
            
        NN = NN + 1
    
    



if drz_or_phot == '3':

    #### IMPORT ####    
    from PIL import Image
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord 
    from astropy.stats import sigma_clipped_stats
    from astropy.stats import SigmaClip
    from photutils.aperture import aperture_photometry
    from photutils.aperture import CircularAperture
    from photutils.aperture import CircularAnnulus
    from photutils.background import Background2D, MedianBackground
    from photutils.utils import calc_total_error
    import skimage
    
    folders = ['EPOCH1/if2901010','EPOCH1/if2901020'] 
    filterns = ['F555W','F814W']
    
    ee_corr = np.loadtxt('encircled_correction.txt')

    NN = 0
    for names in folders: #first 555W then 814W
        os.chdir('/home/ashley/Desktop/FBOTS_2023/HST/'+names)   
        FILE = 'output_sci.fits'
        
        filtern = filterns[NN]
        print(names)
        
        #### THE DATA ####  
        # Loading in image data and header information        
        ##### 010=F555W, 020=F814W ######
        hstimageB = fits.open(FILE)
        hduB = hstimageB #[1]
        dataB = hduB[0].data 
        hdulB = hstimageB[0].header
        w = WCS(hdulB)
        data =  dataB 
        
        if NN == 0:
            VMIN,VMAX = -0.01,0.06
            coord_finch = (1111,1126)  #x,y
            coord_bkg = (1453,1084) #x,y
            coord_ref =  (1495,967)  #x,y 
        elif NN == 1:
            VMIN,VMAX = -0.01,0.06
            coord_finch = (1112,3215)  #x,y
            coord_bkg = (1542,3172) #x,y
            coord_ref =  (1495,3057)
            
        plt.figure(1)
        w = WCS(hstimageB[0].header)
        ax = plt.subplot(1,2,NN+1,projection=w)
        sky1 = SkyCoord(152.0193851, 21.0758833, unit='deg', frame='icrs') 
        sky2 = SkyCoord(152.0122628, 21.0758833, unit='deg', frame='icrs') 
        sky3 = SkyCoord(152.0122555, 21.0705569, unit='deg', frame='icrs') 
        sky4 = SkyCoord(152.0193775, 21.0705569, unit='deg', frame='icrs') 
        X1, Y1 = w.world_to_pixel(sky1)
        X2, Y2 = w.world_to_pixel(sky2)
        X3, Y3 = w.world_to_pixel(sky3)
        X4, Y4 = w.world_to_pixel(sky4)
        ax.imshow(data, origin='lower', cmap='Greys', vmin=VMIN,vmax=VMAX) 
        plt.xlim([int(X1)+120,int(X2)-120]);plt.ylim([int(Y3)+120,int(Y1)-120])
        ax.set_xlabel('Right Ascension',fontsize=14)
        if NN == 0:
            ax.set_ylabel('Declination',fontsize=14)
        elif NN == 1:
            ax.set_ylabel(' ',fontsize=14)
        theta = np.linspace(0,2*np.pi,100)
        plt.plot([coord_finch[0],coord_finch[0]],[coord_finch[1]+15,coord_finch[1]+30],'-r',linewidth='2')
        plt.plot([coord_finch[0]+15,coord_finch[0]+30],[coord_finch[1],coord_finch[1]],'-r',linewidth='2')
        ax.set_title(filterns[NN],fontsize=16)
        #r = 10
        #ax.text(filterns[NN],fontsize=16)
        #ax.text(filterns[NN],fontsize=16)
        # xcirc,ycirc = r*np.sin(theta),r*np.cos(theta)
        # plt.plot(xcirc+coord_finch[0],ycirc+coord_finch[1],'-r',linewidth=1)
        
        if NN == 0 or NN == 1:
            import scipy.ndimage
            plt.figure(1)
            sky1 = SkyCoord(152.0193851, 21.0758833, unit='deg', frame='icrs') 
            sky2 = SkyCoord(152.0122628, 21.0758833, unit='deg', frame='icrs') 
            sky3 = SkyCoord(152.0122555, 21.0705569, unit='deg', frame='icrs') 
            sky4 = SkyCoord(152.0193775, 21.0705569, unit='deg', frame='icrs') 
            X1, Y1 = w.world_to_pixel(sky1)
            X2, Y2 = w.world_to_pixel(sky2)
            X3, Y3 = w.world_to_pixel(sky3)
            X4, Y4 = w.world_to_pixel(sky4)
            VMINB,VMAXB = -0.007+0.003,0.007+0.003
            blurred = scipy.ndimage.gaussian_filter(dataB,sigma=1.5)
            
            axins = ax.inset_axes([0.65, 0.6, 0.35, 0.35])
            axins.imshow(blurred, origin='lower', cmap='Greys', vmin=VMINB, vmax=VMAXB,interpolation='bilinear') #from ds9 99.5%
            if NN == 0:
                x1, x2, y1, y2 = 1030, 1180, 1030+25, 1180+25
            elif NN == 1:
                x1, x2, y1, y2 = 1030, 1180, 3140+3, 3290+3
            axins.set_xlim(x1, x2)
            axins.set_ylim(y1, y2)
            axins.set_xticklabels([])
            axins.set_yticklabels([])
            axins.set_xticks([])
            axins.set_yticks([])
            #ax.indicate_inset_zoom(axins, edgecolor="black")
            # plt.xlim([int(X1)+120,int(X2)-120]);plt.ylim([int(Y3)+120,int(Y1)-120])
            # ax.set_xlabel('Right Ascension',fontsize=14)
            # ax.set_ylabel('Declination',fontsize=14)
            # plt.plot([coord_finch[0],coord_finch[0]],[coord_finch[1]+15,coord_finch[1]+30],'-r',linewidth='2')
            # plt.plot([coord_finch[0]+15,coord_finch[0]+30],[coord_finch[1],coord_finch[1]],'-r',linewidth='2')
            plt.show()
                

        plt.figure(2)
        #interpolate onto a bigger pixel scale, then select central 160x160. Mathches ref images exactly.
        newscale = 160 #new size of image, full width
        boxsize = 40 # size of initial cutout, width is 2*boxsize
        nrtot = 18 #size of pixel area, in final units, used to cutout bit for stats
        #160, 40,18 works well
        print('Scale increase ',newscale/(boxsize*2))
        
        data_finch = data[coord_finch[1]-boxsize:coord_finch[1]+boxsize,coord_finch[0]-boxsize:coord_finch[0]+boxsize]
        data_finch = skimage.transform.resize(data_finch, (newscale,newscale), order=3, mode='reflect', cval=0, clip=True, preserve_range=False, anti_aliasing=None, anti_aliasing_sigma=None)
        data_finch = data_finch[int(newscale*0.1):int(newscale*0.9),int(newscale*0.1):int(newscale*0.9)]
        
        data_bkg = data[coord_bkg[1]-boxsize:coord_bkg[1]+boxsize,coord_bkg[0]-boxsize:coord_bkg[0]+boxsize]
        data_bkg = skimage.transform.resize(data_bkg, (newscale,newscale), order=3, mode='reflect', cval=0, clip=True, preserve_range=False, anti_aliasing=None, anti_aliasing_sigma=None)
        data_bkg = data_bkg[int(newscale*0.1):int(newscale*0.9),int(newscale*0.1):int(newscale*0.9)]


        after_sub = np.empty((15,15,21))
        xshifts,yshifts,heights = [-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7],[-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7],np.linspace(0.9,1.1,21) #[0],[0] #
        x = 0
        for xshift in xshifts:
            y = 0
            for yshift in yshifts:
                h = 0
                for height in heights:
                    data_ref = data[coord_ref[1]-boxsize:coord_ref[1]+boxsize,coord_ref[0]-boxsize:coord_ref[0]+boxsize]
                    data_ref = skimage.transform.resize(data_ref, (newscale,newscale), order=3, mode='reflect', cval=0, clip=True, preserve_range=False, anti_aliasing=None, anti_aliasing_sigma=None)
                    data_ref = data_ref[int(newscale*0.1+yshift):int(newscale*0.9+yshift),int(newscale*0.1+xshift):int(newscale*0.9+xshift)]
                    
                    total_ref = np.sum(data_ref[int(newscale*0.4 - nrtot):int(newscale*0.4 + nrtot),int(newscale*0.4 - nrtot):int(newscale*0.4 + nrtot)]) 
                    total_finch = np.sum(data_finch[int(newscale*0.4 - nrtot):int(newscale*0.4 + nrtot),int(newscale*0.4 - nrtot):int(newscale*0.4 + nrtot)])
                    ratio = total_ref/total_finch
                    
                    scaled_ref = height*(data_ref)/ratio
                    
                    finch_resid = (data_finch) - scaled_ref 
                    
                    std = np.std(finch_resid)
                    
                    after_sub[x,y,h] = std #want to minimise std.
                    
                    h = h + 1
                y = y + 1
            x = x + 1
        
        indexes = np.where(after_sub == np.min(after_sub))
        yshift,xshift,height = yshifts[indexes[1][0]],xshifts[indexes[0][0]],heights[indexes[2][0]]
        data_ref = data[coord_ref[1]-boxsize:coord_ref[1]+boxsize,coord_ref[0]-boxsize:coord_ref[0]+boxsize]
        data_ref = skimage.transform.resize(data_ref, (newscale,newscale), order=3, mode='reflect', cval=0, clip=True, preserve_range=False, anti_aliasing=None, anti_aliasing_sigma=None)
        data_ref = data_ref[int(newscale*0.1+yshift):int(newscale*0.9+yshift),int(newscale*0.1+xshift):int(newscale*0.9+xshift)]

        plt.subplot(2,3,1+NN*3)
        plt.ylabel(filterns[NN],fontsize=16)
        plt.imshow(data_finch, origin='lower', cmap='Greys', vmin=VMIN,vmax=VMAX, interpolation='None')
        # plt.title('BKG')
        plt.xticks([]);plt.yticks([])
        
        plt.subplot(2,3,2+NN*3)
        #plt.title('Ref. star')
        plt.imshow(data_ref, origin='lower', cmap='Greys', vmin=VMIN,vmax=VMAX, interpolation='None') 

        total_ref = np.sum(data_ref[int(newscale*0.4 - nrtot):int(newscale*0.4 + nrtot),int(newscale*0.4 - nrtot):int(newscale*0.4 + nrtot)]) # -data_bkg[int(newscale*0.4 - nrtot):int(newscale*0.9 + nrtot),int(newscale*0.4 - nrtot):int(newscale*0.9 + nrtot)])
        total_finch = np.sum(data_finch[int(newscale*0.4 - nrtot):int(newscale*0.4 + nrtot),int(newscale*0.4 - nrtot):int(newscale*0.4 + nrtot)]) # -data_bkg[int(newscale*0.4 - nrtot):int(newscale*0.9 + nrtot),int(newscale*0.4 - nrtot):int(newscale*0.9 + nrtot)]) #central area ONLY
        ratio = total_ref/total_finch

        scaled_ref = height*(data_ref)/ratio 

        finch_resid = (data_finch) - scaled_ref 
        
        plt.xticks([]);plt.yticks([])
        
        plt.subplot(2,3,3+NN*3)
        plt.imshow(finch_resid, origin='lower', cmap='Greys', vmin=VMIN,vmax=VMAX, interpolation='None') 
        plt.xticks([]);plt.yticks([])

        ###########################################################################################################################
        #Photometry here to confirm non-detecion of any residual flux
        #With various aperture sizes, calculate MEAN +/- STDERRMEAN.
        positions = [(newscale*0.4, newscale*0.4)]
        inner,outer=1.5*(newscale/(boxsize*2)),4.0  #accounts for new pixel scale. 
        aper=0.4
        eenum = 3+NN #3 or 4, f555w or f814w
        aperture_correction = ee_corr[5,eenum] #2 is 0.2, 5 is 0.4. eenum is F555W or F814W
        pixscale = 0.025   #new scale accounted for by inner,outer.
        R = aper/pixscale
        aperture = CircularAperture(positions, r=R)
        annulus_aperture = CircularAnnulus(positions, r_in=inner*R, r_out=outer*R)
        apers = [aperture, annulus_aperture]
        bkg_estimator = MedianBackground()   
        BOX = 50
        PHOTSIG = 3.0
        bkg = Background2D(finch_resid, (BOX, BOX), filter_size=(3, 3),bkg_estimator=bkg_estimator, sigma_clip = SigmaClip(sigma=PHOTSIG, sigma_lower=None, sigma_upper=None, maxiters=5, cenfunc='median', stdfunc='std'), edge_method='pad')
        #Image is blurred 3x3 and the bkg for each pixel determined by 25x25 box centered on it.
        meanexpt = hdulB['EXPTIME']
        gain = 1.5
        effective_gain = meanexpt/gain
        error = calc_total_error(finch_resid - bkg.background, bkg.background_rms, effective_gain)

        # Do the source photometry
        phot_table = aperture_photometry(finch_resid, apers, error=error)
        for col in phot_table.colnames:
            phot_table[col].info.format = '%.8g'
            
        N = int(np.pi*(R**2))

        annulus_masks = annulus_aperture.to_mask(method='center')
        annulus_data = annulus_masks[0].multiply(finch_resid)
        mask = annulus_masks[0].data
        annulus_data_1d = annulus_data[mask > 0]
        mean_sigclip, median_sigclip, std_sigclip = sigma_clipped_stats(annulus_data_1d,sigma=PHOTSIG)
        bkg_sum = median_sigclip * aperture.area  
        
        aacorr = (outer**2) / (inner**2) #Correction to error in bkg sum due to size difference between aperture and annulus.
        bkg_sum_err = phot_table['aperture_sum_err_1']/aacorr  
        aper_sum_err = phot_table['aperture_sum_err_0'] 
        tot_aper_err = np.sqrt(bkg_sum_err**2 + aper_sum_err**2)
        
        final_sum = (phot_table['aperture_sum_0'] - bkg_sum)
        phot_table['residual_aperture_sum'] = final_sum
        phot_table['residual_aperture_sum'].info.format = '%.8g'
        source_flux = phot_table['residual_aperture_sum']
        
        SNR = source_flux/(np.sqrt(N)*std_sigclip)
        
        print('Aperture photometry (ECF corrected) with BKG estimator used: ')
        if SNR >= 3:
            STMAG_source = -2.5*np.log10((source_flux/aperture_correction)*hdulB['PHOTFLAM']) - 21.10
            ABMAG_source = STMAG_source - 5*np.log10(hdulB['PHOTPLAM']) + 18.692
            print('Counterpart magnitude = '+str(ABMAG_source))
            
            #Error up 
            STMAG_source_eu = -2.5*np.log10((np.abs(source_flux-tot_aper_err)/aperture_correction)*hdulB['PHOTFLAM']) - 21.10
            ABMAG_source_eu = (STMAG_source_eu - 5*np.log10(hdulB['PHOTPLAM']) + 18.692) - ABMAG_source
            #Error dn /aperture_correction
            STMAG_source_ed = -2.5*np.log10((np.abs(source_flux+tot_aper_err)/aperture_correction)*hdulB['PHOTFLAM']) - 21.10
            ABMAG_source_ed = ABMAG_source - (STMAG_source_ed - 5*np.log10(hdulB['PHOTPLAM']) + 18.692)
            
            print('Err_up = ',ABMAG_source_eu[0])
            print('Err_dn = ',ABMAG_source_ed[0])
            
        elif SNR < 3:
            #(NOT enc energy corrected, no need for limits?)
            STMAG_source = -2.5*np.log10((3*np.sqrt(N)*std_sigclip*hdulB['PHOTFLAM'])/aperture_correction) - 21.10
            ABMAG_source = STMAG_source - 5*np.log10(hdulB['PHOTPLAM']) + 18.692
            print('3 sigma limit at location = '+str(ABMAG_source))
            
        # theta = np.linspace(0,2*np.pi,100)
        # xcirc,ycirc = 2*R*np.sin(theta),2*R*np.cos(theta)
        # plt.plot(xcirc+newscale*0.4,ycirc+newscale*0.4,'-r',linewidth=1)
        # xcirc,ycirc = 4*R*np.sin(theta),4*R*np.cos(theta)
        # plt.plot(xcirc+newscale*0.4,ycirc+newscale*0.4,'-r',linewidth=1)
        

        
        NN = NN + 1
    






if drz_or_phot == '4':
    
    import statmorph
    import photutils
    import scipy.ndimage as ndi
    from photutils.background import Background2D, MedianBackground
    from astropy.stats import SigmaClip
    from pixel_resample_uvis import resample
    
    BOX,PHOTSIG = 50,3
    
    folders = ['EPOCH1/if2901010','EPOCH1/if2901020']
    
    morph_err = input('Enter 1 to run morph_error code too, ')

    NN = 0  #0 is 555, 1 is 814
    for names in [folders[NN]]: #first 555W then 814W
        os.chdir('/home/ashley/Desktop/FBOTS_2023/HST/'+names)
        print(names)
        
        if morph_err == '1':
            from drizzlepac import astrodrizzle
            #################################################################################
            # Method follows Lyman+ 17, Chrimes+ 19
            # 1. Create new _flt images by adding randomly from their ERR extensions
            # 2. Re-drizzle new _flt images
            # 3. Re-measure morphoparams each time. 
            # 4. Error on r50 etc is std of the output distribution.
            
            spiral_xc,spiral_yc,spiral_r50 = [],[],[]
            sat_xc,sat_yc,sat_r50 = [],[],[]
            
            for trial in np.linspace(0,99,100):
                
                if NN == 0:
                    FILE1 = fits.open('if2901p8q_flc.fits')
                    FILE2 = fits.open('if2901p9q_flc.fits')
                    FILE3 = fits.open('if2901paq_flc.fits')
                elif NN == 1:
                    FILE1 = fits.open('if2901pbq_flc.fits')
                    FILE2 = fits.open('if2901pdq_flc.fits')
                    FILE3 = fits.open('if2901pfq_flc.fits')
                key = FILE2[1].header['WCSNAME']
            
                print('On number ',trial)
                temp1 = resample(FILE1)
                temp2 = resample(FILE2)    
                temp3 = resample(FILE3)       
                
                temp1.writeto('temp1_flc.fits', clobber='True')
                temp2.writeto('temp2_flc.fits', clobber='True')
                temp3.writeto('temp3_flc.fits', clobber='True')
                
                filelist = ['temp1_flc.fits','temp2_flc.fits','temp3_flc.fits']
                astrodrizzle.AstroDrizzle(filelist, output='temp_drc.fits', wcskey=key, driz_sep_bits='64,512', driz_sep_rot=0.0, combine_type='minmed', final_wcs=True, final_pixfrac=0.8, final_scale=0.025, driz_cr_corr=True, driz_cr_scale="1.2 0.7",  driz_cr_snr="3.5 3.0", final_wht_type='EXP', final_rot=0.0)
            
                hstimageB = fits.open('temp_drc_sci.fits')
                data = hstimageB[0].data 
                
                
                if NN == 0:
                    image = data[770:1190,850:1350]
                    psf = data[966-30:966+30,1495-30:1495+30]
                elif NN == 1:
                    image = data[2890:3270,850:1350] 
                    psf = data[3057-30:3057+30,1496-30:1496+30]
    
                threshold = photutils.detect_threshold(image, 1.0)
                npixels = 5 # minimum number of connected pixels
                segm = photutils.detect_sources(image, threshold, npixels)
                
    
                for IN in [0,1]: #spiral, satellite
                    label = np.where(segm.areas == np.flip(np.sort(segm.areas))[IN])[0][0] + 1   
                    segmap = segm.data == label
                    segmap_float = ndi.uniform_filter(np.float64(segmap), size=10) #10
                    segmap = np.int64(segmap_float > 0.5) #in Sextractor this is the 'FLUX_RADIUS'
                    source_morphs = statmorph.source_morphology(image, segmap, gain=1.4, psf=psf) 
                    morph = source_morphs[0]
                    
                    if IN == 0:
                        spiral_xc.append(morph.xc_centroid)
                        spiral_yc.append(morph.yc_centroid)
                        spiral_r50.append(morph.r50)
                    elif IN == 1:
                        sat_xc.append(morph.xc_centroid)
                        sat_yc.append(morph.yc_centroid)
                        sat_r50.append(morph.r50)
            
            np.savetxt('xc_sat.txt',np.array(sat_xc))
            np.savetxt('yc_sat.txt',np.array(sat_yc))
            np.savetxt('r50_sat.txt',np.array(sat_r50))
            np.savetxt('xc_spiral.txt',np.array(spiral_xc))
            np.savetxt('yc_spiral.txt',np.array(spiral_yc))
            np.savetxt('r50_spiral.txt',np.array(spiral_r50))
            #################################################################################
    
    
        # And now do the same, but for the original image #
        # Loading in image data and header information        
        ##### 010=F555W, 020=F814W ######
        FILE = 'output_sci.fits'
        #WMFILE = 'output_wht.fits'
        hstimageB = fits.open(FILE)
        hduB = hstimageB #[1]
        dataB = hduB[0].data 
        hdulB = hstimageB[0].header
        LB = len(dataB[0,:]) 
        w = WCS(hdulB)
        meanexpt = hdulB['EXPTIME']
        data =  dataB 
        
        if NN == 0:
            image = data[770:1190,850:1350] #[480:530,480:530]
            psf = data[966-30:966+30,1495-30:1495+30]
            xfin,yfin=261,356
        elif NN == 1:
            image = data[2890:3270,850:1350] #[480:530,480:530]
            psf = data[3057-30:3057+30,1496-30:1496+30]
            xfin,yfin=261,325
    
        threshold = photutils.detect_threshold(image, 1.0)
        npixels = 5 # minimum number of connected pixels
        segm = photutils.detect_sources(image, threshold, npixels)
            
        for IN in [0,1]:  #spiral, satellite
                    
            label = np.where(segm.areas == np.flip(np.sort(segm.areas))[IN])[0][0] + 1    
            segmap = segm.data == label
         
            segmap_float = ndi.uniform_filter(np.float64(segmap), size=10) 
            segmap = np.int64(segmap_float > 0.5) #in Sextractor this is the 'FLUX_RADIUS'
            
            plt.figure(IN+1)
            plt.subplot(1,2,1)
            plt.imshow(segmap, origin='lower', cmap='gray')
            plt.subplot(1,2,2)
            VMIN,VMAX = -0.01,0.06
            plt.imshow(image, origin='lower', cmap='gray', vmin=VMIN,vmax=VMAX, interpolation='None')
            plt.plot(xfin,yfin,'oc')
            
            source_morphs = statmorph.source_morphology(image, segmap, gain=1.4, psf=psf) 
    
            morph = source_morphs[0]
    
            print('SOURCE NUMBER ',IN)
            print('xc_centroid =', morph.xc_centroid)
            print('yc_centroid =', morph.yc_centroid)
            rarc = np.sqrt( (morph.xc_centroid - xfin)**2 + (morph.yc_centroid - yfin)**2 ) * 0.025
            print('Separation in arcsec is ',rarc)
            print('Finch pos err is 0.002 arcsec.')
            print('Separation in kpc is ',rarc*angscale)
            print('r20 =', morph.r20)
            print('r50 =', morph.r50)
            print('r50 in arcsec is ',morph.r50*0.025)
            print('Host-normed offset is ',rarc/(morph.r50 * 0.025))
            print('r80 =', morph.r80)

            theta = np.linspace(0,2*np.pi,100)
            plt.plot(morph.xc_centroid+morph.r50*np.sin(theta),morph.yc_centroid+morph.r50*np.cos(theta),'-r')
            
            if IN == 0:
                r50_spiral = np.loadtxt('r50_spiral.txt')
                xc_spiral = np.loadtxt('xc_spiral.txt')
                yc_spiral = np.loadtxt('yc_spiral.txt')
                mean_r50_spiral,mean_xc_spiral,mean_yc_spiral = np.mean(r50_spiral),np.mean(xc_spiral),np.mean(yc_spiral)
                std_r50_spiral,std_xc_spiral,std_yc_spiral = np.std(r50_spiral),np.std(xc_spiral),np.std(yc_spiral)
                print('Spiral r50+/-err = ',mean_r50_spiral,std_r50_spiral)
                print('Spiral xc+/-err = ',mean_xc_spiral,std_xc_spiral)
                print('Spiral yc+/-err = ',mean_yc_spiral,std_yc_spiral)
                plt.figure(3)
                plt.subplot(2,3,1)
                plt.hist(xc_spiral,bins=30)
                plt.ylabel('Spiral')
                plt.xlabel('xc')
                plt.subplot(2,3,2)
                plt.hist(yc_spiral,bins=30)
                plt.xlabel('yc')
                plt.subplot(2,3,3)
                plt.hist(r50_spiral,bins=30)
                plt.xlabel('r50')
                
                #Offset in pixels
                rpix_spiral = np.sqrt( (mean_xc_spiral - xfin)**2 + (mean_yc_spiral - yfin)**2 )
                rpix_spiral_err1 = np.sqrt(std_xc_spiral**2 + std_yc_spiral**2)  #galaxy
                rpix_spiral_err2 = 0.0011884466757776967/0.025     #finch, from FWHM/SNR*2.35, sse CHOICE == 2 above
                rpix_spiral_err = np.sqrt(rpix_spiral_err1**2 + rpix_spiral_err2**2)
                
                #Offset in arcsec
                rarc_spiral = rpix_spiral*0.025
                rarc_spiral_err1 = rpix_spiral_err1*0.025
                rarc_spiral_err2 = rpix_spiral_err2*0.025
                rarc_spiral_err = np.sqrt(rarc_spiral_err1**2 + rarc_spiral_err2**2)
                print('   spiral offset in arsec, ',rarc_spiral,rarc_spiral_err)
                
                #Offset in kpc
                rkpc_spiral = rarc_spiral*angscale
                rkpc_spiral_err = rarc_spiral_err*angscale
                print('   spiral offset in kpc, ',rkpc_spiral,rkpc_spiral_err)
                
                #Host-normed offset
                print('   R50 in arcsec, ',mean_r50_spiral*0.025,std_r50_spiral*0.025)
                print('   R50 in kpc, ',mean_r50_spiral*0.025*angscale,std_r50_spiral*0.025*angscale)
                rnorm_spiral = rpix_spiral / mean_r50_spiral
                rnorm_spiral_err = rnorm_spiral * np.sqrt( (rpix_spiral_err/rpix_spiral)**2 + (std_r50_spiral/mean_r50_spiral)**2 )
                print('   Rnorm offset in kpc, ',rnorm_spiral,rnorm_spiral_err)


            elif IN == 1:
                r50_sat = np.loadtxt('r50_sat.txt')
                xc_sat = np.loadtxt('xc_sat.txt')
                yc_sat = np.loadtxt('yc_sat.txt')   
                cond = (xc_sat>0) & (yc_sat>0)
                print('Making sure we pick the right source...redriz means sometimes bject [1] is different.')
                mean_r50_sat,mean_xc_sat,mean_yc_sat = np.mean(r50_sat[cond]),np.mean(xc_sat[cond]),np.mean(yc_sat[cond])            
                std_r50_sat,std_xc_sat,std_yc_sat = np.std(r50_sat[cond]),np.std(xc_sat[cond]),np.std(yc_sat[cond])      
                print('sat r50+/-err = ',mean_r50_sat,std_r50_sat)
                print('sat xc+/-err = ',mean_xc_sat,std_xc_sat)
                print('sat yc+/-err = ',mean_yc_sat,std_yc_sat)
                plt.figure(3)
                plt.subplot(2,3,4)
                plt.hist(xc_sat[cond],bins=30)
                plt.ylabel('Satellite')
                plt.xlabel('xc')
                plt.subplot(2,3,5)
                plt.hist(yc_sat[cond],bins=30)
                plt.xlabel('yc')
                plt.subplot(2,3,6)
                plt.hist(r50_sat[cond],bins=30)
                plt.xlabel('r50')
            
                #Offset in pixels
                rpix_sat = np.sqrt( (mean_xc_sat - xfin)**2 + (mean_yc_sat - yfin)**2 )
                rpix_sat_err1 = np.sqrt(std_xc_sat**2 + std_yc_sat**2)  #galaxy
                rpix_sat_err2 = 0.0017064299207156235/0.025      #finch, from FWHM/SNR*2.35
                rpix_sat_err = np.sqrt(rpix_sat_err1**2 + rpix_sat_err2**2)
                
                #Offset in arcsec
                rarc_sat = rpix_sat*0.025
                rarc_sat_err1 = rpix_sat_err1*0.025
                rarc_sat_err2 = rpix_sat_err2*0.025
                rarc_sat_err = np.sqrt(rarc_sat_err1**2 + rarc_sat_err2**2)
                print('   Sat offset in arsec, ',rarc_sat,rarc_sat_err)
                
                #Offset in kpc
                rkpc_sat = rarc_sat*angscale
                rkpc_sat_err = rarc_sat_err*angscale
                print('   Sat offset in kpc, ',rkpc_sat,rkpc_sat_err)
                
                #Host-normed offset
                print('   R50 in arcsec, ',mean_r50_sat*0.025,std_r50_sat*0.025)
                print('   R50 in kpc, ',mean_r50_sat*0.025*angscale,std_r50_sat*0.025*angscale)
                rnorm_sat = rpix_sat / mean_r50_sat
                rnorm_sat_err = rnorm_sat * np.sqrt( (rpix_sat_err/rpix_sat)**2 + (std_r50_sat/mean_r50_sat)**2 )
                print('   Rnorm offset in kpc, ',rnorm_sat,rnorm_sat_err)


            
if drz_or_phot == '5':
    
    #Pchance calculation
    
    for NN in [0,1]:
        if NN == 0:
            m = 18.94 
            Re = 1.18 #in arcsec
            Ro = 4.35 
            err_pos_finch = 0.0012
        elif NN == 1:
            m = 22.61 #18.94 #22.61 SDSS r-band
            Re = 0.39 #arcsec
            Ro = 1.43
            err_pos_finch = 0.0017
            
        #R depends on Rh and pos_err
        R1 = 2*Re
        R2 = 3*err_pos_finch
        R3 = np.sqrt( Ro**2 + 4*(Re**2) )
        
        R = np.max([R1,R2,R3])
        
        sigma = (1 / (0.33*np.log(10))) * 10**(0.33*(m-24) - 2.44) #Berger2010 version
        
        eta = np.pi*(R**2)*sigma
        
        Pchance = 1 - np.exp(-eta)
            
        print('Pchance % is ',Pchance*100)
        
    
    
    
    
    
    

if drz_or_phot == '6':
    ######## Comparison data
    hostnormedoffset_compare = np.loadtxt('comparison_normed_data.txt',skiprows=2)
    # 0 = FRBs - IR+UV - Bhandari 2022
    # 1 = LGRBs - IR/UVIS - Blanchard et al. 2016
    # 2 = LGRBs - IR - Lyman et al. 2017
    # 3 = SGRBs - UV - Fong & Berger 2013 (inc. Fong et al. 2010)
    # 4 = SGRBs - IR - Fong & Berger 2013 (inc. Fong et al. 2010)
    # 5 = SLSNe
    # 6 = KK2012
    # 7 = Ca-rich SNe, De et al. 2020 
    # 8 = SGRBs - optical - Fong+ 2022

    
    plt.figure(42)
    plt.subplot(1,2,2)
    ## Comparison samples
    lgrbs = np.append(hostnormedoffset_compare[:,1][hostnormedoffset_compare[:,1]<999],hostnormedoffset_compare[:,2][hostnormedoffset_compare[:,2]<999])
    N,bins,patches = plt.hist(lgrbs,histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='c',linestyle='--')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],'--c',linewidth=2,label='LGRBs')
    
    N,bins,patches = plt.hist(hostnormedoffset_compare[:,5][hostnormedoffset_compare[:,5]<999],histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='y',linestyle='-')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],'-y',linewidth=2,label='SLSNe')
  
    N,bins,patches = plt.hist(hostnormedoffset_compare[:,6][hostnormedoffset_compare[:,6]<999],histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='m',linestyle='--')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],'--m',linewidth=2,label='CCSNe')

    N,bins,patches = plt.hist(hostnormedoffset_compare[:,0][hostnormedoffset_compare[:,0]<999],histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='b',linestyle='-')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],'-b',linewidth=2,label='FRBs') #now Bhandari 2022, was Mannings 21 before
    
    N,bins,patches = plt.hist(hostnormedoffset_compare[:,7][hostnormedoffset_compare[:,7]<999],histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='#ff9900',linestyle='--')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],linestyle='--',color='#ff9900',linewidth=2,label='Ca-rich SNe')
    

    N,bins,patches = plt.hist(hostnormedoffset_compare[:,8][hostnormedoffset_compare[:,8]<999],histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='r',linestyle='-')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],'-r',linewidth=2,label='SGRBs')
   

    plt.ylim([0,1])
    plt.xlim([0.02,25])
    plt.xscale('log') 
    plt.xlabel(r'Host-normalised offset $r_{n}$',fontsize=16)  
    plt.ylabel('Cumulative Fraction',fontsize=16) 
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend()    
    
    #spiral
    plt.plot([3.67,3.67],[0,1],'-k',linewidth=3) #555
    plt.plot([4.25,4.25],[0,1],'--k',linewidth=3) #814

    #satellite
    plt.plot([3.98,3.98],[0,1],'-k',linewidth=1) #555
    plt.plot([4.66,4.66],[0,1],'--k',linewidth=1) #814
    
    
    

    plt.figure(42)
    plt.subplot(1,2,1)
    ######## Comparison data
    offset_compare = np.loadtxt('comparison_offset_data.txt',skiprows=2)
    # 0 = FRBs - IR+UV - Bhnadari+22
    # 1 = LGRBs - IR/UVIS - Blanchard et al. 2016
    # 2 = LGRBs - IR - Lyman et al. 2017
    # 3 = SGRBs - IR/UVIS - Fong+22
    # 4 = Type Ia SNe - Wang 2013
    # 5 = All CCSNe - Schulze 2020
    # 6 = SLSNe L15, S20
    # 7 = Ca-rich SNe, De et al. 2020
    # 8 = SGRBs - IR/UVIS - Fong & Berger 2013 (inc. Fong et al. 2010)

    #FBOT offsets
    offset_fbots = np.array([1.7,1.9,0.26,1.19]) #in kpc
    N,bins,patches = plt.hist(offset_fbots,histtype='step',density=True,cumulative=True,bins=1000,color='k',linestyle='-',linewidth=3)
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],color='k',linewidth=3,linestyle='-',label='LFBOTs')


    lgrbs = np.append(offset_compare[:,1][offset_compare[:,1]<999],offset_compare[:,2][offset_compare[:,2]<999])
    N,bins,patches = plt.hist(lgrbs,histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='c',linestyle='--')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],'--c',linewidth=2,label='LGRBs')
    
    N,bins,patches = plt.hist(offset_compare[:,6][offset_compare[:,6]<999],histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='y',linestyle='-')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],'-y',linewidth=2,label='SLSNe')

    N,bins,patches = plt.hist(offset_compare[:,5][offset_compare[:,5]<999],histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='m',linestyle='--')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],'--m',linewidth=2,label='CCSNe')
    
    N,bins,patches = plt.hist(offset_compare[:,0][offset_compare[:,0]<999],histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='b',linestyle='-')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],'-b',linewidth=2,label='FRBs')
    
    N,bins,patches = plt.hist(offset_compare[:,4][offset_compare[:,4]<999],histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='g',linestyle='-')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],'-g',linewidth=2,label='SNe Ia')
    
    N,bins,patches = plt.hist(offset_compare[:,7][offset_compare[:,7]<999],histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='#ff9900',linestyle='--')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],linestyle='--',color='#ff9900',linewidth=2,label='Ca-rich SNe')
    
    N,bins,patches = plt.hist(offset_compare[:,3][offset_compare[:,3]<999],histtype='step',density=True,cumulative=True,bins=1000,linewidth=2,color='r',linestyle='-')
    patches[0].set_xy(patches[0].get_xy()[:-1])
    plt.plot([0,1],[-1,-2],'-r',linewidth=2,label='SGRBs')


    plt.ylim([0,1])
    plt.xlim([0.02,100])
    plt.xscale('log')
    #plt.yticks([])
    #plt.xticks([])  
    plt.xlabel(r'Offset $\delta$r [kpc]',fontsize=16)  
    plt.ylabel('Cumulative Fraction',fontsize=16)  
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()    
    
    #spiral
    plt.plot([16.51,16.51],[0,1],'-k',linewidth=3) #555
    plt.plot([16.55,16.55],[0,1],'--k',linewidth=3) #814

    #satellite
    plt.plot([5.35,5.35],[0,1],'-k',linewidth=1) #555
    plt.plot([5.34,5.34],[0,1],'--k',linewidth=1) #814




    
plt.show()
    
    
    
