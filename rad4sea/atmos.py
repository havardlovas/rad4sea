# A module for simple manipulation of radiometric data 
import spectral as sp
import numpy as np
from Py6S import *

from datetime import datetime
import time


## Based on openhsi/atmos.py

import copy
from tqdm import tqdm
from multiprocessing.dummy import Pool

def _sixs_run_one_wavelength(s, wv:float) -> float:
    """Runs one instance of 6SV for one wavelength wv"""
    s.outputs = None
    a = copy.deepcopy(s)
    # TODO: Pass the wavelength filter
    a.wavelength = Wavelength(wv)
    a.run()
    return SixSHelpers.Wavelengths.recursive_getattr(a.outputs, "pixel_radiance")

# From openhsi/atmos
def run_wavelengths(wavelengths:np.array, s, n_threads:int = 8) -> np.array:
    """Modified version of SixSHelpers.Wavelengths.run_wavelengths that has a progress bar.
    This implementation uses threading (through Python's multiprocessing API)."""
    
    try:
        with Pool(n_threads) as p, tqdm(total=len(wavelengths)) as pbar:
            res = [p.apply_async( _sixs_run_one_wavelength, args=(s, wavelengths[i],), 
                    callback=lambda _: pbar.update(1)) for i in range(len(wavelengths))]
            results = [r.get() for r in res]
    except IndexError:
        print('Multiprocessing failed for Py6S, running single processing with helper func from py6s')
        
        wv, results = SixSHelpers.Wavelengths.run_wavelengths(s, wavelengths, output_name='pixel_radiance')

                
    
    return np.array(results)


# From ISTUTOR
def srf(x, mu, sigma):
    """
    Convoluting a gaussian waveband of central wavelength mu and std sigma"""
    u = (x-mu)/abs(sigma)
    y = (1.0/(np.sqrt(2.0 * np.pi)*abs(sigma)))*np.exp(-u*u/2.0)
    return y/y.sum()

class OceanRad():

    """
    A class representing ocean radiometry data.

    Can be initialized from radiance or reflectance parameters.
    """
    def __init__(self, data_source, parameters):
        """
        Initializes the OceanRad object.

        Args:
            data_source (str): Source of the data ("radiance" or "reflectance").
            parameters (dict): Parameters for the data source.
        """

        self.data_source = data_source
        self.parameters = parameters
        
        # To support different class instantiations (currenly only radiance)
        if data_source == "radiance":
            self.radiance_multiplier = parameters["radiance_multiplier"]
            self.radiance_spy_object = parameters["radiance_spy_object"]
            self.timestamp_unix = parameters["timestamp_unix"]
            self.mask_nodata = parameters["mask_nodata"]
            
            self.altitude_msl = parameters["altitude_msl"]
            self.lon = parameters["lon"]
            self.lat = parameters["lat"]
        
        
        # Read out wavelengths as these be needed for atmos correction
        self.wl = np.array(self.radiance_spy_object.metadata["wavelength"]).astype(np.float64)
        self.fwhm = np.array(self.radiance_spy_object.metadata["fwhm"]).astype(np.float64)
        
        self.E_d_W_m2 = self.get_downwelling_irr()
        
        
    @classmethod
    def from_radiance(cls, rad_parameters):
        """
        Creates an OceanRad object from radiance parameters.

        Args:
            radiance_parameters (dict): Parameters specific to radiance data.

        Returns:
            OceanRad: An OceanRad object initialized with radiance data.
        """

        return cls(data_source="radiance", parameters=rad_parameters)

    @classmethod
    def from_reflectance(cls, refl_parameters):
        """
        Creates an OceanRad object from reflectance parameters.

        Args:
            reflectance_parameters (dict): Parameters specific to reflectance data.

        Returns:
            OceanRad: An OceanRad object initialized with reflectance data.
        """

        return cls(data_source="reflectance", parameters=refl_parameters)
    
    def get_downwelling_irr(self, AEROPROFILE = AeroProfile.Maritime):
        s = SixS()
        # Use gmtime to convert to GMT time structure
        gm_time = time.gmtime(self.timestamp_unix)

        # Use strftime to format the string representation with desired format
        time_format_desired = "%Y-%m-%d %H:%M:%S"
        z_time_str = time.strftime(time_format_desired, gm_time)
        z_time = datetime.strptime(z_time_str,time_format_desired) 
        
        #Viewing and sun geometry
        s.geometry = Geometry.User()
        s.geometry.day = z_time.day
        s.geometry.month = z_time.month
        
        # Position, direction of sensor
        lat = self.lat # latitude in degrees
        lon = self.lon # longitude in degrees
        
        # Stand at zero now
        zen = 0
        azi = 0
        alt = self.altitude_msl / 10**3 # km
        
        date_str = z_time_str.split()[0]
        #print(date_str)
        s.atmos_profile = AtmosProfile.FromLatitudeAndDate(lat, date_str)

        #Viewing and sun geometry
        s.geometry = Geometry.User()
        s.geometry.day = z_time.day
        s.geometry.month = z_time.month
        dt_str = f"{z_time.year}-{z_time.month:02d}-{z_time.day:02d} {z_time.hour:02d}:{z_time.minute:02d}:{z_time.second:02d}"
        s.geometry.from_time_and_location(lat, lon, dt_str, zen, azi)
        
        print(dt_str)
        print(date_str)
        print(f"{lon}, {lat}, {alt}")


        #Altitude
        s.altitudes = Altitudes()
        s.altitudes.set_sensor_custom_altitude(alt) # km
        s.altitudes.set_target_sea_level()
        s.ground_reflectance = GroundReflectance.HomogeneousLambertian(1) # Ground reflectance for spectralon panel


        # Aerosol
        s.aero_profile = AeroProfile.PredefinedType(AEROPROFILE)
        
        
        wl_sim_step = 0.5 # Step nanometers
        wl_sim_start = 390 # Hardcoded lower limit nanometers
        # Simulate and ensure to simulate for higher wavelengths than sensor
        wl_sim_stop  = 1e2 * np.round(self.wl.max()/1e2) + 100
        
        wl_sim_arr_nm = np.arange(wl_sim_start, wl_sim_stop, wl_sim_step)

        wl_sim_arr_um = wl_sim_arr_nm/1e3 # convert to μm for Py6S
        
        # SImulate radiances
        radiances_simulated = run_wavelengths(wl_sim_arr_um, s) # units of (W/m^2/sr/μm)
        radiance_multip_py6s = (1/1e3) # Convert to (W/m^2/sr/nm)
        
        # Calculate spectral responses for all channels
        radiances_resampled = []

        cleaned_wl_sim_arr_um = wl_sim_arr_um[~np.isnan(radiances_simulated)]
        cleaned_radiance = radiances_simulated[~np.isnan(radiances_simulated)]

        for wavelength_nm, width_nm in zip(self.wl, self.fwhm):
            sigma_nm = width_nm/2.3548 # Relationship between standard deviation and FWHM
            
            # Calculating convolution with bands 
            radiance_resampled = np.sum(srf(cleaned_wl_sim_arr_um*1e3, wavelength_nm, sigma_nm) * cleaned_radiance)
            radiances_resampled.append(radiance_resampled)
        
        
        self.wl_sim = cleaned_wl_sim_arr_um*1e3 # Map to nanometers
        self.rad_sim = cleaned_radiance*radiance_multip_py6s
        
        self.rad_sim_resamp = np.array(radiances_resampled)*radiance_multip_py6s
        
        E_d_W_m2 = self.rad_sim_resamp*np.pi # Lambertian
        
        return E_d_W_m2
        
        
        
        
        
        
        
        
    def as_remote_sensing_reflectance(self, is_above_sea = True, do_simple_division = True):
        
        if do_simple_division:
            # Step 1 calculate
            self.E_d_W_m2 = self.get_downwelling_irr()
            self.L_u_W_m2_sr_nm = self.radiance_spy_object[self.mask_nodata, :]*self.radiance_multiplier
            self.Rrs = (self.L_u_W_m2_sr_nm / self.E_d_W_m2).astype(np.float32)

    
    def extract_cube_to_spectrum_list(self, rows, cols):
        """Puts the radiance data (or similar) into a list based on mask_nodata"""
        w_im, h_im, n_wl = self.radiance_spy_object.shape
        item_bytes = 4
        size_datacube = w_im*h_im*n_wl*item_bytes / (1024**3)
        size_spectra = rows.size*n_wl*item_bytes / (1024**3)
        print(f'The datacube has {w_im*h_im} cells')
        print(f'The datacube is {size_datacube} GB')
        print(f'The datacube has {size_spectra} GB of actual data')
        
        self.spectrum_list = np.zeros((rows.size, n_wl), dtype = np.float32)
        # Iterate through bandwize and put data in chunks.
        # Potentially, you can do band batches, e.g. 10 bands at a time to lower the speed loss from 
        for band_nr in range(n_wl):
            band_im = self.radiance_spy_object[:, :, band_nr]
            self.spectrum_list[:, band_nr] = band_im[self.mask_nodata]
            
    
    def write_reflectance(self):
        """Puts the radiance data (or similar) into a list based on mask_nodata"""
        w_im, h_im, n_wl = self.radiance_spy_object.shape
        item_bytes = 4
        size_datacube = w_im*h_im*n_wl*item_bytes / (1024**3)
        size_spectra = self.mask_nodata.size*n_wl*item_bytes / (1024**3)
        #print(f'The datacube has {w_im*h_im} cells')
        #print(f'The datacube is {size_datacube} GB')
        #print(f'The datacube has {size_spectra} GB of actual data')

        
        self.spectrum_list = np.zeros((self.mask_nodata.size, n_wl), dtype = np.float32)
        
        mm = self.radiance_spy_object.open_memmap(writable = True)
        # Iterate through bandwize and put data in chunks.
        # Potentially, you can do band batches, e.g. 10 bands at a time to lower the speed loss from 
        for band_nr in range(n_wl):
            # Iterate bands and modify one band image at a time
            band_im_refl = self.radiance_spy_object[:, :, band_nr] / self.E_d_W_m2[band_nr]
            
            band_im_refl *= self.radiance_multiplier # Scale to get Rrs in [1/sr]
            
            mm[self.mask_nodata.squeeze(), band_nr] = band_im_refl[self.mask_nodata]
        
        del mm
            
        
        
        
            
        
        
        
"""        
        
        
anc_file_path = "/home/notebook/cogs/remoy_202208310800_ntnu_hyperspectral_74m/processed/Output/GIS/AncillaryData/remoy_202208310800_ntnu_hyperspectral_74m_transectnr_5_chunknr_0_north_east.hdr"
# The written ancillary data lis of band names is contained within *.hdr file
anc_image_object = sp.io.envi.open(anc_file_path)
anc_band_list = anc_image_object.metadata['band names']
anc_nodata = float(anc_image_object.metadata['data ignore value'])

# The pixel number in gridded form can e.g. be plotted
unix_time_grid = anc_image_object[:,:, anc_band_list.index('unix_time_grid')]
altitude_grid_msl = anc_image_object[:,:, anc_band_list.index('hsi_alts_msl')]

# Create a generic mask
mask_nodata = np.zeros(unix_time_grid.shape, dtype = bool)
mask_nodata[unix_time_grid != anc_nodata] = True

# Find the lon-lat-
x_ecef = anc_image_object[:,:, anc_band_list.index('position_ecef_0')][mask_nodata].mean()
y_ecef = anc_image_object[:,:, anc_band_list.index('position_ecef_1')][mask_nodata].mean()
z_ecef = anc_image_object[:,:, anc_band_list.index('position_ecef_2')][mask_nodata].mean()
import pymap3d as pm
lat, lon, _ = pm.ecef2geodetic(x_ecef,y_ecef, z_ecef)



timestamp_unix = unix_time_grid[mask_nodata].mean()

altitude_msl = altitude_grid_msl[mask_nodata].mean()

# Example usage
radiance_multiplier = (1 / 1000) #(mW/cm^2*sr*um)*1000.0000 ->(mW/cm^2*sr*um)
radiance_multiplier *= (1e-3 / 1e-4) #(mW/cm^2*sr*um) -> (W/m^2*sr*um)
radiance_multiplier *= (1 / 1e3) # (W/m^2*sr*um) -> (W/m^2*sr*nm)

# The written ancillary data list of band names is contained within *.hdr file

spectral_file_path = "/home/notebook/cogs/remoy_202208310800_ntnu_hyperspectral_74m/processed/Output/GIS/HSIDatacubes/remoy_202208310800_ntnu_hyperspectral_74m_transectnr_5_chunknr_0_north_east.hdr"
radiance_spy_object = sp.io.envi.open(spectral_file_path)

radiance_data = {"radiance_spy_object": radiance_spy_object, 
                 "radiance_multiplier": radiance_multiplier,
                 "timestamp_unix": timestamp_unix,
                 "altitude_msl": altitude_msl,
                 "lon": lon,
                 "lat": lat,
                 "mask_nodata": mask_nodata}"""



