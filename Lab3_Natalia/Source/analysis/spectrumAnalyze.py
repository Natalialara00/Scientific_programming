import numpy as np
import pandas as pd


"""
The Spectrum class allows storing a spectrum with and without etching together with its 
wavelength and refractive index. It contains a method to identify the first minimum 
in the range 1.16 to 1.23 µm in both spectra.
"""

class Spectrum:
  def __init__(self, Wavelength, Transmission, TransmissionEtching, RI):
    self.Wavelength =np.array(Wavelength)  # O(n) 
    self.Transmission = np.array(Transmission)  # O(n) 
    self.TransmissionEtching = np.array(TransmissionEtching)  # O(n) 
    self.RI = RI  # O(1) 

  def firstMinimumPeak (self):
    interval = (self.Wavelength >= 1.16) & (self.Wavelength <= 1.23)  # O(n) 
    # Filter values in the interval 1.16 and 1.23
    Wavelength_interval = self.Wavelength[interval]  # O(n) 
    Transmission_interval  = self.Transmission[interval]  # O(n) 
    TransmissionEtching_interval = self.TransmissionEtching[interval]  # O(n) 
   # Find the position of the minimum values in the interval. Extract the wavelength value.
    minTransmission = Wavelength_interval[np.argmin(Transmission_interval)]  # O(n) 
    minTransmissionEtching = Wavelength_interval[np.argmin(TransmissionEtching_interval)]  # O(n) 

    return minTransmission, minTransmissionEtching  # O(1) 

"""
AnalyzeSpectrum groups several spectra (Spectrum) to analyse the behaviour of the first minimum 
as a function of refractive index. It extracts the minimum transmission values for fibres of different 
diameters and associates them with their corresponding RI. 
"""

class AnalyzeSpectrum:
  def __init__(self, spectra):
      self.spectra = spectra # Different spectrum   # O(1) 

  def landslide(self):
    self.minPeak125 = []  # O(1) 
    self.minPeak25 = []  # O(1) 
    self.RI_values = []  # O(1) 

    for spectrum in self.spectra:
      min_125, min_25 = spectrum.firstMinimumPeak()  # O(n) 
      if min_125 is not None and min_25 is not None:  # O(1) 
        self.minPeak125.append(min_125)  # O(1) 
        self.minPeak25.append(min_25)  # O(1) 
        self.RI_values.append(float(spectrum.RI))  # O(1) 

    return {'RI': self.RI_values, 'min_125': self.minPeak125, 'min_25': self.minPeak25}  # O(1)  

"""
the AnalyzeSpectrum sensitivity daughter class that calculates the spectral sensitivity as the derivative of the first minimum 
with respect to the refractive index (RI) for two different fibre optic diameters.
"""

class Sensitivity(AnalyzeSpectrum):

  def __init__(self, spectra):
      super().__init__(spectra) # To load Landslide method   # O(1) 

  def gradient(self):
    self.landslide()  # O(n) 

    self.sensitivity125 = np.gradient(self.minPeak125, self.RI_values)  # O(n) 
    self.sensitivity25 = np.gradient(self.minPeak25, self.RI_values)  # O(n) 

    return {'RI': self.RI_values, 'sensitivity 125': self.sensitivity125, 'sensitivity 25': self.sensitivity25}  # O(1) 
  

"""
Class designed to calculate covariance matrices between columns of optical data with and without treatment 
(etching and unetching). It allows to analyse the relationship and variability between diameters and RI, being useful 
to identify patterns through representations such as heat maps.
"""

class AnalyzeCovariance:
   
  def __init__(self, Data, unetching, etching):
   self.Data = Data  # O(1) 
   self.unetching = unetching  # O(1) 
   self.etching = etching  # O(1) 

  def unetchingCov(self):
   return self.Data[self.unetching].cov()  # O(n²)
  
  def etchingCov(self):
   return self.Data[self.etching].cov()  # O(n²)

  def RICov(self):
   Covariance = []  # O(1)
   for unetchingCov, etchingCov in zip (self.unetching, self.etching):  # O(n)
     MMFdiameter = self.Data[[unetchingCov, etchingCov]].cov()  # O(n²)
     Covariance.append(MMFdiameter)  # O(1)
   return Covariance  # O(1)

    # This code has a computational time complexity of O(n²)