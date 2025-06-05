import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_prominences, peak_widths
from scipy.stats import f_oneway, shapiro, wilcoxon, pearsonr
from pathlib import Path


"""
The FindPeaks class allows detecting and analysing relevant peaks in transmission spectra. 
It uses the detectPeaks method to identify the peaks in the spectrum and extract their coordinates. 
With analyzePeaks, the most significant peaks are selected according to their prominence (depth) and their 
spectral width, discarding those that could be mistaken for noise. For each relevant peak 
key features are extracted such as: fibre diameter, refractive index, wavelength, 
transmission level, spectral width and prominence.
"""

class FindPeaks:
    def __init__(self, Wavelength, Transmission, RI, MMFDiameter, results):  # O(n)
        self.Wavelength = np.array(Wavelength)  # O(n)
        self.Transmission = np.array(Transmission)  # O(n)
        self.RI = RI  # O(1)
        self.MMFDiameter = MMFDiameter  # O(1)
        self.results = Path(results)  # O(1)
        self.results.mkdir(parents=True, exist_ok=True)  # O(1)
        # Initialize attributes for peaks and results
        self.peaks = []  # O(1)
        self.Peak_Wavelengths = []  # O(1)
        self.Peak_Transmissions = []  # O(1)
        self.Peak_spectral_width = []  # O(1)
        self.peak_prominences = []  # Store most prominent peaks   # O(1)
        self.Peak_Data = []  # Save the most relevant peaks  # O(1)

    def detectPeaks(self):
        # Find peaks
        peaks, _ = find_peaks(-self.Transmission, threshold=0.0001, distance=50)  # O(n)
        self.peaks = peaks  # O(1)
        self.Peak_Wavelengths = self.Wavelength[self.peaks]  # O(n)
        self.Peak_Transmissions = self.Transmission[self.peaks]  # O(n)

        return self  # O(1)

    def analyzePeaks(self):

        wavelength_step = np.mean(np.diff(self.Wavelength))  # O(n)
        # Calculate prominences and spectral widths for the peaks
        prominences = peak_prominences(-self.Transmission, self.peaks)[0]  # O(n)
        spectral_widths = peak_widths(-self.Transmission, self.peaks, rel_height=0.5)[0] * wavelength_step  # O(n)

        # Calculate a score based on prominence to spectral width ratio
        relevant_peak = np.argsort(prominences / spectral_widths)[-min(3, len(prominences)):][::-1]  # O(n)

        # Save the most relevant peaks
        for i in relevant_peak:  # O(1)
            self.Peak_Data.append({
                'MMF Diameter': self.MMFDiameter,
                'surrounding environment': self.RI,
                'Wavelength': float(self.Peak_Wavelengths[i]),
                'Transmission': float(self.Peak_Transmissions[i]),
                'Spectral width': float(spectral_widths[i]),
                'Prominence': float(prominences[i])})  # O(1)
            
        Features_Data = pd.DataFrame(self.Peak_Data)  ## O(n)
        
        return self.Peak_Data  # O(1)


"""
The StaticalAnalysis class performs statistical analysis on the features extracted from the peaks. 
It includes summary statistics (mean, median, variance, extreme values) using summary_statistics. 
It also allows to assess the normality of the data with Shapiro, and to compare groups using ANOVA and Wilcoxon. 
The get_paired_data method matches data by refractive index environment, allowing Wilcoxon 
to be applied comparing specific conditions (e.g. different fibre diameters).
"""

class StaticalAnalysis:
    def __init__(self, Peak_Features_Data, result_save):  
        self.Data = Peak_Features_Data  # O(1)
        self.result_save = Path(result_save)   # O(1)
        self.result_save.mkdir(parents=True, exist_ok=True)  # O(1)
       

    def summary_statistics(self):
        sumary = self.Data.describe(include='all')  # O(n)
        sumary.to_csv(self.result_save / "summary_statistics.csv")  # O(n)
        return sumary  # O(n)
    
    def get_paired_data(self, column='Wavelength', diameter1='125 µm', diameter2='25 µm'):
        MMF_125_Wavelength = self.Data[self.Data['MMF Diameter'] == diameter1][['surrounding environment', column]]  # O(n)
        MMF_25_Wavelength = self.Data[self.Data['MMF Diameter'] == diameter2][['surrounding environment', column]]  # O(n)

        MMF_125_Wavelength = MMF_125_Wavelength.rename(columns={column: f'{column}_125'})  # O(n)
        MMF_25_Wavelength = MMF_25_Wavelength.rename(columns={column: f'{column}_25'})  # O(n)

        paired = pd.merge(MMF_125_Wavelength, MMF_25_Wavelength, on='surrounding environment').dropna()  # O(n Log n)
        return paired  # O(1)


    def shapiro_test(self, column='Wavelength', group_by='MMF Diameter', alpha=0.05):
        results = []  # O(1)
        for Diameter, Shapiro_Data in self.Data.groupby(group_by):  # O(n)
            if len(Shapiro_Data[column]) < 3:  # O(1)
                results.append({'Diameter': Diameter, 'p-value': np.nan, 'Result': 'few samples'})  # O(1)
                continue
            stat, p = shapiro(Shapiro_Data[column])  # O(n)
            result = 'Normal' if p > alpha else 'Not normal'  # O(1)
            results.append({'Diameter': Diameter, 'p-value': p, 'Result': result})  # O(1)
        
        Shapiro_results = pd.DataFrame(results)  # O(n)
        Shapiro_results.to_csv(self.result_save / 'Shapiro_results.csv', index=False)  # O(n)
        return Shapiro_results  # O(1)
    
    def ANOVA_Test(self, column='Wavelength', group_by='MMF Diameter', alpha=0.05): 
        groups = [group[column].values for _, group in self.Data.groupby(group_by)]  # O(n)
        h_stat, p_value  = f_oneway(*groups)  # O(n)
        ANOVA_results = 'Significant difference' if p_value < alpha else 'No significant difference'  # O(1)
        ANOVA_Data = pd.DataFrame([{'H-statistic': h_stat,'p-value': p_value,'Result': ANOVA_results}])  # O(1)
        ANOVA_Data.to_csv(self.result_save / 'anova_results.csv', index=False)  # O(1)
        return ANOVA_Data  # O(1)
    
    def Wilcoxon_test(self, Diameter1='125 µm', Diameter2='25 µm', column='Wavelength', alpha=0.05):
        paired_data = self.get_paired_data(column=column, diameter1=Diameter1, diameter2=Diameter2)  # O(n)
        if len(paired_data) < 3:  # O(1)
            return pd.DataFrame([{'Error': 'Insufficient paired samples'}])  # O(1)

        stat, p_value = wilcoxon( paired_data[f'{column}_125'], paired_data[f'{column}_25'],
            alternative='greater')

        result = 'Significant' if p_value < alpha else 'Not significant'  # O(1)

        Wilcoxon_Data = pd.DataFrame([{'Diameter1': Diameter1, 'Diameter2': Diameter2, 'n samples': len(paired_data),
            'W': stat, 'p-value': p_value, 'Result': result}])
        Wilcoxon_Data.to_csv(self.result_save / 'Wilcoxon_results.csv', index=False)  # O(1)
        return Wilcoxon_Data  # O(1)
    
        # This code has a computational time complexity of O(n Log n)