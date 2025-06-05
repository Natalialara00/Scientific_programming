import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path

from Source.analysis.spectrumAnalyze import Spectrum, AnalyzeSpectrum, Sensitivity, AnalyzeCovariance
from Source.analysis.featuresANDstaticalanalyze import FindPeaks, StaticalAnalysis


class Visualizer:
    def __init__(self, save_figure, result_save):
        self.save_figure = Path(save_figure)  # O(1)
        self.result_save = Path(result_save)  # O(1)

        self.save_figure.mkdir(parents=True, exist_ok=True)  # O(1)
        self.result_save.mkdir(parents=True, exist_ok=True)  # O(1)

    def plot_transmission_spectra(self, Data, unetching, etching, RI_values, colors):
        fig, ax = plt.subplots(2, 3, figsize=(12, 6))  # O(1)
        ax = ax.flatten()  # O(1)
        for i in range(len(RI_values)):  # O(n
            ax[i].plot(Data['Wavelength'], Data[unetching[i]], '--', color=colors[i+2], label=f'RI {RI_values[i]} - 125 µm')
            ax[i].plot(Data['Wavelength'], Data[etching[i]], '-', color=colors[i], label=f'RI {RI_values[i]} - 25 µm')
            ax[i].set_title(f"RI = {RI_values[i]}")
            ax[i].set_xlabel('Wavelength (µm)')
            ax[i].set_ylabel('Transmission (dB)')
            ax[i].legend()
        #ax[-1].set_axis_off()
        plt.tight_layout()  # O(1)
        file_figure = self.save_figure / "1_spectra_by_RI.png"  # O(1)
        plt.savefig(file_figure, dpi=300)  # O(1)
        plt.close()  # O(1)

    def plotLandslide(self, Data, unetching, etching, RI_values):
        spectra = [Spectrum(Data['Wavelength'], Data[unetching[i]], Data[etching[i]], RI_values[i]) for i in range(len(RI_values))]  # O(n)
        analyzer = AnalyzeSpectrum(spectra)  # O(1)
        result_landslide = analyzer.landslide()  # O(1)

        plt.figure(figsize=(8, 5))  # O(1)
        plt.plot(result_landslide['RI'], result_landslide['min_125'], 'o--', label='125 µm')  # O(1)
        plt.plot(result_landslide['RI'], result_landslide['min_25'], '*-', label='25 µm')  # O(1)
        plt.title('First peak movement')  # O(1)
        plt.xlabel('Refractive Index (RI)')  # O(1)
        plt.ylabel('Wavelength ($\lambda$ = ($\mu m$))')  # O(1)
        plt.legend()  # O(1)
        plt.grid(True)  # O(1)
        plt.tight_layout()  # O(1)
        file_figure = self.save_figure / "2_first_peak_landslide.png"  # O(1)
        plt.savefig(file_figure, dpi=300)  # O(1)
        plt.close()  # O(1)

    def plotSensitivity(self, Data, unetching, etching, RI_values):
        spectra = [Spectrum(Data['Wavelength'], Data[unetching[i]], Data[etching[i]], RI_values[i]) for i in range(len(RI_values))]  # O(n)
        analyzerSensitivity = Sensitivity(spectra)  # O(1)
        result_Sensitivity = analyzerSensitivity.gradient()  # O(n)

        plt.figure(figsize=(8, 5))  # O(1)
        plt.plot(result_Sensitivity['RI'], result_Sensitivity['sensitivity 125'], '--o', label='125 µm')  # O(1)
        plt.plot(result_Sensitivity['RI'], result_Sensitivity['sensitivity 25'], '-s', label='25 µm')  # O(1)
        plt.title('Sensitivity in relation to the first movement of the peak')  # O(1)
        plt.xlabel('Refractive Index (RI)')  # O(1)
        plt.ylabel('Sensivity ($\Delta \mu m$/$\Delta RI$)')  # O(1)
        plt.legend()  # O(1)
        plt.grid(True)  # O(1)
        plt.tight_layout()  # O(1)
        file_figure = self.save_figure / "3_sensitivity.png"  # O(1)
        plt.savefig(file_figure, dpi=300)  # O(1)
        plt.close()  # O(1)

    def plotHeatmaps(self, Data, unetching, etching, RI_values):
        analyzerCovariance = AnalyzeCovariance(Data, unetching, etching)

        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        sns.heatmap(analyzerCovariance.unetchingCov(), annot=True, cmap='viridis', square=True, ax=ax[0])  # O(1)
        ax[0].set_title('Heatmap MMF unetching')  # O(1)

        sns.heatmap(analyzerCovariance.etchingCov(), annot=True, cmap='viridis', square=True, ax=ax[1])  # O(1)
        ax[1].set_title('Heatmap MMF etching')  # O(1)

        plt.tight_layout()  # O(1)
        file_figure = self.save_figure / "4_Heatmap_MMFDiameter.png"  # O(1)
        plt.savefig(file_figure, dpi=300)  # O(1)
        plt.close()  # O(1)

        fig, ax = plt.subplots(1, len(RI_values), figsize=(15, 6))  # O(1)
        fig.suptitle('Heatmap for every RI unetching and etching')  # O(1)

        color_map = sns.color_palette("viridis", as_cmap=True)  # O(1)

        for i in range(len(RI_values)):  # O(n)
            MMFdiameter_cov = Data[[unetching[i], etching[i]]].cov()
            sns.heatmap(MMFdiameter_cov, annot=True, cmap='viridis', square=True, ax=ax[i])
            ax[i].set_title(f'RI = {RI_values[i]}')

        plt.tight_layout(rect=[0, 0, 1, 0.95])   # O(1)
        file_figure = self.save_figure / "4b_Heatmap_all_RI.png"  # O(1)
        plt.savefig(file_figure, dpi=300)  # O(1)
        plt.close()  # O(1)

    def plot_relevant_peaks(self, Data, RI, unetching, etching, colors):
        relevant_peak_data = []  # O(1)
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 10))  # O(1)

        for i, (ri, col_unetch, col_etch) in enumerate(zip(RI, unetching, etching)):  # O(n)
            
            analyzer_125 = FindPeaks(Data['Wavelength'], Data[col_unetch], ri, "125 µm", results=self.save_figure)
            analyzer_125.detectPeaks()
            features_125 = analyzer_125.analyzePeaks()
            relevant_peak_data.extend(features_125)

            
            ax[0].plot(Data['Wavelength'], Data[col_unetch], label=f'RI {ri} - 125 µm', color=colors[i])

            # Plot relevant peaks
            peaks_125 = pd.DataFrame(features_125)
            if not peaks_125.empty:
                ax[0].plot(peaks_125['Wavelength'], peaks_125['Transmission'], 'o', color=colors[i])

            
            analyzer_25 = FindPeaks(Data['Wavelength'], Data[col_etch], ri, "25 µm", results=self.save_figure)
            analyzer_25.detectPeaks()
            features_25 = analyzer_25.analyzePeaks()
            relevant_peak_data.extend(features_25)

            ax[1].plot(Data['Wavelength'], Data[col_etch], label=f'RI {ri} - 25 µm', color=colors[i])
            peaks_25 = pd.DataFrame(features_25)
            if not peaks_25.empty:
                ax[1].plot(peaks_25['Wavelength'], peaks_25['Transmission'], 'o', color=colors[i])

        
        ax[0].set_xlabel('Wavelength ($\lambda$)')  # O(1)
        ax[0].set_ylabel("Transmission (dB)")  # O(1)
        ax[0].set_title('Transmission spectrum for 125$\mu m$')  # O(1)
        ax[0].legend()  # O(1)

        ax[1].set_xlabel('Wavelength ($\lambda$)')  # O(1)
        ax[1].set_ylabel("Transmission (dB)")  # O(1)
        ax[1].set_title('Transmission spectrum for 25$\mu m$')  # O(1)
        ax[1].legend()  # O(1)

        fig.tight_layout()  # O(1)
        file_figure = self.save_figure / "5_relevant_peaks_combined.png"  # O(1)
        plt.savefig(file_figure, dpi=300)  # O(1)
        plt.close()  # O(1)

        return pd.DataFrame(relevant_peak_data)  # O(n)

    def plot_histograms(self, Features_Data):
        features = ['Wavelength', 'Transmission', 'Spectral width','Prominence']  # O(1)
        plt.figure(figsize=(12, 10))  # O(1)
        for i, feature in enumerate(features):  # O(1)
            plt.subplot(2, 2, i + 1)  # O(1)
            sns.histplot(Features_Data[feature], kde=True)  # O(n)
            plt.title(f'Distribution of {feature}')  # O(1)
        plt.tight_layout()  # O(1)
        file_figure = self.save_figure / "6_histogram.png"  # O(1)
        plt.savefig(file_figure, dpi=300)  # O(1)
        plt.close()  # O(1)

    def plot_KDE(self, features_df):
        features = ['Wavelength', 'Transmission', 'Spectral width','Prominence']  # O(1)
        plt.figure(figsize=(12, 10))  # O(1)
        for i, feature in enumerate(features):  # O(1)
            plt.subplot(2, 2, i + 1)  # O(1)
            sns.kdeplot(data=features_df, x=feature, hue='MMF Diameter', fill=True, alpha=0.4)  # O(n)
            plt.title(f'Density of {feature} by MMF Diameter')  # O(1)
        plt.tight_layout()  # O(1)
        file_figure = self.save_figure / "7_density_by_group.png"  # O(1)
        plt.savefig(file_figure, dpi=300)  # O(1)
        plt.close()  # O(1)


   # This code has a computational time complexity of O(n)

