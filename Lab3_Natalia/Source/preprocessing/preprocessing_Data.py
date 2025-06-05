import pandas as pd
from pathlib import Path

"""
The Preprocessing class is defined, which includes methods to read data from a CSV file, 
filter by a specific range of wavelengths (1.04 to 1.43 µm), normalise specific 
columns using z-score and save the resulting DataFrame already processed.
"""

class Preprocessing:                                        
    def __init__(self, min_range=1.04, max_range=1.43):   
            self.min_range = min_range  # O(1)    
            self.max_range = max_range  # O(1)        

    def read_data(self, path):
        return pd.read_csv(path)  # O(n)    

    def wavelength_range_Data(self, Data, col = 'Wavelength'):
        return Data[(Data[col] >= self.min_range) & (Data[col] <= self.max_range)]  # O(n)    

    def normalize_RI(self, Data, columns):
        Data[columns] = Data[columns].transform(lambda x: ((x-x.mean())/x.std()))  # O(n)  
        return Data

    #def pre_Data_SP(path):
    #    unetching = ['RI_Water', 'RI_B', 'RI_C', 'RI_D', 'RI_E', 'RI_F']
    #    Data = read_data(path)
    #    Data = wavelength_range_Data(Data, col = 'Wavelength', min_range = 1.20, max_range = 1.27)
    #    Data = normalize_RI(Data, unetching)
    #    return Data

    def pre_Data_1104(self, path, save_path=None):
        unetching = ['RI_Water', 'RI_B', 'RI_C', 'RI_D', 'RI_E', 'RI_F']  # O(1)
        etching = ['RI_Water_etching', 'RI_B_etching', 'RI_C_etching', 'RI_D_etching', 'RI_E_etching', 'RI_F_etching']  # O(1)
        Data = self.read_data(path)  # O(n)
        Data = self.wavelength_range_Data(Data, col = 'Wavelength')  # O(n)
        Data = self.normalize_RI(Data, unetching + etching)  # O(n)

        # Save the dataset in the range min_range = 1.04 and max_range = 1.43.
        if save_path is None:  # O(1)
            save_path = Path("Data")/"processed"/"Data_processed.csv"   # O(1)

        save_path.parent.mkdir(parents=True, exist_ok=True)
        Data.to_csv(save_path, index=False)  # O(n)   
        print(f'DataFrame saved: {save_path}')  # O(1)

        return Data  # O(1)

    # This code has a computational time complexity of O(n)