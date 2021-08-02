import pandas as pd
from pandas import DataFrame
import numpy as np
from scipy import signal
from numpy import trapz
from typing import Optional, List, Dict, Tuple

import plotly.graph_objects as go
import plotly.express as px
import plotly

from peak import Peak

class Chromatogram:
    def __init__(self, filename:str, info_rows:int=11):
        self.filename:str = filename
        self.intensity:DataFrame = pd.read_table(filename, delimiter="\t", skiprows=info_rows+1, header=None)
        self.intensity = self.intensity.fillna(0)
        self.peaks:Dict[int, Peak] = {}
        self.ignored_peaks:Dict[int, Peak] = {}
        self.ratios:Dict[int, float] = {}
        self.calibration:Dict[int, float] = {}
        self.base_chromatogram = None
        
        self.info = {}

        with open(filename) as f:
            for i in range(info_rows):
                line = f.readline().strip().split("\t")
                self.info[line[0].lower()] = line[1]
                if line[1].endswith("msec"):
                    self.info[line[0].lower()] = float(line[1][:-4])/1000
                if line[1].endswith("nm"):
                    self.info[line[0].lower()] = float(line[1][:-2])
                if line[1].endswith("min"):
                    self.info[line[0].lower()] = float(line[1][:-3])
                if line[0].endswith("POINTS"):
                    self.info[line[0].lower()] = int(line[1])
                if line[1]=="true" or line[1]=="false":
                    self.info[line[0].lower()] = eval(line[1].capitalize())
                
            
    def set_base_chromatogram(self, base_chromatogram):
        self.base_chromatogram = base_chromatogram
            
    def _optional(self, value, default):
        """Helper function that sets the optional value with given or default.

        Parameters
        ----------
        value : Any
            the argument of the function that is checked against None
        default : Any
            A default value to return if the given argument is None

        Returns
        -------
        Any
            given value, or default value if None is given
        """
        if value is None:
            return default
        else:
            return value
            
    def calculate_peaks(self, wl_index:Optional[int]=None, wl:Optional[int]=None, 
                        height_diff:float=0.3, min_height_diff:float=0.01):
        """Function that scans the chromatogram intensity at given index or wavelength to identify peaks.
        Creates a class atribute "peaks" as a dictionary with the following information:
        - "wl_index": The Wl index at which the peaks were found;
        - "peaks": centers of the peaks;
        - "left" and "right": local minima to the left and right side of the peak.

        Parameters
        ----------
        wl_index : int, optional
            if given, will scan the intensity at this index. 0 by default
        wl : int, optinal
            if given, calculates the index to scan instead
        height_diff: float
            determines the minimum height of peaks to detect, related to the highest intensity

        Returns
        -------
        None
        """
        wl_index = self._optional(wl_index, 0)
        if wl is not None:
            wl_index = int((wl-self.start_wl)/self.wl_interval)
            
        if self.base_chromatogram is None:
        #
            peaks, properties = signal.find_peaks(self.intensity[wl_index],
                                                  prominence=self.intensity[wl_index].max()*height_diff)
            mins = signal.argrelextrema(np.array(self.intensity[wl_index]), comparator=np.less_equal, order=10)[0]
            right_sides = np.searchsorted(mins, peaks, side="right")
            for i in range(len(peaks)):
                self.peaks[i] = Peak(wl_index, list(self.intensity[wl_index].index[peaks])[i],
                                list(self.intensity[wl_index].index[mins[right_sides-1]])[i],
                                list(self.intensity[wl_index].index[mins[right_sides]])[i])
        else:
            peaks, properties = signal.find_peaks(self.intensity[wl_index],
                                                  prominence=self.intensity[wl_index].max()*min_height_diff)
                
            mins = signal.argrelextrema(np.array(self.intensity[wl_index]), comparator=np.less_equal, order=10)[0]
            right_sides = np.searchsorted(mins, peaks, side="right")
            
            assigned_peaks = []
            for i in range(len(peaks)):
                peak = Peak(wl_index, list(self.intensity[wl_index].index[peaks])[i],
                                list(self.intensity[wl_index].index[mins[right_sides-1]])[i],
                                list(self.intensity[wl_index].index[mins[right_sides]])[i])
                for j in self.base_chromatogram.peaks.keys():
                    if j not in assigned_peaks:
                        if peak.intersects(self.base_chromatogram.peaks[j]):
                            self.peaks[j] = peak
                            assigned_peaks.append(j)
                            break
        
    def ignore_peak(self, peak_index):
        """Function that removes the peak from the list of peaks.

        Parameters
        ----------
        peak_index : int
            Removes peak information (center, left, and right) at the given position (index)

        Returns
        -------
        None
        """
        self.ignored_peaks[peak_index] = self.peaks.pop(peak_index)
            
    def fit_all_peaks(self):
        """Function that fits a baseline and a Gaussian into every peak.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        for ind, peak in self.peaks.items():
            peak.fit_baseline(self.intensity)
            peak.fit_gaussian(self.intensity)
        
    def cut_time(self, before:Optional[int]=None, after:Optional[int]=None, unit:str="s"):
        """Function that cuts the rows of the intensity table before and/or after a given point.

        Parameters
        ----------
        before : int, optional
            if given, remove any points before this one
        after : int, optional
            if given, remove any points after this one
        unit : str, "s" or "p"
            if "s", will consider time in seconds, thus recalculating the position with the Chromatogram's time interval
            otherwise, considers only indices in the table

        Returns
        -------
        None
        """
        
        if before and unit=="s":
            before = int(before/self.info["time_interval"])
        if not before:
            before = 0
        if after and unit == "s":
            after = int(after/self.info["time_interval"])
        if not after:
            after = self.intensity.shape[0]
        self.intensity = self.intensity.loc[range(before, after)]
        
    def calibrate(self, calibration_values:Dict[int,float]):
        for k, v in calibration_values.items():
            if type(k) is tuple:
                coefficient = v/100./sum([self.ratios[i] for i in k])
                for i in k:
                    self.calibration[i] = coefficient
            else:
                coefficient = v/100./self.ratios[k]
                self.calibration[k] = coefficient
                
    def calculate_calibrated_ratios(self, external_calibration_values:Optional[Dict[int,float]]=None):
        calibrated_ratios = []
        if external_calibration_values is not None:
            ratio_sum = 0
            for ind, ratio in self.ratios.items():
                self.ratios[ind] = ratio*external_calibration_values[ind]
                ratio_sum += ratio*external_calibration_values[ind]
            if ratio_sum>1:
                for ind in self.ratios.keys():
                    self.ratios[ind] = self.ratios[ind]/ratio_sum
        else:
            
            for ind, coef in self.calibration.items():
                self.ratios[ind] = self.ratios[ind]*coef
    
    def calculate_areas(self, cut:bool=False)->float:
        """Function that calculates the area of all peaks using their Gaussian functions. Uses numpy.trapz.

        Parameters
        ----------
        cut: bool
            if given, will cut the Gaussian at +-3*sigma. False by default
        Returns
        -------
        List(float)
            the areas of all the peaks
        """
        areas = []
        for ind, peak in self.peaks.items():
            peak.calculate_area(self.intensity, self.info["time_interval"], cut)
            areas.append(peak.area)
        return areas
    
    def calculate_ratios(self, cut:bool=False)->List[float]:
        """Function that calculates the area ratios of all peaks at a given WL index using the Gaussian function. 

        Parameters
        ----------
        cut: bool
            if given, will cut the Gaussian at +-3*sigma. False by default
        Returns
        -------
        List[float]
            the ratios of areas of all peaks
        """
        areas = self.calculate_areas()
        for i in self.peaks.keys():
            self.ratios[i] = self.peaks[i].area/sum(areas)
        return self.ratios
    
    def calculate_er(self, enantiomers:Tuple[int])->List[float]:
        return [self.peaks[i].area/sum([self.peaks[j].area for j in enantiomers]) for i in enantiomers]

    def calculate_conversion(self, reactants_peak:int)->float:
        return (sum([self.peaks[i].area for i in len(self.peaks)])-self.peaks[reactants_peak].area)/sum([self.peaks[i].area for i in len(self.peaks)])

    def plot(self, wl_index:Optional[int]=None, peaks:bool=False, gaussians:bool=False, cut:bool=False):
        def hex_to_rgb(hex_color: str) -> tuple:
            hex_color = hex_color.lstrip("#")
            if len(hex_color) == 3:
                hex_color = hex_color * 2
            return int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
        if len(self.peaks)>0:
            wl_index = self._optional(wl_index, self.peaks[0].wl_index)
        fig = go.Figure()

        colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33']

        fig.add_trace(go.Scatter(x=self.intensity[wl_index].index, y=self.intensity[wl_index], 
                                 name="chromatogram",
                                 line_shape='linear',
                                 line_color="black",
                                 hovertemplate=None
                                    ))
        if peaks:
            for i, p in self.peaks.items():
                fig.add_trace(go.Scatter(x=[p.center], y=[self.intensity.loc[p.center,wl_index]], 
                                     name="Peak "+str(i),mode='markers', 
                                     marker=dict(
                                        symbol="triangle-down",
                                        size=15,
                                        color=colors[i] #set color equal to a variable
                                     ),
                                     customdata=np.dstack((["Peak "+str(i)],
                                                           [p.center*self.info["time_interval"]],
                                                           [p.area],
                                                           [self.ratios[i]])),
                                     hovertemplate='<b>%{customdata[0][0]}</b><br>center: %{customdata[0][1]:.1f}s<br>area: %{customdata[0][2]:f}<br>ratio: %{customdata[0][3]:.1%}'))
                fig.add_trace(go.Scatter(x=[p.left], y=[self.intensity.loc[p.left,wl_index]], 
                                     name="Peak "+str(i),mode='markers', 
                                     marker=dict(
                                        symbol="triangle-up",
                                        size=15,
                                        color=colors[i] #set color equal to a variable
                                     )))
                fig.add_trace(go.Scatter(x=[p.right], y=[self.intensity.loc[p.right,wl_index]], 
                                     name="Peak "+str(i),mode='markers', 
                                     marker=dict(
                                        symbol="triangle-up",
                                        size=15,
                                        color=colors[i] #set color equal to a variable
                                     )))
        if gaussians:
            for i, p in self.peaks.items():
                x = p.gaussian.data.index
                baseline = p.baseline.eval(x=x)
                approx = list(p.gaussian.eval()+baseline)

                fig.add_trace(go.Scatter(x=list(x)+list(x)[::-1], y=approx+list(baseline), 
                                     name="",
                                     line_shape='linear',
                                     fill='toself',
                                     fillcolor=f"rgba{(*hex_to_rgb(colors[i]), 0.2)}",
                                     line = dict(color=colors[i], width=2, dash='dash')
                                     ))

        fig.update_layout(
                plot_bgcolor='rgba(0,0,0,0)',
                xaxis_title="Retention time (s)",
                yaxis_title="Intensity",
                showlegend=False,
                autosize=False,
                width=800,
                height=500,
                margin=dict(l=10, r=10, t=50, b=10),
                paper_bgcolor="white"
        )
        fig.update_xaxes(showline=True, gridwidth=0.5, gridcolor='lightgrey', 
                            zeroline=True, zerolinewidth=0.5, zerolinecolor='black')
        fig.update_yaxes(showline=True, gridwidth=0.5, gridcolor='lightgrey',
                            zeroline=True, zerolinewidth=0.5, zerolinecolor='black')
        fig.update_layout(
            xaxis = dict(
                tickmode = 'array',
                tickvals = np.arange(np.ceil(self.intensity[wl_index].index[0]/500)*500,
                                       np.ceil(self.intensity[wl_index].index[-1]/500)*500,500),
                ticktext = [i*self.info["time_interval"] for i in np.arange(np.ceil(self.intensity[wl_index].index[0]/500)*500,
                                       np.ceil(self.intensity[wl_index].index[-1]/500)*500,500)]
            )
        )
        fig.update_layout(title_text=self.filename+"<br>Wavelength: "+str(int(self.info["start_wl"]+wl_index*self.info["wl_interval"]))+"nm", title_x=0.5)

        return fig