# -*- coding: utf-8 -*-
#
#  Copyright 2021 Pavel Sidorov <pavel.o.sidorov@gmail.com>
#  This file is part of HPLCAnalysis project.
#
#  HPLCAnalysis is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.

import pandas as pd
from pandas import DataFrame
import numpy as np
from scipy import signal
from numpy import trapz
from typing import Optional, List, Dict, Tuple

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
                
            
    def set_base_chromatogram(self, base_chromatogram:"Chromatogram"):
        """Helper function that sets the base chromatogram. This is used for peak detection and calibration.

        Parameters
        ----------
        base_chromatogram : Chromatogram
            the chromatogram to consider as the base for this instance

        Returns
        -------
        None
        """
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
        min_height_diff: float
            determines the minimum height of peaks to detect, in the case when the base chromatogram is set

        Returns
        -------
        None
        """
        wl_index = self._optional(wl_index, 0)
        if wl is not None:
            wl_index = int((wl-self.info["start_wl"])/self.info["wl_interval"])
            
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
             # The detection of peaks in case where base chromatogram is set
            maxes, properties = signal.find_peaks(self.intensity[wl_index],
                                                  prominence=self.intensity[wl_index].max()*min_height_diff)
                
            mins = signal.argrelextrema(np.array(self.intensity[wl_index]), comparator=np.less_equal, order=10)[0]
            right_sides = np.searchsorted(mins, maxes, side="right")
            
            peaks = [Peak(wl_index, list(self.intensity[wl_index].index[maxes])[i],
                                list(self.intensity[wl_index].index[mins[right_sides-1]])[i],
                                list(self.intensity[wl_index].index[mins[right_sides]])[i]) for i in range(len(maxes))]
            
            # First, the loop goes through all peaks in the base chromatogram.
            # Every peak that can be uniquely assigned, is assigned.
            assigned_peaks = []
            for j in self.base_chromatogram.peaks.keys():
                intersections = [i for i in range(len(peaks)) if peaks[i].intersects(self.base_chromatogram.peaks[j])]
                if j not in assigned_peaks and len(intersections)==1:
                    self.peaks[j] = peaks[intersections[0]]
                    assigned_peaks.append(j)
            # Next, if any unassigned peaks are still left in this chromatogram, 
            # they will be assigned to the peaks of the base chromatogram that are left       
            for j in self.base_chromatogram.peaks.keys():
                intersections = [i for i in range(len(peaks)) if peaks[i].intersects(self.base_chromatogram.peaks[j])]
                if len(intersections)>1:
                    for i in intersections:
                        if i not in assigned_peaks:
                            self.peaks[j] = peaks[i]
                            assigned_peaks.append(j)
        
    def ignore_peak(self, peak_index:int):
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
        """Function that caluculated the calibration coefficients given the values to which it should calibrate.

        Parameters
        ----------
        calibration_values : Dict[int,float]
            The values of ratios to which the actual ratios should be calibrated.
            If the sum of several peaks' ratios should be calibrated to a specific value,
            the key in the dictionary can be given as a tuple, e.g., (index1, index2,...):value.

        Returns
        -------
        None
        """
        for k, v in calibration_values.items():
            if type(k) is tuple:
                coefficient = v/100./sum([self.ratios[i] for i in k])
                for i in k:
                    self.calibration[i] = coefficient
            else:
                coefficient = v/100./self.ratios[k]
                self.calibration[k] = coefficient
                
    def calculate_calibrated_ratios(self, external_calibration_values:Optional[Dict[int,float]]=None):
        """Function recalculates the ratios either with own calibration coefficients, or with external ones.

        Parameters
        ----------
        external_calibration_values : Dict[int,float], optional
            The external calibration coefficients. If none are given, will use its own, gotten after calibration.

        Returns
        -------
        None
        """
        calibrated_ratios = []
        if external_calibration_values is not None:
            ratio_sum = 0
            for ind, ratio in self.ratios.items():
                self.ratios[ind] = ratio*external_calibration_values[ind]
                ratio_sum += ratio*external_calibration_values[ind]
            if ratio_sum!=1:
                for ind in self.ratios.keys():
                    self.ratios[ind] = self.ratios[ind]/ratio_sum
        else:
            ratio_sum = 0
            for ind, coef in self.calibration.items():
                self.ratios[ind] = self.ratios[ind]*coef
                ratio_sum += self.ratios[ind]
            if ratio_sum!=1:
                for ind in self.ratios.keys():
                    self.ratios[ind] = self.ratios[ind]/ratio_sum
    
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
        """Function that calculates the enantiomeric ratios based on the peak area ratios given as arguments.
        er(X) = ratio(X)/sum(ratio(enantiomers))

        Parameters
        ----------
        enantiomers: Tuple[int]
            indices of peaks
        Returns
        -------
        List[float]
            the er value
        """
        return [self.ratios[i]/sum([self.ratios[j] for j in enantiomers]) for i in enantiomers]
    
    def calculate_conversion(self, to_exclude:int)->float:
        """Function that calculates the conversion by excluding certain peaks (e.g. excluding substrate)
        conversion(X) = (sum_ratios - ratio(X))/sum_ratios
        It assumes the sum of all ratios is 1.

        Parameters
        ----------
        to exclude: int or Tuple[int]
            indices of peaks which should be exluded from calculation
        Returns
        -------
        float
            the conversion value
        """
        if type(to_exclude)==tuple:
            sum_to_exclude = 0
            for e in to_exclude:
                if e in self.peaks.keys():
                    sum_to_exclude += self.ratios[to_exclude]
            return 1. - sum_to_exclude
        else:
            if to_exclude not in self.peaks.keys():
                return 1.
            return 1.-self.ratios[to_exclude]
