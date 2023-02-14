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

from pandas import DataFrame
from numpy import trapz
from typing import Optional
from lmfit.models import GaussianModel, LinearModel, SkewedGaussianModel


class Peak:
    def __init__(self, wl_index: int, center: int, left: int, right: int):
        self.wl_index: int = wl_index
        self.center: int = center
        self.left: int = left
        self.right: int = right
        self.gaussian_left: int = left
        self.gaussian_right: int = right
        self.gaussian = None
        self.baseline = None
        self.area: float = 0
        self.ignored = False
        
    def fit_baseline(self, intensity: DataFrame):
        lmodel = LinearModel()
        x = list(intensity.loc[self.left:self.right, self.wl_index].index)
        y = list(intensity.loc[self.left:self.right, self.wl_index])
        self.baseline = lmodel.fit([y[0], y[-1]], x=[x[0],x[-1]])
        
    def fit_gaussian(self, intensity: DataFrame):
        # gmodel = GaussianModel()
        # params = gmodel.make_params(center=self.center,
        #                            amplitude=intensity.loc[self.center,self.wl_index], 
        #                            sigma=5)
        gmodel = SkewedGaussianModel()
        params = gmodel.make_params(center=self.center, 
                                    amplitude=intensity.loc[self.center,self.wl_index], 
                                    sigma=5, gamma=0.1)
        x = intensity.loc[self.left:self.right, self.wl_index].index
        y = intensity.loc[self.left:self.right, self.wl_index] - self.baseline.eval(x=x)
        self.gaussian = gmodel.fit(y, params, x=x)
        self.gaussian_left = int(self.gaussian.params["center"].value - 3*self.gaussian.params["sigma"].value)
        self.gaussian_right = int(self.gaussian.params["center"].value + 3*self.gaussian.params["sigma"].value)

    def calculate_area(self,  time_interval: float = 1, cut: bool = False, intensity: Optional[DataFrame] = None):
        if cut:
            # x = self.gaussian.data.index
            # gauss_left = list(x).index(int(self.gaussian.params["center"].value -
            #                               3*self.gaussian.params["sigma"].value))
            # gauss_right = list(x).index(int(self.gaussian.params["center"].value +
            #                               3*self.gaussian.params["sigma"].value))
            approx = self.gaussian.eval()[self.left:self.right]
        else:
            approx = self.gaussian.eval()
        self.area = trapz(approx, dx=time_interval)
    
    def intersects(self, other) -> bool:
        if other.gaussian is not None:
            return bool(set(range(self.left, self.right)).intersection(set(range(other.gaussian_left,
                                                                                 other.gaussian_right))))
        return bool(set(range(self.left, self.right)).intersection(set(range(other.left, other.right))))
    
    def __str__(self):
        return "Peak: center {}, right {}, left {}".format(self.center, self.right, self.left)
