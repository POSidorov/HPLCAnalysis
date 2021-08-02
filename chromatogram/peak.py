import pandas as pd
from pandas import DataFrame
import numpy as np
from numpy import trapz
from typing import Optional, List, Dict, Tuple

from lmfit.models import GaussianModel, ExponentialGaussianModel, LinearModel
from lmfit.model import ModelResult


class Peak:
    def __init__(self, wl_index:int, center:int, left:int, right:int):
        self.wl_index:int = wl_index
        self.center:int = center
        self.left:int = left
        self.right:int = right
        self.gaussian_left:int = left
        self.gaussian_right:int = right
        self.gaussian:ModelResult = None
        self.baseline:ModelResult = None
        self.area:float = 0
        self.ignored = False
        
    def fit_baseline(self, intensity:DataFrame):
        lmodel = LinearModel()
        x = intensity.loc[self.left:self.right,self.wl_index].index
        y = intensity.loc[self.left:self.right,self.wl_index]
        self.baseline = lmodel.fit([list(y)[0]]+[list(y)[-1]], x=[list(x)[0]]+[list(x)[-1]])
        
    def fit_gaussian(self, intensity:DataFrame):
        #gmodel = GaussianModel()
        #params = gmodel.make_params(center=self.center, 
        #                            amplitude=intensity.loc[self.center,self.wl_index], 
        #                            sigma=5)
        gmodel = SkewedGaussianModel()
        params = gmodel.make_params(center=self.center, 
                                    amplitude=intensity.loc[self.center,self.wl_index], 
                                    sigma=5, gamma=0.1)
        x = intensity.loc[self.left:self.right,self.wl_index].index
        y = intensity.loc[self.left:self.right,self.wl_index] - self.baseline.eval(x = x)
        self.gaussian = gmodel.fit(y, params, x=x)
        self.gaussian_left = int(self.gaussian.params["center"].value - 3*self.gaussian.params["sigma"].value)
        self.gaussian_right = int(self.gaussian.params["center"].value + 3*self.gaussian.params["sigma"].value)
            
        
    def calculate_area(self, intensity:DataFrame, time_interval:float=1, cut:bool=False):
        x = self.gaussian.data.index
        if cut:
            gauss_left = list(x).index(int(self.gaussian.params["center"].value - 3*self.gaussian.params["sigma"].value))
            gauss_right = list(x).index(int(self.gaussian.params["center"].value + 3*self.gaussian.params["sigma"].value))
            approx = self.gaussian.eval()[left:right] 
        else:
            approx = self.gaussian.eval()
        self.area = trapz(approx, dx=time_interval)
    
    def intersects(self, other)->bool:
        if other.gaussian is not None:
            return bool(set(range(self.left, self.right)).intersection(set(range(other.gaussian_left, other.gaussian_right))))
        return bool(set(range(self.left, self.right)).intersection(set(range(other.left, other.right))))
    
    def __str__(self):
        return "Peak: center {}, right {}, left {}".format(self.center, self.right, self.left)