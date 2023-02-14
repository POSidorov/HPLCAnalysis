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

import plotly.graph_objects as go
from typing import Optional
from .chromatogram import Chromatogram
import numpy as np


def plot(ch: Chromatogram, wl_index: Optional[int] = None, peaks: bool = False, gaussians: bool = False):
    """Function that creates a figure for the given chromatogram, with information on peaks if required.

        Parameters
        ----------
        ch : Chromatogram
            chromatogram to visualize
        wl_index : int, optional
            if given, visualizes the plot at this specific wavelength. Otherwise, uses the WL of the peaks.
        peaks: bool
            indicate the position of peaks (center, left and right) on the plot
        gaussians: bool
            plot the fitted gaussian for every peak
        cut: bool
            visualize the gaussians with limits put as +-3*sigma
        
        Returns
        -------
        Plotly figure
    """

    # helper function to transform a color from hex to RGB
    def hex_to_rgb(hex_color: str) -> tuple:
        hex_color = hex_color.lstrip("#")
        if len(hex_color) == 3:
            hex_color = hex_color * 2
        return int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    
    if len(ch.peaks) > 0:
        wl_index = optional(wl_index, ch.peaks[0].wl_index)

    fig = go.Figure()
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33']

    fig.add_trace(go.Scatter(x=ch.intensity[wl_index].index, y=ch.intensity[wl_index], 
                             name="chromatogram",
                             line_shape='linear',
                             line_color="black",
                             hovertemplate=None
                             ))
    if peaks:
        # if "peaks" is set to True, will show the position of peaks as well as corresponding minima,
        # and will show the peak number and position (in s) on hover
        for i, p in ch.peaks.items():
            fig.add_trace(go.Scatter(x=[p.center], y=[ch.intensity.loc[p.center,wl_index]], 
                                     name="Peak "+str(i),mode='markers',
                                     marker=dict(
                                     symbol="triangle-down",
                                     size=15,
                                     color=colors[i]  # set color equal to a variable
                                     ),
                                     customdata=np.dstack((["Peak "+str(i)],
                                                       [p.center*ch.info["time_interval"]]
                                                       )),
                                     hovertemplate='<b>%{customdata[0][0]}</b><br>center: %{customdata[0][1]:.1f}s'))
            fig.add_trace(go.Scatter(x=[p.left], y=[ch.intensity.loc[p.left, wl_index]],
                                     name="Peak "+str(i), mode='markers',
                                     marker=dict(
                                     symbol="triangle-up",
                                     size=15,
                                     color=colors[i]  # set color equal to a variable
                                     )))
            fig.add_trace(go.Scatter(x=[p.right], y=[ch.intensity.loc[p.right, wl_index]],
                                     name="Peak "+str(i), mode='markers',
                                     marker=dict(
                                     symbol="triangle-up",
                                     size=15,
                                     color=colors[i]  # set color equal to a variable
                                     )))
    elif gaussians:
        # if "peaks" is set to True, will show the position of peaks as well as corresponding minima,
        # and will show the peak number, position (in s), area and ratio on hover; in addition, will 
        # show the fitted gaussians and baselines for each detected peak.
        for i, p in ch.peaks.items():
            fig.add_trace(go.Scatter(x=[p.center], y=[ch.intensity.loc[p.center,wl_index]], 
                                     name="Peak "+str(i), mode='markers',
                                     marker=dict(
                                     symbol="triangle-down",
                                     size=15,
                                     color=colors[i]  # set color equal to a variable
                                     ),
                                     customdata=np.dstack((["Peak "+str(i)],
                                                       [p.center*ch.info["time_interval"]],
                                                       [p.area],
                                                       [ch.ratios[i]])),
                                     hovertemplate='<b>%{customdata[0][0]}</b><br>center: %{customdata[0][1]:.1f}s<br>area: %{customdata[0][2]:f}<br>ratio: %{customdata[0][3]:.1%}'))
            fig.add_trace(go.Scatter(x=[p.left], y=[ch.intensity.loc[p.left,wl_index]], 
                                     name="Peak "+str(i),mode='markers',
                                     marker=dict(
                                     symbol="triangle-up",
                                     size=15,
                                     color=colors[i]  # set color equal to a variable
                                     )))
            fig.add_trace(go.Scatter(x=[p.right], y=[ch.intensity.loc[p.right,wl_index]], 
                                     name="Peak "+str(i), mode='markers',
                                     marker=dict(
                                     symbol="triangle-up",
                                     size=15,
                                     color=colors[i]  # set color equal to a variable
                                     )))

            x = np.linspace(p.left, p.right, len(p.gaussian.data))
            baseline = p.baseline.eval(x=x)
            approx = list(p.gaussian.eval()+baseline)
            fig.add_trace(go.Scatter(x=list(x)+list(x)[::-1], y=approx+list(baseline)[::1], 
                                     name="",
                                     line_shape='linear',
                                     fill='toself',
                                     fillcolor=f"rgba{(*hex_to_rgb(colors[i]), 0.2)}",
                                     line=dict(color=colors[i], width=2, dash='dash')
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
    # the X axis is rescaled to be shown in seconds
    fig.update_layout(
        xaxis=dict(tickmode='array',
                   tickvals=np.arange(np.ceil(ch.intensity[wl_index].index[0]/500)*500,
                                      np.ceil(ch.intensity[wl_index].index[-1]/500)*500, 500),
                   ticktext=[i*ch.info["time_interval"] for i in
                             np.arange(np.ceil(ch.intensity[wl_index].index[0]/500)*500,
                             np.ceil(ch.intensity[wl_index].index[-1]/500)*500, 500)]
        )
    )
    fig.update_layout(title_text=ch.filename+"<br>Wavelength: " +
                      str(int(ch.info["start_wl"]+wl_index*ch.info["wl_interval"]))+"nm", title_x=0.5)
    return fig


def optional(value, default):
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
