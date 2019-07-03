#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 15:36:37 2011
@ author:                  Sat Kumar Tomer
@ author's webpage:        http://civil.iisc.ernet.in/~satkumar/
@ author's email id:       satkumartomer@gmail.com
@ author's website:        www.ambhas.com
A libray with Python functions for calculations of
micrometeorological parameters and some miscellaneous
utilities.
functions:
    pc_bias : percentage bias
    apb :     absolute percent bias
    rmse :    root mean square error
    mae :     mean absolute error
    bias :    bias
    NS :      Nash-Sutcliffe Coefficient
    L:        likelihood estimation
    correlation: correlation

"""

# import required modules
import numpy as np
from random import randrange
import matplotlib.pyplot as plt
from scipy.stats import kendalltau


def filter_nan(s, o):
    """
    this functions removed the data  from simulated and observed data
    whereever the observed data contains nan

    this is used by all other functions, otherwise they will produce nan as
    output
    """
    if np.sum(~np.isnan(s * o)) >= 1:
        data = np.array([s.flatten(), o.flatten()])
        data = np.transpose(data)
        data = data[~np.isnan(data).any(1)]
        s = data[:, 0]
        o = data[:, 1]
    return s, o


def kendalltau_nan(s, o):
    """
    kendall's tau
    input:
        s: simulated
        o: observed
    output:
        kendall's tau
    """
    s, o = filter_nan(s, o)
    return kendalltau(s, o)[0]


def pc_bias(s, o):
    """
    Percent Bias
    input:
        s: simulated
        o: observed
    output:
        pc_bias: percent bias
    """
    s, o = filter_nan(s, o)
    return 100.0 * sum(s - o) / sum(o)


def apb(s, o):
    """
    Absolute Percent Bias
    input:
        s: simulated
        o: observed
    output:
        apb_bias: absolute percent bias
    """
    s, o = filter_nan(s, o)
    return 100.0 * sum(abs(s - o)) / sum(o)


def rmse(s, o):
    """
    Root Mean Squared Error
    input:
        s: simulated
        o: observed
    output:
        rmses: root mean squared error
    """
    s, o = filter_nan(s, o)
    return np.sqrt(np.mean((s - o) ** 2))


def mae(s, o):
    """
    Mean Absolute Error
    input:
        s: simulated
        o: observed
    output:
        maes: mean absolute error
    """
    s, o = filter_nan(s, o)
    return np.mean(abs(s - o))


def bias(s, o):
    """
    Bias
    input:
        s: simulated
        o: observed
    output:
        bias: bias
    """
    s, o = filter_nan(s, o)
    return np.mean(s - o)


def NS(s, o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
        s: simulated
        o: observed
    output:
        ns: Nash Sutcliffe efficient coefficient
    """
    s, o = filter_nan(s, o)
    return 1 - sum((s - o) ** 2) / sum((o - np.mean(o)) ** 2)


def L(s, o, N=5):
    """
    Likelihood
    input:
        s: simulated
        o: observed
    output:
        L: likelihood
    """
    s, o = filter_nan(s, o)
    return np.exp(-N * sum((s - o) ** 2) / sum((o - np.mean(o)) ** 2))


def correlation(s, o):
    """
    correlation coefficient
    input:
        s: simulated
        o: observed
    output:
        correlation: correlation coefficient
    """
    s, o = filter_nan(s, o)
    if s.size == 0:
        corr = np.NaN
    else:
        corr = np.corrcoef(o, s)[0, 1]

    return corr


def index_agreement(s, o):
    """
	index of agreement

	Willmott (1981, 1982)
	input:
        s: simulated
        o: observed
    output:
        ia: index of agreement
    """
    s, o = filter_nan(s, o)
    ia = 1 - (np.sum((o - s) ** 2)) / (np.sum(
        (np.abs(s - np.mean(o)) + np.abs(o - np.mean(o))) ** 2))
    return ia


def agreement_coefficient(s, o):
    """
    agreement coefficient

    An Agreement Coefficient for Image Comparison by Lei Ji and Kevin Gallo
    input:
        s: simulated
        o: observed
    output:
        ac: agreement coefficient
    """
    s, o = filter_nan(s, o)
    sbar = np.mean(s)
    obar = np.mean(o)
    ac = 1 - (np.sum((s - o) ** 2)) / np.sum(
        (np.abs(sbar - obar) + np.abs(s - sbar)) * (
        np.abs(sbar - obar) + np.abs(o - obar)))
    return ac


def KGE(s, o):
    """
    Kling-Gupta Efficiency
    input:
        s: simulated
        o: observed
    output:
        kge: Kling-Gupta Efficiency
        cc: correlation
        alpha: ratio of the standard deviation
        beta: ratio of the mean
    """
    s, o = filter_nan(s, o)
    cc = correlation(s, o)
    alpha = np.std(s) / np.std(o)
    beta = np.sum(s) / np.sum(o)
    kge = 1 - np.sqrt((cc - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
    return kge, cc, alpha, beta


def assimilation_eff(assimilated, simulated, observed):
    """
    Assimilation efficiency (Aubert et al., 2003)
    Input:
        assimilated: assimilated flow
        simulated: simulated flow
        observed: observed flow
    Output:
        Eff
    """
    s, o = filter_nan(simulated, observed)
    a, o = filter_nan(assimilated, observed)

    Eff = 100 * (1 - np.sum((a - o) ** 2) / np.sum((s - o) ** 2))
    return Eff