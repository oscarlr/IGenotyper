#!/bin/env python
import pandas
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

def plot_histogram(vals,title,plotfn):
    values = pandas.Series(vals)
    fig   = pyplot.figure()
    axes  = values.hist(color='gray', bins=100)
    fig   = pyplot.gcf()
    axes.set_xlabel(title)
    axes.set_ylabel('Count')
    axes.xaxis.grid(False)
    axes.yaxis.grid(False)
    _, max_ = pyplot.ylim()
    fig.savefig(plotfn, format='png')

def plot_barplot(xvals,yvals,xlab,ylab,plotfn):
    fig   = pyplot.figure()
    axes = fig.add_axes([0,0,1,1])
    axes.bar(xvals,yvals)
    fig   = pyplot.gcf()
    axes.set_xlabel(title)
    axes.set_ylabel('Count')
    axes.xaxis.grid(False)
    axes.yaxis.grid(False)
    _, max_ = pyplot.ylim()
    fig.savefig(plotfn, format='png')
    
