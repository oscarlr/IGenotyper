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
    #axes.set_title(title)
    axes.set_xlabel(title)
    axes.set_ylabel('Count')
    axes.xaxis.grid(False)
    axes.yaxis.grid(False)
    #axes.axvline(values.mean(), color='k', linestyle='dashed', linewidth=1)
    _, max_ = pyplot.ylim()
    # axes.text(values.mean() + values.mean()/10,
    #           max_ - max_/10,
    #           'Mean: {:.2f}'.format(values.mean()))
    # axes.text(values.mean() + values.mean()/10,
    #           max_ - (max_/10)*2,
    #           'Std: {:.2f}'.format(values.std()))  
    # axes.text(values.mean() + values.mean()/10,
    #           max_ - (max_/10)*3,
    #           '.5 q: {:.2f}'.format(values.quantile(.5)))  
    # axes.text(values.mean() + values.mean()/10,
    #           max_ - (max_/10)*4,
    #           '.75 q: {:.2f}'.format(values.quantile(.75)))
    fig.savefig(plotfn, format='png')
