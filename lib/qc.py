#!/usr/bin/env python

import os, sys, copy, pprint
import urllib, urllib2, json
import flotplot

API_URL = 'http://api.metagenomics.anl.gov/'

def obj_from_url(url):
    try:
        req = urllib2.Request(url, headers={'Accept': 'application/json'})
        res = urllib2.urlopen(req)
    except urllib2.HTTPError, error:
        print "ERROR: ", error.read()
        return None
    if not res:
        print "ERROR: no results returned"
        return None
    obj = json.loads(res.read())
    if not obj:
        print "ERROR: return structure not valid json format"
        return None
    return obj

class QC:
    def __init__(self, mgid):
        self.mgid = mgid
        self.drisee = Drisee(mgid)
        self.kmer = Kmer(mgid)
        self.bp_histo = NucleoProfile(mgid)

class Drisee:
    def __init__(self, mgid):
        self.mgid = mgid
        self.data = self.get_drisee()
        self.error = self.data['errors'] if self.data and ('error' in self.data) else None
        self.profile = self.data['count_profile'] if self.data and ('count_profile' in self.data) else None
        
    def get_drisee(self):
        return obj_from_url(API_URL+'drisee/'+self.mgid)

    def plot_drisee(self):
        if not (self.data and ('percent_profile' in self.data)):
            return None
        x = self.data['percent_profile']['rows']
        l = self.data['percent_profile']['columns']
        yA = map(lambda y: y[0], self.data['percent_profile']['data'])
        yT = map(lambda y: y[1], self.data['percent_profile']['data'])
        yC = map(lambda y: y[2], self.data['percent_profile']['data'])
        yG = map(lambda y: y[3], self.data['percent_profile']['data'])
        yN = map(lambda y: y[4], self.data['percent_profile']['data'])
        yX = map(lambda y: y[5], self.data['percent_profile']['data'])
        yTot = map(lambda y: y[6], self.data['percent_profile']['data'])
        plt = flotplot.Plot()
        plt.legendloc ='nw'
        return plt.plot_figure([x,x,x,x,x,x,x],[yA,yT,yC,yG,yN,yX,yTot],label=l)

class NucleoProfile:
    def __init__(self, mgid):
        self.mgid = mgid
        self.data = self.get_bp_profile()
        self.profile = self.data['counts'] if self.data and ('counts' in self.data) else None

    def get_bp_profile(self):
        return obj_from_url(API_URL+'bp_histogram/'+self.mgid)

    def plot_bps:
        if not (self.data and ('percents' in self.data)):
            return None
        x = self.data['percents']['rows']
        l = self.data['percents']['columns']
        yA = map(lambda y: y[0], self.data['percents']['data'])
        yT = map(lambda y: y[1], self.data['percents']['data'])
        yC = map(lambda y: y[2], self.data['percents']['data'])
        yG = map(lambda y: y[3], self.data['percents']['data'])
        yN = map(lambda y: y[4], self.data['percents']['data'])
        plt = flotplot.Plot()
        plt.legendloc ='nw'
        return plt.plot_figure([x,x,x,x,x],[yA,yT,yC,yG,yN],label=l)

class Kmer:
    def __init__(self, mgid):
        self.mgid = mgid
        self.data = self.get_kmer()

    def get_kmer(self):
        return obj_from_url(API_URL+'kmer/'+self.mgid)

    def plot_kmer_abundance(self):
        return None

    def plot_kmer_ranked(self):
        return None

    def plot_kmer_spectrum(self):
        return None
    
