#!/usr/bin/env python

import math
import metagenome
from ipyTools import *

class QC:
    def __init__(self, aID):
        self.metagenome = metagenome.Metagenome(aID, True, True)
        self.drisee     = Drisee(self.metagenome)
        self.kmer       = Kmer(self.metagenome)
        self.bp_histo   = NucleoProfile(self.metagenome)

class Drisee:
    def __init__(self, mgObj=None, mgData=None):
        data = None
        if mgObj:
            data = self._get_drisee(mgObj)
        if mgData:
            data = mgData
        self.summary = data['summary'] if has_profile('summary', data) else None
        self.count   = data['counts'] if has_profile('counts', data) else None
        self.percent = data['percents'] if has_profile('percents', data) else None
        
    def _get_drisee(self, mgObj):
        try:
            return mgObj.stats.qc.drisee
        except:
            return None

    def dump(self, filename=None, type='count'):
        if not filename:
            filename = 'drisee_'+type+'_'+random_str()+'.txt'
        profile = None
        if (type == 'count') and self.count:
            profile = self.count
        elif (type == 'percent') and self.percent:
            profile = self.percent
        else:
            return None
        fhdl = open(filename, 'w')
        fhdl.write("#\t"+"\t".join(profile['columns'])+"\n")
        for row in profile['data']:
            fhdl.write("\t".join(map(str, row))+"\n")
        fhdl.close()
        return filename

    def plot(self):
        if not self.percent:
            return None
        l = self.percent['columns']
        x  = map(lambda y: y[0], self.percent['data'])
        yA = map(lambda y: y[1], self.percent['data'])
        yT = map(lambda y: y[2], self.percent['data'])
        yC = map(lambda y: y[3], self.percent['data'])
        yG = map(lambda y: y[4], self.percent['data'])
        yN = map(lambda y: y[5], self.percent['data'])
        yX = map(lambda y: y[6], self.percent['data'])
        yTot = map(lambda y: y[7], self.percent['data'])
        FL_PLOT.legendloc = 'nw'
        FL_PLOT.plot_figure([x,x,x,x,x,x,x],[yA,yT,yC,yG,yN,yX,yTot],label=l)

class NucleoProfile:
    def __init__(self, mgObj):
        data = self._get_bp_profile(mgObj)
        self.count   = data['counts'] if has_profile('counts', data) else None
        self.percent = data['percents'] if has_profile('percents', data) else None

    def _get_bp_profile(self, mgObj):
        try:
            return mgObj.stats.qc.bp_profile
        except:
            return None

    def plot(self):
        if not self.percent:
            return None
        l = self.percent['columns']
        x  = map(lambda y: y[0], self.percent['data'])
        yA = map(lambda y: y[1], self.percent['data'])
        yT = map(lambda y: y[2], self.percent['data'])
        yC = map(lambda y: y[3], self.percent['data'])
        yG = map(lambda y: y[4], self.percent['data'])
        yN = map(lambda y: y[5], self.percent['data'])
        FL_PLOT.legendloc = 'se'
        FL_PLOT.plot_figure([x,x,x,x,x],[yA,yT,yC,yG,yN],label=l)

class Kmer:
    def __init__(self, mgObj):
        self.profile = self._get_kmer(mgObj)

    def _get_kmer(self, mgObj):
        try:
            return mgObj.stats.qc.kmer['15_mer']
        except:
            try:
                return mgObj.stats.qc.kmer['6_mer']
            except:
                return None

    def plot_abundance(self):
        if not (self.profile and ('data' in self.profile)):
            return None
        x = map(lambda z: math.log(z[3], 10), self.profile['data'])
        y = map(lambda z: math.log(z[0], 10), self.profile['data'])
        FL_PLOT.legendloc = 'sw'
        FL_PLOT.plot_figure(x,y,label='kmer rank abundance')

    def plot_ranked(self):
        if not (self.profile and ('data' in self.profile)):
            return None
        x = map(lambda z: math.log(z[3], 10), self.profile['data'])
        y = map(lambda z: 1 - (1.0 * z[5]), self.profile['data'])
        FL_PLOT.legendloc = 'sw'
        FL_PLOT.plot_figure(x,y,label='ranked kmer consumed')

    def plot_spectrum(self):
        if not (self.profile and ('data' in self.profile)):
            return None
        x = map(lambda z: math.log(z[0], 10), self.profile['data'])
        y = map(lambda z: math.log(z[1], 10), self.profile['data'])
        FL_PLOT.legendloc = 'sw'
        FL_PLOT.plot_figure(x,y,label='kmer spectrum')

def has_profile(profile, data):
    if data and (profile in data) and ('data' in data[profile]) and (len(data[profile]['data']) > 0):
        return True
    else:
        return False

def merge_drisee_profile(qc_set, profile='count'):
    if profile == 'count':
        columns = qc_set[0].drisee.count['columns']
        colMax  = len(columns)
        rowMax  = max([ len(x.drisee.count['rows']) for x in qc_set if x.drisee.count ])
        rows    = map(lambda x: x+1, range(rowMax))
        mMatrix = [[0 for i in range(colMax)] for j in range(rowMax)]
        for qc in qc_set:
            if not qc.drisee.count:
                continue
            curLen = len(qc.drisee.count['rows'])
            for r in range(rowMax):
                if r == curLen:
                    break
                for c in range(colMax):
                    mMatrix[r][c] += qc.drisee.count['data'][r][c]
        mData = {'count_profile': {'rows': rows, 'columns': columns, 'data': mMatrix}}
        return Drisee(mgData=mData)
    elif profile == 'percent':
        columns = qc_set[0].drisee.percent['columns']
        colMax  = len(columns)
        rowMax  = max([ len(x.drisee.percent['rows']) for x in qc_set if x.drisee.percent ])
        rows    = map(lambda x: x+51, range(rowMax))
        rowNums = [0 for i in range(rowMax)]
        mMatrix = [[0 for i in range(colMax)] for j in range(rowMax)]
        for qc in qc_set:
            if not qc.drisee.percent:
                continue
            curLen = len(qc.drisee.percent['rows'])
            for r in range(rowMax):
                if r == curLen:
                    break
                rowNums[r] += 1
                for c in range(colMax):
                    mMatrix[r][c] += qc.drisee.percent['data'][r][c]
        for r in range(rowMax):
            for c in range(colMax):
                mMatrix[r][c] = mMatrix[r][c] / rowNums[r]
        mData = {'percent_profile': {'rows': rows, 'columns': columns, 'data': mMatrix}}
        return Drisee(mgData=mData)
    else:
        return None
