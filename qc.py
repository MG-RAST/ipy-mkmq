#!/usr/bin/env python

import math
import flotplot
from ipyTools import *

PLOT = flotplot.Plot()

class QC:
    def __init__(self, aID):
        self.ID       = aID
        self.drisee   = Drisee(self.ID)
        self.kmer     = Kmer(self.ID)
        self.bp_histo = NucleoProfile(self.ID)

class Drisee:
    def __init__(self, aID=None, data=None):
        if aID:
            self.ID   = aID
            self.data = self.get_drisee()
        if data:
            self.data = data
        self.error   = self.data['errors'] if self.data and ('errors' in self.data) else None
        self.count   = self.data['count_profile'] if self.data and ('count_profile' in self.data) else None
        self.percent = self.data['percent_profile'] if self.data and ('percent_profile' in self.data) else None
        if self.count and (not self.percent):
            self.percent = self.count_to_percent()
        
    def get_drisee(self):
        return obj_from_url(API_URL+'drisee/'+self.ID)

    def count_to_percent(self):
        if not self.count:
            return None
        rows = []
        data = []
        for i in range(self.count['rows']):
            if i < 50:
                continue
            cRow  = self.count['data'][i]
            total = sum(cRow)
            percs = map(lambda x: 100 * ((x * 1.0) / total), cRow)
            rows.append(i)
            data.append([ percs[6:11], sum(percs[6:11]) ])
        return {'rows': rows, 'columns': ['A','T','C','G','N','InDel','Total'], 'data': data}

    def plot(self):
        if not self.percent:
            return None
        x = self.percent['rows']
        l = self.percent['columns']
        yA = map(lambda y: y[0], self.percent['data'])
        yT = map(lambda y: y[1], self.percent['data'])
        yC = map(lambda y: y[2], self.percent['data'])
        yG = map(lambda y: y[3], self.percent['data'])
        yN = map(lambda y: y[4], self.percent['data'])
        yX = map(lambda y: y[5], self.percent['data'])
        yTot = map(lambda y: y[6], self.percent['data'])
        PLOT.legendloc = 'nw'
        PLOT.plot_figure([x,x,x,x,x,x,x],[yA,yT,yC,yG,yN,yX,yTot],label=l)

class NucleoProfile:
    def __init__(self, aID):
        self.ID      = aID
        self.data    = self.get_bp_profile()
        self.profile = self.data['counts'] if self.data and ('counts' in self.data) else None

    def get_bp_profile(self):
        return obj_from_url(API_URL+'bp_histogram/'+self.ID)

    def plot(self):
        if not (self.data and ('percents' in self.data)):
            return None
        x = self.data['percents']['rows']
        l = self.data['percents']['columns']
        yA = map(lambda y: y[0], self.data['percents']['data'])
        yT = map(lambda y: y[1], self.data['percents']['data'])
        yC = map(lambda y: y[2], self.data['percents']['data'])
        yG = map(lambda y: y[3], self.data['percents']['data'])
        yN = map(lambda y: y[4], self.data['percents']['data'])
        PLOT.legendloc = 'se'
        PLOT.plot_figure([x,x,x,x,x],[yA,yT,yC,yG,yN],label=l)

class Kmer:
    def __init__(self, aID):
        self.ID   = aID
        self.data = self.get_kmer()

    def get_kmer(self):
        return obj_from_url(API_URL+'kmer/'+self.ID)

    def plot_abundance(self):
        if not (self.data and ('profile' in self.data)):
            return None
        x = map(lambda z: math.log(z[3], 10), self.data['profile'])
        y = map(lambda z: math.log(z[0], 10), self.data['profile'])
        PLOT.legendloc = 'sw'
        PLOT.plot_figure(x,y,label='kmer rank abundance')

    def plot_ranked(self):
        if not (self.data and ('profile' in self.data)):
            return None
        x = map(lambda z: math.log(z[3], 10), self.data['profile'])
        y = map(lambda z: 1 - (1.0 * z[5]), self.data['profile'])
        PLOT.legendloc = 'sw'
        PLOT.plot_figure(x,y,label='ranked kmer consumed')

    def plot_spectrum(self):
        if not (self.data and ('profile' in self.data)):
            return None
        x = map(lambda z: math.log(z[0], 10), self.data['profile'])
        y = map(lambda z: math.log(z[1], 10), self.data['profile'])
        PLOT.legendloc = 'sw'
        PLOT.plot_figure(x,y,label='kmer spectrum')

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
        return Drisee(data=mData)
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
                    mMatrix[r][c] += qc.drisee.count['data'][r][c]
        for r in range(rowMax):
            for c in range(colMax):
                mMatrix[r][c] = mMatrix[r][c] / rowNums[r]
        mData = {'percent_profile': {'rows': rows, 'columns': columns, 'data': mMatrix}}
        return Drisee(data=mData)
    else:
        return None
