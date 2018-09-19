from __future__ import division, print_function

import numpy as np
import csv
import config
import math
import numpy.lib.recfunctions as rf
from Visit import Visit, PROP_DD
from lsst.sims.speedObservatory.telescope import Telescope
from lsst.sims.speedObservatory import sky
from config import NORTH, SOUTHEAST

def _loadDDFs(filename="ddFields_new.txt"):
    allDDFs = {}
    with open(filename, "r") as ddfFile:
        reader = csv.DictReader(ddfFile)
        for row in reader:
            name = row["name"]
            ra = float(row["ra"])
            dec = float(row["dec"])
            uFrac = float(row["u"])
            gFrac = float(row["g"])
            rFrac = float(row["r"])
            iFrac = float(row["i"])
            zFrac = float(row["z"])
            yFrac = float(row["y"])
            assert(math.isclose(uFrac + gFrac + rFrac + iFrac + zFrac + yFrac, 1))
            allDDFs[name] = DDF(name, ra, dec, [],
                                uFrac, gFrac, rFrac, iFrac, zFrac, yFrac)
    return allDDFs


class DDScheduler:
    def __init__(self, nightNum, direction, leftOutFilter):
        self.nightNum = nightNum
        self.direction = direction
        # for the record, this is hacky. Things are starting to get complicated
        # enough to justify a Telescope object with state *gasp*
        self.leftOutFilter = leftOutFilter

        self.allFields = _loadDDFs("ddFields_new.txt")

    def telTimeTonight(self):
        # 40 minutes per night for 3/4 of the year
        # TODO magic numbers
        tonightsFields = self._whichDDsTonight()
        if len(tonightsFields) == 0:
            return 0
        else:
            return config.DDTimePerNight

    def _whichDDsTonight(self):
        # divide deep fields into two groups of 2 based on RA, and COSMOS

        # TODO magic numbers
        # for days -120-60, switch between the four non-COSMOS fields
        #if self.nightNum % 365 > 200 and self.nightNum % 365 <= 201: #45
        if 170 < self.nightNum % 365 <= 170+180:
            # switch off between the two pairs every other night
            if self.direction == SOUTHEAST:
                return [self.allFields["unselected"], self.allFields["ELAIS-S1"]]
            else:
                return [self.allFields["XMM-LSS"], self.allFields["Chandra"]]
        elif 20 < self.nightNum % 365 and self.nightNum % 365 <= 20+90:
            # doing COSMOS tonight
            return [self.allFields["COSMOS"]]
        else:
            return []

    def desiredStartTime(self):
        nightStart = sky.nightStart(config.surveyStartTime, self.nightNum)

        tonightsFields = self._whichDDsTonight()
        if len(tonightsFields) == 0:
            # no DD planned for tonight
            return 0

        assert(self.telTimeTonight() > 0)

        # handle the case that the DDfs in the group cluster around
        # ra of 0
        ras = np.array([field.ra for field in tonightsFields])
        if ras.max() - ras.min() < np.pi:
            ra = ras.mean()
        else:
            ra = ((ras + np.pi) % (2*np.pi)).mean() - np.pi

        if self.direction == NORTH:
            # fudge factor to prevent Chandra from hitting the zenith avoid zone
            ra -= np.radians(10)

        # figure out when ra crosses the meridian, and return a desired start
        # equal to that time minus half the DD observing time
        # TODO hack including the telescope's
        raDiff = ra - sky.raOfMeridian(nightStart)
        if raDiff < -1 * np.pi:
            raDiff += 2*np.pi
        if raDiff < 0:
            # we're already too late -- better grab the DD asap
            #print("too late")
            return nightStart - 1
        else:
            # figure out how long until the field reaches the meridian
            timeToTransit = raDiff * 24*3600 / (2*np.pi)
            #print("time to transit", timeToTransit / 60)
            return nightStart + timeToTransit - self.telTimeTonight() / 2

    def schedule(self, startFilter, endFilter):
        tonightsDDs = self._whichDDsTonight()
        visitTime = config.DDExpTime + config.DDOverheadTime
        nVisitsPerField = self.telTimeTonight() / len(tonightsDDs) / visitTime
        for field in tonightsDDs:
            # zero out the filter that we don't have in the changer tonight
            filterFracs = np.array(field.filterFracs)
            filterFracs[Telescope.filterId[self.leftOutFilter]] = 0
            # and renormalize
            filterFracs /= filterFracs.sum()

            nVisitsPerFilter = np.round(filterFracs * nVisitsPerField).astype(np.int)

            # figure out which order to put the filters in, to avoid
            # changing filters when transitioning between WFD and DD
            filterOrder = []
            # abutting filters are the ones that transition to/from WFD
            nonAbuttingFilters = set(Telescope.filters)
            if startFilter is not None:
                assert(startFilter != self.leftOutFilter)
                filterOrder.append(startFilter)
                nonAbuttingFilters.remove(startFilter)

            if endFilter is not None and endFilter != startFilter:
                assert(endFilter != self.leftOutFilter)
                nonAbuttingFilters.remove(endFilter)

            nonAbuttingFilters.remove(self.leftOutFilter)

            for f in list(nonAbuttingFilters):
                filterOrder.append(f)

            if endFilter is not None and endFilter != startFilter:
                filterOrder.append(endFilter)

            for f in filterOrder:
                nVisits = nVisitsPerFilter[Telescope.filterId[f]]
                for i in range(nVisits):
                    visit = Visit(PROP_DD, None, field.ra, field.dec, 0,
                                  config.DDExpTime, filter=f)
                    yield visit

"""
class DDF:
    # DDF class
    # input: csv file with (Ra,Dec)+number of visits
    # output: numpy array with following dtype:
    # Ra,Dec,Ng,Nr,Ni,Nz,Ny,Nvisits,exptime (total: exptime+overheads+telescope included)
    def __init__(self, config_file,telescope):
        r=[]
        names=None
        with open(config_file, "r") as ddFile:
            fields = csv.DictReader(ddFile)
            if names is None:
                names=fields.fieldnames
            
            for field in fields:
                ra=[]
                for key in fields.fieldnames:
                    #print(key,field[key])
                    what=float(field[key])
                    if key[0]=='N':
                        what=int(field[key])
                    ra.append(what)
                r.append(tuple(ra))

        #print(r,names)
        self.DDF=np.rec.fromrecords(r,names=names)
        self.telescope=telescope
        self.exptime()
        
    def exptime(self):
        
        r=[]
        for dd in self.DDF:
            Nvisits=0 #for each of the visit, exp time and overhead 
            teltime=0 #filter change, readout and settletime of the telescope
            for name in dd.dtype.names:
                if name[0] == 'N':
                    Nvisits+=dd[name]
                    teltime+=self.telescope.filterChangeTime+self.telescope.readoutTime+self.telescope.settleTime
            r.append((config.DDExpTime,Nvisits,teltime+Nvisits*(config.DDExpTime+config.visitOverheadTime)))

        myarr=np.rec.fromrecords(r,names=['DDExpTime','Nvisits','TotExpTime'])
        self.DDF=rf.append_fields(self.DDF,['DDExpTime','Nvisits','TotExpTime'],[myarr['DDExpTime'],myarr['Nvisits'],myarr['TotExpTime']])
 """       

class Sequence:
    # A sequence is a list of visits plus years of observation for this sequence
    # Input parameters :
    # seq_filters : sequence of filters (list)
    # seq_exptimes : exposure time for each filter (list)
    # date_obs : year where this sequence should be executed
    def __init__(self, seq_filters, seq_visits, year_obs, telescope):
        self.filters = seq_filters
        self.visits = seq_visits
        self.years = year_obs

        teltime = (telescope.filterChangeTime + telescope.readoutTime +
                   telescope.settleTime)
        self.nVisits = np.sum(self.visits)
        self.totExpTime = ((self.nVisits *
                            (config.DDExpTime + config.visitOverheadTime))
                           + len(self.filters) * teltime)
        
    def schedule(self, ra, dec):
        for i in range(len(self.filters)):
            for j in range(self.visits[i]):
                visit = Visit(PROP_DD, None, ra, dec, 0, config.DDExpTime,
                              self.filters[i])
                yield visit

class DDF:
    # DDF class
    # Input : ra, dec and sequences of observations
    # Setters and getters :
    #    - prevNight : last night this field was observed
    #    - last_sequence : last sequence observed
    #    - current_sequence: the sequence to be observed
    #    - time_obs: time at which this field should be observed
    #    - ha : hour angle of this field at time_obs
    
    def __init__(self, name, ra, dec, sequences,
                 uFrac, gFrac, rFrac, iFrac, zFrac, yFrac):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.sequences = sequences

        self.filterFracs = [uFrac, gFrac, rFrac, iFrac, zFrac, yFrac]

        self._last_sequence = -1
        self._current_sequence = -1
        self._prevNight = -1
        self._time_obs = -1.
        self._ha = -999.

    @property
    def nseq(self):
        return len(self.sequences)
    
    @property
    def prevNight(self):
        return self._prevNight

    @prevNight.setter
    def prevNight(self, value):
        self._prevNight = value
        
    @property
    def last_sequence(self):
        return self._last_sequence

    @last_sequence.setter
    def last_sequence(self, value):
        self._last_sequence = value

    @property
    def current_sequence(self):
        return self._current_sequence

    @current_sequence.setter
    def current_sequence(self, value):
        self._current_sequence = value

    @property
    def time_obs(self):
        return self._time_obs

    @time_obs.setter
    def time_obs(self, value):
        self._time_obs = value

    @property
    def ha(self):
        return self._ha

    @ha.setter
    def ha(self, value):
        self._ha = value 

    def __repr__(self):
        return "RA: " + str(self.ra)+ \
            ", DEC: " + str(self.dec) + \
            ", prevNight: " + str(self._prevNight) + \
            ", last sequence: " + str(self._last_sequence) + \
            ", current sequence: " + str(self._current_sequence)+ \
            ", time_obs: " + str(self._time_obs) + \
            ", ha: " + str(self._ha)


"""
def loadDDFs(filename, telescope):
    #    To load the DDF at the beginning of the survey (in NightScheduler.py)
    #      Input: filename : name of the file
    #             telescope: telescope model
    #      Output: a list of DD (class above)

    fo = open(filename, "r")
    print("DD  file: ", fo.name)
    first_line = fo.readlines(1)[0]
    names = first_line.strip().split(' ')[1:]
    params = np.loadtxt(filename,dtype={'names': tuple(names),'formats': tuple(['S20']*len(names))})

    allDDFs = []
    for par in params:
        ra = float(par['ra'])
        dec = float(par['dec'])
        print(len(par))
        nseq = int((len(par) - 2) / 3)
        sequences = []
        for i in range(nseq):
            filters = par[3 * i + 2].decode().split(',')
            nExp = [int(s) for s in par[3 * i + 3].decode().split(',')]
            years = int(par[3 * i + 4])
            print(filters, nExp, years)
            sequences.append(Sequence(filters, nExp, years, telescope))
        allDDFs.append(DDF(ra, dec, sequences))
        
    return allDDFs
"""
