import numpy as np
import csv
import config
import numpy.lib.recfunctions as rf
from Visit import Visit, PROP_DD

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
    #seq_exptimes : exposure time for each filter (list)
    # date_obs : year where this sequence should be executed
    def __init__(self,seq_filters,seq_visits,year_obs,telescope):
        self.filters=seq_filters
        self.visits=seq_visits
        self.years=year_obs

        teltime=telescope.filterChangeTime+telescope.readoutTime+telescope.settleTime
        self.Nvisits=np.sum(self.visits)
        self.TotExpTime=self.Nvisits*(config.DDExpTime+config.visitOverheadTime)+len(self.filters)*teltime
        
    def Visits(self,ra,dec):
        for i in range(len(self.filters)):
            for j in range(self.visits[i]):
                visit=Visit(PROP_DD, None, ra,dec, 0,config.DDExpTime ,self.filters[i])
                yield visit

class DD:
    # DD class
    # Input : ra, dec and sequences of observations
    # Setters and getters :
    #    - last_night : last night this field was observed
    #    - last_sequence : last sequence observed
    #    - current_sequence: the sequence to be observed
    #   - time_obs: time at which this field should be observed
    #  - ha : hour angle of this field at time_obs
    
    def __init__(self,ra,dec,sequences):
        self.ra=ra
        self.dec=dec
        self.sequences = sequences

        self._last_sequence=-1
        self._current_sequence=-1
        self._last_night=-1
        self._time_obs=-1.
        self._ha=-999.

    @property
    def nseq(self):
        return len(self.sequences)
    
    @property
    def last_night(self):
        return self._last_night

    @last_night.setter
    def last_night(self, value):
        self._last_night= value
        
    @property
    def last_sequence(self):
        return self._last_sequence

    @last_sequence.setter
    def last_sequence(self, value):
        self._last_sequence= value

    @property
    def current_sequence(self):
        return self._current_sequence

    @current_sequence.setter
    def current_sequence(self, value):
        self._current_sequence= value

    @property
    def time_obs(self):
        return self._time_obs

    @time_obs.setter
    def time_obs(self, value):
        self._time_obs= value

    @property
    def ha(self):
        return self._ha

    @ha.setter
    def ha(self, value):
        self._ha= value 

    def __repr__(self):
        return "RA: " + str(self.ra)+\
            ", DEC: " + str(self.dec) + \
            ", last_night: " + str(self._last_night) + \
            ", last sequence: " + str(self._last_sequence) +\
            ", current sequence: " + str(self._current_sequence)+\
            ", time_obs: " + str(self._time_obs) +\
            ", ha: " + str(self._ha)

        
def Load_DD(filename,telescope):
    """ To load the DDF at the beginning of the survey (in NightScheduler.py)
          Input: filename : name of the file
                      telescope: telescope model
          Output: a list of DD (class above)
    """
    fo = open(filename, "r")
    print("DD  file: ", fo.name)
    first_line=fo.readlines(1)[0]
    names=first_line.strip().split(' ')[1:]
    #print(names)
    params=np.loadtxt(filename,dtype={'names': tuple(names),'formats': tuple(['S20']*len(names))})
    #print(params)
    list_dd=[]
    for par in params:
        ra=float(par['ra'])
        dec=float(par['dec'])
        print(len(par))
        nseq=int((len(par)-2)/3)
        sequences=[]
        for i in range(nseq):
            filters=par[3*i+2].decode().split(',')
            Nexp=[int(s) for s in par[3*i+3].decode().split(',')]
            years=int(par[3*i+4])
            print(filters,Nexp,years)
            sequences.append(Sequence(filters,Nexp,years,telescope))
        list_dd.append(DD(ra,dec,sequences))
        
    return list_dd
