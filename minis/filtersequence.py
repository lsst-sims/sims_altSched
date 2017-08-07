import Scheduler
from lsst.sims.speedObservatory import Telescope
from lsst.sims.speedObservatory import sky
from Config import NORTH, SOUTHEAST
import Config

# keep track of which filter we should start out in when looking
# in the N, S, and E
_startFilters =   {NORTH:     Telescope.filters[1], 
                   SOUTHEAST: Telescope.filters[1]}
_leftOutFilters = {NORTH:     Telescope.filters[0], 
                   SOUTHEAST: Telescope.filters[0]}

def _fId(f):
    # get the id of a filter such as 'u' or 'r'
    return Telescope.filterId[f]

def _inc(f):
    # helper to increment mod len(Telescope.filters)
    fId = _fId(f)
    return Telescope.filters[(fId + 1) % len(Telescope.filters)]

def _attemptFilter(direction, attempt, backup):
    # helper to choose a filter based on which are available tonight
    if attempt == _leftOutFilters[direction]:
        return backup
    else:
        return attempt

def notifyNightCompleted(direction):
    # set up startFilterId for next night
    # it needs to be incremented by 2, with the leftout filter skipped over
    # if necessary 
    # TODO explain better
    for i in range(2):
        nextFilter = _inc(_startFilters[direction])
        if _startFilters[direction] == _leftOutFilters[direction]:
            _startFilters[direction] = nextFilter
            nextFilter = _inc(nextFilter)
        _startFilters[direction] = nextFilter

    # also update the left-out filter
    _leftOutFilters[direction] = _inc(_leftOutFilters[direction])

def _avoidFilters(time):
    moonPhase = sky.phaseOfMoon(time)
    moonRa, moonDec = sky.radecOfMoon(time)
    moonAlt, _ = sky.radec2altaz(moonRa, moonDec, time)

    uMoonUp = moonAlt > Config.moonUMaxAlt
    uMoonBright = moonPhase > Config.moonUMaxPhase
    gMoonUp = moonAlt > Config.moonGMaxAlt
    gMoonBright = moonPhase > Config.moonGMaxPhase

    avoidFilters = set()
    if uMoonUp and uMoonBright:
        avoidFilters.add("u")
    if gMoonUp and gMoonBright:
        avoidFilters.add("g")

    return avoidFilters

def twoFChanges(nightNum, direction, nScans):
    nightStartTime = sky.nightStart(Config.surveyStartTime, nightNum)
    nightEndTime   = sky.nightEnd(  Config.surveyStartTime, nightNum)

    # build a list `fs` of filters to use for each of `nScans` scans
    f = _startFilters[direction]
    fs = []

    # at the beginning and end of each night, do the y observations
    bookendFilter = _attemptFilter(direction, "y", "z")
    fs.append(bookendFilter)
    fs.append(bookendFilter)

    if f == _leftOutFilters[direction]:
        f = _inc(f)
    if f == "y":
        f = _attemptFilter(direction, "r", "i")
    if f in _avoidFilters(nightStartTime) or f in _avoidFilters(nightEndTime):
        f = _attemptFilter(direction, "z", "y")

    for i in range(int((nScans-4))):
        fs.append(f)

    # put on the bookend id for the last 2 scans of the night
    fs.append(bookendFilter)
    fs.append(bookendFilter)

    return fs

def maxFChanges(nightNum, direction,  nScans):
    # keep track of what time we expect it to be for each scan
    timePerScan = sky.nightLength(Config.surveyStartTime, nightNum) / nScans
    time = sky.twilStart(Config.surveyStartTime, nightNum)

    f = _startFilters[direction]
    fs = []

    # at the beginning and end of each night, do the y observations
    bookendFilter = _attemptFilter(direction, "y", "z")
    fs.append(bookendFilter)
    fs.append(bookendFilter)

    # loop through the non-bookend scans. Divide by 4 since we schedule
    # filters in sets of four scans
    assert(nScans % 4 == 0)
    for i in range(int((nScans-4)/4)):
        if f == _leftOutFilters[direction]:
            f = _inc(f)

        if f in _avoidFilters(time):
            # replace u/g observations when the moon is up with z observations
            # unless z is not in the filter changer, in which case observe
            # in y
            fs.append(_attemptFilter(direction, "z", "y"))
        elif f == "y":
            # replace any non-bookend scans that would be y with r (or i)
            # since we care about r (and i) a lot and y not that much
            fs.append(_attemptFilter(direction, "r", "i"))
        else:
            # just use the next filter in line
            fs.append(f)

        # repeat the same filter for the next three scans
        fs.append(fs[-1])
        fs.append(fs[-1])
        fs.append(fs[-1])

        f = _inc(f)
        time += timePerScan * 4

    # put on the bookend id for the last 2 scans of the night
    fs.append(bookendFilter)
    fs.append(bookendFilter)

    return fs