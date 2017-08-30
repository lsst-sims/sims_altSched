# constants
PROP_WFD = 0
PROP_DD = 1

class VisitPair:
    """ A class representing a visit/revisit pair """
    def __init__(self, ra, dec, rotation):
        self.ra = ra
        self.dec = dec
        self.rotation = rotation

        # don't assign the exposure time now -- instead leave it up
        # to the scheduler at the time of exposure
        expTime = -1
        self.visit1 = Visit(PROP_WFD, self, ra, dec, rotation, expTime)
        self.visit2 = Visit(PROP_WFD, self, ra, dec, rotation, expTime)

    def __repr__(self):
        return "RA: " + str(self.ra) + \
               ", DEC: " + str(self.dec) + \
               ", ROT: " + str(self.rotation)

class Visit:
    """ A class representing a single visit """
    # TODO could move args into an ExpParams class since there are so many...
    def __init__(self, prop, visitPair, ra, dec, rotation, expTime, filter=None):
        self.prop = prop
        # this can be None if there is no assigned visitPair, e.g. if this
        # is a DD visit
        self.visitPair = visitPair

        self.ra = ra
        self.dec = dec
        self.rotation = rotation
        self.expTime = expTime
        self.filter = filter

        self.isComplete = False
        self.timeOfCompletion = 0

    def __repr__(self):
        return "RA: " + str(self.ra) + ", DEC: " + str(self.dec) + \
                ", ROT: " + str(self.rotation) + ", expTime: " + \
                str(self.expTime) + ", isComplete?: " + str(self.isComplete)
