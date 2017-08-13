PROP_WFD = 0
PROP_DD = 1

class VisitPair:
    def __init__(self, ra, dec, rotation):

        self.ra = ra
        self.dec = dec
        self.rotation = rotation

        # these are generated at the time of the exposure since we don't
        # know what exposure time to use until then

        self.visit1 = Visit(PROP_WFD, self, ra, dec, rotation, -1)
        self.visit2 = Visit(PROP_WFD, self, ra, dec, rotation, -1)

    def __repr__(self):
        return "RA: " + str(self.ra) + \
               ", DEC: " + str(self.dec) + \
               ", ROT: " + str(self.rotation)

class Visit:
    # TODO ExpParams?
    def __init__(self, prop, visitPair, ra, dec, rotation, expTime, filter=None):
        self.prop = prop
        # if prop is WFD, this visit will be part of a visitpair
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
