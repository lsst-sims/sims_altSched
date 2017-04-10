class VisitPair:
    def __init__(self, miniSurvey, ra, dec):
        self.miniSurvey = miniSurvey

        self.ra = ra
        self.dec = dec

        # these are generated at the time of the exposure since we don't
        # know what exposure time to use until then

        self.visit1 = Visit(self, ra, dec, miniSurvey.rotation, -1)
        self.visit2 = Visit(self, ra, dec, miniSurvey.rotation, -1)

    def __repr__(self):
        return "RA: " + str(self.ra) + ", DEC: " + str(self.dec)

class Visit:
    # TODO ExpParams?
    def __init__(self, visitPair, ra, dec, rotation, expTime):
        self.visitPair = visitPair

        self.ra = ra
        self.dec = dec
        self.rotation = rotation
        self.expTime = expTime

        self.isComplete = False
        self.timeOfCompletion = 0

    def __repr__(self):
        return "RA: " + str(self.ra) + ", DEC: " + str(self.dec) + \
                ", rot: " + str(self.rotation) + ", expTime: " + \
                str(self.expTime) + ", isComplete?: " + str(self.isComplete)
