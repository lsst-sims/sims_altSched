from __future__ import division
import numpy as np

class RotationGenerator:
    """ Class to generate rotator angles """
    def __init__(self, initialValue=-np.pi/2, initialDirection=1):
        self.value = initialValue
        self.direction = initialDirection

    def rotations(self, telescope):
        """ Generator that yields rotations

        The rotations bounce between `minimum` and `maximum`
        """
        minimum = telescope.minRotation
        maximum = telescope.maxRotation

        # don't use a random rotation for each mini because then we'll
        # get long rotational slew times
        # instead, periodically increase and decrease the rotation
        # of each survey between -pi/2 and pi/2 (like pong)

        # TODO is it desirable to never use the same exact rotation twice
        # to maximize the number of angles we use over the whole survey?

        # the very first mini should have rotation -pi/2
        # each subsequent survey adds rotationIncrement until the rotation
        # reaches +pi/2, at which point the rotation decrements with each
        # new mini until it reaches -pi/2. Then the cycle restarts

        # TODO parametrize
        increment = np.pi / 10
 
        while True:
            if self.direction == 1:
                # we've been going up
                if self.value + increment <= maximum:
                    # keep going
                    self.value += increment
                else:
                    # turn around
                    self.value -= increment
                    self.direction = -1
            else:
                # we've been going down
                if self.value - increment >= minimum:
                    # keep going
                    self.value -= increment
                else:
                    # turn around
                    self.value += increment
                    self.direction = 1
            yield self.value

