"""
The transient object to match against
"""
from .mpc import MPC

class Transient:

    def __init__(self, ra:float, dec:float):

        self.ra = ra
        self.dec = dec

    def minor_planet_match(self, discovery_mjd, filter_radius=25.0, **kwargs):
        print("Loading moving object catalog...")
        mpc = MPC()
        print("Moving object catalog loaded...")
        print("Finding moving objects nearby your object...")
        return mpc.minor_planet_filter(self.ra, self.dec, discovery_mjd, filter_radius, **kwargs)
