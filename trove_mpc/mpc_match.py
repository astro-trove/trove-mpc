"""
The match object that the MPC check returns
"""

from dataclasses import dataclass
from typing import Optional
from astropy.coordinates import SkyCoord

@dataclass
class Match(object):
    target_coord:SkyCoord
    match_coord:SkyCoord
    distance:float
    match_name:Optional[str] = ""

    def todict(self):
        return dict(
            transient_coord = dict(
                ra = self.target_coord.ra.deg,
                dec = self.target_coord.dec.deg
            ),
            match_coord = dict(
                ra = self.match_coord.ra.deg,
                dec = self.match_coord.dec.deg
            ),
            distance = self.distance,
            match_name = self.match_name
        )
            
        
