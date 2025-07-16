"""
The MPC class
"""
from .mpc_match import Match
import os
import time
import kete
import pandas as pd
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u

class MPC:

    def __init__(self):
        self.data = self.get_mpc_data()
    
    def get_mpc_data(self, file_path:str="MPCORB.csv"):
        """Download the MPCORB data if not already downloaded or outdated."""
        url = "https://minorplanetcenter.net/Extended_Files/mpcorb_extended.json.gz"
        if not os.path.exists(file_path) or (time.time() - os.path.getmtime(file_path) > 24 * 60 * 60):
            print("Downloading MPCORB.csv...")
            mpcorb_df = kete.mpc.fetch_known_orbit_data(url=url)
            mpcorb_df.to_csv(file_path)
            return mpcorb_df
        return pd.read_csv(file_path)

    def _kete_state_ra_dec(self, state, sun2earth):
        """
        Private method to wrap computing the coordinates of a single state
        """
        obj_earth_pos = state.pos - sun2earth 
        state_vec_equitorial = obj_earth_pos.change_frame(kete.vector.Frames.Equatorial)
        return state_vec_equitorial.ra, state_vec_equitorial.dec
    
    def minor_planet_filter(
            self,
            s_ra:float,
            s_dec:float,
            obsmjd:float,
            filter_radius:float,
            first_cut_radius:float = 10
    ) -> Match|None:
        """
        Match to the minor planet catalog

        Args:
            orbit_catalog (pd.DataFrame) : The MPCORB catalog downloaded by kete.mpc.fetch_known_orbit_data
            s_ra (float) : The RA of the target (or "source") that we want to crossmatch with the MPC
            s_dec (float) : The Dec of the target (or "source") that we want to crossmatch with the MPC
            obsmjd (float) : The MJD to calculate the ephemeris at
            filter_radius (float) : the filter radius in arcsec
        Return:
            Either a Match object if there is match or None if there is no match
        """

        orbit_catalog = self.data
        
        # some variables we need for computing the RA/Dec
        jd = kete.Time.from_mjd(obsmjd).jd
        target_skycoord = SkyCoord(s_ra, s_dec, unit="deg")
        sun2earth = kete.spice.get_state("Earth", jd).pos

        # convert the MPCORB orbits table to kete State objects
        objs = kete.mpc.table_to_states(orbit_catalog)

        # Propagate the objects
        ### first use the a 2 body approximation
        states = kete.propagate_two_body(objs, jd)

        ras, decs = np.array(
            [self._kete_state_ra_dec(state, sun2earth) for state in states]
        ).T
        coord_is_na = np.isnan(ras) + np.isnan(decs)
        mpc_skycoords_approx = SkyCoord(ras, decs, unit="deg")[~coord_is_na]
        idxs = np.where(target_skycoord.separation(mpc_skycoords_approx) < first_cut_radius*u.deg)[0]

        ### now use the full n-body solution for all of these objects
        objs_win_first_cut = np.array(objs)[idxs]
        print(f"Found {len(objs_win_first_cut)} w/in {first_cut_radius}deg, running full n-body on those...")

        states = kete.propagate_n_body(objs_win_first_cut, jd)
        ras, decs = np.array([self._kete_state_ra_dec(state, sun2earth) for state in states]).T
        coord_is_na = np.isnan(ras) + np.isnan(decs)
        mpc_skycoords = SkyCoord(ras, decs, unit="deg")[~coord_is_na]

        # find the separation
        mpc_catalog_idx, sep, _ = target_skycoord.match_to_catalog_sky(mpc_skycoords)

        diff = np.array([c1.separation(c2).deg for c1,c2 in zip(mpc_skycoords, mpc_skycoords_approx[idxs])])
        print(f"The fraction of objects within {first_cut_radius}deg is {len(diff[diff < first_cut_radius])/len(diff) * 100}")
        if sep.arcsec < filter_radius:
            return Match(
                target_coord=target_skycoord,
                match_coord=mpc_skycoords[mpc_catalog_idx],
                distance=sep.arcsec,
                match_name=objs_win_first_cut[mpc_catalog_idx].desig
            )
        return
