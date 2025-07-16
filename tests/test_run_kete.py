"""
Some tests of the functions for MPC checking
"""
import time
import json
import pandas as pd
import numpy as np
from astropy.time import Time
import trove_mpc

def main():

    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--ffr", help="First filter radius, in deg", default=1)
    args = p.parse_args()
    
    test_objs = pd.read_csv("test_objs_mpc.csv", sep=",")
    test_objs.columns = test_objs.columns.str.strip()
    test_objs.name = test_objs.name.str.strip()

    tns_data = pd.read_csv("tns_public_objects.csv", skiprows=1)
    tns_data["full_name"] = tns_data.name_prefix + tns_data.name

    all_data = test_objs.merge(tns_data, left_on="name", right_on="full_name")

    n_test = 1_000
    all_times = []
    correct_matches = 0
    result_dict = []
    for i, row in all_data.iterrows():
        if i >= n_test: break
        
        print(row.full_name, "ra=", row.ra_x,
              "dec=", row.dec,
              "Discovery MJD=", Time(row.discoverydate, format="iso").mjd,
            )

        t = trove_mpc.Transient(row.ra_x, row.dec)
        
        start = time.time()
        filt_radius = 100
        the_match = t.minor_planet_match(
            Time(row.minor_planet_date.split("+")[0], format="iso").mjd,
            filter_radius=filt_radius, #arcsec
            first_cut_radius=float(args.ffr) #deg
        )
        dt = time.time() - start
        all_times.append(dt)
        
        if the_match is not None:
            print(f"Match found: {the_match.match_name} @ {the_match.distance[0]}'', Match on SAGUARO: {row.minor_planet_match}")
            is_corr = the_match.match_name in row.minor_planet_match or row.minor_planet_match in the_match.match_name
            if is_corr:
                correct_matches += 1

            d = the_match.todict()
            d["dt"] = dt
            d["original_match"] = str(row.minor_planet_match)
            d["same_match"] = is_corr
            d["distance"] = d["distance"][0]
            
            result_dict.append(d)
            
        else:
            print(f"NO MATCH FOUND WITHIN {filt_radius}''")
        print(f"Completed in {dt}s")
        print()
        
    print(f"On average, the MPC match took {np.mean(all_times)}s")
    print(f"There were {correct_matches} correct matches out of {min(n_test, len(all_data))}")

    print(result_dict)
    
    with open("test_result.json", "w") as j:
        json.dump(result_dict, j, indent=4)
    
if __name__ == "__main__":
    main()
