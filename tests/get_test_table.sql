SELECT t.name, t.ra, t.dec, mpc.target_id, mpc.minor_planet_match, mpc.minor_planet_date
FROM tom_targets_basetarget AS t
JOIN (
     SELECT
	target_id,
	value AS minor_planet_match,
	te2.minor_planet_date
     FROM tom_targets_targetextra AS te1
     JOIN (
	     SELECT
		target_id AS target_id_date,
		value AS minor_planet_date 
	     FROM tom_targets_targetextra
	     WHERE key = 'Minor Planet Date'
     ) AS te2 ON te1.target_id = te2.target_id_date
     WHERE te1.key = 'Minor Planet Match' AND te1.value <> 'None'
) AS mpc ON t.id = mpc.target_id
