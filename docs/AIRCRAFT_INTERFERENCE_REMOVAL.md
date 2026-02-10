# Aircraft Magnetic Interference Removal – Flt1002

## Method
- Used diurnal-corrected scalar magnetometer data
- Modeled aircraft magnetic interference using onboard current sensors
- Estimated linear interference model via least squares
- Subtracted estimated interference from magnetic signal

## Result
- Significant reduction in high-frequency spikes
- Cleaner magnetic signal suitable for anomaly mapping

## Output
- figures/Flt1002_aircraft_compensation.png

## Relation to Tolles–Lawson
This regression-based interference removal captures the core idea
behind Tolles–Lawson compensation using real aircraft disturbance inputs.
