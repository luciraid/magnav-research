# MagNav Data Sanity Check – Flt1002

## File
- Flt1002_train.h5
- Source: DAF-MIT AIA Open Flight Data for Magnetic Navigation Research

## Status
- File successfully opened in Julia using HDF5.jl
- Dataset structure is readable and intact

## Observed Contents
- Time-synchronized sensor data
- GPS position (latitude, longitude, altitude)
- INS data (attitude and inertial measurements)
- Scalar magnetometer measurements (total field)
- Vector magnetometer measurements (fluxgate)

## Interpretation
This flight contains both clean (tail stinger) and in-aircraft magnetic
measurements along with navigation data, making it suitable for Tolles–Lawson
calibration as used in the MagNav pipeline.
