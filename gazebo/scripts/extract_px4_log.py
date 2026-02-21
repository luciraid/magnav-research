#!/usr/bin/env python3
"""
Extract sensor data from PX4 .ulg log files for MagNav processing.

Usage:
    python3 extract_px4_log.py /path/to/log.ulg output.csv
"""

import sys
import pandas as pd
from pyulog import ULog

def extract_log(ulg_path, output_csv):
    """Extract IMU, GPS, magnetometer data from ULog."""
    
    log = ULog(ulg_path)
    
    # Extract datasets
    data = {}
    
    # Time (microseconds)
    data['timestamp'] = log.get_dataset('vehicle_gps_position').data['timestamp']
    
    # GPS
    gps = log.get_dataset('vehicle_gps_position').data
    data['lat'] = gps['lat'] / 1e7  # degrees
    data['lon'] = gps['lon'] / 1e7
    data['alt'] = gps['alt'] / 1e3  # meters
    
    # IMU (sync to GPS timestamps via interpolation)
    imu = log.get_dataset('sensor_combined').data
    data['accel_x'] = pd.Series(imu['accelerometer_m_s2[0]']).reindex(data['timestamp'], method='nearest')
    data['accel_y'] = pd.Series(imu['accelerometer_m_s2[1]']).reindex(data['timestamp'], method='nearest')
    data['accel_z'] = pd.Series(imu['accelerometer_m_s2[2]']).reindex(data['timestamp'], method='nearest')
    data['gyro_x'] = pd.Series(imu['gyro_rad[0]']).reindex(data['timestamp'], method='nearest')
    data['gyro_y'] = pd.Series(imu['gyro_rad[1]']).reindex(data['timestamp'], method='nearest')
    data['gyro_z'] = pd.Series(imu['gyro_rad[2]']).reindex(data['timestamp'], method='nearest')
    
    # Magnetometer
    mag = log.get_dataset('vehicle_magnetometer').data
    data['mag_x'] = pd.Series(mag['magnetometer_ga[0]']).reindex(data['timestamp'], method='nearest')
    data['mag_y'] = pd.Series(mag['magnetometer_ga[1]']).reindex(data['timestamp'], method='nearest')
    data['mag_z'] = pd.Series(mag['magnetometer_ga[2]']).reindex(data['timestamp'], method='nearest')
    
    # Attitude (heading)
    att = log.get_dataset('vehicle_attitude').data
    data['roll'] = pd.Series(att['roll']).reindex(data['timestamp'], method='nearest')
    data['pitch'] = pd.Series(att['pitch']).reindex(data['timestamp'], method='nearest')
    data['yaw'] = pd.Series(att['yaw']).reindex(data['timestamp'], method='nearest')
    
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    print(f"✓ Saved {len(df)} samples to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 extract_px4_log.py input.ulg output.csv")
        sys.exit(1)
    
    extract_log(sys.argv[1], sys.argv[2])
