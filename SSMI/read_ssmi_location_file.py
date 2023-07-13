import numpy as np

def read_location_file(sat_name='SSMI',ksat=13,beamwidth=70,footprint_type='target',freq=0):
    if footprint_type == 'target':
        path = f"L:access/resampling/{sat_name}/f{ksat:02d}/target_gains/locs/"
        file_name = f"{path}find_target_gain_f{ksat:02d}_circular.loc"
    elif footprint_type == 'source':
        path = f"L:access/resampling/{sat_name}/f{ksat:02d}/source_gains/locs/"
        file_name = f"{path}find_source_gain_f{ksat:02d}_freq_{freq:02d}_locs.txt"
    scan = []
    fov = []
    lats = []
    lons = []
    azim = []
    with open(file_name,'r') as f:
        for line in f:
            x = line.split()
            if footprint_type == 'target':
                scan.append(int(x[0]))
                fov.append(int(x[1]))
                lats.append(float(x[2]))
                lons.append(float(x[3]))
            elif footprint_type == 'source':
                scan.append(int(x[1]))
                fov.append(int(x[2]))
                lats.append(float(x[3]))
                lons.append(float(x[4]))
                azim.append(float(x[5]))
           
    scan = np.array(scan)
    fov = np.array(fov)
    lats = np.array(lats)
    lons = np.array(lons)

    if footprint_type=='target':
        scan = np.reshape(scan,(4,253))
        fov = np.reshape(fov,(4,253))
        lats = np.reshape(lats,(4,253))
        lons = np.reshape(lons,(4,253))
    elif footprint_type=='source':
        if freq < 3:
            scan = np.reshape(scan,(29,64))
            fov = np.reshape(fov,(29,64))
            lats = np.reshape(lats,(29,64))
            lons = np.reshape(lons,(29,64))
            azim = np.reshape(azim,(29,64))
        else:
            scan = np.reshape(scan,(57,128))
            fov = np.reshape(fov,(57,128))
            lats = np.reshape(lats,(57,128))
            lons = np.reshape(lons,(57,128))
            azim = np.reshape(azim,(57,128))
    else:
        raise ValueError(f'Invalid footprint type: {footprint_type}')

    if footprint_type == 'target':
        return lats,lons
    else:
        return lats,lons,azim