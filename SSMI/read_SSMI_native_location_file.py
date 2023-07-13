# reads SSMIS location file



import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def read_SSMI_native_location_file(ksat,orbit_num):
    ssmis_loc_path  = Path('M:/job_access/resampling/SSMI_Locations')
    filename = ssmis_loc_path / f'ssmi_footprint_location_file_{ksat:02d}_{orbit_num:05d}.dat'

    print(f'filename = {filename}')
    with open(filename,mode='rb') as file:

        num_cels = int(np.fromfile(file,dtype=np.int32,count=1,offset=0))
        num_scans = int(np.fromfile(file,dtype=np.int32,count=1,offset=0))
        print(num_cels,num_scans)
        count_to_read = num_cels*num_scans
    
        lons = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        lons19 = np.reshape(lons,(num_scans,num_cels))
        lats = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        lats19 = np.reshape(lats,(num_scans,num_cels))
        tht = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        tht19 = np.reshape(tht,(num_scans,num_cels))
        phi = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        phi19 = np.reshape(phi,(num_scans,num_cels))
        rng = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        rng19 = np.reshape(rng,(num_scans,num_cels))

        num_cels = int(np.fromfile(file,dtype=np.int32,count=1,offset=0))
        num_scans = int(np.fromfile(file,dtype=np.int32,count=1,offset=0))
        print(num_cels,num_scans)
        count_to_read = num_cels*num_scans
    
        lons = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        lons22 = np.reshape(lons,(num_scans,num_cels))
        lats = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        lats22 = np.reshape(lats,(num_scans,num_cels))
        tht = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        tht22 = np.reshape(tht,(num_scans,num_cels))
        phi = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        phi22 = np.reshape(phi,(num_scans,num_cels))
        rng = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        rng22 = np.reshape(rng,(num_scans,num_cels))


        num_cels = int(np.fromfile(file,dtype=np.int32,count=1,offset=0))
        num_scans = int(np.fromfile(file,dtype=np.int32,count=1,offset=0))
        print(num_cels,num_scans)
        count_to_read = num_cels*num_scans

        lons = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        lons37 = np.reshape(lons,(num_scans,num_cels))
        lats = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        lats37 = np.reshape(lats,(num_scans,num_cels))
        tht = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        tht37 = np.reshape(tht,(num_scans,num_cels))
        phi = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        phi37 = np.reshape(phi,(num_scans,num_cels))
        rng = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        rng37 = np.reshape(rng,(num_scans,num_cels))

        num_cels = int(np.fromfile(file,dtype=np.int32,count=1,offset=0))
        num_scans = int(np.fromfile(file,dtype=np.int32,count=1,offset=0))
        print(num_cels,num_scans)
        count_to_read = num_cels*num_scans

        lons = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        lons85 = np.reshape(lons,(num_scans,num_cels))
        lats = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        lats85 = np.reshape(lats,(num_scans,num_cels))
        tht = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        tht85 = np.reshape(tht,(num_scans,num_cels))
        phi = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        phi85 = np.reshape(phi,(num_scans,num_cels))
        rng = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        rng85 = np.reshape(rng,(num_scans,num_cels))

        return {'lats19' : lats19,
                'lons19' : lons19,
                'tht19'  : tht19,
                'phi19'  : phi19,
                'rng19'  : rng19,

                'lats22' : lats22,
                'lons22' : lons22,
                'tht22'  : tht22,
                'phi22'  : phi22,
                'rng22'  : rng22,

                'lats37' : lats37,
                'lons37' : lons37,
                'tht37'  : tht37,
                'phi37'  : phi37,
                'rng37'  : rng37,

                'lats85' : lats85,
                'lons85' : lons85,
                'tht85'  : tht85,
                'phi85'  : phi85,
                'rng85'  : rng85,
                }




if __name__ == '__main__':

    ksat = 13
    orbit_num = 5000

    locs = read_SSMI_native_location_file(ksat,orbit_num)

    lats19 = locs['lats19']
    lons19 = locs['lons19']
    tht19 = locs['tht19']
    phi19 = locs['phi19']
    rng19 = locs['rng19']

    tht37 = locs['tht37']
    phi37 = locs['phi37']
    rng37 = locs['rng37']

    print(f'19 GHz range {np.mean(rng19[rng19 > 100000.0])}')
    print(f'19 GHz tht {np.mean(tht19[tht19 > 10.0])}')

    print(f'37 GHz range {np.mean(rng37[rng37 > 100000.0])}')
    print(f'37 GHz tht {np.mean(tht37[tht37 > 10.0])}')

    lats37 = locs['lats37']
    lons37 = locs['lons37']

    scan0 = 456

    print(lats19[scan0,0:5])
    print(lons19[scan0,0:5])
    print(tht19[scan0,0:5])
    print(0.5*(lons19[scan0,0]+lons19[scan0,2]))
    print((0.5*(lons19[scan0,0]+lons19[scan0,2]) - lons19[scan0,1])*111.0)
    print((0.5*(lats19[scan0,0]+lats19[scan0,2]) - lats19[scan0,1])*111.0)

    print((lons19[scan0,0]-lons19[scan0,1])*111.)
    print((lats19[scan0,0]-lats19[scan0,1])*111.)

    print
    print((0.5*(lons19[scan0,0]+lons19[scan0+2,0]) - lons19[scan0+1,0])*111.0)
    print((0.5*(lats19[scan0,0]+lats19[scan0+2,0]) - lats19[scan0+1,0])*111.0)

    fig,ax = plt.subplots()
    ax.plot(lons37[scan0,:],lats37[scan0,:])

    fig,ax = plt.subplots()
    lon_err = (0.5*(lons19[scan0,0:-2]+lons19[scan0,2:]) - lons19[scan0,1:-1])*111.0
    lat_err = (0.5*(lats19[scan0,0:-2]+lats19[scan0,2:]) - lats19[scan0,1:-1])*111.0
    ax.plot(lon_err*0.25)
    ax.plot(lat_err*0.25)
    fig,ax = plt.subplots()
    lon_err2 = (0.5*(lons19[scan0,0:-4]+lons19[scan0,4:]) - lons19[scan0,2:-2])*111.0
    ax.plot(lon_err2)
    ax.plot(lon_err)
    ax.plot(lon_err*4.)
    
    
    plt.show()
    print