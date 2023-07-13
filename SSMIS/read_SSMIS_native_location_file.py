# reads SSMIS location file



import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def read_SSMIS_native_location_file(ksat,orbit_num):
    ssmis_loc_path  = Path('M:/job_access/resampling/SSMIS_Locations')
    filename = ssmis_loc_path / f'footprint_location_file_{ksat:02d}_{orbit_num:05d}.dat'

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
        lons89 = np.reshape(lons,(num_scans,num_cels))
        lats = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        lats89 = np.reshape(lats,(num_scans,num_cels))
        tht = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        tht89 = np.reshape(tht,(num_scans,num_cels))
        phi = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        phi89 = np.reshape(phi,(num_scans,num_cels))
        rng = np.fromfile(file,dtype=np.float32,count=count_to_read,offset=0)
        rng89 = np.reshape(rng,(num_scans,num_cels))

        lats22 = np.copy(lats19)
        lons22 = np.copy(lons19)
        tht22 = np.copy(tht19)
        phi22 = np.copy(phi19)
        rng22 = np.copy(rng19)

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

                'lats89' : lats89,
                'lons89' : lons89,
                'tht89'  : tht89,
                'phi89'  : phi89,
                'rng89'  : rng89,
                }




if __name__ == '__main__':

    ksat = 18
    orbit_num = 10001

    locs = read_SSMIS_native_location_file(ksat,orbit_num)

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

    print(lats19[920,0:5])
    print(lons19[920,0:5])
    print(tht19[920,0:5])
    print(0.5*(lons19[920,0]+lons19[920,2]))
    print((0.5*(lons19[920,0]+lons19[920,2]) - lons19[920,1])*111.0)
    print((0.5*(lats19[920,0]+lats19[920,2]) - lats19[920,1])*111.0)

    print((lons19[920,0]-lons19[920,1])*111.)
    print((lats19[920,0]-lats19[920,1])*111.)

    print
    print((0.5*(lons19[920,0]+lons19[922,0]) - lons19[921,0])*111.0)
    print((0.5*(lats19[920,0]+lats19[922,0]) - lats19[921,0])*111.0)

    print
    print(lons37[0,0:5])


    fig,ax = plt.subplots()
    ax.plot(lats19[:,45])

    fig,ax = plt.subplots()
    ax.plot(lons19[:,45])

    fig,ax = plt.subplots()
    ax.plot(lats37[:,45])

    fig,ax = plt.subplots()
    ax.plot(lons37[:,45])


    fig,ax = plt.subplots()
    ax.plot(lons37[920,:],lats37[920,:])

    fig,ax = plt.subplots()
    lon_err = (0.5*(lons19[920,0:-2]+lons19[920,2:]) - lons19[920,1:-1])*111.0
    lat_err = (0.5*(lats19[920,0:-2]+lats19[920,2:]) - lats19[920,1:-1])*111.0
    ax.plot(lon_err*0.25)
    ax.plot(lat_err*0.25)
    fig,ax = plt.subplots()
    lon_err2 = (0.5*(lons19[920,0:-4]+lons19[920,4:]) - lons19[920,2:-2])*111.0
    ax.plot(lon_err2)
    ax.plot(lon_err)
    ax.plot(lon_err*4.)
    
    
    plt.show()
    print