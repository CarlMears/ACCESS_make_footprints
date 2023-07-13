
'''
	module l1c_content
	implicit none

 	integer(4), parameter :: maxscan=3600  !1 orbit plus 5% overlap	
 	integer(4), parameter :: maxcel=128 
 	integer(4), parameter :: maxchn=8 

c     l1c file content
	integer(4) ksat_l1c,iorbit_l1c,numscan
	
	integer(4) iqual_flag(maxscan)
	
    real(8) scan_time(maxscan),orbit(maxscan),zang(maxscan),scpos(3,maxscan),scvel(3,maxscan)
	real(4) scloc(3,maxscan),therm(7,maxscan)

	real(4) cellat(maxcel,maxscan),cellon(maxcel,maxscan),celtht(maxcel,maxscan)
	real(4) celphi(maxcel,maxscan),celsun(maxcel,maxscan),celrfi(maxcel,maxscan)
 	
 	real(4) ta_nat(maxchn,maxcel,maxscan)
 	real(4) ta_rsp(maxchn,maxcel,maxscan)
 	real(4) tb_nat(maxchn,maxcel,maxscan)
 	real(4) tb_rsp(maxchn,maxcel,maxscan)
	end module l1c_content

subroutine read_l1c_file(ksat,iorbit,filename_l1c, start_time)
      use l1c_content
      implicit none
      
      integer(4) ksat,iorbit
	  character(120) filename_l1c
      real(8) start_time

      call openbig(3,filename_l1c,'old')
      read(3) ksat_l1c,iorbit_l1c,numscan
      read(3) iqual_flag
      read(3) scan_time,orbit,zang,scpos,scvel,scloc,therm
	  read(3) cellat,cellon,celtht,celphi,celsun,celrfi
      read(3) ta_nat,ta_rsp,tb_nat,tb_rsp
	close(3)
	
	if(numscan.le.0 .or. numscan.gt.maxscan) stop 'error in ingesting numscan, pgm stopped'
	
	if(ksat.ne.ksat_l1c .or. iorbit.ne.iorbit_l1c) stop 'error in ingesting ksat or iorbit, pgm stopped'
	
	start_time=scan_time(1)
	
	return
	end
    '''
from pathlib import Path 
import matplotlib.pyplot as plt

def read_ssmis_l1c_file(*,
                        ksat:int,
                        orbit:int,
                        base_path:Path = Path('S:/ssmi/F13/L1C/'),
                        debug=False,
                        verbose=False):
    
    import numpy as np
    import xarray as xr

    maxscan=3600  
    maxcel=128 
    maxchn=8 

    orbit_lower = int((orbit-1)/5000)*5000+1
    orbit_upper = orbit_lower + 4999



    path = base_path /f'F{ksat:02d}' / 'L1C' / f'r{orbit_lower:06d}_{orbit_upper:06d}'
    filename = path / f'f{ksat:02d}_r{orbit:06d}.dat'
    if verbose:
        print(f'filename = {filename}')
    with open(filename,'rb') as f:
        temp = np.fromfile(f,dtype=np.int32,count=3)
        ksat_l1c = temp[0]
        iorbit_l1c = temp[1]
        numscan = temp[2]

        iqual_flag = np.fromfile(f,dtype=np.int32,count=maxscan)
        if debug:
            fig0,ax0 = plt.subplots(1,1)
            ax0.plot(iqual_flag)
            ax0.set_title('iqual_flag')
        bad = np.where(iqual_flag != 0)
        scan_time = np.fromfile(f,dtype=np.float64,count=maxscan) 
        
        orbit = np.fromfile(f,dtype=np.float64,count=maxscan)
        
        bad2 = np.where(np.any([(orbit < 1.0),iqual_flag != 0],axis=0))
        scan_time[bad2] = np.nan
        orbit[bad2] = np.nan
        if debug:
            fig1,ax1 = plt.subplots(1,1)
            ax1.plot(scan_time)
            ax1.set_title('scan_time')
        if debug:
            fig2,ax2 = plt.subplots(1,1)
            ax2.plot(orbit)
            ax2.set_title('orbit')
        zang = np.fromfile(f,dtype=np.float64,count=maxscan)
        zang[bad2] = np.nan
        if debug:
            fig3,ax3 = plt.subplots(1,1)
            ax3.plot(zang)
            ax3.set_title('zang')

        scpos = np.fromfile(f,dtype=np.float64,count=3*maxscan).reshape((maxscan,3))
        scvel = np.fromfile(f,dtype=np.float64,count=3*maxscan).reshape((maxscan,3))
        scloc = np.fromfile(f,dtype=np.float32,count=3*maxscan).reshape((maxscan,3))
        therm = np.fromfile(f,dtype=np.float32,count=7*maxscan).reshape((maxscan,7))

        cellat = np.fromfile(f,dtype=np.float32,count=maxcel*maxscan).reshape((maxscan,maxcel))
        cellon = np.fromfile(f,dtype=np.float32,count=maxcel*maxscan).reshape((maxscan,maxcel))
        celtht = np.fromfile(f,dtype=np.float32,count=maxcel*maxscan).reshape((maxscan,maxcel))
        celphi = np.fromfile(f,dtype=np.float32,count=maxcel*maxscan).reshape((maxscan,maxcel))
        celsun = np.fromfile(f,dtype=np.float32,count=maxcel*maxscan).reshape((maxscan,maxcel))
        celrfi = np.fromfile(f,dtype=np.float32,count=maxcel*maxscan).reshape((maxscan,maxcel))

        ta_nat = np.fromfile(f,dtype=np.float32,count=maxchn*maxcel*maxscan).reshape((maxscan,maxcel,maxchn))
        ta_rsp = np.fromfile(f,dtype=np.float32,count=maxchn*maxcel*maxscan).reshape((maxscan,maxcel,maxchn))
        tb_nat = np.fromfile(f,dtype=np.float32,count=maxchn*maxcel*maxscan).reshape((maxscan,maxcel,maxchn))
        tb_rsp = np.fromfile(f,dtype=np.float32,count=maxchn*maxcel*maxscan).reshape((maxscan,maxcel,maxchn))

    ds = xr.Dataset(
        data_vars={ 'scan_time':(['scan'],scan_time),
                    'orbit':(['scan'],orbit),
                    'zang':(['scan'],zang),
                    'scpos':(['scan','xyz'],scpos),
                    'scvel':(['scan','xyz'],scvel),
                    'scloc':(['scan','xyz'],scloc),
                    'therm':(['scan','therm_num'],therm),
                    'cellat':(['scan','cell'],cellat),
                    'cellon':(['scan','cell'],cellon),
                    'celtht':(['scan','cell'],celtht),
                    'celphi':(['scan','cell'],celphi),
                    'celsun':(['scan','cell'],celsun),
                    'celrfi':(['scan','cell'],celrfi),
                    'ta_nat':(['scan','cell','channel'],ta_nat),
                    'ta_rsp':(['scan','cell','channel'],ta_rsp),
                    'tb_nat':(['scan','cell','channel'],tb_nat),
                    'tb_rsp':(['scan','cell','channel'],tb_rsp)
                },
        
        coords={'scan':np.arange(1,maxscan+1),
                'cell':np.arange(1,maxcel+1),
                'channel':np.arange(1,maxchn+1),
                'xyz':np.arange(1,4),
                'therm_num':np.arange(1,8)},
        attrs={ 'ksat':ksat_l1c,
                'iorbit':iorbit_l1c,
                'numscan':numscan,
                'maxscan':maxscan,
                'maxcel':maxcel,
                'maxchn':maxchn}
    )

    return ds


if __name__ == '__main__':

    ksat=13
    orbit = 5005


    ds = read_ssmis_l1c_file(ksat=ksat,
                            orbit=orbit,
                            base_path = Path('S:\ssmi'),
                            debug=True,
                            verbose=True)

    test_file = Path(f'M:/job_access/python/resample_wts/SSMI/test_nc/f13_r{orbit:06d}.nc')
    print(f'Writing {test_file}')
    ds.to_netcdf(test_file)
