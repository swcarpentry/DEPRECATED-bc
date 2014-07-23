import os, sys
import datetime
from git import Repo

from netCDF4 import Dataset


def main():
    # Read command line arguments
    script = sys.argv[0]
    inFile = sys.argv[1]
    uVar = sys.argv[2]
    vVar = sys.argv[3]
    outfile_name = sys.argv[4]

    # Read input data 
    uData, vData, input_DATA = read_data(inFile, uVar, vVar)
    
    # Calculate the current speed
    spData = calc_speed(uData, vData)
    
    # Write the output file
    outfile = Dataset(outfile_name, 'w', format='NETCDF4')
    set_global_atts(input_DATA, outfile)
    copy_dimensions(input_DATA, outfile)
    copy_variables(input_DATA, outfile)
    write_speed(input_DATA, outfile, spData) 
    
    outfile.close()
    

def read_data(ifile, uVar, vVar):
    """Read data from ifile corresponding to the U and V variable"""

    input_DATA = Dataset(ifile)
    uData = input_DATA.variables[uVar][:]
    vData = input_DATA.variables[vVar][:]

    return uData, vData, input_DATA


def calc_speed(u, v):
    """Calculate the speed"""

    speed = (u**2 + v**2)**0.5

    return speed


def copy_dimensions(infile, outfile):
    """Copy the dimensions of the infile to the outfile"""
        
    for dimName, dimData in infile.dimensions.iteritems():
        outfile.createDimension(dimName, len(dimData))


def set_global_atts(infile, outfile):
    """Set the global attributes for outfile.
        
    Note that the global attributes are simply copied from
    infile and the history attribute updated accordingly.
        
    """
        
    global_atts = {}
    for att in infile.ncattrs():
        global_atts[att] = eval('infile.'+att)  
        
    new_history = create_history()
    global_atts['history'] = """%s\n%s""" %(new_history,  global_atts['history'])
    outfile.setncatts(global_atts)


def create_history():
    """Create the new entry for the global history file attribute"""

    time_stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    exe = sys.executable
    args = " ".join(sys.argv)
    git_hash = Repo(os.getcwd()).heads[0].commit

    return """%s: %s %s (Git hash: %s)""" %(time_stamp, exe, args, str(git_hash)[0:7])


def copy_variables(infile, outfile):
    """Create variables corresponding to the file dimensions 
    by copying from infile"""
            
    for var_name in ['TIME', 'LATITUDE', 'LONGITUDE']:
        varin = infile.variables[var_name]
        outVar = outfile.createVariable(var_name, varin.datatype, 
                                        varin.dimensions, 
                                        fill_value=varin._FillValue)
        outVar[:] = varin[:]
            
        var_atts = {}
        for att in varin.ncattrs():
            if not att == '_FillValue':
                var_atts[att] = eval('varin.'+att) 
        outVar.setncatts(var_atts)


def write_speed(infile, outfile, spData):
    """Write the current speed data to outfile"""
        
    u = infile.variables['UCUR']   
    spcur = outfile.createVariable('SPCUR', u.datatype, u.dimensions, fill_value=u._FillValue)
    spcur[:,:,:] = spData    
    
    spcur.standard_name = 'sea_water_speed'
    spcur.long_name = 'sea water speed'
    spcur.units = u.units
    spcur.coordinates = u.coordinates
    

main()
