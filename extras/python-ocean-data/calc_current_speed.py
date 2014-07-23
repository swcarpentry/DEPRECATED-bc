import os, sys
import datetime
from git import Repo

from netCDF4 import Dataset


def main():
    script = sys.argv[0]
    inFile = sys.argv[1]
    uVar = sys.argv[2]
    vVar = sys.argv[3]
    outfile_name = sys.argv[4]

    uData, vData, input_DATA = read_data(inFile, uVar, vVar)
    wspData = calc_wsp(uData, vData)
    write_output(wsp_data, input_DATA, outfile_name)


def read_data(ifile, uVar, vVar):
    """Read data from ifile corresponding to the U and V variable"""

    input_DATA = Dataset(acorn_URL)
    uData = input_DATA.variables[uVar][:]
    vData = input_DATA.variables[vVar][:]

    return uData, vData, input_DATA


def calc_wsp(uwnd, vwnd):
    """Calculate the wind speed"""

    wsp = (uwnd**2 + vwnd**2)**0.5

    return wsp


def write_output(spcurData, input_DATA, outfile_name):
    """Write the output file"""

    outfile = Dataset(outfile_name, 'w', format='NETCDF4')

    # Create file dimensions (by copying from input file)
    for dimName, dimData in input_DATA.dimensions.iteritems():
        outfile.createDimension(dimName, len(dimData))

    # Create variables corresponding to those dimensions (by copying)    
    for var_name in ['TIME', 'LATITUDE', 'LONGITUDE']:
        varin = input_DATA.variables[var_name]
        outVar = outfile.createVariable(var_name, varin.datatype, varin.dimensions)
        outVar[:] = varin[:]
    
    # Create the speed variable
    uData = input_DATA.variables['UCUR']   
    spcur = outfile.createVariable(varname, uData.datatype, uData.dimensions)
    spcur[:,:,:] = spcurData    
    
    # Global file attributes
    global_atts = {}
    for att in input_DATA.ncattrs():
        global_atts[att] = input_DATA.att  
    new_history = create_history()
    global_atts['history'] = """%s\n%s""" %(new_history,  global_atts['history'])
    outfile.setncatts(global_atts)
    
    outfile.close()


def create_history():
    """Create the new entry for the global history file attribute"""

    time_stamp = datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    exe = sys.executable
    args = " ".join(sys.argv)
    git_hash = Repo(os.getcwd()).head.commit.hexsha

    return """%s: %s %s (Git hash: %s)""" %(time_stamp, exe, args, git_hash[0:7])


main()
