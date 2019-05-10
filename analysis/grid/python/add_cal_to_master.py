from sys import *
from scipy import *
from numpy import *
from os import *
from os.path import isfile 
import matplotlib
matplotlib.use('Agg')
from MyFunction import *
from fitting_functions import *
import pyfits

#Reading command line

#The master catalogue
infile = sys.argv[1]

#Column base name (e.g. mycalibration)
col_name = sys.argv[2]

#MV: the normalised version of the resolution
#anything else: unnormalised version (R_gal/R_PSF)
res_switch = sys.argv[3]

#Table: fine sampled calibration (100x100 or so). Calibration is based on
#a nearest neighbour search in the table.
#Bin: calibration type in IFC17. You need the surface of doom 
cal_type = sys.argv[4]

#The surface of doom file (direct output of the pipeline in the 2D case)
cal_file = sys.argv[5]

#Output fits file containg the additional columns
outfile = sys.argv[6]

#Reading full FITS catalogue

if(os.path.exists(infile) and os.path.exists(cal_file)):
    fits_data = pyfits.open(infile)
    
    #Getting relevant data from FITS table
    snr = fits_data[1].data['snr_model']
    psf_size = fits_data[1].data['psf_size_in']
    size = fits_data[1].data['size_out']
    e1 = fits_data[1].data['e1']
    e2 = fits_data[1].data['e2']
    e = sqrt(e1*e1 + e2*e2)
    size_ab = size*sqrt((1-e)/(1+e))
    #Creating sqrt(ab) size
    
    if(res_switch == 'MV'):
        resolution = psf_size/((size_ab**2)+psf_size)
    else:
        resolution = psf_size/(size_ab**2)

    #Getting everything necessary for the m calibration

    #easy for table cal
    table_data=genfromtxt(cal_file)
    
    #need tree for table cal
    if(cal_type == "table"):
        x_table=table_data[:,0]
        y_table=table_data[:,1]
        m1_table=table_data[:,2]
        c1_table=table_data[:,3]
        m2_table=table_data[:,4]
        c2_table=table_data[:,5]
        shift_x = min(x_table)
        shift_y = min(y_table)
        x_table = x_table - shift_x
        y_table = y_table - shift_y
        scale_x = max(x_table)
        scale_y = max(y_table)
        x_table=x_table/scale_x
        y_table=y_table/scale_y
        tree_pos = numpy.dstack([x_table,y_table])[0]
        Tree=scipy.spatial.cKDTree(tree_pos,leafsize=128)
    elif(cal_type == "function"):
        m1_params = table_data[:,0]
        m2_params = table_data[:,1]
        #c1_params = table_data[:,2]
        #c2_params = table_data[:,3]

    function = henk_multiplicative

    #Looping through all galaxies
    m1=[]
    m2=[]
    c1=[]
    c2=[]
    for i in range(0,len(snr)):
        if(cal_type=="table"):
            cals=bias_2D_table_query_all(snr[i], resolution[i],m1_table,m2_table,c1_table,c2_table,Tree,scale_x,scale_y,shift_x,shift_y)
        elif(cal_type=="bin"):
            cals=bias_2D_bin_query(snr[i],resolution[i],table_data)
        elif(cal_type=="function"):
            cals=bias_2D_function_query(snr[i],resolution[i],function,function,m1_params,m2_params)


        m1.append(cals[0])
        m2.append(cals[1])
        c1.append(cals[2])
        c2.append(cals[3])

    m1=array(m1)
    m2=array(m2)
    c1=array(c1)
    c2=array(c2)

    nc1=pyfits.Column(name='m1_'+col_name,format='D',array=m1)
    nc2=pyfits.Column(name='m2_'+col_name,format='D',array=m2)
    nc3=pyfits.Column(name='c1_'+col_name,format='D',array=c1)
    nc4=pyfits.Column(name='c2_'+col_name,format='D',array=c2)
    old_cols = fits_data[1].columns
    cols = (nc1,nc2,nc3,nc4)
    new_cols = pyfits.ColDefs(cols)
    new_hdu = pyfits.new_table(old_cols+new_cols)
    new_hdu.writeto(outfile,clobber = True)

else:
   print "Chosen input or calibration file does not exist."
