import numpy 
from numpy import *
import scipy 
from scipy import *


def func_multiplicative1(xy_tuple,*params):
    (sn, res)=xy_tuple
    params=array(params)
    m=(params[0]/sn**params[1])*exp(-1./res**params[2])
    return m.ravel()

def func_m_snonly1(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=params[0]*(sn**params[1])
    return m.ravel()

def func_m_snonly2(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=params[0]+(sn/params[1])**params[2]*(1+(sn/params[1]))**params[3]
    return m.ravel()

def func_m_snonly3(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=params[0]+(sn/params[1])**params[2]*(1+(sn/params[3]))**params[4]
    return m.ravel()

def func_m_snonly4(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=params[0]+sin((sn-params[1])/params[2])*((sn-params[1])**params[3])
    return m.ravel()

def func_m_magonly1(xy_tuple,*params):
    (sn, res)=xy_tuple
    params=array(params)
    m=-params[0]*sn**params[1]
    return m.ravel()

def func_m_magonly2(xy_tuple,*params):
    (sn, res)=xy_tuple
    params=array(params)
    m=-exp((sn/params[0])**params[1])+params[2]
    return m.ravel()

def func_m_magonly3(xy_tuple,*params):
    (sn, res)=xy_tuple
    params=array(params)
    m=-params[0]*(sn/params[1])**params[2]
    return m.ravel()

def func_m_magonly4(xy_tuple,*params):
    (sn, res)=xy_tuple
    params=array(params)
    m=-(sn/params[0])**params[1]+params[2]
    return m.ravel()

def func_m_magonly5(xy_tuple,*params):
    (sn, res)=xy_tuple
    params=array(params)
    m=params[0]*(1/(exp((sn-params[1])/params[2])+1)-1.)
    return m.ravel()

def func_c_magonly1(xy_tuple,*params):
    (sn, res)=xy_tuple
    params=array(params)
    m=params[0]
    return m.ravel()

def func_multiplicative2(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=params[0]/sn*exp(-(params[1]/res)**params[2])
    return m.ravel()

def func_multiplicative3(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=(1/sn)**params[0]*exp(-(params[1]/res)**params[2])
    return m.ravel()

def func_multiplicative4(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=(params[0]/sn)**params[1]*-log(params[2]*res)
    return m.ravel()

def func_multiplicative5(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=params[0]/sn**params[1]*-log(params[2]*res)
    return m.ravel()

def func_multiplicative6(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=params[0]/sn**params[1]*-res**params[2]
    return m.ravel()

def func_multiplicative7(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=params[0]/sn**params[1]*exp(-params[3]/res**params[2])
    return m.ravel()
    
def func_multiplicative8(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=params[0]/log(sn)*-log(params[2]*res**params[1])
    return m.ravel()

def func_multiplicative9(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=params[0]/sn**params[1]*exp(-((sn/res)*params[3])**params[2])
    return m.ravel()

def func_multiplicative10(xy_tuple,*params):
    (sn, res)=xy_tuple
    m=exp(-(params[2]*(sn/res))**params[3])*(params[0]/(sn)**params[1])#(params[0]*log(abs(params[1]/sn))+params[2])#*exp(-(params[1]/res)**params[2])(params[0]/sn**params[1])
    #m=(params[0]/sn**params[1])*exp(-1./res**params[2])

    return m.ravel()

def func_multiplicative11(xy_tuple,*params):
    (sn, res)=xy_tuple
    params=array(params)
    m=params[0]*exp(-params[1]*sn**params[3])*exp(-1/res**params[2])
    return m.ravel()

def func_multiplicative12(xy_tuple,*params):
    (sn, res)=xy_tuple
    params=array(params)
    m=params[0]*exp(-params[1]*sn**(params[3]*res))*exp(-1/res**params[2])
    return m.ravel()


def func_additive1(xy_tuple,*params):
    (sn, res)=xy_tuple
    params=array(params)
    A=(params[0]/sn**params[1])*(res**params[2])
    return A.ravel()
    
def func_additive2(xy_tuple,*params):
    (sn, res)=xy_tuple
    A=(params[0]/sn**params[1])*exp(-1./res**params[2])
    return A.ravel()


def henk_multiplicative(xy_tuple,*params):
    (sn, res)=xy_tuple
    params=array(params)
    f0 = params[0] + params[1]/sn + params[2]/(sn*sn) + params[3]/sqrt(sn);
    f1 = params[4] + params[5]/sn + params[6]/(sn*sn) + params[7]/sqrt(sn);
    f2 = params[8] + params[9]/sn + params[10]/(sn*sn) + params[11]/sqrt(sn);
    f3 = params[12] + params[13]/sn + params[14]/(sn*sn) + params[15]/sqrt(sn);

    m=f0 + f1/res + f2*res + f3*res*res
    return m.ravel()



#dict for the functions and their number of parameters
function_set={'func_multiplicative1' : func_multiplicative1,'func_multiplicative2' : func_multiplicative2,'func_multiplicative3' : func_multiplicative3,'func_multiplicative4' : func_multiplicative4, 'func_multiplicative5' : func_multiplicative5,'func_multiplicative6' : func_multiplicative6,'func_multiplicative7' : func_multiplicative7, 'func_multiplicative8' : func_multiplicative8,'func_multiplicative9' : func_multiplicative9,'func_multiplicative10' : func_multiplicative10,'func_multiplicative11' : func_multiplicative11,'func_multiplicative12' : func_multiplicative12, 'func_additive1' : func_additive1, 'func_additive2' : func_additive2,'func_m_snonly1' : func_m_snonly1,'func_m_magonly1' : func_m_magonly1,'func_c_magonly1' : func_c_magonly1, 'func_m_magonly2' : func_m_magonly2, 'func_m_magonly3' : func_m_magonly3,  'func_m_magonly4' : func_m_magonly4, 'func_m_magonly5' : func_m_magonly5, 'func_m_snonly2' : func_m_snonly2, 'func_m_snonly3' : func_m_snonly3, 'func_m_snonly4' : func_m_snonly4, 'henk_multiplicative' : henk_multiplicative}

function_params={'func_multiplicative1' : 3,'func_multiplicative2' : 3,'func_multiplicative3' : 3,'func_multiplicative4' : 3, 'func_multiplicative5' : 3,'func_multiplicative6' : 3,'func_multiplicative7' : 4, 'func_multiplicative8' : 3,'func_multiplicative9' : 4,'func_multiplicative10' : 4, 'func_multiplicative11' : 4, 'func_multiplicative12' : 4, 'func_additive1' : 3, 'func_additive2' : 3, 'func_m_snonly1' : 2, 'func_m_magonly1' : 2, 'func_c_magonly1' : 1, 'func_m_magonly2' : 3, 'func_m_magonly3' : 3,  'func_m_magonly4' : 3, 'func_m_magonly5' : 3, 'func_m_snonly2' : 4, 'func_m_snonly3' : 5, 'func_m_snonly4' : 4, 'henk_multiplicative' : 16}
