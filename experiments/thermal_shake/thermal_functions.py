import random
import numpy as np

# Thermal functions
def new_name(main,i,n):
    """ Crea un nombre en formato string para un archivo nuevo correspondiente 
    a las N pruebas del análisis de Montecarlo"""
    i=str(i)
    num=i
    if len(i)<len(str(n)):
        dif=len(str(n))-len(i)
        num='0'*dif+i
    return main.split('.')[0]+num+'.'+main.split('.')[1]

def randomX(x,R): # R indica el radio atómico en Angstroms 
    """Teniendo un vector x, se agita una esfera que lo contiene entregando una 
    nueva posición"""
    r=random.uniform(0,R)
    theta=random.uniform(0,2*np.pi)
    phi=random.uniform(0,np.pi)
    new_x=[0,0,0]
    new_x[0]=x[0]+r*np.cos(theta)*np.sin(phi)
    new_x[1]=x[1]+r*np.sin(theta)*np.sin(phi)
    new_x[2]=x[2]+r*np.cos(phi)
    return new_x

def dec3(x):
    """Redondeo de un número float a 3 decimales para correcta escritura en .pqr"""
    valor= str(round(x,3))
    if '.' in valor:
        decimal=valor.split('.')[-1]
        if len(decimal)<3:
            decimal=decimal+'000'
            decimal=decimal[:3]
            return valor.split('.')[0]+'.'+decimal
        else:
            return valor
    else:
        return valor+'.000'

def nombre_atomo(atomo,es_biomolecula=False):
    """ Tomando un string correspondiente al átomo de un archivo pdb (formato C1, C2, etc) 
    se retorna sólo el elemento químico sin el número.
    En caso de que sea biomolecula, puede que algunos atomos se nombren en formato CA (carbono alpha)
    y derivados, en ese caso se retorna solo CHONSP"""
    final=''
    for caracter in atomo:
        if caracter not in "0123456789":
            final+=caracter
    if es_biomolecula:
        for elem in 'CHONSP':
            if elem in final:
                return elem
    return final
