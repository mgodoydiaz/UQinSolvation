# Script para calcular incertidumbre en la molecula 5pti

import os
import numpy as np
import bem_electrostatics
from time import time
import gc

def completar_archivo_entrada(dic_parametros):
    """En base a los parametros de la funcion se completa el archivo de parametros de entrada para el solver externo"""
    arch=open('input_config.prm') # Este archivo siempre se llamara asi
    temp=open('temp.txt','w')
    for linea in arch:
        if '=' in linea:
            par,_=linea.strip().split('=')
            temp.write('{}={}\n'.format(par,dic_parametros[par]))
        else:
            temp.write(linea)
    arch.close()
    temp.close()

    # Se elimina archivo original y renombra archivo temporal
    os.remove('input_config.prm')
    os.rename('temp.txt','input_config.prm')
    return None

# Comienza programa para ejecutar N simulaciones de Montecarlo
START_TIME=time()

from shake_pdb import new_name

#Opciones de Montecarlo y directorios
mainfile='1pgb.pdb'
radios_prueba=[0.1,0.5,1,2,5]
n_tests=100
generador_mesh='nanoshaper'
imprimir_resultados=True
#Diccionario de parametros individuales de cada test
parametros={
        "PDBFILE":mainfile,
        "N_TESTS":n_tests,
        "MESH_GENERATOR":generador_mesh,
        "PARAM1":1,
        "PRINT":imprimir_resultados,
        "OUTPUTFILE":'',
        "TESTDIR":''
        }

for R in radios_prueba:
    i=radios_prueba.index(R)
    path='tests-'+str(i)+'/'
    # Se ajustan parametros en diccionario
    parametros["TESTDIR"]=path
    parametros['OUTPUTFILE']='results_'+path[:-1]+'.csv'
    parametros['PARAM1']=R
    
    #Se actualizan parametros de external_solver
    completar_archivo_entrada(parametros)

    #Se calcula energia de solvatacion con external_solver
    os.system('python external_solver.py')

    gc.collect()

print('Total ellapsed time:',time()-START_TIME,'seconds')