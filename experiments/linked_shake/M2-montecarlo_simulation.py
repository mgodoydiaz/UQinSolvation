# Programa para correr Analisis de Montecarlo en moleculas de Mobley et. al. 2009
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
from shake_functions import agitar_m2,diccionario_ID_nombre,new_name,lista2pqr

# MSMS es archivo ejecutable

os.system('chmod +x bem_electrostatics/ExternalSoftware/MSMS/msms')
#!chmod +x bem_electrostatics/ExternalSoftware/MSMS/msms

#Opciones de Montecarlo y directorios
n_tests=100
path="tests/"
mol2_path="charged_mol2files/"
pqr_path="mobley_test_pqr/"
generador_mesh='msms'
modelo_agitacion='m2'
imprimir_resultados=True
#Diccionario de parametros individuales de cada test
parametros={
        "PQRFILE":'',
        "N_TESTS":n_tests,
        "MOL2FILE":'', 
        "MESH_GENERATOR":generador_mesh,
        "SHAKE_MODEL":modelo_agitacion,
        "PARAM1":0.2,
        "PARAM2":0.2,
        "PRINT":imprimir_resultados,
        "OUTPUTFILE":'',
        "TESTDIR":path
        }
# Se importan nombres de moleculas de interes
from utilities import moleculas

#En caso de pausar ejecucion, se guarda un archivo .txt con el historial de montecarlos realizados para continuar desde el ultimo punto

hist_file=open('historial.txt','a') #si no esta creado, se crea
# Cada linea tendra la forma molecula-p_r-alpha(4 decimales)
hist_file.close()

# Ejecucion programa principal S. Montecarlo
PORC=0
for p_r in [0.05, 0.03,0.01]:
    PORC+=1
    ALPHA=0
    for alpha in [np.pi/4,np.pi/6]:
        ALPHA+=1
        dir_results='R{}-{}/'.format(PORC,ALPHA)
        if dir_results[:-1] not in os.listdir(): os.mkdir(dir_results[:-1])
        for molecula in moleculas:

            # Comprobacion en caso de que ya se haya ejecutado esta prueba
            skip_test=False
            hist_line='{}-{}-{}\n'.format(molecula,p_r,round(alpha,4))
            hist_file=open('historial.txt')
            for line in hist_file:
                if hist_line == line:
                    if imprimir_resultados:
                        print(hist_line.strip()+' SKIPPED')
                    skip_test=True
            hist_file.close()            
            if skip_test: continue
            # Se ajustan parametros en diccionario
            parametros["TESTDIR"]=path
            parametros['OUTPUTFILE']=dir_results+'results_'+molecula+'.csv'
            parametros['PARAM1']=p_r
            parametros['PARAM2']=alpha
            parametros['PQRFILE']=pqr_path+molecula+'.pqr'
            parametros['MOL2FILE']=mol2_path+molecula+'.mol2'
            
            #Se actualizan parametros de external_solver
            completar_archivo_entrada(parametros)

            #Se calcula energia de solvatacion con external_solver
            os.system('python external_solver.py')
            #!python external_solver.py
            if 'error_file.txt' in os.listdir():
                break
            #Se agrega molecula y valores de test a historial
            hist_file=open('historial.txt','a')
            hist_file.write(hist_line)
            hist_file.close()

            # Si es posible se limpia la memoria RAM
            gc.collect()
        if 'error_file.txt' in os.listdir(): break
    if 'error_file.txt' in os.listdir():
        os.system('rm error_file.txt')
        break
