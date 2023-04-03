# Programa para resolver energia de solvatacion de una molecula externamente.
# Para su correcto funcionamiento requiere que el archivo input_config.prm se encuentre en el mismo directorio

def leer_archivo_entrada(prm_file):
    """Esta funcion lee un archivo de parametros de entrada para el solver y retorna un diccionario, de llaves predefinidas.
    El formato de cada linea del archivo de paramatros de entrada `input_config.prm` es PARAMETRO=VALOR"""
    #Diccionario de parametros por defecto
    parametros={
        "PQRFILE":'',
        "N_TESTS":40,
        "MOL2DIR":'', 
        "MESH_GENERATOR":'msms',
        "SHAKE_MODEL":'m2',
        "PARAM1":0.2,
        "PARAM2":0.2,
        "PRINT":False,
        "OUTPUTFILE":''
        }
    
    # A continuacion se lee el archivo de configuracion de entrada para actualizar diccionario
    arch=open(prm_file)
    for line in arch:
        if '=' in line:
            parametro,valor=line.strip().split('=')
            if parametro in ['PARAM1','PARAM2']:
                parametros[parametro]=float(valor)
            elif parametro in ['N_TESTS']:
                parametros[parametro]=int(valor)
            elif parametro=='PRINT':
                parametros[parametro]=bool(valor)
            else:
                parametros[parametro]=valor
    arch.close()
    return parametros

# Comienza programa para calcular energia de solvataciones de molecula informada en input_config.prm del mismo directorio


#Librerias importantes
import os
import numpy as np
import bem_electrostatics
from time import time
from shake_functions import agitar_m2,diccionario_ID_nombre,new_name,lista2pqr
import gc
from scipy.sparse import csr_matrix, SparseEfficiencyWarning
import warnings
warnings.simplefilter("ignore",SparseEfficiencyWarning)

#Se carga diccionario con datos de input_config.prm

parametros=leer_archivo_entrada('input_config.prm')

#Se carga diccionario de valores de epsilon in moleculas de Mobley
from utilities import dict_epsilon

#Se comienza escritura archivo de salida
csv_data=open(parametros['OUTPUTFILE'],'w')
csv_data.write('ITERATION,SOLV. ENERGY,GMRES Iterations,Elapsed time,Number of elements\n')
csv_data.close()

# Se agitan las moleculas utilizando los parametros iniciales
if parametros['TESTDIR'][:-1] not in os.listdir(): os.mkdir(parametros['TESTDIR'][:-1])
diccionario=diccionario_ID_nombre(parametros['MOL2FILE'])

mol_name=parametros['MOL2FILE'].split('/')[-1].split('.')[0]

for i in range(parametros['N_TESTS']):
    lista_pqr=agitar_m2(parametros['PQRFILE'],parametros['MOL2FILE'],diccionario,parametros['PARAM2'],parametros['PARAM1'])
    pqr_test_file=new_name(parametros['TESTDIR']+mol_name+'.pqr',i,parametros['N_TESTS'])
    lista2pqr(lista_pqr,pqr_test_file)

#Se eliminan variables innecesarias
del diccionario
del lista_pqr
gc.collect() 
for i in range(parametros['N_TESTS']):
    pqr_test_file=new_name(parametros['TESTDIR']+mol_name+'.pqr',i,parametros['N_TESTS'])
    try:
        if parametros['PRINT']:
            start_time=time()
        protein=bem_electrostatics.solute(pqr_test_file,mesh_generator=parametros['MESH_GENERATOR'],mesh_density=3)
        protein.gmres_tolerance=1e-5
        protein.gmres_max_iterations=150
        protein.kappa=0
        protein.ep_in = dict_epsilon[mol_name]
        
        protein.calculate_solvation_energy()
        if parametros['PRINT']: 
            ET=time()-start_time
            print('INFO: {vmol} i={0},\t{1} (kcal/mol),{2} GMRES its.,{3} (s),\t{4} mesh elmnts.'.format(i, 
            round(protein.solvation_energy,6),
            protein.solver_iteration_count, 
            round(ET,3),
            protein.mesh.number_of_elements,
            vmol=mol_name))
        csv_data=open(parametros['OUTPUTFILE'],'a')
        csv_data.write('{0},{1},{2},{3},{4}\n'.format(
            i,
            protein.solvation_energy,
            protein.solver_iteration_count,
            round(ET,3),
            protein.mesh.number_of_elements))
        csv_data.close()
    except KeyboardInterrupt:
        print('Interrupcion')
        error_file=open('error_file.txt','w')
        error_file.close()
        break
    except OSError:
        #Algunas veces ocurren errores y mesh_temp queda ocupado, con esto se evita arrastrar el error a las moleculas siguientes
        if 'mesh_temp' in os.listdir():
            os.chdir('mesh_temp')
            for arch in os.listdir():
                os.remove(arch)
            os.chdir('..')
            os.rmdir('mesh_temp')
    except Exception:
        import sys
        exc_type,value,traceback=sys.exc_info()
        print(exc_type.__name__)
        print(traceback)

        #Algunas veces ocurren errores y mesh_temp queda ocupado, con esto se evita arrastrar el error a las moleculas siguientes
        if 'mesh_temp' in os.listdir():
            os.chdir('mesh_temp')
            for arch in os.listdir():
                os.remove(arch)
            os.chdir('..')
            os.rmdir('mesh_temp')


#Se eliminan archivos de prueba
for i in range(parametros['N_TESTS']):
    pqr_test_file=new_name(parametros['TESTDIR']+mol_name+'.pqr',i,parametros['N_TESTS'])
    os.remove(pqr_test_file)
