# Programa para resolver energia de solvatacion de una molecula externamente.
# Para su correcto funcionamiento requiere que el archivo input_config.prm se encuentre en el mismo directorio

def leer_archivo_entrada(prm_file):
    """Esta funcion lee un archivo de parametros de entrada para el solver y retorna un diccionario, de llaves predefinidas.
    El formato de cada linea del archivo de paramatros de entrada `input_config.prm` es PARAMETRO=VALOR"""
    #Diccionario de parametros por defecto
    parametros={
        "PDBFILE":'',
        "N_TESTS":100, 
        "MESH_GENERATOR":'nanoshaper',
        "PARAM1":1,
        "PRINT":False,
        "OUTPUTFILE":'',
        "TESTDIR":''
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
import gc
from scipy.sparse import csr_matrix, SparseEfficiencyWarning
from shake_pdb import new_name
import warnings
warnings.simplefilter("ignore",SparseEfficiencyWarning)
#Se carga diccionario con datos de input_config.prm

parametros=leer_archivo_entrada('input_config.prm')

#Se comienza escritura archivo de salida
csv_data=open(parametros['OUTPUTFILE'],'w')
csv_data.write('ITERATION,SOLV. ENERGY,GMRES Iterations,Elapsed time,Number of elements\n')
csv_data.close()
path=parametros['TESTDIR']
start_time,ET=0,0 # En caso de no necesitar imprimir tiempos por pantalla
for i in range(parametros['N_TESTS']):
    pdb_test_file=path+new_name(parametros['PDBFILE'],i,parametros['N_TESTS'])
    try:
        if parametros['PRINT']:
            start_time=time()
        protein=bem_electrostatics.solute(pdb_test_file,mesh_generator=parametros['MESH_GENERATOR'])
        protein.gmres_tolerance=1e-5
        protein.gmres_max_iterations=150
        protein.calculate_solvation_energy()
        R_res=parametros['PARAM1']
        if parametros['PRINT']: 
            ET=time()-start_time
            print('INFO: R={0}, i:{1}, {2}, {3}, {4} [s],{5}'.format(R_res,i, 
            protein.solvation_energy,
            protein.solver_iteration_count, 
            round(ET,3),
            protein.mesh.number_of_elements))
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

print('')
