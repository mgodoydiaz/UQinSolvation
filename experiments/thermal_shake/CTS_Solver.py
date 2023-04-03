#CTS_Solver Characteristic time solver

#Programa externo para resolver energia de solvatacion de una molecula externamente y saltar fuga de memoria
#Para su correcto funcionamiento requiere que el archivo input_config.prm se encuentre en el mismo directorio

def leer_archivo_entrada(prm_file):
    """Esta funcion lee un archivo de parametros de entrada para el solver y retorna un diccionario, de llaves predefinidas.
    El formato de cada linea del archivo de paramatros de entrada `input_config.prm` es PARAMETRO=VALOR"""
    #Diccionario de parametros por defecto
    parametros = {
        'NAME':'',
        'TESTPATH':'',
        'PQRFILE':'',
        'N_TESTS':100,
        'TIEMPO':0,
        'SIGMA_MOBLEY':0,
        'EP_IN':1,
        'KAPPA':0,
        'DENSIDAD':3,
        'OUTPUT_FILE':''
    }
    
    # A continuacion se lee el archivo de configuracion de entrada para actualizar diccionario
    arch=open(prm_file)
    for line in arch:
        if '=' in line:
            parametro,valor=line.strip().split('=')
            # Se separa transformacion de valores dependiendo del tipo de dato
            if parametro in ['TIEMPO','EP_IN','KAPPA','DENSIDAD']:
                parametros[parametro]=float(valor)
            elif parametro in ['N_TESTS']:
                parametros[parametro]=int(valor)
            else: 
                # Para valores de tipo strings
                parametros[parametro]=valor
    arch.close()
    return parametros

# Comienza programa para calcular energias de solvatacion de molecula informada en input_config.prm del directorio

#Librerias importantes
import os 
import numpy as np
import bem_electrostatics
from time import time
from utilities import masas
from thermal_functions import nombre_atomo, randomX, new_name
import gc
from scipy.sparse import csr_matrix, SparseEfficiencyWarning
import warnings
warnings.simplefilter("ignore",SparseEfficiencyWarning)

#Se carga diccionario con datos de input_config.prm
parametros=leer_archivo_entrada('input_config.prm')

# Se comienza escritura de archivo de salida
output_file = parametros['OUTPUT_FILE']
csv_data=open(output_file,'w')
csv_data.write('ITERATION,SOLV. ENERGY,GMRES Iterations,Elapsed time,Number of elements\n')
csv_data.close()

# Average Thermal Length: definicion de funciones y creacion de diccionario
ATL_dic = {}
kB = 1.3806e-23 # Constante de Boltzmann
t  = parametros['TIEMPO'] # Tiempo característico
T  = 298 # Temperatura en Kelvin
vt = lambda m,T: np.sqrt(kB*t/m) # thermal velocity en funcion de masa y temperatura (298K por defecto)
L  = lambda m,T: vt(m,T)*t*1e10  #1e10 convierte la unidad a Angstroms
for atomo,masa in masas.items():
    ATL_dic[atomo] = L(masa,T)

# Se agitan las moleculas utilizando los parametros iniciales
path = parametros['TESTPATH']

n_tests         = parametros['N_TESTS']
mainfile        = parametros['PQRFILE']
nombre_molecula = parametros['NAME']

# Se cargan datos de archivo pqr
remarks=""""""
atoms=""""""
end=""""""
arch=open(mainfile)
post_atoms=False
for line in arch:
    data=line.split()
    if data[0]!='ATOM' and not post_atoms:
        remarks+=line
    elif data[0]=='ATOM':
        post_atoms=True
        atoms+=line
    else:
        end+=line
arch.close()

#Ahora es conveniente definie la función para crear un archivo aleatorio en base a la agitación térmica

def shake_file(mainfile,i,path,remarks=remarks): 
    file = open( os.path.join(path,new_name(nombre_molecula+'.pqr',i,n_tests)) ,'w')
    file.write(remarks)
    list_atoms=atoms.split('\n')
    for i in range(len(list_atoms)):
        linea=list_atoms[i].split()
        if len(linea)==0: continue
        atomo=nombre_atomo(linea[2],True) # 2 para pqr, -1 para pdb
        _x=list(map(float,linea[5:8])) # 5:8 pqr, 6:9 pdb
        _x=randomX(_x,ATL_dic[atomo])
        x,y,z=list(map(str,_x))
        lineareal=list(list_atoms[i])
        extension= mainfile.split('.')[-1]
        # Si es pdb 
        if extension == 'pdb':
            c=0
            for j in z[::-1]+' ':
                lineareal[53-c]=j
                c+=1
            c=0
            for j in y[::-1]+' ':
                lineareal[45-c]=j
                c+=1
            c=0
            for j in x[::-1]+' ':
                lineareal[37-c]=j
                c+=1
            file.write(''.join(lineareal)+'\n')
        # Si es pqr
        elif extension=='pqr':
            linea[5:8]=list(map(str,_x))
            file.write('\t'.join(linea)+'\n')
    file.write(end)
    file.close()
    return None

solvataciones = [] # Contiene energias de solvatacion 
desviaciones  = [] # Contiene desviaction estandar de las energias de solvatacion

for i in range(n_tests):
    # Se crea nuevo archivo pqr
    test_file=os.path.join(path,new_name(nombre_molecula+'.pqr',i,n_tests))
    shake_file(test_file,i,path)
    
    # Se verifica si la desviacion estandar se ha estancado en un valor, en este caso se termina la ejecucion

    if len(solvataciones)>40: # Cantidad minima de pruebas de montecarlo es de 40
        test_desv=desviaciones[-20:] #Se toma una ventana de las ultima 20 pruebas
        min_desv=min(test_desv)
        max_desv=max(test_desv)
        delta=max_desv-min_desv
        if delta<0.002:
            break

    try:
        start_time=time()
        
        # Carga de archivo pqr 
        molecule=bem_electrostatics.solute(test_file,mesh_generator='msms',mesh_density=parametros['DENSIDAD'])
        
        # Parametros del objeto molecule
        molecule.gmres_tolerance        = 1e-5
        molecule.gmres_max_iterations   = 400
        molecule.ep_in                  = parametros['EP_IN']
        molecule.kappa                  = parametros['KAPPA']

        #Calculo de energia de solvatacion electrostatica
        molecule.calculate_solvation_energy()
        
        #Impresion por pantalla
        ET=time()-start_time
        print('INFO: {0}\ti:{1},\t {2}\t(kcal/mol), {3} it.,\t {4} [s],{5} elem.'.format(output_file.split('.')[0],i, 
        molecule.solvation_energy,
        molecule.solver_iteration_count, 
        round(ET,3),
        molecule.mesh.number_of_elements))

        # Escritura en archivo de salida
        csv_data=open(output_file,'a')
        csv_data.write('{0},{1},{2},{3},{4}\n'.format(
            i,
            molecule.solvation_energy,
            molecule.solver_iteration_count,
            round(ET,3),
            molecule.mesh.number_of_elements))
        csv_data.close()
        
        # Se agrega energia de solvatacion
        if not np.isnan(molecule.solvation_energy):
            solvataciones.append(molecule.solvation_energy)
            desviaciones.append(np.std(solvataciones))

    # Pausa o interrupcion de usuario
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
