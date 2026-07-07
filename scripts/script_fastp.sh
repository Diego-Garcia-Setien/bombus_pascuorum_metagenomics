#!/bin/bash

#Como bamos a trabajar con SLURM, tenemos que indicar las instrucciones

#Con esto indicamos el nombre que tendrá la tarea en la lista de espera
#SBATCH --job-name=bp_fastp_results_tfm

#Con esto le decimos que los errores irán en la carpeta logs/ en un archivo .err
#SBATCH --error=logs/%x-%j.err

#Lo que salga bien irá a la carpeta logs/ en un archivo .out
#SBATCH --output=logs/%x-%j.out

# %x se reemplaza por el nombre del trabajo y %j por el número de ID único que Slurm le asigne a la tarea 

#Las CPUs están divididas en secciones o particiones, le decimos que asigne la partición general, la estandar
#SBATCH --partition=general

#Le asignamos la prioridad al trabajo, regular es por defecto
#SBATCH --qos=regular

#Le decimos al cluster cuandos núcleos de procesador queremos
#SBATCH --cpus-per-task=8

#Le decimos que esos 8 núcleos estén físicamente en la misma máquina
#SBATCH --nodes=1

#Le decimos que solo vamos a lanzar una tarea por nodo
#SBATCH --ntasks-per-node=1

#Le decimos el límite de tiempo que queremos que tarde el script
#SBATCH --time=03:00:00

#Le decimos la memoria RAM para dedicar al trabajo, en este caso 12GB o 12000MB
#SBATCH --mem=12000

#Le decimos que cree un Job Array o matriz de tareas, para que el cluster haga 93 tareas al mismo tiempo
#SBATCH --array=1-93%93

#Creamos la carpeta para guardar los resultados de fastp

mkdir -p ./bombus_pascuorum_metagenomics/data/fastp_results

#Creamos también la carpeta donde se guardaran las lecturas que no superan los filtros

mkdir -p ./bombus_pascuorum_metagenomics/data/fastp_failed

#Vamos a utilizar fastp
#Tenemos que trabajar con los archivos forward y reverse
#Primero buscamos los archivos forward
for fwd in *_1.fq.gz; do
	
	#Asegurar que el archivo existe
	[ -e "$fwd" ] || continue

	#Identificar el archivo reverse correspondiente (_2.fastq.gz)
	rev="${fwd/_1.fq.gz/_2.fq.gz}"

	#Cambiamos el nombre corto y limpio para la muestra resultante
	base=$(basename "$fwd" _1-fq.gz)

	echo "Procesando muestra: $base"
	echo " -> [R1]: $fwd"
	echo " -> [R2]: $rev"

	#Ejecutamos fastp emparejado
	fastp \
		-i "$fwd" \
		-I "$rev" \
		-o "fastp_results/${base}_clean_R1.fq.gz" \
		-O "fastp_results/${base}_clean_R2.fq.gz"\
		--detect_adapter_for_pe \
		--trim_poly_g \
		--trim_poly_x \
		--cut_front 10 \OB
		--cut_tail \
		--n_base_limit \
		--qualified_quality_phred 33 \
		--html "fastp_results/${base}_report.html" \
		--json "fastp_results/${base}_report.json" \
		--failed_out "fastp_failed/${base}_failed_R1.fq.gz" \
		--failed_out_R2 "fastp_failed/${base}_failed_R2.fq.gz" \
		--thread 8
