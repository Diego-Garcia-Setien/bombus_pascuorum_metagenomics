#!/bin/bash

#Primero creamos la carpeta para guardar los resultados de fastp

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
		--cut_front \
		--cut_tail \
		--html "fastp_results/${base}_report.html" \
		--json "fastp_results/${base}_report.json" \
		--failed_out "fastp_failed/${base}_failed_R1.fq.gz" \
		--failed_out_R2 "fastp_failed/${base}_failed_R2.fq.gz" \
		--thread 4
