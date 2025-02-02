#!/bin/bash

PROGRAM="KAntsClustering"

INSTANCE_FOLDER="Instances"
INSTANCE_TENSOR="${INSTANCE_FOLDER}/data30000.bin"

OUTPUT_FOLDER="Result-Tensor"
OUTPUT_FOLDER_TENSOR_8="${OUTPUT_FOLDER}/tensor8/"
OUTPUT_FOLDER_TENSOR_16="${OUTPUT_FOLDER}/tensor16/"
OUTPUT_FOLDER_TENSOR_32="${OUTPUT_FOLDER}/tensor32/"
OUTPUT_FOLDER_TENSOR_64="${OUTPUT_FOLDER}/tensor64/"
OUTPUT_FOLDER_TENSOR_128="${OUTPUT_FOLDER}/tensor128/"
OUTPUT_FOLDER_TENSOR_256="${OUTPUT_FOLDER}/tensor256/"

rm -rf ${OUTPUT_FOLDER}
mkdir -p ${OUTPUT_FOLDER_TENSOR_8} ${OUTPUT_FOLDER_TENSOR_16} ${OUTPUT_FOLDER_TENSOR_32} ${OUTPUT_FOLDER_TENSOR_64} ${OUTPUT_FOLDER_TENSOR_128} ${OUTPUT_FOLDER_TENSOR_256}

echo "Running 8 clusters: ${INSTANCE_TENSOR}"
./${PROGRAM} -b ${INSTANCE_TENSOR} 30000 64 8 30 5000 ${OUTPUT_FOLDER_TENSOR_8} >> "${OUTPUT_FOLDER_TENSOR_8}console.out"

echo "Running 16 clusters: ${INSTANCE_TENSOR}"
./${PROGRAM} -b ${INSTANCE_TENSOR} 30000 64 16 30 5000 ${OUTPUT_FOLDER_TENSOR_16} >> "${OUTPUT_FOLDER_TENSOR_16}console.out"

echo "Running 32 clusters: ${INSTANCE_TENSOR}"
./${PROGRAM} -b ${INSTANCE_TENSOR} 30000 64 32 30 5000 ${OUTPUT_FOLDER_TENSOR_32} >> "${OUTPUT_FOLDER_TENSOR_32}console.out"

echo "Running 64 clusters: ${INSTANCE_TENSOR}"
./${PROGRAM} -b ${INSTANCE_TENSOR} 30000 64 64 30 5000 ${OUTPUT_FOLDER_TENSOR_64} >> "${OUTPUT_FOLDER_TENSOR_64}console.out"

echo "Running 128 clusters: ${INSTANCE_TENSOR}"
./${PROGRAM} -b ${INSTANCE_TENSOR} 30000 64 128 30 5000 ${OUTPUT_FOLDER_TENSOR_128} >> "${OUTPUT_FOLDER_TENSOR_128}console.out"

echo "Running 256 clusters: ${INSTANCE_TENSOR}"
./${PROGRAM} -b ${INSTANCE_TENSOR} 30000 64 256 30 5000 ${OUTPUT_FOLDER_TENSOR_256} >> "${OUTPUT_FOLDER_TENSOR_256}console.out"
