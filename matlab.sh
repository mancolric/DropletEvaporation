#!/bin/sh
#
# Filename: run_unico.sh
#
# ./run_unico.sh
#
echo "Ejecutando matlab.sh en host: $(hostname)$"
echo "Argumentos recibido: $*"
# run matlab in text mode
exec /mnt/cephfs/software/bin/matlab2023a -nojvm -nodisplay -nosplash -r "$*"
