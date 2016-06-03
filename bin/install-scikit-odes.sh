#!/bin/bash
set -e

# This script installs the python bindings to SUNDIALS called scikit ODES.
# It requires SUNDIALS 2.6.2 to be installed.

HERE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
FIDIMAG_DIR="$(dirname "$HERE_DIR")"
LIBS_DIR=${FIDIMAG_DIR}/local

echo "Downloading ODES to "${LIBS_DIR}"."
mkdir -p ${LIBS_DIR}
cd ${LIBS_DIR}

git clone git://github.com/bmcage/odes.git odes
cd odes

# Paths to SUNDIALS are hardcoded in odes/scikits/odes/sundials/setup.py.
# The README recommends to edit that file if you want to change the paths...

sed -i \
    "s|LIB_DIRS_SUNDIALS  = \[base_path\,|LIB_DIRS_SUNDIALS = [base_path, '${LIBS_DIR}\/lib',|g" \
    scikits/odes/sundials/setup.py

python setup.py build
sudo python setup.py install

cd ${HERE_DIR}
