#!/bin/bash
set -e

# This script installs any SUNDIALS version starting with version 2.6.0
# when SUNDIALS moved to a CMake-based installation. Will install locally.
# It may need environment variables to work, like `export CC=gcc` in ARCHER.

SUNDIALS=sundials-2.6.2

# Make sure CMake is installed, since SUNDIALS requires it.
type cmake >/dev/null 2>&1 || { printf "CMake required to build SUNDIALS. You can install it by typing: \nsudo apt install cmake\n"; exit 1;}


HERE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
FIDIMAG_DIR="$(dirname "$HERE_DIR")"

if [ -z $LIBS_DIR ]; then LIBS_DIR=${FIDIMAG_DIR}/local; fi

echo "Will install SUNDIALS to "${LIBS_DIR}" using CC="${CC}"."
mkdir -p ${LIBS_DIR}
cd ${LIBS_DIR}

download_and_cmake_install() {
    # $1 name of the package
    # $2 URL where ${1}.tar.gz can be obtained
    # $3 configure options
    if [ ! -e ${1}.tar.gz ]; then
        echo "Downloading "${1}"."
        wget -q ${2}/${1}.tar.gz
    fi;

    if [ ! -e ${1} ]; then
        tar -xzf ${1}.tar.gz

        echo "Configuring "${1}"."
        mkdir ${1}_build
        cd ${1}_build
        cmake ${3} ../${1}

        echo "Compiling and installing "${1}"."
        {
            make
            make install
        } > /dev/null

        echo "Cleaning up."
        cd ..
        rm -rf ${1}
        rm -rf ${1}_build

        cd ${HERE_DIR}
        echo "Done."
    fi;
}

download_and_cmake_install \
    ${SUNDIALS} \
    http://computation.llnl.gov/projects/sundials-suite-nonlinear-differential-algebraic-equation-solvers/download \
    "-DBUILD_STATIC_LIBS=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX="${LIBS_DIR}" -DEXAMPLES_ENABLE=OFF -DLAPACK_ENABLE=ON -DOPENMP_ENABLE=ON --enable-sse2"
