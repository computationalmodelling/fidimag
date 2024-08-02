#!/bin/bash
set -e

# This script installs any SUNDIALS version starting with version 2.6.0
# when SUNDIALS moved to a CMake-based installation. Will install locally.
# It may need environment variables to work, like `export CC=gcc` in ARCHER.

# Github release from Sundials repository
# https://github.com/LLNL/sundials
SUNDIALS_TAG=v6.6.1
SUNDIALS=sundials-6.6.1

# Make sure CMake is installed, since SUNDIALS requires it.
type cmake >/dev/null 2>&1 || { printf "CMake required to build SUNDIALS. You can install it by typing: \nsudo apt install cmake\n"; exit 1;}


HERE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
FIDIMAG_DIR="$(dirname "$HERE_DIR")"

if [ -z $LIBS_DIR ]; then LIBS_DIR=${FIDIMAG_DIR}/local; fi

echo "Will install SUNDIALS to "${LIBS_DIR}" using CC="${CC}"."
mkdir -p ${LIBS_DIR}
cd ${LIBS_DIR}

download_and_cmake_install() {
    # $1 tag of the package
    # $2 name of the package
    # $3 URL where ${1}.tar.gz can be obtained
    # $4 configure options
    if [ ! -e ${2}.tar.gz ]; then
        echo "Downloading "${3}/${1}/${2}.tar.gz"."
        wget -q -O ${2}.tar.gz ${3}/${1}/${2}.tar.gz
        tar -xzf ${2}.tar.gz
    fi;

    echo "Configuring "${2}"."
    if [ -e ${2}_build ]; then
        rm -r ${2}_build
    fi;
    mkdir -p ${2}_build
    cd ${2}_build
    cmake ${4} ../${2}

    echo "Compiling and installing "${2}"."
    {
        make -j2
        make install
    } > /dev/null

    echo "Cleaning up."
    cd ..
    rm -rf ${2}
    rm -rf ${2}_build

    cd ${HERE_DIR}
    echo "Done."
}

download_and_cmake_install \
    ${SUNDIALS_TAG} \
    ${SUNDIALS} \
    https://github.com/LLNL/sundials/releases/download \
    "-DBUILD_STATIC_LIBS=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX="${LIBS_DIR}" -DEXAMPLES_ENABLE=OFF -DENABLE_LAPACK=ON -DENABLE_OPENMP=ON"
