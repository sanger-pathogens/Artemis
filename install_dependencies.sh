#!/bin/bash

set -x
set -e

start_dir=$(pwd)

EMBOSS_VERSION="6.6.0"

EMBOSS_DOWNLOAD_URL="ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-${EMBOSS_VERSION}.tar.gz"

# Make an install location
if [ ! -d 'build' ]; then
  mkdir build
fi
cd build
build_dir=$(pwd)

# DOWNLOAD ALL THE THINGS
download () {
  url=$1
  download_location=$2

  if [ -e $download_location ]; then
    echo "Skipping download of $url, $download_location already exists"
  else
    echo "Downloading $url to $download_location"
    wget $url -O $download_location
  fi
}

download $EMBOSS_DOWNLOAD_URL "${build_dir}/emboss-${EMBOSS_VERSION}.tgz"

# Build all the things
cd $build_dir

## Emboss
emboss_dir=$(pwd)/EMBOSS-${EMBOSS_VERSION}
if [ ! -d $emboss_dir ]; then
  tar xzf emboss-${EMBOSS_VERSION}.tgz
fi
cd $emboss_dir
if [ -e "${emboss_dir}/build/bin/restrict" ]; then
  echo "Already built Emboss; skipping build"
else
  mkdir build
  ./configure --prefix ${emboss_dir}/build
  make
  make install
fi

export EMBOSS_ROOT=${emboss_dir}/build

cd $start_dir

set +x
set +e
