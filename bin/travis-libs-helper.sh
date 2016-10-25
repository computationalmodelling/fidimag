if [ ! -d $HOME/libraries ]; then
  echo "No cache - building FFTW & Sundials"
  export LIBS_DIR=$HOME/libraries
  bash bin/install-sundials.sh
  bash bin/install-fftw.sh
else
  echo "FFTW & Sundials cache found!"
fi
