sudo bash install-ubuntu-packages.sh
bash install-fftw.sh
bash install-sundials-2.5.sh
sudo pip install cython --upgrade
cd ..
make
echo "export PYTHONPATH=$(dirname $PWD):$PYTHONPATH" >> $HOME/.bashrc
