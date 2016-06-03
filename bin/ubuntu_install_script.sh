sudo bash install-ubuntu-packages.sh
bash install-fftw.sh
bash install-sundials.sh
bash install-scikit-odes.sh
make -C ..
echo "export PYTHONPATH=$(dirname $PWD):\$PYTHONPATH" >> /home/$USER/.bashrc
