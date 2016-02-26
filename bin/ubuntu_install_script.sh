sudo bash install.sh
sudo bash install-ubuntu-packages.sh
sudo pip install cython --upgrade
make
echo "export PYTHONPATH=$(dirname $PWD):$PYTHONPATH" >> /home/$USER/.bashrc
