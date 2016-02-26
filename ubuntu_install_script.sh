sudo bash bin/install.sh
sudo bash bin/install-ubuntu-packages.sh
sudo pip install cython --upgrade
make
echo "PYTHONPATH=/home/$PWD/fidimag.sh" >> /home/$USER/.bashrc
