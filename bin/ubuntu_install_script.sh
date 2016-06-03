bash install.sh
sudo bash install-ubuntu-packages.sh
make -C ..
echo "export PYTHONPATH=$(dirname $PWD):\$PYTHONPATH" >> /home/$USER/.bashrc
