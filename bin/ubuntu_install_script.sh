echo "Super user authentication required to add packages (you will be \
prompted for confirmation)."
sudo bash install-ubuntu-packages.sh
bash install-fftw.sh
bash install-sundials-2.5.sh
sudo pip install cython --upgrade
pushd .. > /dev/null
make

# Adds Fidimag environment variables to profile.d, if they're not already
# there.
FIDIMAG_PROFILE_PATH=/etc/profile.d/fidimag.sh
if [ ! -e "$FIDIMAG_PROFILE_PATH" ]; then
    echo "Super user authentication required to add paths."
    sudo mkdir --parents "$(dirname $FIDIMAG_PROFILE_PATH)"
    sudo bash -c "echo \"export PYTHONPATH=$PWD/:\\\$PYTHONPATH\"\
> $FIDIMAG_PROFILE_PATH"
    sudo bash -c "echo \"export LD_LIBRARY_PATH=$PWD/local/lib:\
\\\$LD_LIBRARY_PATH\" >> $FIDIMAG_PROFILE_PATH"
    sudo chmod 0644 "$FIDIMAG_PROFILE_PATH"
    echo "Path written to $FIDIMAG_PROFILE_PATH."
else
    echo "Path added previously. Skipping."
fi
