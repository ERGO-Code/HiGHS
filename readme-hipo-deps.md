Dependencies install ubuntu / macos 

1. Clone GKLib

git clone https://github.com/KarypisLab/GKlib.git

2. Clone METIS

git clone https://github.com/KarypisLab/GKlib.git

3. Create installs dir 

mkdir installs

4. Install GKlib

cd GKlib
make config prefix=${{runner.workspace}}/installs
make
make install

5. Install METIS

cd METIS
make config prefix=${{runner.workspace}}/installs
make
make install