# SpClust
SpCLUST is a package for divergent nucleotide sequences clustering. Contrarely to traditional clustering methods that focuses on the speed of clustering highly similar sequences, SpCLUST uses a Machine Learning Gaussian Mixture Model and targets the clustering accuracy of divergent sequences with the best possible speed.
The current version of SpCLUST uses Edgar, R.C.'s MUSCLE module (www.drive5.com) for sequences alignment.

# Prerequisite
SpCLUST uses MPI for parallel computation and the executable building for the installation package. Below are some basic instructions for installing MPI on your system.
- For Linux users:
  •	install mpich
  •	install openmpi
  •	install openmpi-devel
  •	echo "export PATH=$PATH:/usr/lib64/openmpi/bin" >> ~/.bashrc
- For Windows users:
  •	Download MS-MPI SDK and Redist installers from Microsoft's website: https://msdn.microsoft.com/en-us/library/bb524831.aspx
  •	Install the downloaded packages

# Installation on Linux
- Get the installation package from the "Linux" folder in our repository: "wget https://github.com/johnymatar/SpCLUST/raw/master/Linux/install.tar.xz"
- Extract the package: "tar -xvf install.tar.xz"
- Run the following commands: "cd spclust", "./configure", "make"
- Run the following command as a sudoer: "make install"
- You can now call the executables from anywhere with the desired arguments
- For serial computation use "spclust with the desired arguments
- For parallel computation use "mpispclust" with the desired arguments
- To use the graphical interface, install mono (run "apt install mono-complete" as a sudoer) and then call "guispclust"

# Usage without installation on Linux
- Get the standalone package from the "Linux" folder in our repository: "wget https://github.com/johnymatar/SpCLUST/raw/master/Linux/standalone.tar.xz"
- Extract the package: "tar -xvf standalone.tar.xz"
- Keep the extracted files together in a same directory and, for each use, browse to that directory from the console: e.g. "cd ~/spclust"
- For serial computation use "./spclust" with the desired arguments
- For parallel computation use "./mpispclust" with the desired arguments
- To use the graphical interface, install mono (run "apt install mono-complete" as a sudoer) and then call "./guispclust"

# Usage without installation on M.S. Windows
- As a prerequisite, install Python 2.7 along with its following libraries: numpy, matplotlib, scikit-learn, scipy, cogent, BioPython
- Download the files from the "Windows" folder in our repository
- Put them together in a same directory, go to that directory from the command line
- Call "spclust" with the desired arguments from there
- For parallel computation call "mpispclust"
- For a graphical interface run "guispclust"

# Integration into Galaxy
Follow the instructions in the README file in the Galaxy folder

# Current version features
- Test for alignment failure
- Parallel computation for the distance matrix using MPI.
- Checking the presence of all required modules and files prior run
- Supports arguments by using -mdist -in -out -alignMode -gapOpen -gapExtend for spclust:
 available mdist are EDNAFULL, BLOSUM62, and PAM250
 available alignMode are: fast, moderate, maxPrecision
- Calls LaplacianAndGMM from python sub-module spclustGMM.py on Windows
- LaplacianAndGMM python sub-module is converted to spclustGMM executable on Linux for easier use
- Runs, via the executable, from anywhere from the original directory without the need of environment variables. Needs all the dependent file to be present in its directory (spclust, mpispclust, muscle executable, and spclustGMM)
- Can enable figures displaying from the python sub-module.
- Compatible for Galaxy integration
