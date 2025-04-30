#!/bin/bash
set -e  # Exit if any command fails

echo "Updating packages..."
sudo apt update && sudo apt upgrade -y

echo "Installing Miniconda..."
cd ~
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh

echo " Initializing Conda..."
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init


echo "Installing Mamba into base environment..."
conda install -y -n base -c conda-forge mamba
echo " Initializing Mamba for shell integration..."
mamba shell init --shell bash --root-prefix="$HOME/miniconda3"
eval "$(mamba shell hook --shell bash)"

echo "Installing unzip and zip"
sudo apt install -y unzip zip

echo " Installing Java (for Nextflow)..."
echo "Installing sdk..."
curl -s https://get.sdkman.io | bash
source "$HOME/.sdkman/bin/sdkman-init.sh"
sdk install java 17.0.9-tem




echo " Installing Nextflow..."
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/

echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
export PATH="$HOME/.local/bin:$PATH"




echo "Configuring Conda channels for bioinformatics..."
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict


echo "Creating ngs environment..."

mamba create -y -n ngs -c conda-forge -c bioconda \
  samtools \
  fastqc \
  bwa \
  trimmomatic \
  sra-tools


echo "Done! >> run source .bashrc >> conda activate ngs >> run biosetup.sh"

