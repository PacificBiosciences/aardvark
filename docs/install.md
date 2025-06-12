# Installing Aardvark
## From conda
The easiest way to install Aardvark is through [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html).
This will download the source code and build Aardvark locally:

```bash
# create a brand new conda environment and install latest Aardvark
conda create -n aardvark -c bioconda aardvark
# OR install latest into current conda environment
conda install aardvark
# OR install a specific version into current conda environment
conda install aardvark=0.7.2
```

## From pre-compiled binary
We provide a pre-compiled binary file for x86_64 Linux distributions.
Use the following instructions to get the most recent version directly from GitHub:

1. Navigate to the [latest release](https://github.com/PacificBiosciences/Aardvark/releases/latest) and download the tarball file (e.g. `aardvark-{version}-x86_64-unknown-linux-gnu.tar.gz`).
2. Decompress the tar file.
3. (Optional) Verify the md5 checksum.
4. Test the binary file by running it with the help option (`-h`).
5. Visit the [User guide](./user_guide.md) for details on running Aardvark.

### Example with v0.7.2
```bash
# modify this to update the version
VERSION="v0.7.2"
# get the release file
wget https://github.com/PacificBiosciences/Aardvark/releases/download/${VERSION}/aardvark-${VERSION}-x86_64-unknown-linux-gnu.tar.gz
# decompress the file into folder ${VERSION}
tar -xzvf aardvark-${VERSION}-x86_64-unknown-linux-gnu.tar.gz
cd aardvark-${VERSION}-x86_64-unknown-linux-gnu
# optional, check the md5 sum
md5sum -c aardvark.md5
# execute help instructions
./aardvark -h
```
