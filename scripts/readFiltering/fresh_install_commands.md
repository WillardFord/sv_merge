# Fresh Install

## Conda 

All commands used while installing sv_merge for the first time on a new computer:

```
conda create -n "sv_merge"
conda activate sv_merge

conda install cmake
conda install compilers # mounted cmake library directories into conda. Not sure if packages are actually useful.
conda install zlib
conda install libcurl # Theoretically installed with cmake but needed to install again here, not sure if version difference or something to do with order of installing compilers and cmake
conda install autoconf # in pursuit of autoreconf
conda install jansson

cd [/path/to/sv_merge]
mkdir build
cd build
# Worked with some warnings
cmake ..
make -j 128
```

Other files downloaded through above commands:
```
[output from installing cmake]
The following NEW packages will be INSTALLED:

  _libgcc_mutex      conda-forge/linux-64::_libgcc_mutex-0.1-conda_forge 
  _openmp_mutex      conda-forge/linux-64::_openmp_mutex-4.5-2_gnu 
  bzip2              conda-forge/linux-64::bzip2-1.0.8-h4bc722e_7 
  c-ares             conda-forge/linux-64::c-ares-1.32.2-h4bc722e_0 
  ca-certificates    conda-forge/linux-64::ca-certificates-2024.7.4-hbcca054_0 
  cmake              conda-forge/linux-64::cmake-3.30.1-hf8c4bd3_0 
  keyutils           conda-forge/linux-64::keyutils-1.6.1-h166bdaf_0 
  krb5               conda-forge/linux-64::krb5-1.21.3-h659f571_0 
  libcurl            conda-forge/linux-64::libcurl-8.8.0-hca28451_1 
  libedit            conda-forge/linux-64::libedit-3.1.20191231-he28a2e2_2 
  libev              conda-forge/linux-64::libev-4.33-hd590300_2 
  libexpat           conda-forge/linux-64::libexpat-2.6.2-h59595ed_0 
  libgcc-ng          conda-forge/linux-64::libgcc-ng-14.1.0-h77fa898_0 
  libgomp            conda-forge/linux-64::libgomp-14.1.0-h77fa898_0 
  libnghttp2         conda-forge/linux-64::libnghttp2-1.58.0-h47da74e_1 
  libssh2            conda-forge/linux-64::libssh2-1.11.0-h0841786_0 
  libstdcxx-ng       conda-forge/linux-64::libstdcxx-ng-14.1.0-hc0a3c3a_0 
  libuv              conda-forge/linux-64::libuv-1.48.0-hd590300_0 
  libzlib            conda-forge/linux-64::libzlib-1.3.1-h4ab18f5_1 
  ncurses            conda-forge/linux-64::ncurses-6.5-h59595ed_0 
  openssl            conda-forge/linux-64::openssl-3.3.1-h4bc722e_2 
  rhash              conda-forge/linux-64::rhash-1.4.4-hd590300_0 
  xz                 conda-forge/linux-64::xz-5.2.6-h166bdaf_0 
  zstd               conda-forge/linux-64::zstd-1.5.6-ha6fb4c9_0 

[output from installing compilers]
The following NEW packages will be INSTALLED:

  _sysroot_linux-64_curr_repodata_hack-3-h69a702a_16
  binutils-2.40-h4852527_7
  binutils_impl_linux-64-2.40-ha1999f0_7
  binutils_linux-64-2.40-hb3c18ed_0
  c-compiler-1.7.0-hd590300_1
  compilers-1.7.0-ha770c72_1
  cxx-compiler-1.7.0-h00ab1b0_1
  fortran-compiler-1.7.0-heb67821_1
  gcc-12.4.0-h236703b_0
  gcc_impl_linux-64-12.4.0-hb2e57f8_0
  gcc_linux-64-12.4.0-h6b7512a_0
  gfortran-12.4.0-h236703b_0
  gfortran_impl_linux-64-12.4.0-hc568b83_0
  gfortran_linux-64-12.4.0-hd748a6a_0
  gxx-12.4.0-h236703b_0
  gxx_impl_linux-64-12.4.0-h557a472_0
  gxx_linux-64-12.4.0-h8489865_0
  kernel-headers_linux-64-3.10.0-h4a8ded7_16
  ld_impl_linux-64-2.40-hf3520f5_7
  libgcc-devel_linux-64-12.4.0-ha4f9413_100
  libgfortran5-14.1.0-hc5f4f2c_0
  libsanitizer-12.4.0-h46f95d5_0
  libstdcxx-devel_linux-64-12.4.0-ha4f9413_100
  sysroot_linux-64-2.17-h4a8ded7_16
  tzdata-2024a-h0c530f3_0
```



## Getting new reads and generating all signature matrices again



Run from somewhere where gcloud works properly
Then from `build` directory
```
../scripts/readFiltering/getFASTQ.sh
```
in scripts/readFiltering dir on my laptop
```
./copyChromToZep.sh
```

On zep again:
```
cd sv_merge/scripts/readFiltering
./splitFastq.sh 
./countKmers.sh # Requires conda environment with fastk, I use a different environment than sv_merge
```

Then generate the characteristic matrices.
Requires a few python packages in a new environment.
```
conda install numpy, matplotlib, joblib

python generateROC.py euclidean
python generateROC.py sketch
```







