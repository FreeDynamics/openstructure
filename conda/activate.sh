PY_VERSION=`python -c 'import sys;print("%d.%d"%(sys.version_info.major,sys.version_info.minor))'`
old_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$CONDA_PREFIX/lib64:$LD_LIBRARY_PATH
export PYTHONPATH=$CONDA_PREFIX/lib64/python$PY_VERSION/site-packages
export OST_ROOT=$CONDA_PREFIX
