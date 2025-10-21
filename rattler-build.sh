# build matrix: 
# python: 3.10, 3.11, 3.12
# libboost: 1.84, 1.86
# openmm: 8.1, 8.2, 8.3


rm -rf output/ # clean previous builds
rm -rf rattl*.log
for python_ver in 3.10 3.11 3.12; do
  for boost_ver in 1.84 1.86; do
    for omm_ver in 8.1 8.2 8.3; do
      rattler-build build --recipe conda/meta.yaml --variant python=${python_ver} --variant libboost=${boost_ver} --variant openmm=${omm_ver} >& rattler-build-${python_ver}-boost${boost_ver}-omm${omm_ver}.log
    done
  done
done

