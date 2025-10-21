set -e -x

mkdir -p build
cd build
build_dir=`pwd`

$BUILD_PREFIX/bin/cmake .. \
    -D CMAKE_BUILD_TYPE=Release \
    -D OPTIMIZE=ON \
    -D CMAKE_INSTALL_PREFIX=$PREFIX \
    -D CMAKE_PREFIX_PATH=$PREFIX \
    -D ENABLE_INFO=off \
    -D ENABLE_GFX=on \
    -D ENABLE_MM=on \
    -D USE_NUMPY=1 \
    -D EIGEN3_INCLUDE_DIR=$BUILD_PREFIX/include/eigen3 \
    -D TIFF_INCLUDE_DIR=$BUILD_PREFIX/include \
    -D TIFF_LIBRARY=$BUILD_PREFIX/lib/libtiff.so \
    -D PNG_PNG_INCLUDE_DIR=$BUILD_PREFIX/include \
    -D PNG_LIBRARY=$BUILD_PREFIX/lib/libpng.so \
    -D ZLIB_INCLUDE_DIR=$BUILD_PREFIX/include \
    -D ZLIB_LIBRARY=$BUILD_PREFIX/lib/libz.so \
    -D OPEN_MM_LIBRARY=$BUILD_PREFIX/lib/libOpenMM.so \
    -D OPEN_MM_INCLUDE_DIR=$BUILD_PREFIX/include \
    -D OPEN_MM_PLUGIN_DIR=$BUILD_PREFIX/lib/plugins \
    -D CMAKE_CXX_FLAGS="-std=c++14 -fpermissive"

$BUILD_PREFIX/bin/make -j$CPU_COUNT
$BUILD_PREFIX/bin/make install

for CHANGE in "activate" "deactivate"
do
    mkdir -p "${PREFIX}/etc/conda/${CHANGE}.d"
    cp "${RECIPE_DIR}/${CHANGE}.sh" "${PREFIX}/etc/conda/${CHANGE}.d/${PKG_NAME}_${CHANGE}.sh"
done