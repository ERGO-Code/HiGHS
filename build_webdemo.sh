script_path=$(realpath $(dirname $0))
build_dir="build_webdemo"

# interesting for Node.js: -sNODERAWFS=1, allows direct os filesystem access
LDFLAGS=$(echo $(cat << EOF
-sINVOKE_RUN=0
-sMODULARIZE=1
-sEXPORTED_RUNTIME_METHODS=['FS','callMain']
-sEXPORT_NAME='HiGHS'
-sSINGLE_FILE
-sALLOW_MEMORY_GROWTH=1
--shell-file ${script_path}/app/highs_webdemo_shell.html
EOF
))

# compare https://stackoverflow.com/a/6481016
physical_cores=$(grep ^cpu\\scores /proc/cpuinfo | uniq | awk '{print $4}')

cd ${script_path} &&
rm -rf ${build_dir} &&
mkdir ${build_dir} &&
cd ${build_dir} &&
LDFLAGS=$LDFLAGS emcmake cmake .. \
-DCMAKE_BUILD_TYPE=Release -DEMSCRIPTEN_HTML=on &&
emmake make -j ${physical_cores}