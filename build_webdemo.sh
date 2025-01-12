# Webdemo
# -------

# HiGHS can directly be compiled into a single HTML file and used
# in a browser. This requires `emscripten` to be installed from their website
# (unfortunately, e.g. `sudo apt install emscripten` in Ubuntu Linux is broken):

#     https://emscripten.org/docs/getting_started/downloads.html

# Then, run

#     sh build_webdemo.sh

# This will create the file `build_webdemo/bin/highs.html`. For fast edit
# iterations run

#     find src app | entr -rs 'make -C build_webdemo highs; echo'

# This will rebuild `highs.html` every time a source file is modified (e.g.
# from Visual Studio Code).

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