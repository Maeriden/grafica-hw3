# echo "[$(date +%T)] Compiling libhwlib.a"
# pushd obj/ > /dev/null
# g++ -c -std=gnu++14 -fPIC -O2 -DYGLTF_NO_IMAGE -DYOBJ_NO_IMAGE -DYSCN_NO_IMAGE ../src/image.cpp ../src/scene.cpp ../src/ext/yocto_scn.cpp ../src/ext/yocto_obj.cpp ../src/ext/yocto_gltf.cpp
# popd > /dev/null

# echo "[$(date +%T)] Archiving libhwlib.a"
# ar qcs bin/libhwlib.a obj/image.o obj/scene.o obj/yocto_scn.o obj/yocto_obj.o obj/yocto_gltf.o

defines="-DYGLTF_NO_IMAGE -DYOBJ_NO_IMAGE -DYSCN_NO_IMAGE"
features=""
sse_flags="-msse -msse2 -msse3 -msse4 -mfpmath=sse"
debug_flags="-g -DENABLE_ASSERT=1"
# profile_flags="-pg -no-pie"
link_flags="-Lbin -lhwlib -pthread"

echo "[$(date +%T)] Building model"
g++ -obin/model -std=gnu++14 ${debug_flags} ${profile_flags} ${sse_flags} ${defines} ${features} src/model.cpp ${link_flags}
