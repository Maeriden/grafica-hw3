./bin/model -o in/empty.obj empty
./bin/raytrace -r 720 -s 3 -o out/empty.png in/empty.obj
# [ -f out/empty.png ] && xdg-open out/empty.png&

./bin/model -o in/simple.obj simple
./bin/raytrace -r 720 -s 3 -o out/simple.png in/simple.obj
# [ -f out/simple.png ] && xdg-open out/simple.png&

./bin/model -o in/displace.obj displace
./bin/raytrace -r 720 -s 3 -o out/displace.png in/displace.obj
# [ -f out/displace.png ] && xdg-open out/displace.png&

./bin/model -o in/instances.obj instances
./bin/raytrace -r 720 -s 3 -o out/instances.png in/instances.obj
# [ -f out/instances.png ] && xdg-open out/instances.png&

./bin/model -o in/normals.obj normals
./bin/raytrace -r 720 -s 3 -o out/normals.png in/normals.obj
# [ -f out/normals.png ] && xdg-open out/normals.png&

./bin/model -o in/subdiv.obj subdiv
./bin/raytrace -r 720 -s 3 -o out/subdiv.png in/subdiv.obj
# [ -f out/subdiv.png ] && xdg-open out/subdiv.png&
