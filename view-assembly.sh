#!/bin/bash

find ./build/src -name "fields.dot" && dot -Tpng ./build/src/fields.dot -o ./build/src/fields.png && display ./build/src/fields.png

find ./build/src -name "graph_residual.dot" && dot -Tpng ./build/src/graph_residual.dot -o ./build/src/graph_residual.png && display ./build/src/graph_residual.png

find ./build/src -name "graph_jacobian.dot" && dot -Tpng ./build/src/graph_jacobian.dot -o ./build/src/graph_jacobian.png && display ./build/src/graph_jacobian.png

find ./build/src -name "error.dot" && dot -Tpng ./build/src/error.dot -o ./build/src/assembly-error.png && display ./build/src/assembly-error.png

