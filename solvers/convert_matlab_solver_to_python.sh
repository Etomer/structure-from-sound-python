#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 path_to_cpp_file"
    exit 1
fi

INPUT_FILE="$1"
BASENAME=$(basename "$INPUT_FILE" .cpp)
#OUTPUT_FILE="${BASENAME}_pybind.cpp"
OUTPUT_FILE="./generated_solvers/${BASENAME}_pybind.cpp"

# Extract the function name (assumes it's on one line as: MatrixXcd name(...))
FUNC_NAME=$(awk '/MatrixXcd/ && /\(/ {
    sub(/.*MatrixXcd[ \t]+/, "")
    sub(/\(.*/, "")
    print
    exit
}' "$INPUT_FILE")

# Strip out mexFunction and generate Pybind11 module
awk '
/void mexFunction/,/}/ { skip=1; next } 
{ if (!skip) print }
' "$INPUT_FILE" > "$OUTPUT_FILE"

# Add Pybind11 includes (with complex support)
if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' '1i\
#include <pybind11/pybind11.h>\
#include <pybind11/eigen.h>\
#include <pybind11/complex.h>
' "$OUTPUT_FILE"
else
    sed -i '1i#include <pybind11/pybind11.h>\n#include <pybind11/eigen.h>\n#include <pybind11/complex.h>' "$OUTPUT_FILE"
fi

# Remove the #include "mex.h" line
sed -i '' '/#include *"mex.h"/d' "$OUTPUT_FILE"

# Append module definition
cat << EOF >> "$OUTPUT_FILE"



namespace py = pybind11;

PYBIND11_MODULE(${BASENAME}, m) {
    m.def("${FUNC_NAME}", &${FUNC_NAME}, "Compute solver result");
}
EOF

echo "Generated: $OUTPUT_FILE"
