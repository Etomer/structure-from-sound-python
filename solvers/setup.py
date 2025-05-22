import sys
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
from pathlib import Path

# --- Parse custom input ---
if len(sys.argv) < 2 or not sys.argv[1].endswith("_pybind.cpp"):
    print("Usage: python setup.py <your_file_pybind.cpp> build_ext --inplace")
    sys.exit(1)

cpp_file = sys.argv[1]
sys.argv.pop(1)  # Remove the filename so setuptools doesn't choke

# --- Extract module name from filename ---
module_name = Path(cpp_file).stem.replace("_pybind", "")
print(module_name)
# --- Build extension ---
ext_modules = [
    Pybind11Extension(
        module_name,
        [cpp_file],
        include_dirs=["/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3"],
    ),
]

setup(
    name=module_name,
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)

