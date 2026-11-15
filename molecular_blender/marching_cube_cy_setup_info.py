import platform
import subprocess

def getExtensionArgs(Utils):
    if platform.system() == "Darwin":
        try:
            brew_prefix = subprocess.check_output(
                ["brew", "--prefix", "libomp"], stderr=subprocess.DEVNULL
            ).decode().strip()
            omp_include = f"{brew_prefix}/include"
            omp_lib = f"{brew_prefix}/lib"
        except Exception:
            omp_include = "/opt/local/include/libomp"
            omp_lib = "/opt/local/lib/libomp"
        return {
            "extra_compile_args": ["-Xpreprocessor", "-fopenmp", f"-I{omp_include}"],
            "extra_link_args": [f"-L{omp_lib}", "-lomp"],
        }
    else:
        return {
            "extra_compile_args": ["-fopenmp"],
            "extra_link_args": ["-fopenmp"],
        }
