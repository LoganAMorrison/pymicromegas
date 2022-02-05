import os
import subprocess
from setuptools import find_packages, setup
from pathlib import Path

SETUP_DIR = Path(__file__).parent.absolute()
PYMICROMEGAS_DIR = SETUP_DIR.joinpath("pymicromegas")
MICROMEGAS_DIR = SETUP_DIR.joinpath("micromegas")
MSSM_DIR = MICROMEGAS_DIR.joinpath("MSSM")


def get_extension_suffix() -> str:
    proc = subprocess.Popen(
        ["python3-config", "--extension-suffix"], stdout=subprocess.PIPE
    )
    return proc.stdout.readline().decode("ascii").strip("\n")


def move_extension(ext_name: str) -> None:
    os.rename(
        MSSM_DIR.joinpath(ext_name),
        PYMICROMEGAS_DIR.joinpath(ext_name),
    )


def build_micromegas(softsusy: bool = False, suspect: bool = True):
    """
    Build the micrOMEGAs library, the pymicromegas extension, the
    micromegas-api extension, the spheno extension and optionally the softsusy
    and suspect extensions.
    """
    # Build micrOMEGAs and CalcHEP
    subprocess.check_call(["make", "-C", "micromegas"])
    # Build MSSM
    subprocess.check_call(["make", "-C", "micromegas/MSSM"])
    # Build pymicromegas, micromegas-api and spheno
    subprocess.check_call(
        ["make", "-C", "micromegas/MSSM", "--file=Makefile.pymicromegas"]
    )
    if suspect:
        # Build suspect extension
        subprocess.check_call(
            ["make", "-C", "micromegas/MSSM", "--file=Makefile.suspect"]
        )
    if softsusy:
        # Build softsusy extension
        subprocess.check_call(
            ["make", "-C", "micromegas/MSSM", "--file=Makefile.suspect"]
        )
    # Move .so's to package directory
    ext_suffix = get_extension_suffix()
    micromegas_ext_name = "micromegas" + ext_suffix
    pymicromegas_ext_name = "pymicromegas" + ext_suffix
    spheno_ext_name = "spheno" + ext_suffix

    move_extension(micromegas_ext_name)
    move_extension(pymicromegas_ext_name)
    move_extension(spheno_ext_name)

    if suspect:
        suspect_ext_name = "suspect" + ext_suffix
        move_extension(suspect_ext_name)
    if softsusy:
        softsusy_ext_name = "softsusy" + ext_suffix
        move_extension(softsusy_ext_name)


if __name__ == "__main__":
    build_micromegas()
    setup(
        name="pymicromegas",
        version="1.0.0",
        author="Logan Morrison",
        author_email="loanmorr@ucsc.edu",
        maintainer="Logan Morrison",
        maintainer_email="loanmorr@ucsc.edu",
        description="Python interface to micrOMEGAs",
        packages=find_packages(),
        package_data={"pymicromegas": ["*.so"]},
        install_requires=["pip"],
        python_requires=">=3",
        zip_safe=False,
        include_package_data=True,
        license="gpl-3.0",
        platforms="MacOS and Linux",
        classifiers=["Programming Language :: Python"],
    )
