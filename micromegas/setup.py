from setuptools import find_packages, setup

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
