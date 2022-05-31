import setuptools
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


def get_requirements(path):
    with open(path, "r") as fh:
        content = fh.read()
    return [
        req
        for req in content.split("\n")
        if req != '' and not req.startswith('#')
    ]


install_requires = get_requirements('requirements.txt')

setuptools.setup(
    name="Quagga",
    version="0.0.3",
    author="Fan Feng, Sean Patrick Moran",
    author_email="fanfeng@umich.edu",
    description="A user-friendly package for calling stripes from contact maps",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dmcbffeng/StripeCaller",
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    include_package_data=True,
    scripts=['Quagga/bin/stripe_caller'],  # call from command line
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
