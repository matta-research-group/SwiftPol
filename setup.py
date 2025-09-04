from setuptools import find_packages, setup

# Package setup
setup(
    name='swiftpol',
    version='0.1.4',
    packages=find_packages(include=['swiftpol', 'swiftpol.*']),
)
