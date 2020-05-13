from setuptools import setup, find_packages

setup(
    name='mrpod',
    version='0.0.1',
    description=('Multi-resolution proper orthogonal decomposition.'),
    long_description='',
    keywords='principle component analysis, wavelet transform',
    author=('Zhiyao Yin, German Aerospace Center (DLR)'),
    #author_email='',
    url='https://mrpod.readthedocs.io',
    maintainer=('Zhiyao Yin'),
    maintainer_email='zhiyao.yin@dlr.de',
    license='BSD 3-Clause',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
        ],
    packages=['mrpod'],
    package_dir={'mrpod': 'mrpod'},
    install_requires=['numpy']
    )