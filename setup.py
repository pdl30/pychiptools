import os
from setuptools import setup, find_packages

setup(name='pychiptools',
      version='0.0.1',
      packages=find_packages(),
      description='pychiptools is a Python module to analyze ChIP-seq NGS data',
      author='Patrick Lombard',
      author_email='ptk.lmb55@gmail.com',
     # packages=['pychiptools'],
      package_data={"pychiptools":['data/*']},
      scripts=['scripts/bdg2bw', 'scripts/pychip_align.py', 'scripts/pychip_diff_bind.py', 'scripts/pychip_motifs.py', 'scripts/pychip_peak_anno.py', 
        'scripts/pychip_peak_call.py', 'scripts/pychip_viz.py'],
      install_requires=['pysam', 'pybedtools'],
      license='GPLv3',
      platforms='any',
      classifiers=[
         'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
         'Development Status :: 3 - Alpha',
         'Programming Language :: Python :: 2.7',
         'Environment :: Console',
      ],
      long_description="""

pychiptools is a Python module to analyze ChIP-seq NGS data

 Contact
=============

If you have any questions or comments about pychiptools, please feel free to contact me via
eMail: ptk.lmb55@gmail.com

""",
    )
