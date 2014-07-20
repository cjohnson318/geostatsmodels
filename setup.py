from distutils.core import setup

setup(name='geostatsmodels',
      version='0.1',
      description='Geostatistics models',
      author='Connor Johnson',
      author_email='connor.labaume.johnson@gmail.com',
      url='https://github.com/cjohnson318/geostatsmodels',
      packages=['geostatsmodels'],
      scripts = ['kriging_example.py','simulation_example.py'],
)
