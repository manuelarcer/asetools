# as suggested in the following
# https://towardsdatascience.com/setuptools-python-571e7d5500f2#:~:text=be%20more%20appropriate.-,The%20setup.,as%20the%20command%20line%20interface.

from setuptools import setup
if __name__ == '__main__':
    setup(
        name='asetools',
        version='0.1.0',
        description='Complement to ASE to perform analysis and manipulate atoms',
        author='Juan M Arce-Ramos',
        author_email='manuel.arceram@gmail.com',
        url='https://github.com/manuelarcer/asetools.git',
        packages=find_packages(),
        install_requires=[
        # List of packages your project depends on
        entry_points={
            'console_scripts': [
                'executable_script1=asetools.bin.executable_script1:main',
                'executable_script2=asetools.bin.executable_script2:main',
                # etc.
            ],
        }
        ...
    )
