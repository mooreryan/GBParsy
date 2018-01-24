from distutils.core import setup, Extension 
 
setup(name = "gbfpy", 
        version = "0.1", 
        description="gbff parser extension module", 
        author = "Tae-ho Lee", 
        ext_modules=[Extension("gbfpy", ["gbfpy.c", "gbfp.c"])] 
) 
