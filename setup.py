#!/usr/bin/env python3
import setuptools
import site
from numpy.distutils.core import setup, Extension

site.ENABLE_USER_SITE = True

setup(
    ext_modules=[
        Extension(name="f2pydemo.pyprod", sources=["src/f2pydemo/prod.f90"]),
        Extension(name="f2pydemo.badprec", sources=["src/f2pydemo/badprec.f90"]),
        Extension(name="f2pydemo.vecmed", sources=["src/f2pydemo/vecmedian.f90"])
    ]
)
