#!/usr/bin/env python

import os
import sys

from distutils.core import setup
from distutils.core import Extension

unroll = 0
unroll_form_error = "setup.py: error: option --unroll-order must have form --unroll-order=INT."
unroll_value_error = "setup.py: error: INT for --unroll-order must be positive."
other_args = list()
for opt in sys.argv:
    if opt.startswith("--unroll-order="):
        unroll = opt.split("--unroll-order=")
        if len(unroll) != 2:
            print(unroll_form_error)
            exit(1)
        try:
            unroll = int(unroll[1])
        except Exception:
            print(unroll_form_error)
            exit(1)
        if unroll < 0:
            print(unroll_value_error)
            exit(1)
    else:
        other_args.append(opt)

sys.argv = other_args

if unroll:

    try:
        source = open(os.path.join("source", "pairwise.c"), "rb")
    except Exception:
        print("setup.py: error: failed to open pairwise.c for reading.")
        exit(1)

    try:
        sink = open(os.path.join("source", "pairwise-unrolled.c"), "wb")
    except Exception:
        print("setup.py: error: failed to open pairwise-unrolled.c for writing.")
        exit(1)

    print("setup.py: notice: manually unrolling time-intensive loop in pairwise.c to order %d." % unroll)

    code = source.read().split("\n")
    source.close()

    expanded = list()

    for line in code:

        stripped = line.lstrip()
        
        if not stripped.startswith("//###UNROLLING_"):
            expanded.append(line)
                
        elif unroll:
        
            if stripped.startswith("//###UNROLLING_DECLARATIONS"):
                indent = line.split("//###UNROLLING_DECLARATIONS")[0]
                
                for i in xrange(unroll):
                    expanded.append(indent + "register size_t b%d;" % i)
                    expanded.append(indent)
                    expanded.append(indent + "register double dx%d;" % i)
                    expanded.append(indent + "register double dy%d;" % i)
                    expanded.append(indent + "register double dz%d;" % i)
                    if i != (unroll - 1):
                        expanded.append(indent)
                
            if stripped.startswith("//###UNROLLING_ITERATIONS"):
                indent = line.split("//###UNROLLING_ITERATIONS")[0]
                
                expanded.append(indent + "for(; j + %d < n; j += %d) {" % (unroll - 1, unroll))
                expanded.append(indent)
                
                indent2 = indent + "    "
                
                for i in xrange(unroll):
                    expanded.append(indent2 + "b%d = (j + %d) * 3;" % (i, i))
                    expanded.append(indent2)
                    expanded.append(indent2 + "dx%d = ax - *(p + b%d + 0);" % (i, i))
                    expanded.append(indent2 + "dy%d = ay - *(p + b%d + 1);" % (i, i))
                    expanded.append(indent2 + "dz%d = az - *(p + b%d + 2);" % (i, i))
                    expanded.append(indent2)
                    expanded.append(indent2 + "*(d + %d) = sqrt((dx%d * dx%d) + (dy%d * dy%d) + (dz%d * dz%d));" \
                                    % (i, i, i, i, i, i, i))
                    expanded.append(indent2)
                
                expanded.append(indent2 + "d += %d;" % unroll)
                expanded.append(indent2)                
                
                expanded.append(indent + "}")

    expanded = "\n".join(expanded)

    try:
        sink.write(expanded)
        sink.flush()
        os.fsync(sink)
        sink.close()
    except Exception:
        print("setup.py: error: failed to write modified source to pairwise-unrolled.c.")
        exit(1)

    print("setup.py: notice: wrote modified source to pairwise-unrolled.c.")
    print("setup.py: notice: passing control to distutils.")

os.environ["CC"] = "gcc"

gcc_args = list()
gcc_args.append("-Ofast")

pw = Extension("pairwise",
               [os.path.join("source", "pairwise%s.c" % ("-unrolled" if unroll else ""))],
               extra_compile_args = gcc_args)

NAME = "pairwise"
VERSION = "1.0"
DESCRIPTION = "Fast multi-threaded calculation of Euclidean pairwise distances."
AUTHOR = "Adam J. Marsden"
EXT_MODULES = [pw]

setup(name = NAME,
      version = VERSION,
      description = DESCRIPTION,
      author = AUTHOR,
      ext_modules = EXT_MODULES)
