project: ISOFT
src_dir: ./src
         ./plotting
		 ./scripts
include: ~/.local/include
media_dir: ./media
page_dir: doc_pages
output_dir: ./doc
summary: Ice Shelf/Ocean Fluid- and Thermodynamics: a suite of tools and
         programs to use for simulating the evolution of ice shelves
		 interacting with the ocean.
project_github: https://github.com/cmacmackin/isoft
display: public
         protected
         private
graph: false
source: false
search: false
mathjax_config: mj-config.js
extra_filetypes: py #
extra_mods: factual_mod:https://cmacmackin.github.io/factual
            f95_lapack:http://www.netlib.org/lapack95/lug95/
			logger_mod:https://cmacmackin.github.io/flogging/
			penf:http://szaghi.github.io/PENF/index.html
			h5lt:https://support.hdfgroup.org/HDF5/doc/HL/RM_H5LT.html
			hdf5:https://support.hdfgroup.org/HDF5/doc/RM/RM_H5Front.html#LowLevelAPIs
			chebyshev_mod:https://cmacmackin.github.io/factual/module/chebyshev_mod.html
author: Chris MacMackin
author_description: I am a graduate student at the University of Oxford, studying
                    the melting and evolution of ice shelves.
website: https://cmacmackin.github.io
github: https://github.com/cmacmackin
linkedin: https://www.linkedin.com/in/christopher-macmackin-a11283173/
email: cmacmackin@gmail.com

DOI LICENSE

ISOFT is a piece of software/suite of tools which I developed while
working on my PhD thesis to simulate the evolution of ice shelves
coupled to ocean plumes. It attempts to provide an object oriented,
extensible framework with which to model glacial flow using modern
Fortran. Though developed for ice shelves, it could in principle be
modified to simulated grounded ice dynamics as well. I've published
the code and documentation to GitHub in the hopes that it might be
useful to others.

As much as practical, ISOFT was kept agnostic as to whether it was
handling a 1-D or a 2-D representation of an ice shelf/plume. However,
the current implementation does explicitly assume a 1-D system in a
number of instances and would thus need to be modified to handle 2-D
problems.  1-D simulations were sufficiently fast that they could be
run in serial. Multithreading could easily be implemented in many
parts of the code where arithmetic is performed on arrays. Indeed,
most of these cases are simple enough that a compiler may be able to
parallelise them automatically. More sophisticated approaches
involving message passing would likely be necessary to make 2-D
simulations practical, but this would be far more difficult to
implement and would likely require substantial refactoring of the
code. In particular, the
[nonlinear solvers](|pages|3-codeDesign/3-solvers.html) would likely
need to be replaced.

## Quick Start



## Documentation



## License

ISOFT is licensed under the GNU General Public License (GPL) v3.0 or
later. The terms are provided in the file `LICENSE`. The Lesser General
Public License (LGPL) would have been used, but the Chebyshev pseudo-spectral
implementation uses the FFTW3 library, which is licensed under the GPL.
