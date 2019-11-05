The LP models used in our comparisons were obtained from the following sources:

* J. Castro [http://www-eio.upc.es/~jcastro/]
  - huge CTA instances
  - integrated refinery problems
  - L1.zip
  - Linf.zip

* J. A. J. Hall [http://www.maths.ed.ac.uk/hall/PublicLP/]

* Kennington collection [http://www.netlib.org/lp/data/kennington/index.html]

* C. Meszaros [http://old.sztaki.hu/~meszaros/public_ftp/lptestset/]
  - misc
  - New
  - problematic
  - stochlp

* MIPLIB 2010 [http://miplib.zib.de/download/miplib2010-complete.tgz]

* H. Mittelmann [http://plato.asu.edu/ftp/lptestset/]
  - fome
  - misc
  - nug
  - pds
  - rail

* Netlib collection [http://www.netlib.org/lp/data/index.html]

Each "testset" is a collection of LP models. Each model has a name, a group
(which is the source of the model) and a 4-digit ID. The dimensions given in the
files <testset>.tbl refer to the constraint matrix after presolving the model
using the Gurobi LP presolver (without running the optimizer).
