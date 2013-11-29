PySBR - Shape-based registration in Python 
==========================================


A tool for shape-based registration in fully automated SPECT imaging quantification protocols.


PySBR is a tool implemented in python that provides a fully automated approach to achieve
SPECT to SPECT shape-based matching. Specifically, the tool comprises the following parts: 
1. A fast intensity-based registration workflow using nipype and ANTS (http://stnava.github.io/ANTs/) to search for an 
affine transformation that initialises non-linear registration;
2. an image-to-atlas registration process, in which we optimise a metric based on shape 
analysis of the image to hierarchically locate a set of landmarks of interest;
3. a resampling module to map the subject’s image to the template


Data
----

We are deciding the final availability of datasets, we will post here the url soon.


References
----------

Papers making use of PySBR should cite:


Papers making use of the template released alongwith PySBR should cite:


Papers making use of the datasets released alongwith PySBR should cite:




Credits
-------

PySBR methodology and implementation: 
Oscar Esteban, Gert Wollny

Template building:
Aida Niñerola-Baizán

Simulated database:
Judith Gallego, Albert Cot

Quantification software: 
Aida Niñerola-Baizán, Berta Martí-Fuster

Clinical application and supervision:
Francisco Lomeña (flomena@clinic.ub.es)

Senior researchers:
Francisco Lomeña (flomena@clinic.ub.es, Xavier Setoain (setoain@clinic.ub.es),
Javier Pavía (jpavia@clinic.ub.es), Domènec Ros (dros@ub.edu),
Andrés Santos (andres@die.upm.es), M.-Jesús Ledesma-Carbayo (mledesma@die.upm.es)


Acknowledgements
----------------

This work was supported in part by Multimodal Imaging tools for Neurological Diseases
(MIND-t) project of Biomedical Research Networking center in Bioengineering, Biomaterials and
Nanomedicine (CIBER-BBN), by Spain’s Ministry of Science and Innovation through SAF2009-
08076, TEC2011-28972-C02-02, IPT-300000-2010-003 and CDTICENIT (AMIT project) and Fondo de
Investigaciones Sanitarias (PI12-00390). B. Martí-Fuster was awarded a PhD fellowship (App Form Call
07-2009) of Institute for Bioengineering of Catalonia (IBEC).


