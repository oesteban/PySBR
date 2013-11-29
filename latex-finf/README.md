
Frontiers in Neuroinformatics - Call for Research Topics (Python in neuroscience II)
====================================================================================


About the call
--------------

Please visit this information here: http://www.frontiersin.org/Neuroinformatics/researchtopics/Python_in_Neuroscience_II/1591 



Author Guidelines
-----------------

  * Visit them at http://www.frontiersin.org/Neuroinformatics/authorguidelines
  * I've downloaded a latex template ("FRONTIERS LateX.zip")



Submitted outline
-----------------

Single-photon emission computed tomography (SPECT) imaging of the dopaminergic system using
specific agents (e.g. 123 I- or 99mTc- based) is of great value to discriminate parkinsonian syndromes
from other movement disorders by quantifying the dopamine transporter (DAT) binding. Quantification
of DAT-scan imaging can facilitate early diagnosis, following up progression and helping assess the
effects of treatment strategies. In clinical routine, this application is generally performed by visual
assessment. Nonetheless, some efforts have been done to automate the process and standardise the
outcome by means of image registration of the subject's study to a template on which a number of
quantification regions are defined. This registration is challenging due to the fact that an structural
reference (e.g. T1 weighted MRI) is generally not available for the subject under study.

DAT-scan images of normal subjects are characterised by presenting a 'comma-shaped' pattern due to
the significant uptake located in the caudate nucleus and putamen. However, in subjects affected by
parkinsonian syndromes, the neuronal degeneration is more pronounced in the putamen and the
appearance becomes 'point-like'. Quantification protocols define different regions of interest along the
striatum structures to measure the uptake and characterise the pattern. Without structural reference, the
application of intensity-based non-linear registration approaches is cumbersome: the registration is not
accurate enough under very strong regularisation conditions, and more free conditions hinder the
quantification turning 'comma-shaped' patterns into 'point-like' and viceversa depending on the choice
of the template.

In this work we present PySBR (Shape Based Registration -SBR- in python), a tool implemented in
python that provides a fully automated approach to achieve this alignment by using shape based image
registration. Specifically, the tool comprises the following parts: 1. A fast intensity-based registration
process to search for an affine transformation that initialises non-linear registration; 2. an
image-to-atlas registration process, in which we optimise a metric based on shape analysis of the image
to hierarchically locate a set of landmarks of interest; 3. a resampling module to map the subject's
image to the template.

In order to evaluate PySBR applicability to the problem at hand, the software was used along with
nipype for testing the non-linear alignment of ~30 simulated data sets and 10 real data sets to a certain
template. Additionally, complementary results were obtained by means of intensity-based registration,
at different levels or regularisation. For the simulated data sets, the alignment obtained on the
DAT-scan was applied to the reference T1 weighted MRI to evaluate the overlap found between
specific ROIs defined in the template and those obtained from the subject. For the real data sets,
cross-comparison with widely used intensity-based registration techniques was provided. The source
code of PySBR, as well as the evaluation framework coded in nipype are made publicly available.


