HDR-VDP-2: A calibrated visual metric for visibility and quality
predictions in all luminance conditions

This directory contains matlab code of the HDR-VDP-2 - a visual
difference predictor for high dynamic range images. This is the
successor of the original HDR-VDP.

Always check for the latest release of the metric at:

http://hdrvdp.sourceforge.net/

The current version and the list of changes can be found in the
ChangeLog.txt.

-----------------------------------------------------------------
To install the metric:

1. Download matlabPyrTools from http://www.cns.nyu.edu/~lcv/software.html
2. Add matlabPyrTools to the matlab path

   It is important to compile and installed MEX files for
   matlabPyrTools (check README.txt and the MEX directory). The metric
   may return different results if matlabPyrTools is used without
   compiled MEX files.

   If you are using matlabPyrTools 1.4, you will need to install a
   patch. Refer to the project web-site
   (http://hdrvdp.sourceforge.net/) for more information.

3. Add the 'hdrvdp-VERSION' directory to the matlab path

-----------------------------------------------------------------
To run the metric:

Check Contents.m and the documentation for hdrvdp.m,
hdrvdp_visualize.m and hdrvdp_pix_per_deg.m


-----------------------------------------------------------------
Citiations:

If you find this metric useful, please cite the paper:

HDR-VDP-2: A calibrated visual metric for visibility and quality predictions in all luminance conditions
Rafał Mantiuk, Kil Joong Kim, Allan G. Rempel and Wolfgang Heidrich.
In: ACM Transactions on Graphics (Proc. of SIGGRAPH'11), 30(4), article no. 40, 2011

AND the current version of the metric, for example HDR-VDP 2.0. Check
ChangeLog.txt for the current version of the metric.

-----------------------------------------------------------------
Contact:

http://groups.google.com/group/hdrvdp
Rafal Mantiuk <mantiuk@gmail.com>

-----------------------------------------------------------------
For more more information, refer to the project web-site:
http://hdrvdp.sourceforge.net/
