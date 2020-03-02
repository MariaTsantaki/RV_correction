RV correction
======

Correct for radial velocity shifts. 

It is run:

```python
python rv_correction.py fname.fits
```
The fits file has to be in a specific format. 

1) a text file (.dat or .txt) with 2 columns for the wavelength and flux.

2) a fits file (.fits) from the ESO archive (ESO SDP 1D spectrum) and also 1D files which include the key words 'CRVAL1' and 'CDELT1' in their header.

3) a table fits format (.spec). This is the format [FASMA](https://github.com/MariaTsantaki/FASMA-synthesis) saves the output synthetic spectra if requested.  


AUTHOR
-------

    M. Tsantaki

LICENCE
-------

It uses the MIT licence.
Copyright © 2020 Maria Tsantaki

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
