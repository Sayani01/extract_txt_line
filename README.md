# extract_txt_line
A project to extract text lines from a document with english/bengali handwriting. This repository includes C++ implementation of the following feature extraction algorithms:

-- deriving connected components
-- deriving rectangular bounding box in the directions 0deg, +/-45deg, 90deg
-- thinning algorithms - FPTA, GHPTA, KWK
-- Gaussian blur with rotation
-- binarization algorithms - NICK, Sauvola, Niblack, Otsu
-- projection profiles in 0deg, +/-45deg, 90deg directions

The project is to identify textlines from a document written in multiple directions in English or Bengali. The codes have been tested on ICDAR dataset as well as some randomly collected handwriting samples and a 98% textline detection accuracy have been achieved.
