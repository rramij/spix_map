# spix_map
Make spectral index and curvature maps from two or more FITS images \
Run on Conda 2 or 3; \
required: numpy version > 1.12.0 \
required: astropy \

Instructions:
1. Edit the 'spix_map_NEW.py' file and provide the inputs. make sure to keep the image and other informations in the same order in the List.

2. To properly align the image pixels, choose a point source (in CARTA or DS9) and note the central pix position in the 'align.dat' file. Again keep the information order same.

3. Finally run: python spix_map_NEW.py

4. That's all!!
