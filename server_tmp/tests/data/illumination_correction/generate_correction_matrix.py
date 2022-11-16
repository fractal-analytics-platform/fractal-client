import numpy as np
from PIL import Image

im1 = np.ones(2160 * 2560, dtype=np.uint16).reshape(2160, 2560)
Image.fromarray(im1).save("illum_corr_matrix.png", optimize=True)
