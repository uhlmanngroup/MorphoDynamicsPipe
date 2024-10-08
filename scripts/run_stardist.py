import skimage
#from datetime import datetime

#30 seconds
#start = datetime.now(); print(start)
from stardist.models import StarDist2D
#print(datetime.now() - start)

model = StarDist2D.from_pretrained('2D_versatile_fluo')

import sys
#this_input = sys.argv[1]
#this_output = sys.argv[2]

this_input = list(snakemake.input)
this_output = list(snakemake.output)

im = skimage.io.imread(this_input[0])
labels, _ = model.predict_instances(skimage.exposure.equalize_adapthist(im))
#    print(output_name)
#    np.save(output_name, labels)
skimage.io.imsave(this_output[0], labels, check_contrast=False)
