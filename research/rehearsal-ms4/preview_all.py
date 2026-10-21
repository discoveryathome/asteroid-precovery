from pathlib import Path
import json,csv,io
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
cases=[c for c in json.loads((R/'exposure-inventory.json').read_text()) if c['identity_status']!='wrong-default-exposure']
fig,axs=plt.subplots(len(cases),2,figsize=(9,4*len(cases)))
for row,c in enumerate(cases):
 data=fits.getdata(R/c['archive_file']);h=fits.getheader(R/c['archive_file']);w=WCS(h)
 text=(R/'data/raw'/f"{c['id']}-horizons.txt").read_text()
 ep=next(csv.reader(io.StringIO(text.split('$$SOE')[1].split('$$EOE')[0].strip())))
 sky=[float(ep[3]),float(ep[4])];px=w.all_world2pix([sky],0)[0];pub=w.all_world2pix([[c['ra_deg'],c['dec_deg']]],0)[0]
 orig=np.asarray(Image.open(R/'data/raw'/c['observer_image']).convert('RGB'))
 axs[row,0].imshow(orig);axs[row,0].set_title(f"Published {c['date']} / {c['published_region']}")
 lo,hi=np.percentile(data,[10,99]);ax=axs[row,1]
 ax.imshow(data,cmap='gray_r',origin='lower',vmin=lo,vmax=hi)
 ax.plot(*pub,'o',mfc='none',mec='tab:orange',ms=12,label='Published position')
 ax.plot(*px,'+',color='tab:blue',ms=12,label='JPL prediction')
 ax.set_title(f"Retrieved {c['archive_region']} / {c['archive_plate_id']}\nHeader DATE-OBS {c['archive_date_obs']}",fontsize=10)
 for ax in axs[row]:ax.axis('off')
axs[0,1].legend(loc='upper left',fontsize=7)
fig.suptitle('2002 MS4: five archival cutouts for follow-up\nReference images: Andrew Lowe / DSS. Independent astrometry and timing review pending.',fontsize=12)
fig.tight_layout(rect=[0,0,1,.97]);fig.savefig(R/'five-exposure-preview.png',dpi=130)
print('Rendered five-exposure-preview.png')
