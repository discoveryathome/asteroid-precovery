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
cases=json.loads((R/'cases.json').read_text())['cases']
fig,axs=plt.subplots(2,3,figsize=(13,8))
preds=[]
for row,(case,region) in enumerate(zip(cases,['XE740','XE741'])):
 path=R/'data/raw'/f'{region}.fits';data=fits.getdata(path);h=fits.getheader(path);w=WCS(h)
 text=(R/'data/raw'/f"{case['id']}-horizons.txt").read_text()
 ep=next(csv.reader(io.StringIO(text.split('$$SOE')[1].split('$$EOE')[0].strip())))
 sky=[float(ep[3]),float(ep[4])];px=w.all_world2pix([sky],0)[0];pub=w.all_world2pix([[case['ra_deg'],case['dec_deg']]],0)[0]
 residual=np.hypot((sky[0]-case['ra_deg'])*np.cos(np.deg2rad(sky[1])),sky[1]-case['dec_deg'])*3600
 preds.append({'case':case['id'],'prediction_time_horizons_UT1':ep[0].strip(),'predicted_radec_deg':sky,'published_position_radec_deg':[case['ra_deg'],case['dec_deg']],'published_minus_prediction_distance_arcsec':float(residual),'not_an_independent_image_measurement':True,'sky_motion_arcsec_hour':[float(ep[5]),float(ep[6])],'expected_motion_during_exposure_arcsec':float(np.hypot(float(ep[5]),float(ep[6]))*case['exposure_minutes']/60),'ellipse_3sigma_arcsec':[float(ep[11]),float(ep[12])],'ellipse_theta_deg':float(ep[13]),'pixel_xy':px.tolist()})
 orig=np.asarray(Image.open(R/'data/raw'/f"observer-{case['date'].replace('-','')}.gif").convert('RGB'))
 axs[row,0].imshow(orig);axs[row,0].set_title(f"Investigator's published image ({case['id']})")
 lo,hi=np.percentile(data,[10,99])
 for col in [1,2]:
  ax=axs[row,col];ax.imshow(data,cmap='gray_r',origin='lower',vmin=lo,vmax=hi)
  ax.plot(*pub,'o',mfc='none',mec='tab:orange',ms=12,label='Published position')
  ax.plot(*px,'+',color='tab:blue',ms=12,label='JPL prediction')
  ax.set_title(region+(' native cutout' if col==1 else ' central 68 arcsec'))
  if col==2:ax.set_xlim(px[0]-20,px[0]+20);ax.set_ylim(px[1]-20,px[1]+20)
 for ax in axs[row]:ax.axis('off')
axs[0,2].legend(loc='upper left',fontsize=7)
fig.suptitle('2002 MS4 archive leads — published position / current-orbit prediction; image measurement pending')
fig.tight_layout();fig.savefig(R/'preview.png',dpi=140)
(R/'prediction-audit.json').write_text(json.dumps(preds,indent=2)+'\n')
print(json.dumps(preds,indent=2))
