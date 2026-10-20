"""Offline reduction of frozen 1999 AN10 precovery images, with explicit review gates."""
import csv
import io
import json
import math
import os
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from scipy.special import erf
from scipy.ndimage import median_filter
from .acquire import ROOT
from .verify import screen, strict_time, check_manifest


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def ephemeris():
    text = (ROOT/'data/raw/horizons-an10-1955.txt').read_text()
    if 'JPL#100' not in text or 'Palomar Mountain' not in text: raise ValueError('Unexpected orbit or site')
    rows = []
    for row in csv.reader(io.StringIO(text.split('$$SOE')[1].split('$$EOE')[0].strip())):
        t=datetime.strptime(row[0].strip(),'%Y-%b-%d %H:%M:%S')
        rows.append([(t-datetime(1970,1,1)).total_seconds(), *[float(row[i]) for i in [3,4,11,12,13]]])
    return np.array(rows)


def prediction(t, eph):
    stamp = (t-datetime(1970,1,1)).total_seconds()
    if not eph[0,0] <= stamp <= eph[-1,0]: raise ValueError('Prediction outside frozen ephemeris')
    a = [float(np.interp(stamp,eph[:,0],eph[:,i])) for i in range(1,6)]
    if min(a[2:4]) <= 0: raise ValueError('Missing or invalid ephemeris uncertainty')
    theta = np.deg2rad(a[4]);u=np.array([np.cos(theta),np.sin(theta)])
    v=np.array([-np.sin(theta),np.cos(theta)])
    cov=(a[2]/3)**2*np.outer(u,u)+(a[3]/3)**2*np.outer(v,v)
    return dict(ra_deg=a[0],dec_deg=a[1],smaa_3sigma_arcsec=a[2],smia_3sigma_arcsec=a[3],
                theta_deg=a[4],covariance_east_north_arcsec2=cov.tolist(),
                time=t.isoformat(),timescale='UT1',frame='ICRF astrometric topocentric',
                orbit='JPL#100',center='675',interpolation='linear at 1-minute spacing')


def offsets(a,b):
    """Local east/north offsets, wrap RA and include cos(dec)."""
    return np.array([((a[0]-b[0]+180)%360-180)*np.cos(np.deg2rad(b[1])),a[1]-b[1]])*3600


def catalog():
    rows=list(csv.DictReader((ROOT/'data/raw/gaia-dr3.csv').open()))
    good=[]
    dt=1955.070-2016.0
    for r in rows:
        if not r['ruwe'] or not 12.5<float(r['phot_g_mean_mag'])<18.5 or float(r['ruwe'])>1.4: continue
        if r['parallax'] and float(r['parallax'])>5: continue
        if np.hypot(float(r['pmra']),float(r['pmdec']))>200: continue
        dec=float(r['dec']);ra=float(r['ra'])
        ra+=float(r['pmra'])*dt/(3.6e6*np.cos(np.deg2rad(dec)))
        dec+=float(r['pmdec'])*dt/3.6e6
        good.append((r['source_id'],ra,dec,float(r['phot_g_mean_mag'])))
    return good


def calibrate(data,w,cat,scale):
    stars=[]
    noise=1.4826*np.median(np.abs(data-np.median(data)))
    for sid,ra,dec,mag in cat:
        x,y=w.all_world2pix([[ra,dec]],0)[0]
        ix,iy=int(round(x)),int(round(y));rad=6
        if min(ix,iy)<rad or ix+rad>=data.shape[1] or iy+rad>=data.shape[0]: continue
        patch=data[iy-rad:iy+rad+1,ix-rad:ix+rad+1].astype(float)
        yy,xx=np.mgrid[iy-rad:iy+rad+1,ix-rad:ix+rad+1]
        bg=np.median(patch);peak=patch.max()-bg
        if peak<6*noise: continue
        def model(p):
            cx,cy,sig,amp,b=p
            return b+amp*np.exp(-((xx-cx)**2+(yy-cy)**2)/(2*sig**2))
        fit=least_squares(lambda p:((model(p)-patch)/noise).ravel(),[x,y,2.,peak,bg],
                          bounds=([x-4,y-4,.6,0,bg-3*noise],[x+4,y+4,5,peak*3,bg+3*noise]),loss='soft_l1')
        cx,cy,sig,amp,b=fit.x
        if sig>4.5 or np.hypot(cx-x,cy-y)>3.5 or amp/noise<6: continue
        stars.append(dict(source_id=sid,gaia_epoch_ra_deg=ra,gaia_epoch_dec_deg=dec,
                          x=float(cx),y=float(cy),expected_x=float(x),expected_y=float(y),sigma_pixel=float(sig)))
    if len(stars)<8: raise ValueError(f'Only {len(stars)} calibration stars')
    pts=np.array([[s['x'],s['y']] for s in stars]);target=np.array([[s['expected_x'],s['expected_y']] for s in stars])
    design=np.column_stack([np.ones(len(pts)),pts]);keep=np.ones(len(pts),dtype=bool)
    for _ in range(4):
        coeff=np.linalg.lstsq(design[keep],target[keep],rcond=None)[0]
        err=np.linalg.norm(design@coeff-target,axis=1)
        cutoff=max(.6,3*np.median(err[keep]));keep=err<cutoff
    if keep.sum()<8: raise ValueError('Insufficient inlier stars')
    coeff=np.linalg.lstsq(design[keep],target[keep],rcond=None)[0]
    err=np.linalg.norm(design@coeff-target,axis=1)
    # Independent held-out prediction errors: leave-one-out of accepted stars.
    loo=[]
    for i in np.where(keep)[0]:
        train=keep.copy();train[i]=False
        c=np.linalg.lstsq(design[train],target[train],rcond=None)[0]
        loo.append((design[i]@c-target[i])*scale)
    rms=float(np.sqrt(np.mean(np.array(loo)**2)))
    for s,k,e in zip(stars,keep,err):s.update(used=bool(k),fit_residual_arcsec=float(e*scale))
    return coeff,dict(catalogue='Gaia DR3 gaiadr3.gaia_source',epoch=2016.0,
                      propagated_epoch=1955.070,propagation='linear proper motion; parallax <5 mas, PM <200 mas/year; perspective acceleration omitted',
                      model='six-parameter affine correction of native DSS WCS',
                      original_plate_solution='GSC1 polynomial coefficients preserved in FITS',
                      inliers=int(keep.sum()),loo_rms_per_axis_arcsec=rms,
                      affine_observed_to_catalogue_pixel=coeff.tolist(),stars=stars)


def measured_sky(xy,w,coeff):
    corrected=np.r_[1.,xy]@coeff
    return w.all_pix2world([corrected],0)[0]


def predicted_pixel(sky,w,coeff):
    pix=w.all_world2pix([sky],0)[0]
    return np.linalg.solve(coeff[1:].T,pix-coeff[0])


def fit_trail(data, center, vector, mask_stars=()):
    length=np.linalg.norm(vector);theta=np.arctan2(vector[1],vector[0])
    # Use an elongated native-pixel window enclosing the trail and blank margins.
    halfx=int(abs(vector[0])/2+20);halfy=int(abs(vector[1])/2+20)
    cx,cy=center; x0=max(0,int(cx)-halfx);x1=min(data.shape[1],int(cx)+halfx+1)
    y0=max(0,int(cy)-halfy);y1=min(data.shape[0],int(cy)+halfy+1)
    patch=data[y0:y1,x0:x1].astype(float);yy,xx=np.mgrid[y0:y1,x0:x1]
    bg=float(np.median(patch));noise=float(1.4826*np.median(np.abs(patch-bg)))
    valid=np.isfinite(patch)
    for star in mask_stars:
        valid &= (xx-star['x'])**2+(yy-star['y'])**2>max(5,3*star['sigma_pixel'])**2
    if valid.sum()<200: raise ValueError('Too few unmasked pixels')
    def model(p):
        x,y,L,ang,sig,amp,b=p
        along=(xx-x)*np.cos(ang)+(yy-y)*np.sin(ang)
        cross=-(xx-x)*np.sin(ang)+(yy-y)*np.cos(ang)
        profile=.5*(erf((along+L/2)/(np.sqrt(2)*sig))-erf((along-L/2)/(np.sqrt(2)*sig)))
        return b+amp*profile*np.exp(-cross**2/(2*sig**2))
    initial=[cx,cy,length,theta,1.5,max(noise, np.percentile(patch,90)-bg),bg]
    lo=[cx-12,cy-12,length*.45,theta-.3,.5,0,bg-3*noise]
    hi=[cx+12,cy+12,length*1.55,theta+.3,5.,max(patch.max()-bg,noise)*3,bg+3*noise]
    fit=least_squares(lambda p:((model(p)-patch)/noise)[valid],initial,bounds=(lo,hi),loss='soft_l1',max_nfev=400)
    p=fit.x;fitted=model(p)
    def point(p):
        x,y,sig,amp,b=p
        return b+amp*np.exp(-((xx-x)**2+(yy-y)**2)/(2*sig**2))
    pf=least_squares(lambda p:((point(p)-patch)/noise)[valid],
                     [cx,cy,2.,max(noise,patch.max()-bg),bg],
                     bounds=([cx-12,cy-12,.5,0,bg-3*noise],[cx+12,cy+12,5.,(patch.max()-bg)*3,bg+3*noise]),loss='soft_l1')
    # Robust loss-based BIC-like discriminator, NOT a calibrated likelihood ratio.
    def deviance(a): return float(np.sum(2*(np.sqrt(1+((a-patch)[valid]/noise)**2)-1)))
    delta=deviance(point(pf.x))-deviance(fitted)-2*np.log(valid.sum())
    template=(fitted-p[6])/max(p[5],1e-12)
    snr=float(p[5]*np.sqrt(np.sum(template[valid]**2))/noise/2)
    # Factor two noise inflation for correlated/compressed photographic pixels.
    # Empirical offset controls are reported separately, no Gaussian p-value claimed.
    cov=np.linalg.pinv(fit.jac.T@fit.jac)*4
    return dict(x=float(p[0]),y=float(p[1]),length_pixel=float(p[2]),angle_rad=float(p[3]),
                sigma_pixel=float(p[4]),amplitude=float(p[5]),background=float(p[6]),noise=noise,
                snr=snr,delta_bic_vs_point=float(delta),length_ratio=float(p[2]/length),
                angle_error_deg=float(abs(np.rad2deg(p[3]-theta))),
                centroid_covariance_pixel2=cov[:2,:2].tolist(),
                fit_converged=bool(fit.success),unmasked_pixels=int(valid.sum()),
                pixel_bounds=[x0,x1,y0,y1]),(patch,fitted,valid)


def plot_evidence(output, entries):
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.cache/matplotlib'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(len(entries),3,figsize=(13,3.5*len(entries)),squeeze=False)
    for axesrow,(m,arrays) in zip(axes,entries):
        image,model,mask=arrays
        lo,hi=np.percentile(image,[5,97])
        for ax,img,label in zip(axesrow,[image,model,image-model],['Observed native pixels','Fitted trail + background','Residual']):
            ax.imshow(img,origin='lower',cmap='gray_r',vmin=lo if label!='Residual' else -3*m['noise'],vmax=hi if label!='Residual' else 3*m['noise'])
            ax.set_title(f"{m['id']}\n{label}",fontsize=10);ax.set_xlabel('Cutout x (pixel)');ax.set_ylabel('Cutout y (pixel)')
        x0,x1,y0,y1=m['pixel_bounds'];axesrow[0].plot(m['x']-x0,m['y']-y0,'r+',ms=10)
    fig.tight_layout();fig.savefig(output/'evidence.png',dpi=160);plt.close(fig)


def run(output=None, verify=True):
    output=Path(output or ROOT/'data/derived');output.mkdir(parents=True,exist_ok=True)
    if verify: check_manifest(ROOT,json.loads((ROOT/'data/manifest.json').read_text()))
    controls=json.loads((ROOT/'data/controls.json').read_text());eph=ephemeris();cat=catalog()
    entries=[];measurements=[];predictions=[];calibrations={};contexts=[]
    for control in controls['positive_controls']:
        data,header=fits.getdata(ROOT/control['file'],header=True);w=WCS(header)
        scale=abs(header['XPIXELSZ']*header['PLTSCALE']/1000)
        coeff,cal=calibrate(data,w,cat,scale);calibrations[control['region']]=cal
        start,end=map(strict_time,[control['start'],control['end']]);mid=start+(end-start)/2
        pred=[prediction(t,eph) for t in [start,mid,end]]
        for p in pred:p['control']=control['id']
        predictions.extend(pred)
        pix=np.array([predicted_pixel([p['ra_deg'],p['dec_deg']],w,coeff) for p in pred])
        m,arrays=fit_trail(data,pix[1],pix[2]-pix[0],cal['stars'])
        sky=measured_sky([m['x'],m['y']],w,coeff);res=offsets(sky,[pred[1]['ra_deg'],pred[1]['dec_deg']])
        # Position covariance: centroid fit + stellar calibration + timing floor.
        local=np.column_stack([offsets(measured_sky([m['x']+1,m['y']],w,coeff),sky),
                               offsets(measured_sky([m['x'],m['y']+1],w,coeff),sky)])
        rate=offsets([pred[2]['ra_deg'],pred[2]['dec_deg']],[pred[0]['ra_deg'],pred[0]['dec_deg']])/(end-start).total_seconds()
        meascov=local@np.array(m['centroid_covariance_pixel2'])@local.T+np.eye(2)*max(cal['loo_rms_per_axis_arcsec'],.5)**2+np.outer(rate,rate)*control['time_sigma_seconds']**2
        total=meascov+np.array(pred[1]['covariance_east_north_arcsec2'])
        nearest=min(np.linalg.norm(offsets(sky,[s[1],s[2]])) for s in cat)
        m.update(id=control['id'],control_kind='documented-positive',region=control['region'],image_file=control['file'],
                 measured_ra_deg=float(sky[0]),measured_dec_deg=float(sky[1]),
                 observation_time=mid.isoformat(),timescale='UT1',observation_time_origin='MPEC trail-end midpoint override',
                 header_date_obs=header.get('DATE-OBS'),exposure_seconds=(end-start).total_seconds(),
                 residual_east_north_arcsec=res.tolist(),residual_arcsec=float(np.linalg.norm(res)),
                 measurement_covariance_east_north_arcsec2=meascov.tolist(),
                 combined_mahalanobis_squared=float(res@np.linalg.solve(total,res)),
                 nearest_stationary_catalogue_arcsec=float(nearest),stationary_match=bool(nearest<4),
                 epoch_matches=True,time_provenance=True,
                 inside_coverage=bool(np.all(pix[:,0]>10)&np.all(pix[:,0]<data.shape[1]-10)&np.all(pix[:,1]>10)&np.all(pix[:,1]<data.shape[0]-10)),
                 calibration_valid=bool(cal['inliers']>=8 and cal['loo_rms_per_axis_arcsec']<2),
                 artifact_flag=False,edge_flag=False)
        m['verdict']=screen(m);measurements.append(m);entries.append((m,arrays))
        contexts.append((data,w,coeff,cal,scale,pix,m))
        # Save an actual FITS measurement cutout; preserve parent WCS in CNPIX.
        x0,x1,y0,y1=m['pixel_bounds'];h=header.copy();h['CNPIX1']+=x0;h['CNPIX2']+=y0
        fits.PrimaryHDU(data[y0:y1,x0:x1],header=h).writeto(output/(control['id']+'.fits'),overwrite=True)
    # Real stationary source, selected deterministically from calibrated stars.
    data,w,coeff,cal,scale,pix,base=contexts[1]
    usable=[s for s in cal['stars'] if s['used'] and 60<s['x']<data.shape[1]-60 and 60<s['y']<data.shape[0]-60]
    star=sorted(usable,key=lambda s:s['source_id'])[0]
    for kind,center,array in [('stationary-gaia-star',np.array([star['x'],star['y']]),data),
                              ('wrong-match-offset-60arcsec',pix[1]+np.array([0,60/scale]),data)]:
        m,arrays=fit_trail(array,center,pix[2]-pix[0])
        sky=measured_sky([m['x'],m['y']],w,coeff)
        nearest=min(np.linalg.norm(offsets(sky,[s[1],s[2]])) for s in cat)
        residual=float(np.linalg.norm(offsets(sky,[predictions[4]['ra_deg'],predictions[4]['dec_deg']])))
        m.update(id=kind,control_kind='real-pixel-negative',image_file='data/raw/XO036.fits',reference_source_id=star['source_id'] if kind=='stationary-gaia-star' else None,residual_arcsec=residual,
                 stationary_match=bool(nearest<4),nearest_stationary_catalogue_arcsec=float(nearest),
                 artifact_flag=False,edge_flag=False,epoch_matches=True,time_provenance=True,inside_coverage=True,calibration_valid=True)
        m['verdict']=screen(m);measurements.append(m);entries.append((m,arrays))
        x0,x1,y0,y1=m['pixel_bounds'];h=fits.getheader(ROOT/'data/raw/XO036.fits');h['CNPIX1']+=x0;h['CNPIX2']+=y0
        fits.PrimaryHDU(data[y0:y1,x0:x1],header=h).writeto(output/(kind+'.fits'),overwrite=True)
    # Synthetic artifact is explicitly distinguished from a real archive defect.
    artifact=np.full_like(data,float(np.median(data)),dtype=float)
    rng=np.random.default_rng(137108);artifact+=rng.normal(0,base['noise'],artifact.shape)
    ix,iy=np.round(pix[1]).astype(int);artifact[iy,ix]+=100*base['noise']
    m,arrays=fit_trail(artifact,pix[1],pix[2]-pix[0]);m.update(id='injected-single-pixel-artifact',control_kind='synthetic-negative',
        residual_arcsec=float(np.hypot(m['x']-pix[1,0],m['y']-pix[1,1])*scale),
        artifact_flag=bool((artifact.max()-np.median(artifact))>30*base['noise']),stationary_match=False,edge_flag=False,
        epoch_matches=True,time_provenance=True,inside_coverage=True,calibration_valid=True)
    m['verdict']=screen(m);measurements.append(m);entries.append((m,arrays))
    fits.PrimaryHDU(artifact).writeto(output/'synthetic-artifact.fits',overwrite=True)
    wrong=fits.getheader(ROOT/'data/raw/irsa-XE035.fits')
    wrong_m=dict(id='wrong-epoch-XE035',control_kind='real-metadata-negative',
                     file='data/raw/irsa-XE035.fits', measurement_performed=False,
                     header_date_obs=wrong['DATE-OBS'],epoch_matches=wrong['DATE-OBS'].startswith('1955-01-26'),time_provenance=False)
    wrong_m['verdict']=screen(wrong_m);measurements.append(wrong_m)
    # Published endpoint consistency is independent of pixel centroid fitting.
    reference=[]
    for c in controls['positive_controls']:
        for time,sky in zip([c['start'],c['end']],c['reference_endpoints_deg']):
            p=prediction(strict_time(time),eph);r=offsets(sky,[p['ra_deg'],p['dec_deg']])
            reference.append(dict(control=c['id'],time=time,published_ra_dec_deg=sky,residual_arcsec=float(np.linalg.norm(r))))
    observed_motion=offsets([measurements[1]['measured_ra_deg'],measurements[1]['measured_dec_deg']],
                            [measurements[0]['measured_ra_deg'],measurements[0]['measured_dec_deg']])
    predicted_motion=offsets([predictions[4]['ra_deg'],predictions[4]['dec_deg']],
                             [predictions[1]['ra_deg'],predictions[1]['dec_deg']])
    motion_error=float(np.linalg.norm(observed_motion-predicted_motion))
    pair=dict(observed_displacement_east_north_arcsec=observed_motion.tolist(),
              predicted_displacement_east_north_arcsec=predicted_motion.tolist(),
              displacement_error_arcsec=motion_error, consistent=motion_error<12,
              independent_exposures=2, stationary_hypothesis_displacement_arcsec=0)
    dump(output/'pair-motion.json',pair)
    metrics=dict(documented_exposures=2,published_endpoints=4,pair_motion_consistent=pair['consistent'],
                 positives_passing=sum(m['verdict']['screen_pass'] for m in measurements[:2]),
                 negatives_rejected=sum(not m['verdict']['screen_pass'] for m in measurements[2:]),
                 historical_arc_start_extension_days=(strict_time(controls['discovery_date'])-strict_time(controls['positive_controls'][0]['start'])).total_seconds()/86400,
                 historical_extension_precision='Discovery date known to day; midnight convention used, so interpret extension approximately.',
                 new_arc_extension_days=0,orbital_uncertainty_improvement=None,
                 orbital_improvement_status='Not estimated: requires consistent orbit refits with/without precovery, same weighting and epoch; current JPL solution already includes these data.',
                 scientific_acceptance=False,mpc_submission=False)
    dump(output/'predictions.json',predictions);dump(output/'calibration.json',calibrations)
    dump(output/'measurements.json',measurements);dump(output/'reference-comparison.json',reference);dump(output/'metrics.json',metrics)
    plot_evidence(output,entries)
    print(json.dumps(metrics,indent=2));return metrics


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);parser.add_argument('--bootstrap',action='store_true')
    args=parser.parse_args();run(args.output,verify=not args.bootstrap)
