import unittest
import numpy as np
from precovery.pipeline import offsets, ephemeris, prediction, fit_trail
from datetime import datetime


class ScientificInvariants(unittest.TestCase):
    def test_ra_wrap_and_polar_projection(self):
        a=offsets([.001,60],[359.999,60])
        np.testing.assert_allclose(a,[3.6,0],atol=1e-7)

    def test_covariance_and_time_bounds(self):
        e=ephemeris();p=prediction(datetime(1955,1,26,10),e)
        eig=np.linalg.eigvalsh(p['covariance_east_north_arcsec2'])
        np.testing.assert_allclose(eig,sorted([(p['smaa_3sigma_arcsec']/3)**2,(p['smia_3sigma_arcsec']/3)**2]))
        with self.assertRaises(ValueError):prediction(datetime(1955,1,25),e)

    def test_known_injected_trail_measured_not_echoed(self):
        from scipy.special import erf
        rng=np.random.default_rng(9);y,x=np.mgrid[:100,:150]
        # True center is deliberately displaced from the predicted initial center.
        signal=20*.5*(erf((x-76+25)/2)-erf((x-76-25)/2))*np.exp(-(y-48)**2/8)
        data=100+signal+rng.normal(0,1,signal.shape)
        m,_=fit_trail(data,[75,50],[50,0])
        self.assertLess(abs(m['x']-76),.3);self.assertLess(abs(m['y']-48),.3)
        self.assertGreater(m['delta_bic_vs_point'],10)


if __name__=='__main__':unittest.main()

class ArchivePixelAgreement(unittest.TestCase):
    def test_byte_range_and_archive_service_cutouts_agree(self):
        from astropy.io import fits
        from precovery.acquire import ROOT
        for region in ['XE036','XO036']:
            a,ha=fits.getdata(ROOT/f'data/raw/{region}.fits',header=True)
            b,hb=fits.getdata(ROOT/f'data/raw/irsa-{region}.fits',header=True)
            dx=int(hb['CNPIX1']-ha['CNPIX1']);dy=int(hb['CNPIX2']-ha['CNPIX2'])
            x0=max(0,dx);y0=max(0,dy);x1=min(a.shape[1],dx+b.shape[1]);y1=min(a.shape[0],dy+b.shape[0])
            self.assertGreater((x1-x0)*(y1-y0),100000)
            np.testing.assert_array_equal(a[y0:y1,x0:x1],b[y0-dy:y1-dy,x0-dx:x1-dx])
