import tempfile
import unittest
from pathlib import Path
from precovery.verify import screen, strict_time, check_manifest


class VerificationTests(unittest.TestCase):
    def good(self):
        return dict(snr=20., residual_arcsec=1., length_ratio=1., angle_error_deg=0.,
                    delta_bic_vs_point=50., fit_converged=True, epoch_matches=True, time_provenance=True,
                    inside_coverage=True, calibration_valid=True, stationary_match=False,
                    artifact_flag=False, edge_flag=False)

    def test_screen_never_accepts_scientifically(self):
        r = screen(self.good())
        self.assertTrue(r['screen_pass'])
        self.assertFalse(r['scientifically_accepted'])
        self.assertFalse(r['submission_allowed'])

    def test_alternatives_override_perfect_residual(self):
        for key in ['stationary_match', 'artifact_flag', 'edge_flag']:
            m = self.good(); m['residual_arcsec'] = 0.; m[key] = True
            self.assertFalse(screen(m)['screen_pass'], key)

    def test_missing_nonfinite_wrong_epoch_and_mismatch(self):
        for key, value in [('snr', float('nan')), ('residual_arcsec', 50.),
                           ('epoch_matches', False), ('time_provenance', False),
                           ('calibration_valid', None), ('fit_converged', False), ('length_ratio', 0.01)]:
            m = self.good(); m[key] = value
            self.assertFalse(screen(m)['screen_pass'], key)

    def test_invalid_archive_time_rejected(self):
        with self.assertRaises(ValueError): strict_time('1955-01-25T09:65:00')

    def test_corruption_rejected(self):
        with tempfile.TemporaryDirectory(dir='.') as d:
            Path(d, 'image').write_text('changed')
            with self.assertRaisesRegex(ValueError, 'Checksum'):
                check_manifest(d, {'files':[{'path':'image','sha256':'0'*64}]})


if __name__ == '__main__': unittest.main()
