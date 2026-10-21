import copy,json,unittest
from review_package import HERE,check
class ReviewTests(unittest.TestCase):
    def setUp(self):self.r=json.loads((HERE/'review-request.json').read_text())
    def test_pending_is_technical_only(self):self.assertFalse(check(self.r)['scientific_acceptance'])
    def test_changed_evidence(self):
        self.r['evidence'][0]['sha256']='0'*64
        with self.assertRaises(ValueError):check(self.r)
    def test_no_anonymous_decision(self):
        self.r['verdict']='supported-for-retrospective-use'
        with self.assertRaises(ValueError):check(self.r)
    def test_no_automatic_acceptance(self):
        self.r['scientifically_accepted']=True
        with self.assertRaises(ValueError):check(self.r)
    def test_no_silent_time_adjustment(self):
        self.r['observations']['replication_start']='1955-01-25T10:05:00'
        with self.assertRaises(ValueError):check(self.r)
if __name__=='__main__':unittest.main()
